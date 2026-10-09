// PACFM 2 core — specificity-aware pathway profiling (TypeScript port of pacfm2/core.py; results must match).
//
// Schemes: naive (sum of all pathway KOs), cutoff (KOs in <= c pathways; c = 1: unique KOs only),
// fractional (KO abundance split equally among its pathways), em (split in proportion to the estimated per-KO
// abundance of each pathway, iterated), fraccov (split in proportion to each pathway's KO coverage ** gamma),
// emcov (em weighted by coverage ** gamma). Global/overview maps (011xx, 012xx) are excluded.
const GLOBAL = /^01[12]\d\d$/;
export function buildKegg(koPathway, hierarchy) {
    const info = new Map(hierarchy.map(([p, a, b, n]) => [p, { levelA: a, levelB: b, name: n }]));
    const pathwayKos = new Map();
    const koPathways = new Map();
    const seen = new Set();
    for (const [ko, p] of koPathway) {
        if (GLOBAL.test(p) || !info.has(p))
            continue;
        const key = ko + "|" + p;
        if (seen.has(key))
            continue;
        seen.add(key);
        (pathwayKos.get(p) ?? pathwayKos.set(p, []).get(p)).push(ko);
        (koPathways.get(ko) ?? koPathways.set(ko, []).get(ko)).push(p);
    }
    return { pathwayKos, koPathways, info };
}
export const nAssoc = (k, ko) => k.koPathways.get(ko)?.length ?? 0;
/** abundances: KO -> value for one sample. Returns pathway -> abundance. */
export function pathwayAbundance(k, abund, scheme, cutoff = 1, maxIter = 200, tol = 1e-8, gamma = 2) {
    const out = new Map();
    const kos = [...abund.keys()].filter(ko => k.koPathways.has(ko));
    if (scheme === "naive" || scheme === "cutoff" || scheme === "fractional") {
        for (const ko of kos) {
            const ps = k.koPathways.get(ko);
            if (scheme === "cutoff" && ps.length > cutoff)
                continue;
            const v = abund.get(ko) / (scheme === "fractional" ? ps.length : 1);
            for (const p of ps)
                out.set(p, (out.get(p) ?? 0) + v);
        }
        return out;
    }
    const paths = [...new Set(kos.flatMap(ko => k.koPathways.get(ko)))];
    const size = new Map(paths.map(p => [p, k.pathwayKos.get(p).length]));
    // KO coverage of each pathway in this sample, raised to gamma (fraccov, emcov); 1 for em
    const prior = new Map(paths.map(p => [p, 0]));
    for (const ko of kos)
        if (abund.get(ko) > 0)
            for (const p of k.koPathways.get(ko))
                prior.set(p, prior.get(p) + 1);
    for (const p of paths)
        prior.set(p, scheme === "em" ? 1 : (prior.get(p) / size.get(p)) ** gamma);
    const allocate = (weight) => {
        const total = new Map(paths.map(p => [p, 0]));
        for (const ko of kos) {
            const ps = k.koPathways.get(ko);
            const den = ps.reduce((s, p) => s + weight(p), 0);
            if (den <= 0)
                continue;
            for (const p of ps)
                total.set(p, total.get(p) + abund.get(ko) * weight(p) / den);
        }
        return total;
    };
    if (scheme === "fraccov")
        return allocate(p => prior.get(p));
    // em, emcov
    let theta = new Map();
    for (const ko of kos)
        for (const p of k.koPathways.get(ko))
            theta.set(p, (theta.get(p) ?? 0) + abund.get(ko));
    for (const p of paths)
        theta.set(p, theta.get(p) / size.get(p));
    let total = new Map();
    for (let it = 0; it < maxIter; it++) {
        const th = theta;
        total = allocate(p => th.get(p) * prior.get(p));
        let diff = 0, mx = 0;
        const next = new Map();
        for (const p of paths) {
            const v = total.get(p) / size.get(p);
            diff = Math.max(diff, Math.abs(v - theta.get(p)));
            mx = Math.max(mx, theta.get(p));
            next.set(p, v);
        }
        theta = next;
        if (diff <= tol * Math.max(1, mx))
            break;
    }
    return total;
}
/** per pathway: coverage (detected/all KOs), uniqueCoverage, uniqueShare of naive abundance */
export function pathwayEvidence(k, abund) {
    const res = new Map();
    for (const [p, kos] of k.pathwayKos) {
        let det = 0, udet = 0, ucount = 0, tot = 0, utot = 0;
        for (const ko of kos) {
            const a = abund.get(ko) ?? 0;
            const uniq = nAssoc(k, ko) === 1;
            if (uniq)
                ucount++;
            if (a > 0) {
                det++;
                if (uniq)
                    udet++;
            }
            tot += a;
            if (uniq)
                utot += a;
        }
        if (tot > 0)
            res.set(p, { coverage: det / kos.length, uniqueCoverage: ucount ? udet / ucount : NaN, uniqueShare: utot / tot });
    }
    return res;
}
/** pathway pairs that share detected KOs, with the shared abundance (for the bundled links) */
export function sharedLinks(k, abund, maxAssoc = 10) {
    const links = new Map();
    for (const [ko, a] of abund) {
        const ps = k.koPathways.get(ko);
        if (!ps || ps.length < 2 || ps.length > maxAssoc || a <= 0)
            continue;
        for (let i = 0; i < ps.length; i++)
            for (let j = i + 1; j < ps.length; j++) {
                const key = ps[i] < ps[j] ? ps[i] + "|" + ps[j] : ps[j] + "|" + ps[i];
                links.set(key, (links.get(key) ?? 0) + a / (ps.length - 1));
            }
    }
    return links;
}
/** pathway pairs sharing at least minShared detected KOs (abundance > 0), for the sheets of the radial plot.
 *  Each pair is oriented [a, b, n, fracA, fracB] so that a is the pathway with the larger share of its KOs in the
 *  overlap (fraction = n / KOs of the pathway); the sheet takes a's colour. At most maxPairs pairs, largest first. */
export function sharedKoSheets(k, abund, minShared = 5, maxPairs = 800) {
    const n = new Map();
    for (const [ko, a] of abund) {
        const ps = k.koPathways.get(ko);
        if (!ps || ps.length < 2 || a <= 0)
            continue;
        for (let i = 0; i < ps.length; i++)
            for (let j = i + 1; j < ps.length; j++) {
                const key = ps[i] < ps[j] ? ps[i] + "|" + ps[j] : ps[j] + "|" + ps[i];
                n.set(key, (n.get(key) ?? 0) + 1);
            }
    }
    const out = [];
    for (const [key, c] of n) {
        if (c < minShared)
            continue;
        const [x, y] = key.split("|");
        const fx = c / k.pathwayKos.get(x).length, fy = c / k.pathwayKos.get(y).length;
        out.push(fx >= fy ? [x, y, c, fx, fy] : [y, x, c, fy, fx]);
    }
    return out.sort((p, q) => q[2] - p[2] || (p[0] + p[1]).localeCompare(q[0] + q[1])).slice(0, maxPairs);
}
// ---------------------------------------------------------------- statistics (group comparison in the browser)
function rank(v) {
    const idx = v.map((x, i) => [x, i]).sort((a, b) => a[0] - b[0]);
    const r = new Array(v.length);
    for (let i = 0; i < idx.length;) {
        let j = i;
        while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0])
            j++;
        for (let t = i; t <= j; t++)
            r[idx[t][1]] = (i + j) / 2 + 1;
        i = j + 1;
    }
    return r;
}
function normSf(z) {
    const t = 1 / (1 + 0.2316419 * Math.abs(z));
    const d = 0.3989423 * Math.exp(-z * z / 2);
    const p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return z > 0 ? p : 1 - p;
}
/** two-sided Mann-Whitney U test, normal approximation with tie correction */
export function mannWhitney(a, b) {
    const n1 = a.length, n2 = b.length;
    if (n1 === 0 || n2 === 0)
        return 1;
    const r = rank([...a, ...b]);
    const r1 = r.slice(0, n1).reduce((s, x) => s + x, 0);
    const u = r1 - n1 * (n1 + 1) / 2;
    const n = n1 + n2;
    const counts = new Map();
    for (const x of [...a, ...b])
        counts.set(x, (counts.get(x) ?? 0) + 1);
    const tie = [...counts.values()].reduce((s, t) => s + t ** 3 - t, 0);
    const sd = Math.sqrt(n1 * n2 / 12 * ((n + 1) - tie / (n * (n - 1))));
    if (sd === 0)
        return 1;
    const z = (Math.abs(u - n1 * n2 / 2) - 0.5) / sd;
    return Math.min(1, 2 * normSf(z));
}
export function benjaminiHochberg(p) {
    const n = p.length, order = p.map((x, i) => [x, i]).sort((a, b) => a[0] - b[0]);
    const q = new Array(n);
    let prev = 1;
    for (let k = n - 1; k >= 0; k--) {
        const [x, i] = order[k];
        prev = Math.min(prev, x * n / (k + 1));
        q[i] = prev;
    }
    return q;
}
