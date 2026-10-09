import { readFileSync, writeFileSync } from "fs";
import { buildKegg, pathwayAbundance } from "./core.mjs";
const ref = JSON.parse(readFileSync(process.argv[2], "utf8"));
const k = buildKegg(ref.kp, ref.h);
const prof = new Map(readFileSync(process.argv[3], "utf8").trim().split("\n").slice(1).map(l => { const [a, b] = l.split("\t"); return [a, Number(b)]; }));
const out = [];
for (const s of ["naive", "fractional", "em", "fraccov", "emcov"]) for (const [p, v] of pathwayAbundance(k, prof, s)) out.push(`${s}\t${p}\t${v}`);
for (const [p, v] of pathwayAbundance(k, prof, "cutoff", 2)) out.push(`cutoff2\t${p}\t${v}`);
writeFileSync(process.argv[4], "scheme\tpathway\tvalue\n" + out.join("\n"));
