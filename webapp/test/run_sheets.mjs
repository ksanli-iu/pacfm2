import { readFileSync, writeFileSync } from "fs";
import { buildKegg, sharedKoSheets } from "./core.mjs";
const ref = JSON.parse(readFileSync(process.argv[2], "utf8"));
const k = buildKegg(ref.kp, ref.h);
const prof = new Map(readFileSync(process.argv[3], "utf8").trim().split("\n").slice(1).map(l => { const [a, b] = l.split("\t"); return [a, Number(b)]; }));
writeFileSync(process.argv[4], JSON.stringify(sharedKoSheets(k, prof, 5)));
