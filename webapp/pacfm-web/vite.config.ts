import { defineConfig } from "vite";

// `npm run build`       standalone static site in dist/ (any static host; kegg_reference.json + demo copied along)
// `npm run build:site`  the lab site's copy: pacfm.js + pacfm.css with fixed names in the Django static folder
//                       (PACFM_SITE_OUT, default the helixome ecobiom backend); the site serves the data files itself
const SITE_OUT = process.env.PACFM_SITE_OUT ?? "dist-site";

export default defineConfig(({ mode }) => mode === "site"
  ? {
      base: "./",
      publicDir: false,
      build: {
        outDir: SITE_OUT, emptyOutDir: true,
        rolldownOptions: {
          input: "src/main.ts",
          output: { entryFileNames: "pacfm.js", assetFileNames: "pacfm[extname]" },
        },
      },
    }
  : { base: "./" });
