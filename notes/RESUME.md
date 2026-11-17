# Resume point (paused 2026-11-16)

**Status 2026-11-17: COMPLETE.**
- **Report:** REPORT.md and REPORT.html are final.
- **Fact-check:** 43 issues were found and all were fixed (`data/derived/report_factcheck.json`).
- **Possible follow-ups:**
  - analyse the archival COS G130M spectrum of QSO 2MASS J07452815+1919525;
  - re-run the DESI analysis with DESI sky fibres if NERSC becomes reachable;
  - request the SWONS raw frames to use as proper-motion epoch 1.

**Update 2026-11-17 (superseded):** waves 3 and 4 are complete (`notes/w3_*`, `notes/w4_*`); the report was drafted and fact-checked.

**Update 2026-11-16 evening:**
- **Wave 2:** completed. Results are in `data/derived/wave2_agent_results.json` and `notes/verify_*.md`.
- **Wave 3:** running as run `wf_307fa6f3-97b`. It has two stages:
  - Resolve: `oii`, `northsouth`, `canon`, `figs`.
  - Judge: `judge_obs`, `judge_phys`, `judge_skeptic`.
  - Outputs go to `notes/w3_*.md`, `figures/w3_figs_*.png` and `data/derived/w3_canon_values.json`.
- **Remaining:** REPORT.md and REPORT.html, then a fact-check.

If interrupted, resume with `Workflow({scriptPath: ".../workflows/scripts/silk-nebula-synthesis-wave3-wf_307fa6f3-97b.js", resumeFromRunId: "wf_307fa6f3-97b"})`. The script lives in the session directory.

## Done
Wave 1 is complete: nine research agents. Outputs:
- notes: `notes/{discovery,paper_review,geometry,context,imaging_optical,imaging_mw,spectra_desi,spectra_sdss_lamost,spectra_ism}.md`
- scripts, `data/raw`, `data/derived` (including `silk_trace.json` v1.0) and 54 figures in `figures/`

The structured agent results are saved in `data/derived/wave1_agent_results.json`.

## Stopped
Wave 2 (run `wf_cc6927c1-6e0`) was stopped before any agent finished, so it has no cached results and all ten agents re-run on resume. It has two parts:
- **Verify**, eight adversarial checks:
  - `paper`: distance units, [S II] density, energetics.
  - `velocity`: LAMOST [O III] +25 versus DESI/SDSS +13 km/s.
  - `desi`: independent re-analysis, including Hg 4358 contamination of [O III] 4363 and fibre-loss calibration.
  - `sdss`: independent re-analysis, including [O II]3727, north/south ratios and precovery MJDs.
  - `hiwham`: HI and WHAM at the revised V_LSR ≈ 0.
  - `xrayfuv`: Monogem X-ray edge null test, FUV/NUV predictions.
  - `precovery`: SDSS imaging-sky reanalysis and a dated precovery table.
  - `distance`: the HD 61717 −17 km/s component, IUE spectra of HD 60848, and a distance synthesis.
- **Model**, two agents:
  - `shockmodel`: complete versus truncated radiative shocks.
  - `altmodel`: conduction front, mixing layer, non-equilibrium cooling and photoionization.

Each writes only `*verify_<tag>*` files.

To resume in this session:
`Workflow({scriptPath: "/Users/mickey/.claude/projects/-Users-mickey-astronomy-silk-nebula/c7d7dd5c-f92f-4bd3-b591-acbd37ba6cef/workflows/scripts/silk-nebula-verify-model-wave2-wf_cc6927c1-6e0.js"})`
If the session-local script is unavailable, re-author wave 2 from the summary above. The wave-1 state-of-knowledge is in the notes.

## Remaining after wave 2
1. Synthesis and judge panel over the hypotheses, plus a completeness critic.
2. `REPORT.md` and a self-contained local `REPORT.html` with figures, opened in the browser, and a README.

## Provisional wave-1 headline findings (unverified)
- **Paper's distance:** the Table 4 distances are in kpc, so 160 pc should read ~1.6 kpc. Kinematic distances are meaningless here in any case.
- **[O III] velocity:** V_helio +12.5 to +14.5 km/s in DESI and SDSS, against the paper's +24 to +27.
- **[O III] brightness:** about 1–4 R, not the paper's 10–13 R.
- **Line ratios:** [O III]/Hβ ≳ 14–20 and [O II]/[O III] ≈ 1.1–1.5. The filament is Balmer-dark and collisionally excited.
- **Red-arm ratios** in the paper match the warm ionized medium (WIM) and also appear in the nebula-free background region P4.
- **Monogem Ring:** the filament lies on the outer edge of Monogem's soft X-rays, with an HI wall just outside. The likely distance is 250–320 pc.
- **Precovery:** SDSS imaging sky levels from 2002–2007, detected at about 4σ.
