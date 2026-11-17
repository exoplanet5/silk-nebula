# Silk Nebula study

Read [REPORT.md](REPORT.md) (English) or [REPORT_zh.md](REPORT_zh.md) (简体中文). The self-contained, bilingual page [REPORT.html](REPORT.html) embeds all figures and has an EN | 中文 switch at the top right. The page remembers your choice, and `REPORT.html#zh` or `#en` opens a given language. Figures carry English labels in both versions; captions are translated. The report covers:
- the discovery history and the characteristics of the [O III] arc;
- archival imaging precovery, back to SDSS in 2002;
- a claim-by-claim review of Wu et al. 2026 (arXiv:2610.13329, LAMOST MRS-N), whose "possible SNR" verdict rests on:
  - a kpc/100-pc unit slip;
  - an air-airglow velocity artefact;
  - a 4–8× brightness overestimate;
- new DESI DR1 and SDSS/BOSS spectroscopy of the filament's own lines;
- a revised hypothesis: an incompletely cooled radiative shock seen edge-on, most plausibly on the outer limb of the Monogem Ring at ~250–320 pc.

## Comment on Wu et al. 2026

[DISPUTE.md](DISPUTE.md) and [DISPUTE.html](DISPUTE.html) form a formal comment on the paper. Its Summary is bilingual (English and 中文, with an EN | 中文 switch in the HTML); the rest is in English. It covers six disputed claims, each with dedicated figures (Figures C1–C6), tables (Tables C1–C7), requested corrections and a verified reference list. The claims are:
- the Table 4 kpc/100-pc unit misreading;
- why no kinematic distance is possible;
- the [O III] wavelength-reference artefact;
- the Hα velocity bias;
- the Table 3 brightness overestimate;
- the remaining diagnostics.

The section drafts are in `notes/dispute/`. Rebuild with:

```bash
~/.venvs/astro313/bin/python scripts/dispute_assemble.py
~/.venvs/astro313/bin/python scripts/make_report_html.py DISPUTE
```

## Where things are

| Path | Contents |
|---|---|
| `notes/` | Per-topic working notes, in four waves: `discovery`…`spectra_ism` (wave 1), `verify_*` (wave 2), `w3_*` (wave 3), `w4_*` (wave 4). The canonical values are in `w3_canon.md`. Later waves supersede earlier ones |
| `figures/` | All figures. The report uses `w3_figs_f1…f7` |
| `scripts/` | Analysis code. Shared helpers: `silk_geometry.py`, `silk_style.py` |
| `data/raw/` | Downloads, one folder per source (**local only, not in git**) |
| `data/derived/` | Measurements, `silk_trace.json`, `w3_canon_values.json`, per-wave agent results, `report_factcheck.json` (**local only, not in git**) |
| `sources/` | Web captures. The paper itself (`sources/arxiv/`) is **local only, not in git**; see below |

## Rebuild the HTML report

```bash
cd ~/astronomy/silk-nebula
~/.venvs/astro313/bin/python scripts/make_report_html.py   # REPORT.md + REPORT_zh.md -> REPORT.html
```

Run the individual analyses with `~/.venvs/astro313/bin/python -I scripts/<name>.py`. Each note lists its own scripts and inputs.

## Access limits encountered (2026-11-16/17)

- NERSC hosts (`data.desi.lbl.gov`, `legacysurvey.org`) were unreachable through the local proxy. As a result:
  - DESI spectra came from NOIRLab SPARCL;
  - DESI sky fibres and per-petal sky models were not used;
  - Legacy Surveys images were not used.
- LAMOST MRS-N nebular fibre spectra are not public.

## Repository and local data

`data/` is kept out of git (see `.gitignore`). It holds ~7 GB of survey downloads and ~0.6 GB of derived products. Each note in `notes/` names the scripts and inputs it used.
- **Re-downloading:** the fetch scripts in `scripts/` (`*_fetch*.py`, `*_query*.py`) re-download the public inputs, from DESI DR1 via SPARCL, SDSS, GALEX, eROSITA DR1, WHAM, GALFA-HI, NSNS and others.
- **Not reproducible from the repository alone:** the results in `data/derived/`. Rebuilding them means re-running the analysis scripts on those downloads.

The paper is also kept out of git, because arXiv's default licence does not allow redistribution. Several scripts expect it at `sources/arxiv/`, so to reproduce the paper-based checks, fetch it there:

```bash
mkdir -p sources/arxiv && cd sources/arxiv
curl -L -o 2610.13329v1.pdf https://arxiv.org/pdf/2610.13329v1
pdftotext -layout 2610.13329v1.pdf 2610.13329v1.txt
```

The page renders and figure crops in `sources/arxiv/figs/` were made by hand from 300-dpi `pdftoppm` renders of the PDF; no script regenerates them. To recreate one (e.g. `fig5b_crop.png`, the crop of the paper's Fig. 5b that Section 1 of the comment digitizes), render the page with `pdftoppm -r 300 -png` and crop it. `scripts/verify_velocity_fig3_vectors.py` reads the Fig. 3 vectors straight from the PDF.

