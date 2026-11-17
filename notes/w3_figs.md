# Wave 3 · composite figures for the final report (tag `figs`)

Written 2026-11-16/17. The figures follow the `scripts/silk_style.py` dark theme: Menlo font, concise title at the upper left, small credit "Illustration: Mickey". All are 200 dpi PNGs, and each was inspected and iterated until labels were clean.

Numbers come from the wave-2 notes (`notes/verify_*.md`) and their data products. Superseded wave-1 values are not used, with two exceptions, both labelled "wave 1" on the figure:
- the LAMOST LRS stellar-residual brightness, which nobody re-tested;
- the paper's MRS-N spectra rescaled by sky lines.

Where this agent computed something new (profiles, a DESI Hα window, SDSS stacks, WHAM/HI heliocentric conversions), it uses the wave-2 code or method unchanged. Those results are marked **derived here**.

Run all scripts with `~/.venvs/astro313/bin/python -I scripts/<name>.py`. Data-prep scripts run before the plotting scripts:
- F1: `w3_figs_f1_data.py`
- F4: `w3_figs_desi_regspec_ha.py` and `w3_figs_sdss_stacks.py`
- F6: `w3_figs_f6_profiles_data.py`

---

## F1 · `figures/w3_figs_f1_overview.png` (overview map)

**Caption.** NSNS DR0.2 [O III] 5007 (asinh stretch, 0–4 R) over RA 108–124°, Dec 8–30°. Overlaid:
- the Silk E (spine), W and W2 strands (dotted) and the Hα patch (red, half-maximum contour) from `silk_trace.json` v1.0;
- the Wu et al. P1–P4 regions and the 2.5°-radius LAMOST MRS-N plate;
- eRASS1 0.2–0.6 keV contours at 6, 7 and 8 cts s⁻¹ deg⁻²;
- the directions to PSR B0656+14 / the Monogem two-SNR centre (11.4° / 12.1° from the southern spine) and to Geminga (17.6°);
- the SDSS/BOSS sky fibres (yellow dots) and DESI DR1 coadds (violet dots) used by the wave-2 spectroscopic verifications.

The southern Silk lies on the outer (eastern) flank of the Monogem soft-X-ray shell. The contours do not single out the spine.

**Sources:**
- `data/raw/nsns_hips/oiii/Norder3` (NSNS DR0.2 HiPS, CC-BY-NC-SA, S. Ziegenbalg), mosaicked to a 1′ CAR grid → `data/derived/w3_figs_f1_maps.npz`.
- `verify_xrayfuv_erass_pixels.npz`: eRASS1 DR1 with the verify_xrayfuv point-source mask; 0.1° bins, smoothed with σ = 0.35° plus 0.15°.
- `silk_trace.json`, `verify_sdss_fibres.csv`, `verify_desi_coadd_meta.npz` (`ok_contrib`).
- Separations: `notes/verify_distance.md` (11.36 ± 0.15°, 12.10 ± 0.12°) and the `silk_trace.json` candidate centres (Geminga 17.56°).

**Caveat.** The X-ray contours are heavily smoothed and are only a guide to the shell geometry. Quantitative edge statistics are in `verify_xrayfuv`.

## F2 · `figures/w3_figs_f2_velocity.png` (velocity reconciliation)

**Caption.** Heliocentric velocities toward the Silk. Paper Table 4 [O III] (P1–P3, +23.8…+26.8 km s⁻¹, red) sits ~10 km s⁻¹ above every other measurement. The −10.1 ± 1.0 km s⁻¹ correction (air airglow re-referencing plus the missing −8.3 km s⁻¹ v_corr) brings it to +13.7…+16.7. The paper's own Fig. 3a, measured against its [N I] airglow, gives +13.6 ± 1.2.

DESI and SDSS [O III], the LAMOST [N II]/[S II] red lines, WHAM Hα and the local HI all cluster at V_helio ≈ +10…+15 km s⁻¹. That is the adopted [O III] value, +13.3 ± 1.0 (V_LSR ≈ +1). LAMOST Hα is offset by +4–7 km s⁻¹, which geocoronal over-subtraction reproduces. The northern [O III] is ~10 km s⁻¹ redder in DESI (+20.9) and in SDSS (+16.4 ± 5.9).

**Sources:**
- Paper Table 4: `sources/arxiv/2610.13329v1.txt`.
- Correction and Fig. 3a: `verify_velocity_summary.json`.
- DESI: `notes/verify_desi.md` rev. 2 (all +11.3 [10.2, 13.9]; Dec < 19.5: +10.6 [10.2, 12.0]; Dec ≥ 19.5: +20.9 [18.5, 23.1]) and the verify_velocity ridge stack (+14.4 ± 1.9).
- SDSS: `verify_sdss_results_base_r0.15.json` (south, Dec < 20: +10.0 ± 1.8 ⊕ 0.8 zero point; north +16.4 ± 5.9).
- Adopted band: `verify_velocity_summary.json` combined value.

**Derived here** (`data/derived/w3_figs_f2_velocity.json`):
- WHAM main component, heliocentric, from the verify_hiwham per-beam Gaussian fits (`verify_hiwham_wham_beams.csv`; V_helio = v_g − v_c,LSR + 13.5):
  - south on-filament: +13.8 ± 0.8 km s⁻¹ (beam std, 6 beams; V_LSR +1.28);
  - north: +10.3 ± 1.3 (16 beams).
- GALFA on-spine HI peak, Dec 12–22, from `verify_hiwham_hi_extra.json`: median V_helio +14.2 (16–84 %: 13.1–15.5); V_LSR +1.8.
- `verify_velocity.md` quotes WHAM as "≈ +12". The per-beam conversion here gives +13.8 for the southern on-filament beams. Both lie inside the cluster.

**Caveat.** The top V_LSR axis uses one −12.2 km s⁻¹ offset. The true offset varies from −13.2 (Dec 12) to −9.3 km s⁻¹ (Dec 28) along the spine. The DESI north/south split is at Dec 19.5 and the SDSS split at Dec 20.

## F3 · `figures/w3_figs_f3_brightness.png` (brightness reconciliation)

**Caption.** [O III] 5007 peak or ridge intensity on a log scale, with the angular scale of each measurement at the right. Every flux-calibrated measurement of the southern ridge lies at about 1.4–2.5 R (green band). That covers DESI fibres after the 0.736 GFA-acceptance correction, SDSS/BOSS total-sky ON−OFF fibres, SDSS imaging sky on 17′ scales, and the NSNS map. The paper's 10–13 R (red band) is 5–8× higher. The QHYCCD press-note value (3.3–6.7 R depending on bandwidth) and the paper's own MRS-N spectra rescaled with sky lines (~3 R, wave 1) lie in between. The northern arc is 0.4–0.8 R.

**Sources:**
- `notes/paper_review.md`: Table 3 → 12.8 / 9.9 / 10.1 R; 8.3–14.9 R for a 50–70 Å width.
- `notes/discovery.md`: QHYCCD.
- `notes/imaging_optical.md`: NSNS on-ridge medians 2.04–2.18 R and 1′ peaks 2.95–3.64 R.
- `notes/verify_precovery.md`: SDSS imaging, PSP 2.31 ± 0.38 R; range 1.65 (all confounders) to 2.55 (Frames).
- `notes/verify_desi.md` and `verify_desi_brightness.json`: south 1.62 R (variants 1.60–1.71, ±8 % sys); E nodes 1.44–2.30; W @ Dec 18 3.88 ± 0.20; north 0.56–0.71.
- `notes/verify_sdss.md`: south 1.39 ± 0.07 ± 0.07 R; W strand 1.93 ± 0.11; P1/P2/P3 1.49/1.78/1.87; per-plate peaks 1.6–2.5 S and 0.4–0.8 N.
- `notes/spectra_sdss_lamost.md` (wave 1, not re-verified): LAMOST LRS 0.88 ± 0.13 R, averaged over |d| < 0.2° and Dec 12–28, so diluted; MRS-N × sky scale 2.8–3.2 R (1.2–8.6).

## F4 · `figures/w3_figs_f4_spectra.png` (spectral evidence)

**Caption.** The [O III]-correlated spectrum of the southern filament, normalised to unit [O III] 5007 flux (energy units, per Å). Top row: DESI DR1 template regression (Dec < 19.5, 16–84 % block-bootstrap band). Bottom row: SDSS/BOSS total-sky ON−OFF stack (Dec < 20).

- **[O II]:** [O II] 3726/29 is as strong as [O III] and resolved in DESI. [O II]/[O III] is 1.1–1.2 (DESI) and 1.8 (SDSS BOSS fibres).
- **Hβ:** Hβ is barely visible. The dashed curve shows Hβ for [O III]/Hβ = 3 (a complete radiative shock) and the dotted curve for = 11. Measured: > 11–18 (DESI, 2σ); ≈ 23 (SDSS, 13–48).
- **4363:** [O III] 4363 appears at the [O III] wavelength (2–3σ), not at the Hg I 4358 street-lamp line.
- **Hα:** Hα is weak against [N II]. In the SDSS Dec < 20 stack it includes Hα structure near the patch; the Dec < 17.6 stack (orange) shows little Hα.

**Sources and method:**
- DESI 3727 / Hβ / 5007 / 4363 windows: `verify_desi_regspec.npz` (wave 2, subset "south").
- DESI Hα + [N II]: `data/derived/w3_figs_desi_regspec_ha.npz`, **derived here** with the unchanged `verify_desi_regspec.contsub`/`regspec` functions and the same base selection (n = 3,852 south).
- SDSS: `data/derived/w3_figs_sdss_stacks.npz`, **derived here** with the unchanged `verify_sdss_figures.stacks` (tag base, f = 0.84/0.60, ON < 0.15°, OFF 0.45–2°).
  - The [O II] panel uses BOSS-only groups, normalised by the [O III] of the same groups.
  - Linear baselines are fitted away from the lines.
  - Normalisation is in `w3_figs_f4_norm.json`. The SDSS line-integral 5007 is 1.55 R, against the fitted 1.39 R of verify_sdss (different estimator).
- Ratios quoted in the panels are the wave-2 values.

## F5 · `figures/w3_figs_f5_diagnostics.png` (excitation diagnostics)

**Caption.** Panel (a) plots [O III]/Hβ against [O II]/[O III]; panel (b) plots [O III] 4363/5007 (right axis: low-density T_e) against [O III]/Hβ. Energy ratios throughout.

**What the figure shows:**
- **Observed:** the southern constraints (DESI box: [O II]/[O III] 1.03–1.24, [O III]/Hβ > 11–18; SDSS point) lie where truncated (incomplete) radiative shocks of 100–160 km s⁻¹ pass. The diamonds mark the best truncation.
- **Truncated-shock tracks** are cumulative ratios against swept column, from the verify_shockmodel NEI code. The arrow shows the ×2.3 R88-calibrated Hβ systematic.
- **Excluded families:**
  - complete shocks (triangles; published range 1.2–4.5);
  - photoionisation at any U (grey curves) for the WD J0736 / 150 kK / 0.14 keV / DIG spectra;
  - conduction fronts (4363/5007 too high);
  - complete cooling flows.
- **North and patch:** the north (SDSS north and DESI-north [O II]) and the Hα patch move toward the low-[O III]/Hβ, high-[O II] corner, where more complete shocks lie.
- **Panel b:** 4363/5007 ≈ 0.06–0.09 needs T_e ≈ 3–6 × 10⁴ K. Photoionisation equilibrium gives ≤ 0.012.

**Sources:**
- `verify_shockmodel_tracks.npz` and `_fit.json`: v = 100–160, X(H⁰) = 0.05, b = 10 tracks, shown from 20 % of the final [O III]; best indices; b = 3 complete runs.
- R88 SE/M: `verify_xrayfuv_uv_predictions.json`; 4363 for SE (0.071) from `notes/verify_xrayfuv.md`.
- `verify_altmodel_photoion.json`, `verify_altmodel_cooling.json` (GS07 isobaric/isochoric Z = 1 stopped at 5, 3, 2, 1 × 10⁴ K, plus the H tail; warm CIE vs T; static and evaporating fronts) and `verify_altmodel_scenarios.json` (PIE + heating, fossil).
- Observed values: `notes/verify_desi.md` rev. 2 and `verify_sdss_results_base_r0.15.json`. SDSS photon ratios are converted to energy (×4861/5007 for [O III]/Hβ).
  - SDSS south: [O III]/Hβ 22.8 (12.4–46.2); [O II]/[O III] 1.84 (95 %: 1.46–2.37).
  - SDSS north: 6.2 (3.1–12.0) and 7.3.
  - Hα patch: 1.73 and 11.8.
  - DESI 4363/5007: 0.06–0.09 with the Gaussian value 0.083 ± 0.029. SDSS: 0.038 ± 0.035.

**Caveats:**
- The verify_shockmodel fit used the withdrawn DESI [O III]/Hβ ≈ 46 as input. The figure shows the revised > 11–18 instead.
- DESI did not measure a northern [O III]/Hβ. The DESI-north box spans the SDSS-north [O III]/Hβ range.
- The T_e axis uses the classical low-density 4363 relation, which is approximate above 3 × 10⁴ K.

## F6 · `figures/w3_figs_f6_profiles.png` (cross-spine profiles, south)

**Caption.** Profiles across the E spine for the southern Silk (spine Dec < 17.6). Positive offset = east, outside Monogem.
- **Top row ([O III], Rayleigh):** the NSNS DR0.2 median profile has a ridge amplitude of 2.36 R (z = 33). DESI DR1 fibre residuals (GFA-corrected, petal-relative) give 1.40 R.
- **Other rows:** each is in units of its own placement-null scatter, with grey bands for the 16–84 % and 2.5–97.5 % null ranges.

Results:
- GALEX FUV and NUV show the same narrow ridge: 172 ± 32 CU (z = 5.2) and 141 ± 33 CU (z = 4.1).
- SHASSA Hα (Dec 12–15.6), GALFA HI (V_LSR −2…+6) and IRIS 100 µm show nothing at the spine.
- eRASS1 0.2–0.6 keV declines from inside to outside. It stays within the band of rotations about PSR B0656+14, so the gradient is generic to the Monogem shell at this radius, and there is no X-ray filament.

**Derived here** (`data/derived/w3_figs_f6_profiles.npz` and `.json`). These are presentation-level estimates; the formal values are in the wave-2 notes.

Method:
- Offsets use `verify_xrayfuv_common.CurveCoords`.
- Ridge amplitude: a Gaussian (FWHM 0.28°) plus a linear baseline, |x| ≤ 1.2° (GALEX 1.0°).
- Nulls exclude pixels within 1.0° (GALEX 0.6°) of the real spine.

Inputs and null placements per band:

| Band | Input | Null placements |
|---|---|---|
| NSNS | `w3_figs_f1_maps.npz`, 2′ blocks | 42 RA/Dec shifts |
| DESI | `verify_desi_linefluxes.npz` × `kR_5008`, with good_fit / ok_contrib | none (standard errors) |
| GALEX | `verify_xrayfuv_galex_measure.prep_band`, per-tile demeaning as in verify_xrayfuv_fig_uv | 20 RA shifts of ±0.8–3° |
| SHASSA | `verify_hiwham_hamaps_part_Sm.npz` (verify_hiwham profile) | 400 nulls from verify_hiwham |
| eRASS1 | `verify_xrayfuv_erass_pixels.npz` | 27 rotations about the PSR, \|θ\| = 20–180° |
| GALFA HI | `imaging_mw_eq_galfa_nhi_m02p06.fits` | 20 usable shifts (the eastern GALFA tiles are missing) |
| IRIS | `imaging_mw_eq_iris100.fits` | 42 shifts |

Comparison with the wave-2 values:
- These numbers agree with wave 2: verify_xrayfuv south FUV contrast +120 ± 26 CU (z = 4.6) and NUV +96 ± 23 (z = 4.2); verify_hiwham SHASSA z = 0.39; verify_hiwham south HI z = +0.6; on-spine X-ray bump −0.38 cts s⁻¹ deg⁻².
- The HI window here (V_LSR −2…+6, wave-1 EQ grid) is wider than the task's "~0…+4". verify_hiwham's per-slice limits (1.5–1.75 × 10¹⁹ cm⁻², 2σ, south) are the formal ones.

## F7 · `figures/w3_figs_f7_timeline.png` (precovery and discovery, 2000–2026)

**Caption.** Top: on-strand [O III] S/N per epoch for each dataset.
- **SDSS imaging g-sky runs.** Filled: Dec 12–20.5. Hollow: the north, a broadband confounder.
- **SDSS/BOSS plates.** Filled: S/N includes the per-plate placebo. Hollow: plates without a placebo test.
- **LAMOST LRS plans.** Stellar residuals.
- **DESI DR1 months.** Partly circular: regressed on a forward model fitted to all epochs. Values above 19 are capped and annotated.

The [O III] filament is present from the first usable epochs:
- earliest ≥ 3σ: SDSS plate 927 on 2002-10-30 (3 fibres, not placebo-tested);
- first imaging detection: run 3560 on 2002-12-31;
- plate 1203 on 2003-01-30, confirmed by two independent analyses;
- first ≥ 5σ: plate 2074 on 2005-03-03.

That is about 21 years before its 2023 discovery. Middle: the epochs of each dataset. Bottom (zoom 2022.75–2027): the independent discoveries and publications by SWONS/Wu et al., NSNS/Ziegenbalg and Fesen/Drechsler/Schaeffer.

**Sources:**
- `verify_sdss_precovery.csv` and `_placebo.csv` (verify_sdss).
- `verify_precovery_table.json`: g PSP gauss028 per run; wave-1-residual per-plate S/N for 927-52577 (5.0σ) and 1263-52708 (3.9σ); LAMOST LRS per plan; DESI months.
- Events from `discovery_timeline.json` and `notes/discovery.md`. Key dates: AstroBin zh47r9 2024-01-14; Cloudy Nights 2024-06-11; QHYCCD 2024-07-12; RNAAS 2025-08-28; Fesen+ arXiv 2026-08-27; Wu+ arXiv 2026-10-09; AJ 172, 298 on 2026-11-04.

**Caveats:**
- The SWONS (2023-11) and NSNS full-array (2024-01) dates are self-reported.
- The DESI bar covers the coadd dates used here (2021-01 to 2022-02), not the full DR1 span.

---

## Files written (all new; no earlier file modified)

**Scripts** (`scripts/`):
- `w3_figs_common.py`
- `w3_figs_f1_data.py`, `w3_figs_f1_overview.py`
- `w3_figs_f2_velocity.py`
- `w3_figs_f3_brightness.py`
- `w3_figs_desi_regspec_ha.py`, `w3_figs_sdss_stacks.py`, `w3_figs_f4_spectra.py`
- `w3_figs_f5_diagnostics.py`
- `w3_figs_f6_profiles_data.py`, `w3_figs_f6_profiles.py`
- `w3_figs_f7_timeline.py`

**Derived** (`data/derived/`):
- `w3_figs_f1_maps.npz`
- `w3_figs_f2_velocity.json`
- `w3_figs_desi_regspec_ha.npz`
- `w3_figs_sdss_stacks.npz`
- `w3_figs_f4_norm.json`
- `w3_figs_f6_profiles.npz`, `w3_figs_f6_profiles.json`

**Figures** (`figures/`):
- `w3_figs_f1_overview.png`
- `w3_figs_f2_velocity.png`
- `w3_figs_f3_brightness.png`
- `w3_figs_f4_spectra.png`
- `w3_figs_f5_diagnostics.png`
- `w3_figs_f6_profiles.png`
- `w3_figs_f7_timeline.png`

No downloads were made. The scripts import the wave-2 modules `verify_desi_regspec`, `verify_desi_common`, `verify_sdss_lib/results/figures`, `verify_xrayfuv_common/galex_measure` and `imaging_optical_common` read-only.
