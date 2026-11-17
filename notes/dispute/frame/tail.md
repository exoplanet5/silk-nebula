## 7. What the paper gets right

- **The [O III] filament is real line emission, and P4 is free of it.** DESI DR1 detects it independently at ~30σ against empirical nulls. SDSS/BOSS sky fibres detect it along its whole length, and SDSS imaging sky levels show it from 2002 onward (`REPORT.md` §2, §4).
- **The red-arm lines (Hα, [N II], [S II]) are foreground/background warm ionized gas** (Sect. 4.2). WHAM shows the same component on and off the filament.
- **The red-arm [N II] and [S II] velocities agree with WHAM,** and the LAMOST stellar heliocentric correction for this plate is correct (Section 3).
- **A shock origin of the [O III] is plausible.** DESI and SDSS spectra give [O III]/Hβ > 11 (2σ), [O II]/[O III] ≈ 1.1–1.2 and He II/[O III] < 0.026. That is a Balmer-dark, collisionally excited spectrum, best matched by an incompletely cooled radiative shock (`REPORT.md` §5). It rests on evidence the paper did not have, since MRS-N has no Hβ.
- **Two conversions are arithmetically correct:** the energy formula (E = 2.01×10⁵¹ erg for the paper's inputs) and the magnitude-to-rayleigh conversion.

## 8. Consolidated requested corrections

Each item is detailed in the "Requested correction" paragraph of the section cited.

1. **Units and distance** (Section 1). In the Table 4 footnotes, change "in units of 100pc" to "in kpc". Withdraw the fiducial ~160 pc from the Abstract and Sects. 3.3, 4.1 and 5, and correct "∼100–190 pc", "approximately 100 pc" and "100–150 pc".
2. **Kinematic distances** (Section 2). State that kinematic distances are not meaningful toward (ℓ, b) ≈ (204°, +19°). Either delete the Table 4 distance columns or label them as formal values without physical meaning. Give the size, radius, energy and age as explicit functions of d.
3. **[O III] velocities** (Section 3). Document the blue-arm wavelength frame. Re-measure the velocities in the vacuum heliocentric frame, or lower them by ≈ 10 km s⁻¹. Withdraw the [O III]/red-line offset, the two-layer-shell discussion and the [O III] distances.
4. **Hα velocities** (Section 4). Flag the +4–7 km s⁻¹ Hα offset and remove the Hα distances. Report the sky-subtraction parameters, and plot per-fibre Hα − [N II] against the sky/nebula ratio.
5. **Surface brightness** (Section 5). Correct or withdraw the Table 3 values and the Sect. 4.1 comparison. Document the photometry: instrument, filter curve, apertures, background model and zero-point units.
6. **Diagnostics and energetics** (Section 6). Restate the red-line ratios as WIM properties, and correct the [S II] densities, T_e/V_t and abundance statements. Correct the age bracket ("+1/3") and R_PDS (17.1 pc), and present E and the age only as scalings with unknown d and n₀. Qualify the [O III] widths and λ4959, and remove the radio argument.
7. **Data.** Release the MRS-N per-fibre spectra of plate MN074209N160100N01 (2024-02-01), red and blue arms, with the subtracted sky templates and wavelength-frame metadata. Our [O III] and Hα corrections are aggregate estimates, and per-fibre spectra are needed to confirm them.
8. **Credit.** Make the AJ acknowledgment consistent with the text: "discovered and named by the SWONS team" and "detected for the first time to its full extent" should acknowledge the independent NSNS (Ziegenbalg 2025, RNAAS, published 14 months earlier) and Schaeffer/Fesen/Drechsler discoveries (`REPORT.md` §1.1).

## 9. Data, code and reproducibility

- **The paper.** The PDF, its layout text and the HTML are in `sources/arxiv/`. Table and figure values were read from the text layer. Fig. 3 spectra were extracted from the PDF vector data (`scripts/verify_velocity_fig3_vectors.py`), and Fig. 5b was digitized from a 300-dpi page crop.
- **Section scripts** (run with `~/.venvs/astro313/bin/python -I scripts/<name>.py`):

  | Section | Scripts | Data |
  |---|---|---|
  | 1 | `dispute_distance_kd.py` (reuses `verify_paper_kd.py`, a re-implementation of `kd`) | `data/derived/dispute_distance_kd.json` |
  | 2 | `dispute_kinematic_kd.py` | `dispute_kinematic_kd.json` |
  | 3 | `dispute_oiii_values.py`, `dispute_oiii_figure.py` | `dispute_oiii_values.json` |
  | 4 | `dispute_halpha_analysis.py`, `dispute_halpha_figure.py` (reuse `verify_velocity_halpha_sim.py`) | `dispute_halpha_results.json` |
  | 5 | `dispute_brightness_*.py` | `dispute_brightness_*.json` |
  | 6 | `dispute_other_fig6.py` | `dispute_other_values.json` |

- **Independent data.** DESI DR1 came via NOIRLab Data Lab/SPARCL. Also used: SDSS DR17/DR18/DR19, WHAM-SS DR1, GALEX, eROSITA DR1 and NSNS DR0.2. Methods and caveats are in `notes/verify_*.md` and `notes/w3_*.md`.
- **Assembly.** This document was assembled by `scripts/dispute_assemble.py` from `notes/dispute/*.md` and rendered by `scripts/make_report_html.py DISPUTE`.
- **Not available to us.** The MRS-N nebular spectra are not public, and neither are the SWONS photometric inputs or DESI sky-fibre spectra (NERSC-hosted). Where these limit a conclusion, the relevant section says so.
