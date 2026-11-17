## 4. The LAMOST Hα velocity offset is plausibly a geocoronal-subtraction artefact

**Claim in the paper.**
- **Table 4 (p. 13).** Hα V_helio is +16.3 ± 1.2, +17.0 ± 1.2, +18.4 ± 1.4 and +19.7 ± 1.4 km s⁻¹ for P1–P4 (V_LSR +4.0…+6.1). Hα kinematic distances are given in both distance columns.
- **Sect. 3.2.3 (p. 10).** P4 "is better characterized by the stronger Hα line, yielding a smaller uncertainty of ∼1.4 km s⁻¹". Fig. 5 maps the V_LSR "measured with Hα".
- **Sect. 2.2 (p. 6).** After the Hα sky subtraction of Zhang et al. (2021), "the radial-velocity precision reaches < 1.5 km s⁻¹".
- **Sects. 3.3 and 4.2 (pp. 11, 16–17).** Hα, [N II] and [S II] are treated as one red-arm velocity system. Sect. 4.2 says only that "the velocity differences among Hα, [N II], and [S II] in this region remain below 20 km s⁻¹" (p. 17). The paper does not discuss why Hα is 4.5–7.0 km s⁻¹ redder than [N II] and [S II] in every region of Table 4.

The paper is right that the red-arm lines trace a foreground/background WIM, the same on and off the filament (Sect. 4.2). This section concerns only the Hα velocity.

**Why it is wrong.** Unless stated otherwise, numbers come from `scripts/dispute_halpha_analysis.py` → `data/derived/dispute_halpha_results.json`, plus `scripts/dispute_halpha_review.py` → `data/derived/dispute_halpha_review.json`. The Hα, [N II], [S II] and WHAM velocities are also compared side by side in [VELOCITY.md](VELOCITY.md) §5–6.

1. **The offset is internal to the paper's own red arm.** From Table 4, with errors added in quadrature:
   - Hα − [N II] = +4.5, +5.0, +6.7 and +6.4 km s⁻¹. The weighted mean is **+5.5 ± 1.1 km s⁻¹ (4.9σ)**, with χ² = 0.7 about the mean for 3 dof. The offset is consistent with a constant.
   - Hα − [S II] = +5.8 ± 1.2 km s⁻¹, while [N II] − [S II] = +0.1 ± 1.4.

   The offset is as large in the nebula-free P4 as it is on the nebula.
2. **WHAM agrees with [N II] and [S II], not with Hα.**
   - **Per region.** In the nearest 1° beams of WHAM-SS DR1 (Haffner et al. 2003), the main component lies at V_helio +11.6…+13.4 (V_LSR −0.3…+1.3; FWHM_int 28–31 km s⁻¹; total intensity 2.9–4.4 R). Hα − WHAM = **+4.8 ± 0.7 km s⁻¹**, whereas [N II] − WHAM = −0.8 ± 1.0 (Table 5).
   - **Across the plate.** The 2335 Hα V_LSR values digitized from Fig. 5 have a median of +5.4 km s⁻¹. The colour bar clips some values, but not enough to move the median. The 25 WHAM beams inside the 2.5° plate have a median of +0.95 (16–84 %: +0.4…+2.0), and WHAM is the same on and off the filament (`notes/verify_hiwham.md`).
   - **No hidden component.** WHAM integrates the same diffuse Hα. A real Hα component 4–5 km s⁻¹ redder than [N II] would therefore appear in it.
3. **Standard rest wavelengths and the wavelength scale do not explain it.**
   - Using 6564.632 Å (NIST Ritz; Kramida et al. 2024) instead of 6564.61 Å moves Hα by only 1.0 km s⁻¹.
   - To produce the offset alone, the rest wavelength would have to be 6564.49 ± 0.02 Å. That is below every H I 2–3 component NIST lists (6564.52–6564.72 Å). The published numbers cannot exclude it. Such an error would also be an artefact, but it would not grow with the sky/nebula ratio (Requested correction 3).
   - The red-arm scale is tied to sky lines at 6287–6553 Å and extrapolated to 6731 Å (Wu et al. 2022, Sect. 3.2; Ren et al. 2021). Yet [S II], the line furthest from that range, agrees with [N II] and with WHAM. A scale error cannot move Hα alone.
4. **On this plate, the sky line is blended with the WIM line and is as bright as it.**
   - **Velocity.** For the exposures of 2024-02-01 (mid-times 13:58–15:12 UT, Xinglong), v_corr = −8.1…−8.9 km s⁻¹ (Astropy Collaboration 2022). Geocoronal Hα is mainly Lyβ-pumped 3p→2s emission. The 2:1-weighted centroid of the NIST 2s–3p lines (6564.5377 and 6564.5844 Å) lies 2.6 km s⁻¹ blueward of 6564.61 Å. The sky line therefore sits at V_helio −10.7…−11.4 km s⁻¹. That is 22–25 km s⁻¹ from the WIM, less than one instrumental FWHM (34.5 km s⁻¹; Zhang et al. 2021).
   - **The sky intensity is predicted, not fitted.** The separation |RV([N II]) − RV(λ6554)| = 20–22 km s⁻¹ is below the 25 km s⁻¹ needed to resolve the two lines (Zhang et al. 2021, Sect. 2.2). The pipeline therefore predicts I(Hα_sky)/I(λ6554) = 0.954 − 0.011 × (−Sunalt) (Wu et al. 2022, Sect. 3.3). At solar altitudes of −51°, −60° and −62°, this gives 0.39, 0.29 and 0.27 (mean 0.32).
   - **The sky is as bright as the nebula.** Table 2 gives the nebular Hα as 0.25–0.37 in the same unit (the OH line the paper calls λ6553). The subtracted sky is therefore **0.86–1.27×** the published line. On the MRS-N validation plates, the sky was "about 40% or even higher" of the nebula. Plates that deviated from the relation kept residuals of up to 7 % (Zhang et al. 2021, Sect. 2.4). Here, the same fractional error in the ratio has 2–3× more leverage.
5. **Within the subtraction model, only over-subtraction reddens the fit.**
   - **Model.** The forward model of `scripts/verify_velocity_halpha_sim.py` adds the WIM line (at the [N II] velocity, with WHAM's width) and the geocoronal line. It subtracts the MRS-N template, with σ = √(σ²₆₅₅₄ + 0.02²) Å, centred at RV(λ6554) + δ. A Gaussian plus a quadratic is then fitted.
   - **Sign test.** With the correct sky intensity or too little subtracted, no δ in −16…+20 km s⁻¹ shifts the centroid by more than +0.4 km s⁻¹. The pipeline placement (δ = 0) gives −0.0…−0.9 km s⁻¹. The literal "6563/6554 × RV(λ6554)" placement gives −0.6…−1.7 km s⁻¹. These ranges cover both normalisations of step 6. A redward shift needs a negative residual on the blue wing (Fig. 4a, b).
6. **How much over-subtraction?**
   - **Self-consistency.** The published Hα flux is itself the output of the biased fit. We therefore iterate the sky/nebula ratio until the model's fitted flux reproduces Table 2. In the adopted fit, the true nebular Hα is 0.35–0.46, and the subtracted sky is 0.69–0.92× the true line.
   - **Centroids alone.** The offsets need **47 %** over-subtraction (≥ 35 % at 1σ; χ² = 0.32 for 2 dof). The fitted Hα would then have an intrinsic FWHM of 13–23 km s⁻¹, against 31.0–32.4 ± 3–4 km s⁻¹ in Table 2.
   - **Centroids and widths (adopted).** With the true width set to WHAM's, the fit gives **30 % (24–33 %)**, with the template 2.5 km s⁻¹ from RV(λ6554), and χ² = 4.25 for 6 dof.
     - It reproduces the widths (27.6–34.3 km s⁻¹) and +3.2…+5.7 km s⁻¹ of the offsets.
     - It leaves **+1.4 ± 1.1 km s⁻¹ (1.2σ)** unexplained, so about three-quarters of the offset is explained.
     - With the template fixed at RV(λ6554), the fit gives 14 % (χ² = 9.1 for 7 dof) and explains about half, leaving +2.9 ± 1.1 km s⁻¹.
     - Taking the Table 2 flux at face value instead gives 16 % and 12 %, explaining 64 % and 51 %.
   - **Required sky ratio.** The true ratio would be 0.25 (0.28 with the template fixed) instead of 0.32. This lies within the 0.18–0.83 of Zhang et al.'s calibration plates, three of whose ten deviate from the relation. The pipeline applies one ratio per plate (Zhang et al. 2021, Sect. 2.1). Evaluated at the first exposure, the relation gives 23 % more sky than the exposure mean, and the epoch used is not stated.
   - **Solar Hα absorption.** The Moon was below the horizon (−20° to −7°), but the field lies at ecliptic latitude −5°. Solar lines carry 21 ± 2 % of the blue-arm continuum in Fig. 3a (`notes/verify_velocity.md`). Solar Hα absorption in such a continuum (Kurucz 2005 atlas) shifts the centroid in the same direction. The shift is +0.6…+2.7 km s⁻¹, instead of −0.5, for 0.1–0.3 Å⁻¹ of continuum relative to the line flux (`verify_velocity_halpha_sim.json`). This could supply the remaining +1.4 km s⁻¹, but the red-arm continuum level is not published.
7. **Status.**
   - **Plausible, not demonstrated.** This is a model only. The MRS-N spectra are unreleased (L.-Y. Zhang et al. 2025, arXiv v1, Sect. VI), and the geocoronal intensity on this night is unknown. With the published widths imposed, the model explains one-half to three-quarters of the offset.
   - **Independent of the model.** Steps 1–2 hold regardless: Hα is +4.8 ± 0.7 km s⁻¹ off WHAM on this plate, and [N II] and [S II] are not.

![Figure 4: geocoronal-subtraction forward model and per-region Hα, [N II], [S II] and WHAM velocities](figures/dispute_halpha_fig4.png)

*Figure 4a–d. (a) Forward model for P1 at the adopted joint fit, which fits centroids and widths with a self-consistent sky/nebula ratio (step 6). It shows: the WIM Hα (V_helio +11.8, the Table 4 [N II] velocity, with intrinsic FWHM 28.7 km s⁻¹ from WHAM); the true geocoronal line at v_corr − 2.6 km s⁻¹ (NIST 2s–3p, 2:1); the MRS-N template, 30 % too bright, 2.5 km s⁻¹ from RV(λ6554); the residual; and the Gaussian + quadratic fit, whose centroid is shifted by +3.8 km s⁻¹. All curves are convolved with the 0.32 Å LSF; the ticks mark the centroids. (b) Fitted centroid shift against over-subtraction (model/true − 1) for P1–P4, at the joint-fit template position. The sky/nebula ratio is made self-consistent at each point; below −30 % the iteration does not converge. The red band is the observed Table 4 Hα − [N II] offset (dark: weighted mean ± 1σ; light: P1–P4 range). The green band and line show the joint fit (1σ range, best value). The dotted line shows the centroid-only fit, whose template position (+2.25 km s⁻¹) is close to that of the curves. (c) Fitted intrinsic Hα FWHM for the same models, with the Table 2 Hα widths (weighted mean ± 1σ) and the WHAM widths used as the true WIM width. (d) V_helio of LAMOST Hα, [N II] and [S II] (Table 4); Hα minus the joint-fit model shift; WHAM-SS DR1 Hα from the nearest beam (error: beam-to-beam scatter within 1°, floor 0.5 km s⁻¹); and the [O III] velocity adopted in Sect. 3. Top axis: V_LSR ≈ V_helio − 12.2 km s⁻¹ (mean correction for P1–P4). Model code: `scripts/verify_velocity_halpha_sim.py` (reused) and `scripts/dispute_halpha_analysis.py`; figure: `scripts/dispute_halpha_figure.py`.*

**Table 5.** Per-region Hα velocities and offsets, in km s⁻¹ (V_helio unless stated). LAMOST values are from Wu et al. Table 4; sky/nebula ratios use Table 2. WHAM values are the main-component velocity of the nearest WHAM-SS DR1 beam (`data/derived/verify_hiwham_wham_beams.csv`), converted with the paper's solar motion (V_LSR − V_helio = −11.9…−12.6). The quoted error is the beam-to-beam scatter within 1°; a floor of 0.5 is applied in the differences. Model rows and distances come from `data/derived/dispute_halpha_results.json`. Distances use our re-implementation of the `kd` (Wenger et al. 2018) traditional solution with the Reid et al. (2014) curve (`scripts/verify_paper_kd.py`).

| Quantity | P1 | P2 | P3 | P4 (no nebula) | Combined |
|---|---|---|---|---|---|
| LAMOST Hα (Table 4) | +16.3 ± 1.2 | +17.0 ± 1.2 | +18.4 ± 1.4 | +19.7 ± 1.4 | — |
| LAMOST [N II] (Table 4) | +11.8 ± 1.6 | +12.0 ± 1.6 | +11.7 ± 1.9 | +13.3 ± 2.3 | — |
| LAMOST [S II] (Table 4) | +11.8 ± 2.5 | +12.0 ± 1.6 | +11.7 ± 1.9 | +12.7 ± 2.1 | — |
| WHAM Hα V_helio (V_LSR; beam offset) | +13.4 ± 0.4 (+1.3; 0.48°) | +11.6 ± 1.1 (−0.3; 0.27°) | +13.4 ± 0.3 (+0.9; 0.35°) | +12.9 ± 0.5 (+0.5; 0.37°) | — |
| Hα − [N II] | +4.5 ± 2.0 | +5.0 ± 2.0 | +6.7 ± 2.4 | +6.4 ± 2.7 | **+5.5 ± 1.1** |
| Hα − [S II] | +4.5 ± 2.8 | +5.0 ± 2.0 | +6.7 ± 2.4 | +7.0 ± 2.5 | +5.8 ± 1.2 |
| Hα − WHAM | +2.9 ± 1.3 | +5.4 ± 1.6 | +5.0 ± 1.5 | +6.8 ± 1.5 | **+4.8 ± 0.7** |
| [N II] − WHAM | −1.6 ± 1.7 | +0.4 ± 1.9 | −1.7 ± 2.0 | +0.4 ± 2.4 | −0.8 ± 1.0 |
| Subtracted sky / nebular Hα: published flux (true flux, joint fit) | 1.06 (0.82) | 0.86 (0.69) | 1.23 (0.92) | 1.27 (0.91) | — |
| Model shift, centroid-only fit (47 %) | +5.4 | +4.5 | +6.0 | +6.5 | χ² 0.32 / 2 dof |
| Model Hα FWHM, centroid-only fit (Table 2: 31.9, 31.0, 32.2, 32.4) | 21.5 | 22.7 | 20.0 | 13.0 | widths not fitted |
| Model shift, joint centroid + width fit (30 %) | +3.8 | +3.2 | +4.3 | +5.7 | χ² 4.25 / 6 dof |
| Model Hα FWHM, joint fit | 32.8 | 30.8 | 34.3 | 27.6 | — |
| Table 4 Hα distance, MC / traditional (as printed) | 0.47 / 0.61 | 0.31 / 0.54 | 0.34 / 0.58 | 0.54 / 0.80 | — |
| Traditional distance (kpc) from Hα / [N II] / WHAM velocity | 0.43 / none / 0.13 | 0.56 / 0.01 / none | 0.57 / none / 0.08 | 0.77 / 0.09 / 0.05 | — |

**Consequence.**
- **The Hα velocities are red-biased.** The Hα V_helio and V_LSR in Table 4 and Fig. 5 carry a +4.5…+6.7 km s⁻¹ systematic error. Neither the quoted 1.2–1.4 km s⁻¹ errors nor the "< 1.5 km s⁻¹" precision covers it. The Hα velocity is no evidence of a separate nebular or DIG component, and it must not be used for distances.
- **The Hα distances come from the offset.** At the local gradient of 8.9–9.8 km s⁻¹ kpc⁻¹, the Hα − [N II] offset alone is worth 0.49–0.69 kpc. The printed Hα distances (0.31–0.80, in kpc per Sect. 1) owe their size to it. With the [N II] or WHAM velocity, the same rotation model gives no solution or 0.01–0.13 kpc, which carries no information here (Sect. 2).
- **Widths, temperatures and line ratios are affected too (model only).** In the adopted fit, the residual changes the Hα widths by −3…+5 km s⁻¹ and the Hα flux by −19…−29 %.
  - The Eq. 2 temperatures shift by −4100…+7400 K. That is comparable to the 3200–7700 K that Eq. 2 gives from the Table 2 widths.
  - [S II]/Hα would fall from 0.78–0.93 to 0.56–0.72, still above the paper's 0.4 criterion. [N II]/Hα would fall from 0.88–1.11 to 0.68–0.86.
- **The DIG conclusion (Sect. 4.2) stands.** The artefact affects fibres on and off the nebula alike.

**Requested correction.**
1. **Flag the offset.** In Table 4, Fig. 5 and Sect. 3.2.3, state that Hα is 4–7 km s⁻¹ redder than [N II], [S II] and WHAM. Withdraw "better characterized by the stronger Hα line".
2. **Remove the Hα distances.** Delete the Hα rows of the Table 4 distance columns, and any Sect. 3.3 red-line distance that uses them.
3. **Refit Hα.**
   - Leave the geocoronal intensity free, or tie the Hα centroid to [N II] (Zhang et al. 2021, Sect. 2.2).
   - Report the solar altitude and the sky ratio that were applied, and the rest wavelength used for Hα.
   - Plot per-fibre Hα − [N II] against the sky/nebula ratio. A subtraction artefact grows with this ratio; a rest-wavelength error does not.
4. **Release the spectra.** Publish the red-arm spectra of P1–P4 with the subtracted sky template.

Suggested erratum: *"The Hα radial velocities in Table 4 and Fig. 5 are 4–7 km s⁻¹ redder than the simultaneously measured [N II] and [S II] lines and than WHAM Hα. On this plate the geocoronal Hα line lies within one resolution element of the nebular line and has comparable strength. The offset is most plausibly due to imperfect subtraction of the geocoronal line. The Hα velocities and Hα-based kinematic distances should not be used."*

**References used in this section.**
- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (2026), The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ 172, 298, doi:10.3847/1538-3881/aeb11c, arXiv:2610.13329 [verified via Crossref and arXiv API]
- Zhang, W., Wu, H., Wu, C.-J., Ren, J.-J., et al. (2021), LAMOST medium-resolution spectral survey of Galactic nebulae (LAMOST-MRS-N): subtraction of geocoronal Hα emission, RAA 21, 280, doi:10.1088/1674-4527/21/11/280, arXiv:2108.08021 [verified via Crossref and arXiv API; quoted text checked in the arXiv PDF]
- Wu, C.-J., Wu, H., Zhang, W., Li, Y., et al. (2022), The Data Processing of the LAMOST Medium-resolution Spectral Survey of Galactic Nebulae (LAMOST MRS-N Pipeline), RAA 22, 075015, doi:10.1088/1674-4527/ac7387, arXiv:2110.04547 [verified via Crossref and arXiv API; Sects. 3.2–3.3 checked in the arXiv PDF]
- Ren, J.-J., Wu, H., Wu, C.-J., Zhang, W., et al. (2021), On the radial velocity calibrations in the LAMOST medium-resolution spectroscopic survey of nebulae, RAA 21, 051, doi:10.1088/1674-4527/21/3/51, arXiv:2008.06236 [verified via Crossref and arXiv API; sky-line list checked in the arXiv PDF]
- Zhang, L.-Y., Wu, C.-J., Fang, X., Zhang, W., et al. (2025), LAMOST Medium-resolution Spectroscopic Survey of the Rosette Nebula, AJ 170, 309, doi:10.3847/1538-3881/ae0f94, arXiv:2510.08950 [verified via Crossref and arXiv API; Sect. VI data statement read in arXiv v1]
- Haffner, L. M., Reynolds, R. J., Tufte, S. L., Madsen, G. J., et al. (2003), The Wisconsin Hα Mapper Northern Sky Survey, ApJS 149, 405, doi:10.1086/378850 (WHAM-SS DR1 data, version 161116) [verified via Crossref]
- Kramida, A., Ralchenko, Yu., Reader, J., & NIST ASD Team (2024), NIST Atomic Spectra Database (ver. 5.12), doi:10.18434/T4W30F, H I lines retrieved 2026-11-16 and re-checked 2026-11-17 [verified via DataCite and the NIST ASD output page]
- Kurucz, R. L. (2005), New atlases for solar flux, irradiance, central intensity, and limb intensity, Mem. Soc. Astron. Ital. Suppl. 8, 189, bibcode 2005MSAIS...8..189K, no DOI; atlas file from kurucz.harvard.edu/sun/fluxatlas2005 [verified via the ADS record and the data-server file header]
- Reid, M. J., Menten, K. M., Brunthaler, A., Zheng, X. W., et al. (2014), Trigonometric Parallaxes of High Mass Star Forming Regions: The Structure and Kinematics of the Milky Way, ApJ 783, 130, doi:10.1088/0004-637X/783/2/130 [verified via Crossref]
- Wenger, T. V., Balser, D. S., Anderson, L. D., & Bania, T. M. (2018), Kinematic Distances: A Monte Carlo Method, ApJ 856, 52, doi:10.3847/1538-4357/aaaec8 [verified via Crossref]
- Astropy Collaboration, Price-Whelan, A. M., Lim, P. L., et al. (2022), The Astropy Project: Sustaining and Growing a Community-oriented Open-source Project and the Latest Major Release (v5.0) of the Core Package, ApJ 935, 167, doi:10.3847/1538-4357/ac7c74 [verified via Crossref]

*Reference verification log: `data/derived/dispute_halpha_refs.json` (`scripts/dispute_halpha_refs.py`); NIST line values: `data/derived/dispute_halpha_review.json`.*
