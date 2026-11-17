## 3. The [O III] velocity offset is a wavelength-reference artefact

**Claim in the paper.**
- Sect. 3.2.3 (p. 10) says the MRS-N spectra "have been corrected for heliocentric motion using the LAMOST official pipeline", so the emission-line RVs "represent heliocentric velocities". Sect. 2.2 (p. 6) quotes a precision of "< 1.5 km s⁻¹" after the secondary wavelength calibration of Ren et al. (2021).
- The paper's Table 4 (p. 13) gives [O III] RVs for P1/P2/P3 of +26.8 ± 4.3, +26.8 ± 4.1 and +23.8 ± 3.3 km s⁻¹ (V_LSR 15.5/14.5/11.2). For [N II]/[S II] it gives +11.7…+13.3.
- The resulting [O III]–red-line offset is "∼5−15 km s⁻¹" (Sect. 3.2.3), or "∼10–16 km s⁻¹" in V_LSR (Sect. 3.3, p. 11). It underpins three things:
  - a "distance discrepancy of 100–150 pc" between blue- and red-arm distances (Sect. 4.2, p. 16);
  - the "different layers of the same expanding shell" interpretation (Sect. 4.2, p. 17);
  - the [O III]-only kinematic distance (Abstract; Sect. 3.3; Sect. 5).

**Why it is wrong.** The red-arm [N II]/[S II] velocities and the LAMOST stellar heliocentric correction are both sound (steps 5–6). The blue-arm [O III] velocities are not. Numbers come from `scripts/dispute_oiii_values.py` → `data/derived/dispute_oiii_values.json`, or from the project files named. A side-by-side table of all [O III] and red-line velocities, with their frames, is in [VELOCITY.md](VELOCITY.md) §5. "The paper's Table 4" is Wu et al.'s table; "Table 4" alone is ours, below.

1. **The paper's plotted wavelength scale is air and topocentric.**
   - The PDF stores Fig. 3 (p. 7) as vector data. The re-extracted Fig. 3a blue-arm spectrum has 3,590 samples at a median step of 0.106 Å, with the axes calibrated on the tick marks to 0.002 Å rms (`verify_velocity_fig3_measure.json`).
   - It contains the [N I] λλ5198, 5200 airglow doublet, which Fig. 3a labels "sky line". The doublet ratio is 1.71; the UVES atlas gives 1.73 (Hanuschik 2003).
   - On the plotted scale the airglow sits **+1.8 ± 0.8 km s⁻¹** from its air rest wavelengths. This uses five sky lines and the mean of the NIST ASD (Kramida et al. 2024) and Hanuschik rest values. Against Hanuschik's [N I] wavelengths alone, the reference the paper cites for this line (p. 6), the offset is +0.4 ± 0.8.
   - Every other frame is excluded. On each alternative scale the airglow would sit at:
     - air, heliocentric: −8.3 km s⁻¹ (excluded at 12σ);
     - air, with the correction applied with the wrong sign: +8.3 km s⁻¹ (7.8σ);
     - vacuum: +75…+83 km s⁻¹ (≥ 87σ).
   - LAMOST MRS 1D products are vacuum and heliocentric (Luo et al. 2015; Ren et al. 2021). This plate's stellar spectra carry `VACUUM = T`, `HELIO = T` and a median `HELIO_RV` of +8.33 km s⁻¹. For the three exposures of 2024-02-01 at Xinglong, astropy (Astropy Collaboration 2022) gives v_corr = −8.32 km s⁻¹ (−8.03…−8.59 over P1–P3 and the exposures).
   - The simplest reading is that the blue arm was re-referenced to air sky-line wavelengths, which also removed the heliocentric correction. Ren et al. (2021) and Wu et al. (2022) describe sky-line recalibration of the red arm only (OH/[O I], 6287–6553 Å), so the blue-arm step is undocumented. Steps 3–4 do not depend on this mechanism.
2. **The paper's Table 4 was measured on this scale.**
   - In Fig. 3a, [O III] λ5007 lies at 5007.245 ± 0.014 Å; the Fig. 3c fit centre is 5007.246 Å. Relative to air rest that is +24.1 ± 0.9 km s⁻¹.
   - The Fig. 3b spectra give +24.9/+28.6/+20.3 km s⁻¹ for P1/P2/P3, within 0.4–1.0σ of the paper's Table 4.
   - Fig. 3b minus the paper's Table 4 averages −1.5 ± 2.4 km s⁻¹ (weighted; χ² = 1.3 for 3 dof about zero, with the paper's quoted errors). Were the table heliocentric and the plot not, the difference would be +10.1; this is excluded at 4.5σ.
3. **The velocity measured against the airglow needs no calibration.**
   - [O III] is measured against the airglow in the same spectrum. Subtract the +1.8 km s⁻¹ sky offset and a +0.4 km s⁻¹ centroid bias from a solar Fraunhofer blend at 5007.22 Å (Kurucz 2005 atlas fit; `verify_velocity_fig3_solar.json`). That gives **V_topo = +21.9 km s⁻¹**; adding v_corr gives **V_helio = +13.6 ± 1.2 km s⁻¹**, whatever produced the plotted scale.
4. **Size of the error.**
   - The paper's Table 4 [O III] values are too high by **10.1 ± 1.0 km s⁻¹**: 8.32 km s⁻¹ of missing v_corr plus the 1.8 km s⁻¹ sky offset.
   - Corrected values are P1 +16.7 ± 4.4, P2 +16.7 ± 4.2 and P3 +13.7 ± 3.5 km s⁻¹.
   - Only Fig. 3a (position not stated) shows the airglow, and the per-fibre MRS-N spectra are not public. The correction is therefore applied in aggregate. Plate zero points (0.2–4 km s⁻¹; Ren et al. 2021) and any per-fibre offsets cannot be checked.
5. **Alternatives are excluded.**
   - *Blue/red-arm offset.* On the paper's plate (MN074209N160100N01), LAMOST DR12 (Cui et al. 2012) per-arm stellar RVs give rv_b0 − rv_r0 = +0.7 ± 0.5 km s⁻¹ (N = 76). On 38 plates in the field the offset ranges from −1.6 to +1.9 km s⁻¹. Against Gaia DR3 (Gaia Collaboration 2023; Katz et al. 2023), rv_r0 − Gaia = +1.1 ± 0.8 and rv_b0 − Gaia = +1.0 ± 2.2 km s⁻¹ (N = 6).
   - *Air/vacuum mix-up.* This would produce 65–102 km s⁻¹ (Morton 2000 conversion).
   - *Solar blend.* It biases the centroid by only +0.4 km s⁻¹, and by at most +2.5 km s⁻¹ even for a purely solar continuum.
6. **Independent spectra agree with the corrected value.**
   - DESI DR1 (DESI Collaboration 2026) is vacuum and barycentric, with the correction applied to the wavelength grid before sky subtraction (desispec; Guy et al. 2023). Science-fibre residuals (flux − Redrock model) give +14.4 ± 1.9 km s⁻¹ from 390 on-ridge fibres (including a ±1.5 zero point) and +11.3 (16–84 %: 10.2–13.9) from a joint regression on both [O III] lines.
   - SDSS sky fibres (Abdurro'uf et al. 2022) give +9.6 ± 2.8 km s⁻¹ (all strands). A second re-analysis gives +12.3 ± 2.5 (all) and +10.0 ± 1.8 ± 0.8 km s⁻¹ at Dec < 20°, which contains P1–P3.
   - DESI and SDSS alone give **+12.9 ± 1.6 km s⁻¹**. The calibration-free Fig. 3a value lies +0.7 ± 2.0 above this. The paper's P1–P3 weighted mean (+25.5 ± 2.2) lies +12.6 ± 2.7 above it (4.6σ), and +2.5 ± 2.8 above it after correction.
   - Per region (Table 4), the paper minus DESI within 0.5° is +14.0 ± 2.4 km s⁻¹ (5.7σ), and +3.9 ± 2.5 (1.5σ) after correction.
   - The weighted mean of DESI, Fig. 3a and SDSS is **+13.3 ± 1.0 km s⁻¹** (χ² = 2.2 for 2 dof). Other estimator choices give +12.4 ± 0.9, hence a ±1 km s⁻¹ systematic. In V_LSR this is +1.2 ± 1.0, or +1.5 with Schönrich et al. (2010).
   - The velocity varies along the filament. DESI gives +9.7…+12.2 km s⁻¹ at Dec < 17.6° and +17…+23 at Dec 17.6–21°, where SDSS gives +10.0 ± 2.2 (`w3_northsouth.md`). P1 and P3 lie in the southern zone; P2 (Dec 17.74°) is at its edge, with a region-matched DESI value of +15.2 ± 2.9 (n = 9). The 10.1 km s⁻¹ error belongs to the wavelength scale and is unaffected.
   - The corrected [O III] agrees with the paper's own [N II] (P1–P4 mean +12.2 ± 0.9) and with WIM Hα from WHAM (Haffner et al. 2003). WHAM gives +11.5…+13.4 in the nearest beams, and V_LSR +1.3 on the southern filament against +1.2 ± 1.0 for [O III] (`verify_hiwham.md`).
   - [O III] − [N II] becomes **+1.1 ± 1.3 km s⁻¹** (+1.4 ± 1.5 from Fig. 3a alone), against +13.7 ± 2.4 from the paper's P1–P3 values.

![Three panels. (a) The [N I] 5198/5200 airglow doublet from the paper's Fig. 3a spectrum, with markers for the doublet positions on air or vacuum and topocentric or heliocentric scales; the data sit on the air-topocentric markers. (b) [O III] 5007 from Fig. 3a and the Fig. 3b P1-P3 spectra, with air and vacuum rest wavelengths and the plotted-scale positions of lines with true heliocentric velocity +13.3 and +25 km/s; the line peaks at the +13.3 position. (c) Velocity ladder: the paper's published Table 4 [O III], the paper's plotted-scale values, corrected values, DESI, SDSS, the paper's red lines, WHAM H-alpha and the paper's H-alpha, with the adopted +13.3 +- 1.0 km/s band.](figures/dispute_oiii_fig3.png)

*Figure 3a–c. The paper's [O III] velocities are on a topocentric, air-referenced scale.*
- *All panels: spectra re-extracted from the vector polylines of arXiv:2610.13329v1, p. 7 (`scripts/verify_velocity_fig3_vectors.py`). Values from `scripts/dispute_oiii_values.py`; plot from `scripts/dispute_oiii_figure.py`. v_corr = −8.32 km s⁻¹ (astropy, Xinglong, 2024-02-01).*
- *(a) [N I] λλ5198, 5200 in Fig. 3a, continuum-subtracted (white), with a two-Gaussian fit (blue). Vertical lines show where the doublet falls on each candidate scale: air-topocentric (green), air-heliocentric (orange dashed), vacuum-topocentric (cyan dotted) and vacuum-heliocentric, the native LAMOST frame (purple dash-dot). Rest wavelengths are from NIST ASD level energies; air values follow Morton (2000).*
- *(b) [O III] λ5007 in Fig. 3a (white; fit in blue) and the Fig. 3b P1–P3 spectra (thin, each normalised by its own fit). Top axis: velocity read off the plotted scale, as in the paper's Table 4. Green and red lines: where lines with true V_helio = +13.3 and +25 km s⁻¹ fall on this scale.*
- *(c) V_helio (bottom axis) and V_LSR ≈ V_helio − 12.2 (top axis; ±0.4 km s⁻¹ over P1–P4).*
  - *Region rows: the paper's Table 4 (red), Fig. 3b as plotted (orange squares), the paper's Table 4 − 10.1 (green), DESI positive-residual science fibres within 0.5° (blue; `verify_desi_velocity.json`).*
  - *Fig. 3a row: the plotted-scale and airglow-referenced values, with the velocity budget.*
  - *DESI row: on-ridge stack and joint regression (16–84 %). SDSS row: all-strand ON−OFF stack and the Dec < 20° stack of the second re-analysis.*
  - *Lower block: the paper's [N II] (purple) and [S II] (×), the WHAM-SS main Hα component in the nearest beam (`spectra_ism_wham_results.json`), and the paper's Hα (grey; treated in Section 4).*
  - *Green band: the adopted +13.3 ± 1.0 km s⁻¹.*

**Table 4.** [O III] velocities by region (km s⁻¹): the paper's Table 4, the paper's Fig. 3b re-measured, corrected, and compared with DESI, the red lines and WHAM.
- "Corrected" = published − 10.1 ± 1.0. The solar-blend bias (+0.4 km s⁻¹ in Fig. 3a; it depends on each spectrum's continuum) is not removed. V_LSR uses our per-region astropy correction with the paper's solar motion (−12.10/−11.87/−12.58 for P1–P3).
- DESI = weighted mean of positive-residual science fibres with S/N > 4 within 0.5°, without the ±1.5 km s⁻¹ zero point. P2 has only 9 fibres.
- WHAM = main Hα component of the nearest 1° beam (0.27–0.48° away; wave-1 fits, fit errors only).
- Differences are published → corrected, with errors added in quadrature. The "P1–P3 weighted" row treats the common ±1.0 correction as independent, so its errors are slightly understated.

| Region | Paper Table 4 V_helio | Paper Fig. 3b, as plotted | Corrected V_helio (V_LSR) | DESI ≤ 0.5° (n) | Paper − DESI: published → corrected | Paper [N II] | [O III] − [N II]: published → corrected | WHAM Hα V_helio |
|---|---|---|---|---|---|---|---|---|
| P1 | +26.8 ± 4.3 | +24.9 ± 2.5 | +16.7 ± 4.4 (+4.6) | +11.9 ± 1.9 (61) | +14.9 ± 4.7 (3.2σ) → +4.8 ± 4.8 (1.0σ) | +11.8 ± 1.6 | +15.0 ± 4.6 → +4.9 ± 4.7 | +13.1 ± 0.2 |
| P2 | +26.8 ± 4.1 | +28.6 ± 1.5 | +16.7 ± 4.2 (+4.8) | +15.2 ± 2.9 (9) | +11.6 ± 5.0 (2.3σ) → +1.5 ± 5.1 (0.3σ) | +12.0 ± 1.6 | +14.8 ± 4.4 → +4.7 ± 4.5 | +11.5 ± 0.1 |
| P3 | +23.8 ± 3.3 | +20.3 ± 1.3 | +13.7 ± 3.5 (+1.1) | +9.2 ± 1.0 (128) | +14.6 ± 3.5 (4.2σ) → +4.5 ± 3.6 (1.2σ) | +11.7 ± 1.9 | +12.1 ± 3.8 → +2.0 ± 3.9 | +13.4 ± 0.4 |
| P4 (background) | — | — | — | — | — | +13.3 ± 2.3 | — | +12.7 ± 0.2 |
| P1–P3 weighted | +25.5 ± 2.2 | — | +15.4 ± 2.3 | — | +14.0 ± 2.4 (5.7σ) → +3.9 ± 2.5 (1.5σ) | — | +13.7 ± 2.4 → +3.7 ± 2.5 | — |
| Fig. 3a (position not given) | — | +24.1 ± 0.9 | **+13.6 ± 1.2** (vs own [N I]) | — | — | +12.2 ± 0.9 (P1–P4 mean) | → +1.4 ± 1.5 | — |
| All data, adopted | — | — | **+13.3 ± 1.0** (+1.2) | — | — | +12.2 ± 0.9 | → +1.1 ± 1.3 | +12.7 (mean) |

**Consequence.**
- The corrected [O III] moves with the paper's own [N II]/[S II] and with the local WIM. The "∼5–15" (V_helio) and "10–16 km s⁻¹" (V_LSR) offsets are therefore not properties of the nebula. This removes the 100–150 pc blue/red distance discrepancy (Sect. 4.2, p. 16) and the motivation for the two-layer-shell and inhomogeneous-bubble scenarios (p. 17).
- The [O III] kinematic distances in the paper's Table 4 and the "fiducial ∼160 pc" rest on V_LSR = 11–16 km s⁻¹. The corrected values are V_LSR +1…+5 km s⁻¹ (P1–P3) and +1.2 ± 1.0 (adopted). At these velocities the kinematic method gives no useful distance (Section 2).
- Velocity no longer separates the [O III] from the DIG/WIM. Any separation must rest on morphology and line ratios.
- The "< 1.5 km s⁻¹" precision (Sect. 2.2) does not apply to the published blue-arm velocities, whose systematic error is about 10 km s⁻¹.

**Requested correction.**
1. State the wavelength frame used for the blue-arm data in the paper's Fig. 3 and Table 4: air or vacuum, topocentric or heliocentric, and the 5198 Å sky-line reference value. Then re-measure [O III] in the vacuum heliocentric frame, or against the per-fibre airglow with v_corr applied.
2. Erratum to the paper's Table 4: lower the [O III] RVs and V_LSR by about 10 km s⁻¹ (our aggregate estimate is 10.1 ± 1.0, giving +16.7/+16.7/+13.7 km s⁻¹ V_helio for P1–P3), or replace them with re-measured values. Withdraw the [O III] Distance¹/Distance² entries.
3. Delete the [O III]/low-ionisation offset statements (Sects. 3.2.3, 3.3, 4.2 and 5) and the two-layer-shell discussion (p. 17). Suggested wording: *"The [O III] RVs originally listed in Table 4 were measured on a topocentric, air-referenced blue-arm wavelength scale and are ≈ 10 km s⁻¹ too high. Once corrected, they agree with the [N II] and [S II] velocities, so the reported [O III]–low-ionisation offset and the [O III]-based kinematic distance are withdrawn."*
4. Release the per-fibre MRS-N blue-arm spectra, at least the 279 [O III] fibres, with their wavelength-frame metadata.

**References used in this section.**
- Abdurro'uf, Accetta, K., Aerts, C., et al. (2022), The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Release of MaNGA, MaStar, and APOGEE-2 Data, ApJS 259, 35, doi:10.3847/1538-4365/ac4414 [verified via Crossref].
- Astropy Collaboration, Price-Whelan, A. M., Lim, P. L., et al. (2022), The Astropy Project: Sustaining and Growing a Community-oriented Open-source Project and the Latest Major Release (v5.0) of the Core Package, ApJ 935, 167, doi:10.3847/1538-4357/ac7c74 [verified via Crossref].
- Cui, X.-Q., Zhao, Y.-H., Chu, Y.-Q., Li, G.-P., et al. (2012), The Large Sky Area Multi-Object Fiber Spectroscopic Telescope (LAMOST), RAA 12, 1197, doi:10.1088/1674-4527/12/9/003 [verified via Crossref].
- DESI Collaboration, Abdul Karim, M., Adame, A. G., et al. (2026), Data Release 1 of the Dark Energy Spectroscopic Instrument, AJ 171, 285, doi:10.3847/1538-3881/ae4c43; arXiv:2503.14745 [verified via Crossref and arXiv; first posted to arXiv in 2025].
- Gaia Collaboration, Vallenari, A., Brown, A. G. A., Prusti, T., et al. (2023), Gaia Data Release 3: Summary of the content and survey properties, A&A 674, A1, doi:10.1051/0004-6361/202243940 [verified via Crossref].
- Guy, J., Bailey, S., Kremin, A., et al. (2023), The Spectroscopic Data Processing Pipeline for the Dark Energy Spectroscopic Instrument, AJ 165, 144, doi:10.3847/1538-3881/acb212 [verified via Crossref; the order of the barycentric correction and sky subtraction was checked in the desispec `extract.py` code, not in the paper].
- Haffner, L. M., Reynolds, R. J., Tufte, S. L., et al. (2003), The Wisconsin Hα Mapper Northern Sky Survey, ApJS 149, 405, doi:10.1086/378850 [verified via Crossref].
- Hanuschik, R. W. (2003), A flux-calibrated, high-resolution atlas of optical sky emission from UVES, A&A 407, 1157, doi:10.1051/0004-6361:20030885 [verified via Crossref; line list from VizieR J/A+A/407/1157].
- Katz, D., Sartoretti, P., Guerrier, A., Panuzzo, P., et al. (2023), Gaia Data Release 3: Properties and validation of the radial velocities, A&A 674, A5, doi:10.1051/0004-6361/202244220 [verified via Crossref].
- Kramida, A., Ralchenko, Yu., Reader, J., & NIST ASD Team (2024), NIST Atomic Spectra Database (ver. 5.12), NIST Standard Reference Database 78, doi:10.18434/T4W30F [verified via DataCite; the N I level energies 19224.464 and 19233.177 cm⁻¹ re-checked in ver. 5.12 on 2026-11-17; the version first consulted was not recorded].
- Kurucz, R. L. (2005), New atlases for solar flux, irradiance, central intensity, and limb intensity, Mem. Soc. Astron. Ital. Suppl. 8, 189, bibcode 2005MSAIS...8..189K, no DOI [verified via the ADS record (web search); atlas file header checked].
- Luo, A.-L., Zhao, Y.-H., Zhao, G., et al. (2015), The first data release (DR1) of the LAMOST regular survey, RAA 15, 1095, doi:10.1088/1674-4527/15/8/002 [verified via Crossref].
- Morton, D. C. (2000), Atomic Data for Resonance Absorption Lines. II. Wavelengths Longward of the Lyman Limit for Heavy Elements, ApJS 130, 403, doi:10.1086/317349 [verified via Crossref].
- Ren, J.-J., Wu, H., Wu, C.-J., et al. (2021), On the radial velocity calibrations in the LAMOST medium-resolution spectroscopic survey of nebulae, RAA 21, 051, doi:10.1088/1674-4527/21/3/51, arXiv:2008.06236 [verified via Crossref and arXiv; the red-band sky-line list, "the MRS-N data reduction will only focus on the red band", "calibrated to vacuum, and also corrected to the heliocentric frame" (of the LAMOST arc-lamp calibration) and the 0.2–4 km s⁻¹ plate deviations were checked in the arXiv PDF].
- Schönrich, R., Binney, J., & Dehnen, W. (2010), Local kinematics and the local standard of rest, MNRAS 403, 1829, doi:10.1111/j.1365-2966.2010.16253.x [verified via Crossref].
- Wu, C.-J., Wu, H., Zhang, W., et al. (2022), The Data Processing of the LAMOST Medium-resolution Spectral Survey of Galactic Nebulae (LAMOST MRS-N Pipeline), RAA 22, 075015, doi:10.1088/1674-4527/ac7387, arXiv:2110.04547 [verified via Crossref and arXiv; "mainly for the red band" was checked in the arXiv PDF].
- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (2026), The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ 172, 298, doi:10.3847/1538-3881/aeb11c, arXiv:2610.13329 [verified via Crossref and arXiv; quotations are from the arXiv v1 text in `sources/arxiv/`].
