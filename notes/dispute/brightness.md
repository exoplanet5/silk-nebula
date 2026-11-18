## 5. Table 3 overestimates the [O III] surface brightness by a factor of ~5–9

**Claim in the paper.** Table 3 (p. 8) lists the [O III]-band surface brightness of P1, P2 and P3 as 23.60 ± 0.028, 23.87 ± 0.035 and 23.85 ± 0.039 mag arcsec⁻². The footnote says these were "derived by interpolation from Gaia and SDSS photometry". The zero points, 11.23–11.32, have no stated units.

Sect. 3.1 (p. 6) describes the method. The team measured "star-free sub-regions within P1–P3 and adjacent nebula-free background regions" and calibrated them on field stars whose Gaia and SDSS magnitudes were transformed to a common system following Jordi et al. (2010). It notes that P1 "appears the faintest in the narrowband image" yet "yields the highest measured surface brightness", and attributes this to stellar light.

Sect. 4.1 (p. 15) converts the values to "∼10–13 Rayleigh — well above the 0.5–5 Rayleigh detection limit of the Fesen et al. (2024) [O III] images". The survey camera (Sect. 2.1, p. 4) is a 50 mm f/0.95 lens with a 6.0 nm FWHM filter in its aperture and ~50″ pixels, under an SQM ≈ 21.5 mag arcsec⁻² sky; a 180 mm f/2.8 reflector was used for follow-up.

**What is right.** The arithmetic is correct. An AB magnitude (Oke & Gunn 1983) in a 60 Å band that contains only λ5007 gives 12.75 / 9.94 / 10.12 R [calc]. The filament is real, and Fesen et al. (2026) imaged it.

**Why it is wrong.**

1. **Conversion choices cannot rescue the values** (Fig. 5a). We varied three choices:
   - band width, 50–70 Å;
   - magnitude system, ±0.1 mag (Blanton & Roweis 2007, Table 1: AB − Vega = −0.08 mag in SDSS g, +0.02 mag in V);
   - λ4959 leakage, from none to full.

   Over these choices Table 3 spans 5.7–16.3 R [calc]. Even a 30 Å band gives 5.0–6.4 R. The paper does not say which instrument produced Table 3. With the 5 nm filter that QHYCCD (2024) reports for the f/2.8 follow-up, the values are 8.3–10.6 R.
2. **Independent measurements give 1.4–2.5 R** (Table 6, Fig. 5b). Like Table 3, these are excesses over the local background; the SDSS P-region means are absolute, but the off-strand level is only 0.07 R.
   - DESI DR1 (DESI Collaboration 2026): 1.62 R (southern-ridge median; calibration uncertain by about 20 %, see below).
   - SDSS/BOSS sky fibres (DR17, DR19; Abdurro'uf et al. 2022; SDSS Collaboration 2026): 1.49 / 1.78 / 1.87 R at P1 / P2 / P3.
   - SDSS imaging sky levels (DR18; Almeida et al. 2023): 2.31 ± 0.38 R on 17′ scales (3.9σ against coherent nulls).
   - NSNS DR0.2 maps (Ziegenbalg 2025; calibration not re-verified), measured here: 2.04 / 2.06 / 2.18 R, at a 51.5″ pixel scale.

   The DESI calibration is less certain than the ±8 % of the measured fibre acceptances. DR1 has no guide-camera fibre acceptance for the two survey-validation tiles 80644 and 80683, which supply 85 % of the fibres within 0.5° of P1, all of those near P3 and 46 % of the usable fibres at Dec < 19.5°; programme medians were used instead (`notes/paper_appendix_desi.md`). The DESI brightness near P1 and P3, and the 1.62 R southern median, are therefore uncertain by about 20 %. The factor to Table 3 is unaffected: Table 3 is 6.2 and 6.6 times the DESI model sky brightness at P1 and P3 (2.07 ± 0.06 and 1.54 ± 0.04 R).

   Over the same regions Table 3 is **4.6–8.6× higher** [calc]. That is 1.7–2.3 mag, or 43–83 times its quoted errors. Without the star-affected P1 the factor is 4.6–6.2×.
3. **Angular scale does not explain the gap.** The Table 3 subregion averages are 2.6–4.3× above even the brightest single measurements [calc]: the brightest SDSS fibre (3.65 R), the brightest DESI node (3.88 R) and the NSNS 1′ peaks (2.95–3.64 R).
4. **The team's own earlier value is fainter.** QHYCCD (2024) reports the SWONS result as 24.3 mag arcsec⁻² for the *brightest* region, calibrated on V-band stars. That is 3.3–6.7 R for a 30–60 Å band, and 0.43–0.70 mag fainter than the Table 3 *region* values [calc].
5. **No identified physical multiplicative factor is large enough** (Fig. 5d). As I_true/I_Table3 [calc]:
   - AB vs Vega: 0.91–1.10.
   - Band width: 0.83–1.17.
   - λ4959 leakage: 0.75–1.00.
   - Filter in the lens aperture (Sect. 2.1): 1.01–1.27 at field angles 0–5°, and 0.93–2.8 at 10–15°. These come from a generic 60 Å 3-cavity filter model, not the real SWONS curve.
   - Filter in an f/2.8 beam: 0.97–1.08.
   - An untuned filter in the f/0.95 converging beam gives 2.0–4.0, which would make Table 3 an *under*estimate.

   As a deliberately generous bound we stacked the separate lowering extremes: 50 Å, +0.1 mag and full λ4959. No single filter satisfies these together, and the filter model already contains the width and λ4959 terms. We added ~0.3 mag of P1 starlight. The result, 5.5–5.8 R, is still **2.6–3.7×** above the independent values. Matching them needs a factor of 0.12–0.22.

   A common factor would also keep the regions in order, but they are not in order.
   - Table 3 ranks P1 brightest, by 0.25–0.27 mag, which Sect. 3.1 attributes to starlight.
   - SDSS ranks P1 faintest, by 0.29–0.38 R (errors ≤ 0.09 R).
   - NSNS finds the three regions equal within 0.14 R.
   - The paper's own image shows P1 faintest.

   So at least the P1 excess, 2.8–3.1 R, is additive, as Sect. 3.1 itself suggests.
6. **An additive error of the required size is plausible** (Fig. 5c). The excess is 7.9–11.3 R [calc]. Two sources could supply it.
   - **Background or flat-field residual.** The excess is 6–20 % of an in-band sky of 21.1–22.0 mag arcsec⁻², or 9–13 % at the paper's SQM reading of 21.5. SDSS/BOSS sky spectra give 21.1–21.5 AB at 5070 Å; the 22.0 limit allows for a darker narrow-band continuum. A residual this size is large for nearby apertures. However, the paper notes that dome and twilight flats "introduced significant vignetting" in this f/0.95 system (Sect. 2.1).
   - **Mismatched starlight.** Gaia DR3 stars with G > 14 within 0.25° of P1–P3 add up to 3.0–4.4 R-equivalent (Gaia Collaboration 2023; calc). Even with a perfectly star-free background aperture, that is only 27–56 % of the excess, so brighter stars or their halos would also be needed. P1 has 15.5–20.6 R-equivalent more starlight than P2/P3, and its extra 2.8–3.1 R is 13–20 % of that surplus.

   Without the frames and aperture definitions we cannot say which additive term dominates. Nor can we exclude a bookkeeping error: aperture area, pixel scale, exposure normalisation or the unitless zero points. The identified calibration factors, however, cannot account for the excess.

![Four-panel check of Table 3 against independent [O III] data](figures/dispute_brightness_5.png)

*Figure 5. [O III] λ5007 surface brightness of the Silk Nebula: Wu et al. (2026) Table 3 compared with independent data.*
- *(a) Conversion of an [O III]-band AB surface brightness to I(5007) for 30–70 Å effective widths, with the band containing only the line. The white line is 60 Å, the stated SWONS FWHM. Grey shading is ±0.1 mag (AB/Vega) at 60 Å. Red symbols are Table 3: P1 circle, P2 square, P3 triangle. Orange is the QHYCCD-note 24.3 mag arcsec⁻² (bar 30–60 Å, diamond 50 Å). Green is the independent 1.4–2.5 R.*
- *(b) All measurements on one scale, with the angular scale of each at right. Symbol shape gives the region (P1/P2/P3) or a ridge-wide value (diamond; the orange diamond is the QHYCCD value at 50 Å). Filled symbols are region means or medians. Open symbols are peaks or alternative estimators: the star is the brightest SDSS fibre, and the DESI open square is the W-strand node at Dec 18° near P2. The NSNS P1 and P2 symbols overlap. The hatched band is the 0.5–5 R range quoted in Sect. 4.1.*
- *(c) Additive background error, as a percentage of the in-band sky, needed to raise I_true to the Table 3 mean of 10.9 R: (I_T3 − I_true)/I_sky. The box spans the measured I_true and plausible skies, 7–17 % for the mean and 6–20 % region by region (P1–P3 matched to NSNS or SDSS). The dashed line is the paper's SQM reading.*
- *(d) Multiplicative systematics as I_true/I_Table3. Filter cases use a generic 3-cavity, 60 Å FWHM interference-filter model (n\* 1.8–2.05, CWL 5007–5030 Å), not the measured SWONS curve. "At stop" is the filter in the lens aperture (Sect. 2.1), where each field point sees a collimated beam. The green band is the required factor. The stack combines the separate width, magnitude-system and λ4959 lowering extremes with ~0.3 mag of P1 starlight.*

*Data: the paper's Sects. 2.1, 3.1, 4.1 and Table 3; DESI DR1; SDSS DR17/DR19 sky fibres and DR18 imaging-sky levels; NSNS DR0.2; Gaia DR3. Computed by `scripts/dispute_brightness_calc.py` and plotted by `scripts/dispute_brightness_fig.py`.*

**Table 6.** [O III] λ5007 brightness of the southern Silk Nebula by method. The last column is the Table 3 value divided by the measurement, region-matched where possible. Values for P1 / P2 / P3 are separated by "/". Sources: the paper; `data/derived/dispute_brightness_calc.json`; the project notes `verify_desi.md`, `paper_appendix_desi.md`, `verify_sdss.md`, `verify_precovery.md`, `imaging_optical.md` and `spectra_sdss_lamost.md`.

| Measurement | Method | Scale | I(5007) (R) | Table 3 / value |
|---|---|---|---|---|
| Wu+26 Table 3 | SWONS narrow-band photometry, Gaia/SDSS zero points; converted here (AB, 60 Å, line only) | star-free subregions, ~50″ px | 12.75 / 9.94 / 10.12 (all variants 5.7–16.3) | 1 |
| Table 3, all lowering extremes stacked | 50 Å, +0.1 mag, full λ4959, P1 −0.3 mag, applied together | — | 5.50 / 5.66 / 5.76 | 2.3 / 1.8 / 1.8 |
| QHYCCD note (SWONS team result, 2024) | 24.3 mag arcsec⁻², V-calibrated | brightest region | 3.3–6.7 (30–60 Å) | 1.5–3.8 |
| LAMOST MRS-N relative × sky-line scale | Table 2 ratios × [N I] 5198 / OH 6553 sky-line scale, APO or UVES (Hanuschik 2003) ratio | fibres | 2.8–3.2 (APO) or 14–16 (UVES) | indeterminate; not used |
| NSNS DR0.2 | calibrated narrow-band survey, on-ridge minus off-ridge (not re-verified) | ≤ 0.25° of ridge, 51.5″ px | 2.04 / 2.06 / 2.18 (1′ peaks 2.95 / 3.64 / 3.37) | 6.3 / 4.8 / 4.6 (peaks 4.3 / 2.7 / 3.0) |
| SDSS imaging sky, g (19 runs, 2002–07) | ridge-template fit to field sky levels; 3.9σ vs coherent nulls | 17′ | 2.31 ± 0.38 (2.39, 2.55 other estimators) | 3.9–5.5 |
| SDSS/BOSS sky fibres | total sky, fibre-acceptance corrected, mean within 0.3° | 2–3″ fibres | 1.49 / 1.78 / 1.87 (±0.05–0.09; ±0.2 zero point); brightest 3.65 | 8.6 / 5.6 / 5.4 |
| SDSS/BOSS south ON−OFF | strand vs 0.45–2° background | r_on 0.15° | 1.39 ± 0.10 | 7.2–9.2 |
| DESI DR1 south ridge | flux-calibrated, GFA fibre acceptance (programme medians for two tiles) | 1.5″ fibres | 1.62 (±20 % sys; nodes 1.4–2.3; W node 3.88) | 6.1–7.9 |
| Fesen et al. 2024 limits (as quoted) | — | — | 0.5–5 | — |

**Consequence.**
- **The Sect. 4.1 comparison fails.** At 1.4–2.5 R the Silk Nebula lies *inside* the quoted 0.5–5 R range of Fesen et al. (2024), near its logarithmic midpoint of 1.6 R. It is not "well above" that range. Its appearance in Fesen et al. (2026) is consistent with this brightness and does not require a bright filament.
- **Downstream quantities.** The paper's energy uses an adopted n₀ = 1.0 ± 0.3 cm⁻³, not Table 3. Anything scaled from Table 3 would be inflated [calc]: emission measure and [O III] luminosity by 4.6–8.6×, an emission-measure density by 2.2–2.9×, and the Sect. 4.1 energy (E ∝ n₀^1.16) by 6–12× if n₀ were tied to I. The "tenuous" of Sect. 3.1 is, if anything, strengthened.
- **Under the paper's own shock interpretation, published shock models make Table 3 inconsistent with the GALEX FUV.**
  - Regressing GALEX FUV (Martin et al. 2005; Morrissey et al. 2007) on NSNS gives 92 ± 10 CU (photons s⁻¹ cm⁻² sr⁻¹ Å⁻¹) per NSNS rayleigh in the south (`verify_xrayfuv_galex_measure.json`).
  - If Table 3 were right, that would be 15–20 CU per true rayleigh [calc]. This is 2.9–3.8× below (2.6× at +1σ) the lowest non-resonant FUV floor of the tabulated shock models at E(B−V) = 0.03.
  - That floor is 57 CU R⁻¹ for the incomplete shocks of Raymond et al. (1988). For the complete shocks of Raymond (1979), Cox & Raymond (1985) and Dopita & Sutherland (1996) it is 85–277 CU R⁻¹ (up to 281 with the magnetised Dopita & Sutherland runs).
  - Resonance scattering cannot remove non-resonant lines. Closing the gap would need E(B−V) ≳ 0.17–0.21 (a model-only extrapolation that ignores the [O III] attenuation, so the true requirement is higher). Schlegel et al. (1998) give 0.03–0.04.

**Requested correction.**
1. Publish an erratum stating that the Table 3 values are 4.6–8.6× brighter than independent measurements of the same regions. Either re-reduce the data, or label the values uncalibrated and withdraw the conversion to rayleighs.
2. Document the photometry:
   - the instrument and filter (50 mm f/0.95 with 6 nm, or 180 mm f/2.8 with 5 nm), the measured filter curve and the beam geometry;
   - the magnitude system and the effective band width;
   - the positions and sizes of the source and background apertures, and the star masking;
   - the flat-field and background model;
   - the units of the zero points 11.23–11.32.
3. Replace the Sect. 4.1 sentence with, for example: "Flux-calibrated DESI and SDSS spectra give an [O III] λ5007 surface brightness of ≈1.4–1.9 R on the southern filament, and SDSS and NSNS imaging give ≈2.0–2.6 R. This lies within the 0.5–5 R range of the Fesen et al. (2024) [O III] images, consistent with the filament's appearance in Fesen et al. (2026)."

**References used in this section.**
- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (2026), The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ 172, 298, doi:10.3847/1538-3881/aeb11c; arXiv:2610.13329 [verified via Crossref and arXiv]
- Fesen, R. A., Drechsler, M., Strottner, X., Falls, B., et al. (2024), Deep Optical Emission-line Images of Nine Known and Three New Galactic Supernova Remnants, ApJS 272, 36, doi:10.3847/1538-4365/ad410a [verified via Crossref]
- Fesen, R. A., Drechsler, M., Falls, B., Schaeffer, T., et al. (2026), Optical Discovery of New and Candidate Galactic Supernova Remnants Plus Optical Imaging of the Monogem Supernova Remnant, arXiv:2608.26666 [verified via arXiv]
- Ziegenbalg, S. (2025), Northern Sky Narrowband Survey DR0.2: First [O III] and [S II] Results, RNAAS 9, 227, doi:10.3847/2515-5172/adfec7 [verified via Crossref]
- DESI Collaboration, Abdul Karim, M., et al. (2026), Data Release 1 of the Dark Energy Spectroscopic Instrument, AJ 171, 285, doi:10.3847/1538-3881/ae4c43; arXiv:2503.14745 [verified via Crossref and arXiv; first posted to arXiv in 2025]
- Almeida, A., Anderson, S. F., Argudo-Fernández, M., et al. (2023), The Eighteenth Data Release of the Sloan Digital Sky Surveys: Targeting and First Spectra from SDSS-V, ApJS 267, 44, doi:10.3847/1538-4365/acda98 [verified via Crossref]
- SDSS Collaboration, Aghakhanloo, M., Aird, J., et al. (2026), The Nineteenth Data Release of the Sloan Digital Sky Survey, ApJS 285, 9, doi:10.3847/1538-4365/ae4697; arXiv:2507.07093 [verified via Crossref and arXiv; first posted to arXiv in 2025; the journal author list is alphabetical after the collaboration]
- Abdurro'uf, Accetta, K., Aerts, C., et al. (2022), The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Release of MaNGA, MaStar, and APOGEE-2 Data, ApJS 259, 35, doi:10.3847/1538-4365/ac4414 [verified via Crossref]
- Gaia Collaboration, Vallenari, A., Brown, A. G. A., et al. (2023), Gaia Data Release 3: Summary of the content and survey properties, A&A 674, A1, doi:10.1051/0004-6361/202243940 [verified via Crossref]
- Jordi, C., Gebran, M., Carrasco, J. M., et al. (2010), Gaia broad band photometry, A&A 523, A48, doi:10.1051/0004-6361/201015441 [verified via Crossref]
- Oke, J. B., & Gunn, J. E. (1983), Secondary standard stars for absolute spectrophotometry, ApJ 266, 713, doi:10.1086/160817 [verified via Crossref]
- Blanton, M. R., & Roweis, S. (2007), K-Corrections and Filter Transformations in the Ultraviolet, Optical, and Near-Infrared, AJ 133, 734, doi:10.1086/510127 [verified via Crossref; Table 1 values read from the arXiv PDF, astro-ph/0606170]
- Martin, D. C., Fanson, J., Schiminovich, D., et al. (2005), The Galaxy Evolution Explorer: A Space Ultraviolet Survey Mission, ApJ 619, L1, doi:10.1086/426387 [verified via Crossref]
- Morrissey, P., Conrow, T., Barlow, T. A., et al. (2007), The Calibration and Data Products of GALEX, ApJS 173, 682, doi:10.1086/520512 [verified via Crossref]
- Raymond, J. C. (1979), Shock waves in the interstellar medium, ApJS 39, 1, doi:10.1086/190562 [verified via Crossref]
- Cox, D. P., & Raymond, J. C. (1985), Preionization-dependent families of radiative shock waves, ApJ 298, 651, doi:10.1086/163649 [verified via Crossref]
- Raymond, J. C., Hester, J. J., Cox, D., et al. (1988), Spatial and spectral interpretation of a bright filament in the Cygnus Loop, ApJ 324, 869, doi:10.1086/165945 [verified via Crossref]
- Dopita, M. A., & Sutherland, R. S. (1996), Spectral Signatures of Fast Shocks. I. Low-Density Model Grid, ApJS 102, 161, doi:10.1086/192255 [verified via Crossref]
- Schlegel, D. J., Finkbeiner, D. P., & Davis, M. (1998), Maps of Dust Infrared Emission for Use in Estimation of Reddening and Cosmic Microwave Background Radiation Foregrounds, ApJ 500, 525, doi:10.1086/305772 [verified via Crossref]
- Hanuschik, R. W. (2003), A flux-calibrated, high-resolution atlas of optical sky emission from UVES, A&A 407, 1157, doi:10.1051/0004-6361:20030885 [verified via Crossref]
- QHYCCD (2024), Silk Nebula Discovered by Enthusiast Astronomers with QHYCCD Cameras, web note (2024-07-12), https://www.qhyccd.com/silk-nebula-discovered-by-enthusiast-astronomers-with-qhyccd-cameras/ [verified via the archived copy `sources/web/discovery_qhyccd_pr.md`, accessed 2026-11-16; not peer-reviewed]

*Provenance: [calc] = `scripts/dispute_brightness_calc.py` → `data/derived/dispute_brightness_calc.json`. The magnitude-to-rayleigh conversion and the filter model are imported unchanged from `scripts/verify_paper_energetics.py`. The Crossref and arXiv records are in `data/derived/dispute_brightness_refs.json` (`scripts/dispute_brightness_refs.py`).*
