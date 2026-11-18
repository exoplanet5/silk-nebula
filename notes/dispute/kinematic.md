## 2. Why no kinematic distance is possible toward the Silk Nebula

**Claim in the paper.** Sect. 3.3 (p. 10) derives kinematic distances with the `kd` package (Wenger et al. 2018) and the Reid et al. (2014) rotation curve. It excludes objects near the Galactic anticentre (GAC) "owing to the negligible contribution of circular rotation to VLSR along these sightlines" and concedes that the Silk "lies near its boundary, which introduces significant uncertainties into the derived kinematic distance". Even so, it adopts "a fiducial distance of ∼160 pc (the approximate median of the six [O III]-based kinematic distance estimates; see Table 4)" (Sect. 3.3, p. 11). The abstract and Sect. 5 (p. 17) repeat it, and Sects. 3.3 and 4.1 (pp. 11, 15–16) scale size, energy and age from it. For P1 [O III], Table 4 (p. 13) lists RV = 26.8 ± 4.3 km s⁻¹, V_LSR = +15.5 km s⁻¹ (no error printed; we adopt the RV error), a Monte Carlo distance of 1.63 (+1.01/−0.84) and a traditional distance of 1.79. Section 1 deals with the units. Here we show that, in either unit, the velocity carries no useful distance information. The paper flags the uncertainty; our point is that it is total.

**Why it is wrong.** Unless another source is named, every number below comes from `scripts/dispute_kinematic_kd.py` → `data/derived/dispute_kinematic_kd.json`, which re-implements the paper's rotation model from the `kd` source (Persic et al. 1996 curve, Reid et al. 2014 A5 parameters, `kd` projection, IAU standard solar motion). This model reproduces Table 4 (Section 1).

1. **Rotation produces almost no velocity here** (Fig. 2a).
   - Across P1–P4 the mean slope dV_LSR/dd over 0–2 kpc is **7.6–9.1 km s⁻¹ kpc⁻¹** (`kd` form, ×cos b projection or flat curve); as d → 0 it is 9.1–10.7.
   - Along the whole [O III] spine it is 4.3–9.8, lowest at the northern end (ℓ = 191.4°).
   - Wenger et al. (2018) define the GAC as ℓ = 180°, where "LSR velocities due to circular motion tend towards zero". Their uncertainty maps mask 160° < ℓ < 200° (captions of their Figs. 5, 8, 19). P1–P4 lie only 2.5–5.5° outside that band. The spine crosses ℓ = 200° at Dec +20.75°: **8.1° of its 16.9° lies inside the band**, 5.9° of it within the paper's own Fig. 1 (Dec ≤ +26.3°; `notes/geometry.md`). Read as a 20° radius about (180°, 0°), the zone misses the spine (24.5–32.8° away), but the gradient is equally small.

2. **Ordinary peculiar velocities equal the whole rotational signal.**
   - Giant molecular clouds within 3 kpc have a one-dimensional cloud-to-cloud dispersion of 7.8 (+0.6/−0.5) km s⁻¹, only ~20 % less without small-scale streaming (Stark & Brand 1989). Anderson et al. (2012) adopt 7 km s⁻¹ streaming motions, which "may be 5 to 10 km s⁻¹". "The vast majority" of masers have peculiar motions ≲ 10 km s⁻¹ (Reid et al. 2019); super-bubbles can accelerate clouds to ≈ 20 km s⁻¹ (Reid et al. 2014).
   - At the slopes of step 1, σ_pec = 5 / 7 / 10 km s⁻¹ means **Δd = 0.55–0.66 / 0.77–0.92 / 1.10–1.32 kpc**.
   - Rotation contributes only 1.4–1.7 km s⁻¹ at 160 pc and 2.6–3.1 km s⁻¹ at 300 pc toward P1–P4 (≤ 3.4 anywhere on the spine). A nearby object's V_LSR is essentially its peculiar velocity.

3. **Gas at a known small distance gets absurd kinematic distances** (Fig. 2b).
   - In the Redfield & Linsky (2008) model the P1 sightline crosses the LIC and Gem clouds, both within 15 pc (rotation there: 0.15 km s⁻¹). Their projected V_helio = +19.1 ± 1.1 and +33.9 ± 1.1 km s⁻¹ (calculator output, `data/raw/spectra_ism_lism/dynlism_P1.html`) are V_LSR = +7.0 and +21.8 km s⁻¹.
   - The paper's curve puts them at **0.74 and 2.70 kpc**, 50 and 180 times their maximum distance.
   - The LIC velocity is observed: Ca II K at V_helio +19.1 and +19.4 km s⁻¹ toward 8 Cnc (63 pc, 4.7° from the Silk; `notes/spectra_ism.md`). Gas within 63 pc thus gets a formal distance ≈ 12 times too large. The Gem velocity is a model projection, not detected in the available Ca II/Na I spectra.

4. **The posterior is set by σ_pec, not by the data** (Fig. 2c, d; Table 3).
   - For the paper's +15.5 ± 4.3 km s⁻¹, the 16–84 % range widens from 1.28–2.51 kpc (σ_pec = 0) to 1.01–3.31 (σ_pec = 7) and 0.94–3.85 kpc (σ_pec = 10). For the corrected +1.2 ± 1.0 km s⁻¹ (Section 3) it widens from 0.06–0.23 to 0.18–1.26 and 0.26–1.90 kpc. At σ_pec = 7–10 every interval spans a factor of 3.3–7.3.
   - For +1.2 at σ_pec = 7, 44 % of Monte Carlo draws fall below V_rot(0) = 0 and have no solution. The median is 0.12 kpc if they are kept at d = 0 and 0.60 kpc in the posterior; `kd` discards them, biasing it upward.
   - The `kd` Monte Carlo behind Table 4 adds no peculiar-velocity term: Wenger et al. (2018) argue that the Reid et al. (2014) parameter errors "inherently" include streaming. That ensemble argument does not cover the random velocity of one nebula.

5. **The velocity field along the filament does not follow rotation.**
   - At a common distance of 0.3–1.6 kpc, rotation predicts the northern spine (Dec +24°) to be 0.8–4.1 km s⁻¹ *lower* than P1.
   - DESI measures it *higher*: ΔV_LSR ≈ +7…+11 (± ~2.5) km s⁻¹, i.e. +4…+9 km s⁻¹ heliocentric at 1.8–3.7σ across the six estimators (`notes/w3_northsouth.md`). SDSS gives a consistent but insignificant +6.7 ± 6.5 km s⁻¹, and the surveys disagree at Dec +17.6…+21°, so this is a ~2–3σ result.
   - Read kinematically, one continuous arc would lie at 0–0.12 kpc in the south (DESI −2.8…−0.2 has no solution; corrected P1 +1.2 gives 0.12) and at 0.86–1.68 kpc in the north.

6. **The paper's own text implies the same.** Sect. 3.3 (p. 10) lists "non-circular motions" as a systematic. Sect. 4.2 (p. 17) proposes that the "∼5–15 km s⁻¹" [O III]/low-ionization offset arises in "different layers of the same expanding shell" (expanding at "∼140 km s⁻¹") or in an inhomogeneously expanding bubble. Section 3 shows the offset is an artefact, but on the paper's own reading either mechanism is a non-circular motion as large as the entire rotational V_LSR at 1.6 kpc toward P1 (14.1 km s⁻¹). Sect. 5 (p. 17) allows "certain systematic uncertainties" yet still adopts the kinematic value.

![Four-panel figure: (a) map of the rotation-velocity gradient over l 180-230, b 0-30 with the Silk spine, P1-P4 and the anticentre exclusion boundaries; (b) rotation-only V_LSR versus distance toward P1 and the northern spine with the paper's, corrected and northern [O III] velocities, and the LIC/Gem local clouds with their formal kinematic distances; (c) and (d) distance posteriors for V_LSR = +15.5 and +1.2 km/s for sigma_pec = 0, 3, 7 and 10 km/s.](figures/dispute_kinematic_fig2.png)

*Figure 2a–d. Kinematic distances toward the Silk Nebula.*
- *Model and code (all panels): the Persic et al. (1996) universal curve with the Reid et al. (2014) A5 parameters as tabulated in `kd`; the `kd` projection; the IAU standard solar motion. Code: `scripts/dispute_kinematic_kd.py`. Numbers: `data/derived/dispute_kinematic_kd.json`.*
- *(a) Mean circular-rotation gradient dV_LSR/dd over 0–2 kpc for ℓ = 180–230°, b = 0–30°, with contours in km s⁻¹ kpc⁻¹. The orange line is the project's [O III] spine (`data/derived/silk_trace.json`), and the white dots are the paper's P1–P4 (Table 1). The dashed line is ℓ = 200°, the edge of the 160° < ℓ < 200° band that Wenger et al. (2018) mask (their Figs. 5, 8 and 19). The dotted line is a 20° radius about (ℓ, b) = (180°, 0°).*
- *(b) Rotation-only V_LSR(d) toward P1 (solid line; grey band ±7 km s⁻¹) and toward the northern spine at Dec +24° (dashed line). The horizontal bands are the paper's P1 [O III] V_LSR (Table 4), the corrected [O III] V_LSR (Section 3; `verify_velocity_summary.json`) and the northern DESI [O III] V_LSR range (`notes/w3_northsouth.md`). The diamonds are the LIC and Gem clouds of the Redfield & Linsky (2008) model, both within 15 pc. Dotted lines lead to their formal (traditional) kinematic distances (open circles).*
- *(c, d) Posterior P(d | V_LSR) toward P1 for σ_pec = 0, 3, 7 and 10 km s⁻¹. The prior is flat on 0–6 kpc and the likelihood Gaussian with σ² = σ_V² + σ_pec². The rotation-curve parameters are marginalised over 800 draws, and the curves are normalised to unit peak. Bars mark 16–84 % intervals and ticks the medians. The orange bar in (c) is the paper's Table 4 Monte Carlo value for P1 [O III], 1.63 (+1.01/−0.84), read as kpc.*

**Table 3.** Kinematic-distance ranges toward P1 (ℓ = 203.8°, b = +19.5°) for each velocity input and peculiar-velocity dispersion σ_pec. "Posterior" uses a flat prior on 0–6 kpc and marginalises over the Reid et al. (2014) parameters. Moving the prior edge to 4 or 10 kpc changes the +15.5 rows at σ_pec = 7 to 1.92 [0.96–2.99] or 2.05 [1.02–3.37], and at σ_pec = 10 to 1.95 [0.84–3.14] or 2.28 [0.97–4.11]. Every other row changes by ≤ 0.05 kpc. "`kd`-style MC" is Wenger et al. (2018) Method C plus a Gaussian σ_pec term; velocities below V_rot(0) are kept at d = 0, whereas `kd` discards them. The lower block lists single-valued (traditional) distances for gas whose location is known independently, and for the two ends of the [O III] spine. All values are computed in `scripts/dispute_kinematic_kd.py`.

| V_LSR input (km s⁻¹; source) | σ_pec (km s⁻¹) | Posterior median [16–84 %] (kpc) | P(d < 0.3 kpc) | P(d > 1 kpc) | `kd`-style MC median [16–84 %] (kpc) |
|---|---|---|---|---|---|
| +15.5 ± 4.3 (Table 4, P1 [O III]; ± = RV error) | 0 | 1.85 [1.28–2.51] | 0.001 | 0.94 | 1.77 [1.22–2.42] |
| | 3 | 1.89 [1.20–2.69] | 0.005 | 0.90 | 1.77 [1.11–2.55] |
| | 7 | 2.03 [1.01–3.31] | 0.025 | 0.84 | 1.78 [0.77–3.01] |
| | 10 | 2.20 [0.94–3.85] | 0.037 | 0.83 | 1.78 [0.48–3.46] |
| *Paper, Table 4 MC (P1 [O III])* | *none* | *1.63 [0.79–2.64] (mode)* | — | — | — |
| +1.2 ± 1.0 (corrected [O III], Section 3) | 0 | 0.14 [0.06–0.23] | 0.95 | 0.00 | 0.12 [0.02–0.22] |
| | 3 | 0.28 [0.09–0.56] | 0.53 | 0.01 | 0.12 [0.00–0.44] |
| | 7 | 0.60 [0.18–1.26] | 0.26 | 0.26 | 0.12 [0.00–0.88] |
| | 10 | 0.87 [0.26–1.90] | 0.18 | 0.44 | 0.12 [0.00–1.25] |
| *Traditional distance, σ_pec = 0* | | | | | |
| LIC cloud, +7.0 (true d < 15 pc; seen at ≤ 63 pc) | 0 | 0.74 | — | — | — |
| Gem cloud, +21.8 (true d < 15 pc) | 0 | 2.70 | — | — | — |
| North [O III] at Dec +24°, +5.7…+10.4 (DESI) | 0 | 0.86–1.68 | — | — | — |
| South [O III] at P1, +1.2 (corrected) / −2.8…−0.2 (DESI) | 0 | 0.12 / no solution | — | — | — |

**Consequence.**
- The paper has no measured distance. Neither "∼160 pc" nor its kpc reading (Section 1) is supported, nor is anything scaled from it:
  - the ∼28 pc extent and the ∼30 pc radius;
  - E = (2.0 ± 0.7)×10⁵¹ erg and t = (6.8 ± 1.8)×10⁴ yr (Sect. 4.1);
  - the "distance discrepancy of 100–150 pc" between the blue- and red-arm lines (Sect. 4.2, p. 16).
- The kinematics cannot distinguish 160 pc, ~280 pc and 1.6 kpc; the corrected velocity fits gas at tens of parsecs equally well.
- What can constrain the distance (`REPORT.md` §5.5–6; `notes/w4_distance.md` §3–4):
  - **High-ion (C IV, Si IV) absorption toward hot stars of known parallax at 120–800 pc.** Optical Ca II/Na I cannot isolate the Silk, because its corrected velocity equals that of the local gas.
  - **A two-epoch [O III] proper motion plus a shock velocity**: d (kpc) = v_s/(4.74 μ), μ in mas yr⁻¹.
  - **A physical tie to the Monogem rim.** At Dec < 17.6° the spine lies 12.10 ± 0.12° from the Knies et al. (2024) two-remnant centre, implying 252–318 pc (68 %) only if the association holds and the shell is centred at the parallax distance of PSR B0656+14. The southern spine is straighter than the Monogem circle (radius of curvature ≈ 67°) and only tangent to it.

**Requested correction.**
1. Withdraw the kinematic distance as a result in the abstract, Sect. 3.3, Sects. 4.1–4.2 and Sect. 5. State instead:

   > "Toward (ℓ, b) ≈ (204°, +19°) the circular-rotation gradient is only ≈ 8–9 km s⁻¹ kpc⁻¹ (≈ 4 km s⁻¹ kpc⁻¹ at the northern end, which lies inside the 160° < ℓ < 200° anticentre zone). This is comparable to typical interstellar peculiar velocities of 5–10 km s⁻¹. Kinematic distances are therefore not meaningful for the Silk Nebula, and its distance is not constrained by the present data."

2. In Table 4, either delete both distance columns or label them "kpc; formal values without peculiar motions; not physically meaningful".
3. Give size, radius, energy and age as explicit functions of d. For example, the spine length is ≈ 0.30 d (16.9° of arc).
4. Delete the "100–150 pc" blue/red distance discrepancy (Sect. 4.2).

**References used in this section.**
- Anderson, L. D., Bania, T. M., Balser, D. S., & Rood, R. T. (2012), The Green Bank Telescope H II Region Discovery Survey. III. Kinematic Distances, ApJ 754, 62, doi:10.1088/0004-637X/754/1/62, arXiv:1205.4228 [verified via Crossref; the "7 km s⁻¹ … may be 5 to 10 km s⁻¹" wording checked in the arXiv full text].
- Knies, J. R., Sasaki, M., Becker, W., Liu, T., et al. (2024), A new understanding of the Gemini-Monoceros X-ray enhancement from discoveries with eROSITA, A&A 688, A90, doi:10.1051/0004-6361/202348834, arXiv:2401.17289 [verified via Crossref and arXiv; the two-remnant centre is taken from the project's `notes/w3_canon.md` and `notes/context.md`].
- Persic, M., Salucci, P., & Stel, F. (1996), The universal rotation curve of spiral galaxies — I. The dark matter connection, MNRAS 281, 27, doi:10.1093/mnras/278.1.27 [verified via Crossref].
- Redfield, S., & Linsky, J. L. (2008), The Structure of the Local Interstellar Medium. IV. Dynamics, Morphology, Physical Properties, and Implications of Cloud-Cloud Interactions, ApJ 673, 283, doi:10.1086/524002, arXiv:0709.4480 [verified via Crossref and the arXiv abstract ("15 warm clouds located within 15 pc of the Sun"); projected velocities from the authors' kinematic calculator, output saved in `data/raw/spectra_ism_lism/`].
- Reid, M. J., Menten, K. M., Brunthaler, A., Zheng, X. W., et al. (2014), Trigonometric Parallaxes of High Mass Star Forming Regions: The Structure and Kinematics of the Milky Way, ApJ 783, 130, doi:10.1088/0004-637X/783/2/130, arXiv:1401.5377 [verified via Crossref; the "≈ 20 km s⁻¹" super-bubble statement (citing Sato et al. 2008) checked in the arXiv full text].
- Reid, M. J., Menten, K. M., Brunthaler, A., Zheng, X. W., et al. (2019), Trigonometric Parallaxes of High-mass Star-forming Regions: Our View of the Milky Way, ApJ 885, 131, doi:10.3847/1538-4357/ab4a11, arXiv:1910.03357 [verified via Crossref; "the vast majority of sources have moderate peculiar motions of ≲ 10 km s⁻¹" (Sect. 5) checked in the arXiv full text].
- Stark, A. A., & Brand, J. (1989), Kinematics of Molecular Clouds. II. New Data on Nearby Giant Molecular Clouds, ApJ 339, 763, doi:10.1086/167334 [verified via Crossref and the ADS scan of the abstract: σ_v = 7.8 (+0.6/−0.5) km s⁻¹ within 3 kpc, ~20 % smaller without small-scale streaming].
- Wenger, T. V., Balser, D. S., Anderson, L. D., & Bania, T. M. (2018), Kinematic Distances: A Monte Carlo Method, ApJ 856, 52, doi:10.3847/1538-4357/aaaec8, arXiv:1802.04203 [verified via Crossref; the GC/GAC definition, the "within 15° of the GC and 20° of the GAC" exclusion, the 160° < ℓ < 200° mask (Figs. 5, 8, 19) and "we do not add any additional streaming motion uncertainty" checked in the arXiv full text].
- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (2026), The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ 172, 298, doi:10.3847/1538-3881/aeb11c, arXiv:2610.13329 [verified via Crossref and arXiv; quotations from the arXiv v1 text in `sources/arxiv/`].
