# Wave 2 verification: the [O III] velocity discrepancy (tag `velocity`)

Date: 2026-11-16. Velocities are heliocentric (V_helio) unless marked V_LSR. V_LSR − V_helio = −12.19 km/s
toward the nebula for the standard solar motion the paper uses (20 km/s toward 18h03m50.29s +30°00′16.8″). For
Schönrich+2010 the value is −11.88. Labels: **measured** = fitted here; **derived** = computed from stated
inputs; **literature** = quoted from a paper or a teammate's file.

## 1. Verdict

**The paper's [O III] velocity offset is a calibration artefact, not a property of the nebula.**

The decisive test uses the paper's own spectrum. Fig. 3 of arXiv:2610.13329v1 is stored in the PDF as vector
data: 3,590 samples at 0.106 Å, the native MRS sampling. That spectrum contains the [N I] 5198/5200 airglow
doublet, which is at topocentric rest. Measuring [O III] against it needs no wavelength calibration at all.
- **Measured:** V_topo([O III]) = +21.9 km/s. Adding v_corr = −8.32 km/s gives **V_helio = +13.6 ± 1.2 km/s**.
- **The published scale is topocentric and in air.** On the paper's plotted wavelength scale, the airglow lines
  sit at +1.8 ± 0.8 km/s of their air rest wavelengths. LAMOST delivers vacuum, heliocentric spectra
  (`VACUUM=T`, `HELIO=T`). So the paper's processing re-referenced the blue arm to air sky wavelengths, which
  undid the −8.3 km/s heliocentric correction.
- **On that scale, the paper's [O III] lands where Table 4 puts it.** The fit centre in Fig. 3c is 5007.246 Å
  (+24.1 km/s). In Fig. 3b, P1/P2/P3 are at +24.9/+28.6/+20.3 km/s, close to Table 4's +26.8/+26.8/+23.8.
- **Derived correction.** Table 4 [O III] minus (8.32 + 1.8) = **−10.1 ± 1.0 km/s**. P1, P2 and P3 become
  +16.7, +16.7 and +13.7 km/s (± 3–4). In V_LSR that is +4.5, +4.5 and +1.5 km/s.

**Best [O III] velocity:** V_helio = **+13.3 ± 1.0 km/s**, i.e. **V_LSR = +1.2 ± 1.0 km/s** (+1.5 with
Schönrich). This is the weighted mean of:
- DESI, my own stack: +14.4 ± 1.9 (including ±1.5 zero point)
- the paper's Fig. 3a against its own airglow: +13.6 ± 1.2
- SDSS, my own ON−OFF: +9.6 ± 2.8

χ² = 2.2 for 2 degrees of freedom. On top of this there is real structure between strands at the ±3–5 km/s
level: DESI gives E +13, W +18.5.

Once the [O III] is corrected, it moves with the WIM:
- [O III] − [N II] (paper red arm, +12.2) = +1.1 ± 1.0 km/s
- WHAM Hα is at V_helio ≈ +12

The kinematic argument for a separate [O III] component, and the kinematic distance, therefore disappear.

## 2. What I checked and how

| # | Check | Method (own code) | Script |
|---|---|---|---|
| A | LAMOST blue vs red arm zero point on the paper's plate | LAMOST DR12 v2.0 `med_catalogue` per-arm stellar RVs (`rv_b0`, `rv_r0`, `rv_b1`, `rv_r1`). Covers 76 stars on MN074209N160100N01, 38 plates in the Silk box and 26 plates from 2024-01-26 to 02-07. Gaia DR3 RV comparison. Header `HELIO_RV` vs astropy. | `scripts/verify_velocity_lamost_arms.py` |
| B | Forensics of possible sky-line reference mistakes | Rest wavelengths from NIST levels and Hanuschik's UVES atlas. 288 combinations of measured feature × reference wavelength × [O III] rest × frame handling. Simulated single-Gaussian fits to the [N I] doublet at R = 7500. Ren+2021's red-arm RV quadratics extrapolated to 5008 Å. | `scripts/verify_velocity_forensics.py` |
| C | **The paper's own spectra** | PyMuPDF extraction of Fig. 3 polylines from `sources/arxiv/2610.13329v1.pdf` p. 7, with the axes calibrated on the vector tick marks (rms 0.002 Å). Gaussian fits to [O III] 4959/5007, [N I] 5198/5200 and three weaker UVES-atlas sky lines. | `scripts/verify_velocity_fig3_vectors.py` |
| D | Sky zero point, solar (Fraunhofer) content of the continuum, [O III] refit | The continuum-normalised Fig. 3a spectrum is regressed on the Kurucz 2005 solar flux atlas (LSF-convolved) over a Doppler grid. [O III] is then refit with that solar component. | `scripts/verify_velocity_fig3_solar.py` |
| E | Hα offset from geocoronal subtraction | Forward model at the MRS-N LSF. WIM Hα at +12.2 km/s, FWHM 27 km/s. Geocoronal line at v_corr with its fine-structure offset. MRS-N template placed at RV(OH 6554), or literally at λ6554 × 6563/6554, with free intensity. Optional solar Hα absorption. Fitted with Gaussian + quadratic. | `scripts/verify_velocity_halpha_sim.py` |
| F | DESI and SDSS frames and velocities, recomputed | DESI: stacks of SPARCL residuals (flux − Redrock model), vacuum and barycentric. SDSS: per-group ON−OFF total sky (flux + sky) with my own [O I] 5577 zero point. Frame-sign test of [O I] 5577 against `HELIO_RV`. | `scripts/verify_velocity_desi_sdss.py` |
| G | Summary figure | | `scripts/verify_velocity_figure.py` → `figures/verify_velocity_summary.png` |

## 3. Results: claim → verdict → my value

| Claim (wave 1 / paper) | Verdict | My value / reason |
|---|---|---|
| Paper Table 4 [O III] +23.8…+26.8 is V_helio | **REFUTED** | It is a topocentric velocity on an air scale tied to the airglow lines (+1.8 km/s sky offset). The heliocentric correction was not applied. Correct by −10.1 ± 1.0 km/s, giving +16.7/+16.7/+13.7. |
| Paper [O III] − red-line offset of ~10–16 km/s is real | **REFUTED** | The calibration-free Fig. 3a value gives [O III] − [N II] = +1.4 ± 1.5 km/s. The offset is the missing −8.3 km/s v_corr plus +1.8 km/s of sky-reference offset. |
| Wave-1 hypothesis: blue-arm arc/secondary zero-point problem | **REVISED** | No arc-lamp problem. MN-plate stars give `rv_b0 − rv_r0` = +0.7 ± 0.5 km/s (N = 76). Other plates in the box give −1.6…+1.9; the 2024-01/02 plates +0.4…+0.9. The official zero points b1−b0 = −0.75 and r1−r0 = 0.0 leave the arms equal within 0.1 km/s. |
| Wave-1 hypothesis: heliocentric slip of 2\|v_corr\| = 16.2 km/s | **REVISED** | The slip is \|v_corr\| = 8.3 km/s (correction absent, not sign-flipped) plus 1.8 km/s of sky-reference offset. A sign flip would give ≈ +18, which the Fig. 3 data exclude. |
| LAMOST helio correction itself (stars) | CONFIRMED correct | MN-plate stars − Gaia DR3 (6 stars): `rv_b0` +1.0 ± 2.2, `rv_r0` +1.1 ± 0.8 km/s. A ±8 or ±16 km/s error is excluded. Header `HELIO_RV` = +7.3…+9.3 (median 8.33) = −v_corr (astropy −8.25…−8.38). |
| DESI [O III] +14.5 km/s (wave 1) | **CONFIRMED** | Own on-ridge stack, 390 fibres with predicted residual > 0.5 R: **+14.4 ± 1.2** (bootstrap). Residual > 1 R: +13.9 ± 1.2. South (Dec < 17.6): +13.1 ± 1.2. Negative-residual fibres, i.e. the petal sky mean: +12.6 ± 1.1, consistent. |
| SDSS sky fibres +12.5 ± 2.6 (wave 1) | CONFIRMED within errors | Own reimplementation: all strands +9.6 ± 2.8 (45 groups). South +7.5 ± 3.4. Reasonable selection variants give +6.8…+10.5, so my method sits about 3 km/s below wave 1. Both are far below +25. |
| DESI spectra vacuum and barycentric; sky residuals keep the nebula's velocity | CONFIRMED | desispec `extract.py` applies the barycentric factor to the wavelength grid before sky subtraction, so sky model, nebula and target share one frame. Positive- and negative-residual stacks agree (14.4 vs 12.6 ± 1.6). |
| SDSS spec-lite vacuum heliocentric, `HELIO_RV` = −v_corr | CONFIRMED | [O I] 5577 velocity against `HELIO_RV` over 288 groups: slope −1.08. Scatter is 12.4 km/s under the −v_corr convention and 42.0 under the opposite one. |
| Geocoronal Hα at V_helio ≈ −8.2 km/s on 2024-02-01 | CONFIRMED, refined | v_corr = −8.25…−8.38 km/s. Lyβ-pumped 3p→2s geocoronal emission sits −3.3 km/s from recombination-weighted Hα (NIST levels). Its effective position on the Hα velocity scale is therefore −11.6 km/s. |
| LAMOST Hα 4–7 km/s redder than WHAM and [N II]/[S II] | CONFIRMED (observed); explained as an artefact (PLAUSIBLE) | A misplaced template alone gives the **wrong sign** (−0.3…−1.4 km/s for both template positions). Over-subtraction does fit: subtracting 1/0.75 = 1.33× the true geocoronal intensity reproduces all four regions (χ² = 0.4 for 2 dof). It also matches the larger offsets in the Hα-faint regions P3/P4. Solar Hα absorption in the continuum adds +0.6…+2.7 km/s for C_sun/I = 0.1–0.3 Å⁻¹. |
| "5198 Å" line used for blue-arm wavelength calibration (task premise) | REVISED | The paper uses 5198 Å only for relative flux. Ren+2021 (arXiv:2008.06236) calibrates the **red arm only**, with seven OH/[O I] lines at 6287–6553 Å. Wu+2022 says MRS-N products are "mainly for the red band". Fig. 3a nevertheless shows the blue arm was re-referenced to the airglow in **air**, with [N I] 5198 at Hanuschik's 5197.928 Å (+1.3 km/s) and NIST 5197.902 (+2.8). |
| Rest wavelengths | measured / derived | [N I]: 5199.349/5201.705 Å vacuum, 5197.902/5200.257 Å air (NIST levels). Hanuschik's air values are 5197.928/5200.286 (+1.5 km/s) with a doublet ratio of 1.73. OH 6553.617 Å air = 6555.4275 Å vacuum (Ren). Hanuschik's 6553.626 Å is a blend (+0.4 km/s). |
| Sky-line mistakes that give +11…+14 | derived | With the resolved strong [N I] component, +8…+18 km/s arises only from frame errors: +8.3 (sky forced to rest in heliocentric spectra) or +16.6 (sign flip). Hanuschik's wavelength adds +1.5. Air/vacuum mismatches give −66…−86 or +68…+85 km/s, which is absurd. A doublet-ratio mismatch (ratio 1.4 vs 1.9) moves a moment centroid by ≤ 10 km/s, but simulated single-Gaussian fits lock onto one component (< 1 km/s) or onto the unresolved pair (+69 km/s). Extrapolating Ren's red-arm quadratics to 5008 Å would err by ±13 km/s rms (−15…+29 km/s for 2018-11 to 2019-01 plates). That is possible in principle, but the Fig. 3 data show it was not what happened. |

## 4. New findings

1. **The paper's Fig. 3 is machine-readable spectral data.** The vector polylines give the blue-arm spectrum
   (3,590 px, 4933–5372 Å), the P1–P4 spectra (603 px each, 4953–5020 Å) and the Gaussian-fit centres.
   Measured from it:
   - [N I] doublet ratio 1.71 (UVES: 1.73)
   - LSF σ = 0.271 Å, so R ≈ 8,100
   - [O III] 5007 σ_obs = 0.44 Å, i.e. an intrinsic **FWHM ≈ 50 km/s** in this example (paper 43–49; DESI ≤ 33)
   - [O III] 4959/5007 velocities on the plotted scale of +31.9/+24.1 km/s. The 4959 line is about 2σ redder.
   - The paper's vertical [O III] markers sit at the observed positions (4959.40, 5007.30 Å), not at rest
     wavelengths.
2. **There is a solar Fraunhofer blend at 5007.22 Å (air) with residual flux 0.20.** On 2024-02-01 it falls
   exactly on the topocentric [O III] position (5007.21 Å). Fe I 5005.71 and 5006.12 sit on the blue wing.
   - In Fig. 3a, 21 ± 2 % of the continuum carries solar absorption lines; Mg b is 13 % deep against 51–63 % for
     the convolved Sun. The continuum is 0.88 Å⁻¹ relative to the [O III] flux.
   - The resulting [O III] centroid bias is −0.4 km/s here and ≤ 2.5 km/s even for a purely solar continuum.
     It is minor for velocities but does remove a few percent of the [O III] flux in any analysis that ignores it.
   - [O III] 4959 is likewise flanked by Fe I 4957.3/4957.6 Å.
3. **Corrected paper velocities fit the DESI strand structure.** P2 corrected to +16.7 (+18.5 on the Fig. 3b
   scale) matches DESI's W strand (+17.5…+18.5). The DESI positive-residual fibres within 0.3° of P1 give
   +7.5 ± 2.2, against the corrected Table 4 P1 value of +16.7 ± 4.4. That is a 1.9σ tension, possibly real
   E/W structure.
4. **Ren+2021's red-arm sky-line RV calibration** carried −4 km/s plate offsets in 2018 (Sc lamp). On the paper's
   plate, the LAMOST stellar RVs agree with Gaia to about 1 km/s in both arms.

## 5. Limitations

- **Fig. 3a is one spectrum,** whose location and possible co-addition are unknown. The identification "Table 4
  was measured on the same scale" rests on three facts:
  - the plotted-scale [O III] velocities of Figs 3a/3b/3c (+20…+29) match Table 4
  - the plotted-scale airglow sits at rest
  - LAMOST products are vacuum and heliocentric

  Per-fibre MRS-N spectra are still not public (arXiv:2510.08950), so the correction is applied in aggregate
  (−10.1 ± 1.0 km/s). Per-fibre zero points (Ren+2021: ±0.3–1.5 km/s; spectrographs ±3 km/s) are not available.
- **The red arm cannot be inspected in the same way**, because the paper has no red-arm spectrum figure. Its
  correctness is inferred from the agreement of [N II]/[S II] with WHAM, not shown directly.
- **[N I] rest wavelengths:** NIST levels and Hanuschik's measurements differ by 1.5 km/s. I adopt the mean
  (±0.8 km/s).
- **SDSS:** my independent stack lies 2–5 km/s below wave 1, depending on group selection and zero point
  (bootstrap 2.5–3.7 km/s). The SDSS 4959 fits were unstable and are not used.
- **DESI absolute zero point:** DESI stars vs Gaia show +1…+3 km/s (wave-1 `rvcheck`), and two backup tiles show
  −19 km/s. I carry ±1.5 km/s. The DESI P2/P3 discs have too few fibres (P3: 3).
- **The Hα simulation is a model, not a measurement.** It assumes the MRS-N template description in Zhang+2021
  and Wu+2022: template at RV(OH 6554), σ = √(σ6554² + 0.02²), intensity from the solar-altitude relation (0.33
  × OH 6554 at sun altitude −51…−62°). How the true 2024 geocoronal intensity compares with that 2017–2020
  calibration is not known, and solar-cycle effects could have either sign. The red-arm continuum level is not
  known, so the solar-Hα term is a range.

## 6. Files

**Scripts:**
- `scripts/verify_velocity_lamost_arms.py`
- `scripts/verify_velocity_forensics.py`
- `scripts/verify_velocity_fig3_vectors.py`
- `scripts/verify_velocity_fig3_solar.py`
- `scripts/verify_velocity_halpha_sim.py`
- `scripts/verify_velocity_desi_sdss.py`
- `scripts/verify_velocity_figure.py`

**Derived data:**
- `data/derived/verify_velocity_lamost_arms.json`
- `data/derived/verify_velocity_forensics.json`
- `data/derived/verify_velocity_fig3_measure.json`
- `data/derived/verify_velocity_fig3_spectra.npz`
- `data/derived/verify_velocity_fig3_solar.json`
- `data/derived/verify_velocity_halpha_sim.json`
- `data/derived/verify_velocity_desi_sdss.json`
- `data/derived/verify_velocity_desi_sdss_stacks.npz`
- `data/derived/verify_velocity_summary.json`

**Raw downloads** (all on 2026-11-16, logged in each folder's `download_log.txt`):
- `data/raw/verify_velocity_literature/`: arXiv PDFs of Ren+2021 (2008.06236), Wu+2022 (2110.04547),
  Zhang+2021 (2108.08021), Wu+2021 (2007.05240) and Li+2022 (2205.14821)
- `data/raw/verify_velocity_lamost/`: three LAMOST DR12 v2.0 SQL queries of per-arm RVs
- `data/raw/verify_velocity_solar/`: Kurucz 2005 solar flux atlas, fetched with HTTP range requests at
  494–526 nm and 650–675 nm

**Figure:** `figures/verify_velocity_summary.png`
