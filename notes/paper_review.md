# Critical review of Wu et al. 2026 (arXiv:2610.13329): what the LAMOST MRS-N data do and do not establish

*Topic tag: `paper_review` · prepared 2026-11-16 · all numbers reproducible with*
`~/.venvs/astro313/bin/python -I scripts/paper_checks.py` → `data/derived/paper_checks.json`
*and* `scripts/paper_review_figures.py` → `figures/paper_review_*.png`.

Paper reviewed: Y.-F. Wu, S.-Y. Jia, C.-J. Wu, C. Zhang et al., "The Silk Nebula: A Newly Discovered Giant [O III]
Emission Filament", arXiv:2610.13329v1 (https://arxiv.org/abs/2610.13329, 2026-10-09), AJ doi:10.3847/1538-3881/aeb11c.
Local copies are `sources/arxiv/2610.13329v1.pdf`, `sources/arxiv/2610.13329v1.txt` and `sources/arxiv_2610.13329v1_html.md`.
Rendered pages are in `sources/arxiv/figs/pages/page-NN.png` (110 dpi). High-resolution crops are `sources/arxiv/figs/fig8_crop.png`,
`fig5b_crop.png` and `fig3c_crop.png`, and embedded images were extracted to `sources/arxiv/figs/images/`.
Section numbers below follow the PDF (3.3, 4.1, …); the HTML capture writes them as III.3, IV.1, ….

Conventions: **measured** = a value in the paper's tables; **reproduced** = recomputed here from the paper's inputs;
**derived** = new quantity computed here; **literature** = an external published value. Velocities are in km s⁻¹,
V_LSR is in the IAU standard LSR unless stated otherwise, and uncertainties are 1σ or 16–84 % unless stated otherwise.

---

## 1. Bottom line

1. **The headline distance is wrong by a factor of 10 because of a unit error.** Table 4's "Distance" columns are labelled
   "in units of 100 pc". They are the native **kpc** output of the `kd` package (Wenger et al. 2018) the paper cites.
   Emulating `kd`'s Reid+2014 rotation curve gives **1.785 kpc** for P1 [O III] (V_LSR = 15.5), against the paper's "1.79".
   At l ≈ 204°, b ≈ 19.5°, Galactic rotation gives only **V_LSR ≈ 1.5 km s⁻¹ at 160 pc** (2.8 at 300 pc) and
   **13–14 km s⁻¹ at 1.6 kpc**. Read literally, the [O III] velocities imply ~1.2–1.9 kpc, not 160 pc.
   Everything scaled from 160 pc (28 pc length, 30 pc radius, E = 2×10⁵¹ erg, age 6.8×10⁴ yr) inherits the error.
2. **No kinematic distance is meaningful here.** The rotational gradient is only 8.9–9.5 km s⁻¹ kpc⁻¹ along this sightline.
   With realistic 5–10 km s⁻¹ non-circular motions, the P1 [O III] Monte-Carlo distance spans 0.5–3.9 kpc (16–84 %).
   The velocity could equally be a 10–14 km s⁻¹ peculiar or expansion motion of a nearby shell.
   Even the corrected ~1.6 kpc is not a distance measurement.
3. **The "shock-excited" verdict rests on red-arm ratios that belong to the foreground/background warm ionized medium (WIM), not the filament.**
   [S II]/Hα, [N II]/Hα, the [S II] doublet, the Hα/[N II] widths, and the N2/N2S2Hα "abundances" are statistically
   identical in the nebula-free region P4 and on the filament (|Δ| ≤ 0.6σ for [S II]/Hα).
   They match faint high-latitude WIM ratios from WHAM ([N II]6583/Hα ≈ 0.5–0.8, [S II]6716/Hα ≳ 0.4; Madsen et al. 2006).
   WHAM-SS shows 2.9–4.4 R of Hα at all four positions with V_LSR = −0.5…+1.5 km s⁻¹ and no excess on the filament.
4. **The only filament-specific spectroscopic facts** are as follows. [O III]λ5007 is detected at S/N ≈ 6 in ~279 fibres along the
   imaged filament, at V_LSR ≈ +11…+16 km s⁻¹ (±3–4 statistical; blue- vs red-arm zero point unverified).
   It is offset by ~10 km s⁻¹ from the WIM lines. The [O III] width (43–49 km s⁻¹) is **not** measured.
   At S/N 6.1, with the paper's FWHM cut, any true width from 5 to 45 km s⁻¹ yields "intrinsic" widths of ~36–47 km s⁻¹.
5. **The [S II] densities (64–150 cm⁻³) are inconsistent with the paper's own Table 2 and with physics.**
   Table 2's ratios (1.03–1.19) give 190–440 cm⁻³ via the paper's own eq. 5. The nebula-free P4 has the *highest* implied density.
   Any n_e ≳ 60 cm⁻³ would need an [O III] path length of ~0.002 pc (filling factor ~10⁻³).
   Filling the observed strand width would make Hα thousands of R.
   The faint [O III] emission requires n_e ~ 1–3 cm⁻³.
6. **What the data *do* constrain (derived here).** The LAMOST Hα tracks WHAM to ±11 % (11.8 ± 1.3 R per relative unit).
   The Hα excess on the brightest [O III] region P3 is then 0.12 ± 0.33 R.
   If the paper's [O III] surface brightness (≈10 R) is right, [O III]5007/Hα > 13 (2σ, photon units) on the filament.
   Collisionally excited [O III] can do this only if **T_e ≳ 1.6×10⁴ K**, even with all oxygen as O²⁺.
   This rules out ordinary WIM/H II photoionization and favours either a shock seen before full recombination or a very hot ionizing source.
   This is the strongest "shock-like" argument available, and the paper does not make it.
7. **Verdict.** "Possible SNR" is still defensible as a hypothesis, but not for the paper's reasons. The paper's
   distance, size, energy, age, density, temperature, abundance, and BPT-based classification are not established by these data.
   The paper's own fallback ("DIG or ionized bubble cannot be excluded") is correct for the red lines.
   The [O III] filament itself is not DIG-like: WIM [O III]/Hα ≈ 0.06, against > 10 here.

---

## 2. Verdict table

| # | Paper claim (section) | Reproduced? | Assessment | Severity |
|---|---|---|---|---|
| 1 | Fiducial distance ~160 pc from [O III] kinematics (Sect. 3.3, Table 4, abstract) | **Yes, numerically.** The `kd`-Reid14 emulation gives 1.785 (paper 1.79) for P1, 1.74 (1.85) for P2 and 1.15 (1.20) for P3, in **kpc**. The MC mode is 1.70 [HPD 1.14–2.36] vs 1.63 (+1.01/−0.84) | Unit error: "100 pc" should be kpc. 160 pc → ~1.6 kpc. At 160 pc rotation predicts V_LSR = 1.5 km s⁻¹, not 11–16 | **Critical** |
| 2 | Kinematic distance usable just outside the anticentre exclusion zone (Sect. 3.3) | Reproduced gradient 8.9–9.5 km s⁻¹ kpc⁻¹ | 7 (10) km s⁻¹ of peculiar motion shifts d by 0.79 (1.12) kpc. P1: 0.76–3.37 kpc (σ_pec = 7). P3: 0.27–2.32 kpc, P(d < 0.35 kpc) = 0.18. No useful constraint | **Critical** |
| 3 | Size 1°×10° → 3 pc × 28 pc, shell radius ~30 pc (Sects. 3.3, 4.1) | Yes: 10° at 160 pc = 28.0 pc | Inherits #1. The traced spine is 16.9° long (geometry v1.0) and the main strand FWHM is 0.28°, not 1°. Radius 30 pc is asserted, not fitted; the geometry fit gives a 21.3° curvature radius (59.5 pc at 160 pc) | Major |
| 4 | V_s = 119 ± 5 km s⁻¹ (Dopita 1979) (Sect. 4.1) | Yes: 119.3 for n_e = 64, n0 = 1 | Error too small: n0 ± 30 % alone gives ±18. With n_e 64–150: V_s = 127–189 (16–84 %). The relation applies to the [S II] zone of a complete radiative shock, but the [S II] here is WIM | Major |
| 5 | E = (2.0 ± 0.7)×10⁵¹ erg, R_PDS = 18 ± 2.8 pc, t = (6.8 ± 1.8)×10⁴ yr (Sect. 4.1) | E = 2.01×10⁵¹ ✓. R_PDS = 17.1 pc (18 needs E = 2.4×10⁵¹). Age 6.75×10⁴ yr only with "+1/3" and R_PDS = 18; the printed formula has "+1" (→ 7.5–8.7×10⁴) | Arithmetic fine. Inputs (d, R_s, n_e) invalid. At 1.6 kpc, E = 2.9×10⁵⁴ erg (n0 = 1) to 3×10⁵³ (n0 = 0.01), which is impossible for one SN. At 288 pc, E = 1.3×10⁵² (n0 = 1) | **Critical** (conclusion) / Minor (formula typo) |
| 6 | V_LSR from V_helio with 20 km s⁻¹ toward 18h03m50.29s +30°00′16.8″ (Sect. 3.2.3) | Correction −11.9 to −12.6 (standard), −11.6 to −12.2 (Schönrich 2010; astropy agrees) | Table 4 V_LSR − V_helio varies line-to-line within a region (−10.0 to −15.2). It deviates from a single correction by −2.8…+2.6 (rms 1.2) | Minor |
| 7 | [S II]/Hα > 0.4 ⇒ shock-excited; [N II]/Hα 0.5–1.5 also SNR-like (Sects. 4.1, 5) | Ratios reproduced: [S II]/Hα 0.81–0.90 on-nebula vs 0.80 in P4 | Background = nebula (Δ ≤ 0.58σ). The values are those of faint WIM ([S II]6716/Hα 0.40–0.47, [N II]6583/Hα 0.70–0.80). The SNR criterion is not discriminating at these surface brightnesses | **Critical** |
| 8 | Sabin+13 / Kopsacheili+20 diagrams place points in the SNR zone; DIG overlap "minor" (Fig. 8) | Inspected Fig. 8 | Grey background points fully overlap the green/red nebula points. The orange "DIG" ellipse (Belfiore+22) is bright extragalactic DIG ([S II]/Hα ≈ 0.3–0.55), not faint Galactic WIM | Major |
| 9 | n_e ≈ 64 cm⁻³ median, ~150 cm⁻³ at centre, lower in background (Sect. 3.4) | Eq. 5 reproduced. Table 2 R gives n_e = 298, 247, 189, 400 cm⁻³ (P1–P4). Sanders+16: 328, 274, 215, 440 | 64–150 cm⁻³ needs R = 1.25–1.36. P4 (background) has the lowest R = highest n_e. 1σ ranges span ~60–1400 cm⁻³ | Major |
| 10 | Physical consistency of n_e with faint emission (not discussed) | Derived | EM ≥ 7.2 pc cm⁻⁶ for 10 R [O III] at 10⁴ K → L = 0.002 pc at n_e = 64. Filling the 0.28° strand (0.77 pc at 160 pc) needs n_e ≈ 3.3 cm⁻³ | Major |
| 11 | T_e ≈ 9000 K (text), 8.7×10³ (Fig. 7), "≈10,000 K" adopted; on-filament 1–2×10⁴ K vs P4 7.5×10³ K (Sect. 3.4) | Eq. 2 coefficient verified (23.55). Table 2 widths give 6780, 6580, 3220, 7670 K (P1–P4), σ_T = 7300–10 600 K, P(T < 0) = 17–38 % | T is unconstrained per region. Hα fine structure biases eq. 2 by +1.2–1.3×10³ K. The on/off-filament difference is not significant | Major |
| 12 | V_t = 29–33 km s⁻¹ shows shock turbulence (eq. 3, Table 4) | Eq. 3 coefficients verified. Table 2 widths give V_t = 26.7, 25.8, 29.9, 26.5 (P1–P4) | Table 4 values (33.1, 32.9, 30.6, 29.9) don't follow from eq. 3. P4 (background) is the same, so this is WIM turbulence | Minor |
| 13 | Intrinsic [O III] FWHM 43–49 km s⁻¹ (Table 2) | MC of the paper's fitting + cut at S/N 6.1 | Passing fits return median 35.8, 36.6, 40.0, 46.9 km s⁻¹ for true widths 5.4, 20, 30, 45 (pass fractions 10, 24, 48, 82 %). The width is not measured | Major |
| 14 | [O III] velocity errors 3–4 km s⁻¹ | Yes: MC σ_v = 3.1–4.3 at S/N 6.1 | Statistical part OK. Blue-arm vs red-arm zero point not demonstrated | Reproduced |
| 15 | Strong [O III]λλ4959,5007; S/N > 3 usability threshold (Sect. 3.2.1, Table 2) | Ratios 2.29 ± 1.33, 2.44 ± 1.28, 2.53 ± 1.17 (z = −0.4…−0.5 vs 2.98) | λ4959 S/N 1.8–2.3 is below the paper's own threshold. It is exactly what 5007 at S/N 6.1 predicts (2.05–2.2). MC ratio 1.7–4.9 (16–84 %). Uninformative, not anomalous | Minor |
| 16 | Near-solar, uniform O/H: N2 → 8.82, N2S2Hα → 8.70; "third-order polynomial" (Sect. 3.5) | Yes: PP04 linear 8.81–8.85, Dopita+16 8.64–8.75 | Eq. 6 is the *linear* PP04 N2 (the cubic gives 9.09–9.19). N2 = −0.10…−0.15 is outside PP04 (−2.5…−0.3) and Marino+13 (−1.6…−0.2) validity. Both indices measure WIM excitation, and P4 gives the same → no information on the filament | Major (interpretation) |
| 17 | [O III] SB 23.6–23.9 mag arcsec⁻² ≈ 10–13 R (Sect. 4.1) | Yes: 12.8, 9.9, 10.1 R (AB, 60 Å FWHM) | ±10 % for Vega vs AB, 8.3–14.9 R for a 50–70 Å effective width. Background subtraction and magnitude system are not stated. P1 "faintest yet brightest" shows ~0.3 mag stellar contamination | Moderate |
| 18 | [O III] vs red-line offset of 5–15 km s⁻¹ from two shell layers; DIG along the line of sight (Sects. 3.3, 4.2) | Offset ≈ 10–16 km s⁻¹ reproduced from Table 4 | Plausible, but the LAMOST Hα (+4…+6) is ~4–5 km s⁻¹ redder than WHAM Hα (−0.5…+1.5), while LAMOST [N II]/[S II] (−1…+1.7) match WHAM. Arm-to-arm and Hα (geocoronal/sky-subtraction) systematics need checking before interpreting the offset | Moderate |
| 19 | No NVSS/TGSS counterpart ⇒ consistent with old SNR (Sect. 4.1) | n/a | Interferometric surveys resolve out degree-scale emission (NVSS largest angular scale ~16′). The non-detection is uninformative | Minor |
| 20 | Detected by Fesen+2026 and "part of the Monogem Ring" while adopting 160 pc (Sect. 4.1) | n/a | Monogem/PSR B0656+14 is at ~288–300 pc, and the corrected kinematic value is ~1.6 kpc. The three are mutually inconsistent and not discussed | Moderate |

---

## 3. What the LAMOST MRS-N data establish, and what they do not

**Established (measured):**
- [O III]λ5007 is detected spectroscopically along the narrow-band filament (279 fibres, region S/N 6.1–6.6). That makes the
  imaging feature real line emission, not an optical artefact. There is no [O III] in P4.
- The [O III] heliocentric velocity is +23.8…+26.8 km s⁻¹ (V_LSR +11…+16), with a statistical precision of 3–4 km s⁻¹
  (confirmed by simulation). The fibre-to-fibre spread in Fig. 5 is ~0–30 km s⁻¹.
- Hα, [N II] and [S II] are present in *all* fibres at nearly equal intensity (relative Hα 0.25–0.37) and equal ratios on and
  off the filament. Their V_LSR (−1…+6) matches the diffuse WIM seen by WHAM.

**Not established:**
- Any distance (kinematic or otherwise). The adopted 160 pc is a units error, and the kinematic method is
  uninformative at |V_LSR| < 16 km s⁻¹ in this direction.
- Shock excitation *of the filament* from red-line ratios, the [S II] density, T_e, V_t, and O/H. All of these describe the WIM.
- The intrinsic [O III] width (S/N too low for R ≈ 7500 with a 40 km s⁻¹ instrumental profile).
- The [O III]/Hβ ratio (no Hβ coverage) or any absolute line ratio between arms (no flux calibration).
- SNR energetics and age.

---

## 4. Detailed checks

### 4.1 LSR conversion (checklist 1)
- The standard solar motion (20 km s⁻¹ toward 18h03m50.29s +30°00′16.8″ J2000) is (l, b) = (56.16°, 22.76°), i.e.
  (U, V, W) = (10.27, 15.32, 7.74) km s⁻¹. These are the IAU values hard-coded in `kd`.
- Correction V_LSR − V_helio: P1 −12.10, P2 −11.87, P3 −12.58, P4 −12.38 km s⁻¹ (standard). With Schönrich et al. 2010
  (11.1, 12.24, 7.25): −11.81, −11.65, −12.19, −12.08; astropy's `LSR` frame gives identical values.
  The two conventions differ by ≤ 0.4 km s⁻¹.
- Table 4 instead implies per-line corrections of −10.0 (P3 [S II]) to −15.2 (P4 Hα) within regions where the true correction
  varies by < 0.7 km s⁻¹. Paper − ours = −2.8…+2.6 km s⁻¹ (rms 1.2). V_LSR and V_helio were probably taken from different
  fibre subsets or statistics. This is a minor issue, but it means Table 4 velocities carry ~1–3 km s⁻¹ of unexplained scatter.
- Observation geometry: at transit at Xinglong (2024-02-01 15:15 UT, altitude 66°) the barycentric correction is
  −8.2 km s⁻¹. Geocoronal Hα therefore lies at V_helio ≈ −8.2 (V_LSR ≈ −20.3) km s⁻¹, only ~25 km s⁻¹ from the WIM Hα
  and within one instrumental FWHM (40 km s⁻¹). Imperfect geocoronal/sky subtraction could shift the LAMOST Hα centroid.

### 4.2 Kinematic distance (checklist 2)
Code: `paper_review_lib.MODELS` (full 3-D projection) and an exact emulation of `kd`'s `reid14_rotcurve.calc_vlsr`
(V = R0 sin l (Θ/R − Θ0/R0), which omits the cos b factor), read from https://github.com/tvwenger/kd (2026-11-16).

Expected V_LSR toward P1 (l = 203.77°, b = 19.49°):

| model | 160 pc | 288 pc | 300 pc | 1.0 kpc | 1.6 kpc | dV/dd (km s⁻¹ kpc⁻¹) |
|---|---|---|---|---|---|---|
| flat, R0 8.15, Θ0 236, std Sun | 1.50 | 2.66 | 2.77 | 8.66 | 13.14 | 8.89 |
| Reid+14 A5 linear (R0 8.34, Θ0 240, −0.2) | 1.59 | 2.75 | 2.86 | 8.77 | 13.29 | 8.91 |
| Reid+19 A5-like linear (R0 8.15, Θ0 236, −0.1)* | 0.08 | 1.25 | 1.35 | 7.26 | 11.77 | 8.92 |
| `kd` Reid+14 URC (emulated) | 1.59 | 2.83 | 2.95 | 9.26 | 14.10 | 9.49 |
| `kd` Reid+19 URC (kd nominal parameters) | 0.51 | 1.73 | 1.85 | 8.09 | 12.87 | 9.39 |

\*Reid+2019 linear-model parameters are quoted from memory (R0 = 8.15 ± 0.15, Θ0 = 236 ± 7, dΘ/dR = −0.1 ± 0.4,
U/V/W⊙ = 10.6/10.7/7.6); the URC parameters are taken from the `kd` source.

- **Units test.** For V_LSR = 15.5 / 14.5 / 11.2 km s⁻¹ (P1/P2/P3 [O III]), the traditional distances are 1.79 / 1.74 / 1.15 kpc
  (`kd`-Reid14 emulation) and 1.74–2.15 kpc across all models. Table 4 gives 1.79 / 1.85 / 1.20. A first-order Oort estimate
  (A sin 2l cos²b ≈ 9.5 km s⁻¹ kpc⁻¹) also gives d ≈ 1.6 kpc for 15 km s⁻¹. **Table 4's distances are in kpc.**
  The text's "~100–190 pc", "fiducial ~160 pc" (median of 1.63, 1.58, 0.98, 1.79, 1.85, 1.20 = 1.605) and the
  "100–150 pc" red/blue discrepancy are all low by a factor of 10.
- **The `kd` MC emulation** (resampled a1, a2, a3, R0 and V; no-solution samples dropped as `kd` does) gives P1 [O III]
  mode 1.70 kpc [68 % HPD 1.14–2.36], P2 1.70 [1.08–2.29] and P3 1.10 [0.74–1.53]. The paper's intervals are ~1.5× wider,
  perhaps from larger velocity errors. For the red lines with V_LSR ≈ −1…0, 35–68 % of samples have no solution and are
  silently discarded by `kd`, so the tabulated red-line distances (0.06–0.28) are truncation artefacts.
- **Realistic MC** (Reid+14 linear with parameter errors, Table 4 velocity errors, plus cloud peculiar velocity σ_pec;
  samples below V(d = 0) kept at d = 0):

| sample | σ_pec = 0 | 5 | 7 | 10 km s⁻¹ |
|---|---|---|---|---|
| P1 [O III] median (16–84 %) kpc | 1.92 (1.19–2.76) | 1.93 (0.94–3.11) | 1.92 (0.76–3.37) | 1.94 (0.47–3.86) |
| P1 P(d < 0.35 kpc) | 0.01 | 0.05 | 0.08 | 0.14 |
| P3 [O III] median (16–84 %) kpc | 1.22 (0.69–1.81) | 1.22 (0.44–2.10) | 1.21 (0.27–2.32) | 1.22 (0.02–2.73) |
| P3 P(d < 0.35 kpc) | 0.05 | 0.13 | 0.18 | 0.24 |

  With Reid+19 the medians are 0.2 kpc larger. **Judgement:** a 160 pc source needs V_rot ≈ 1.5 km s⁻¹. The measured +11…+16
  would then be a 10–14 km s⁻¹ peculiar motion, which is entirely ordinary for a shell or cloud. Conversely, a 1.6 kpc source needs no
  peculiar motion. The data cannot tell these apart, so the kinematic distance carries essentially no information.
  The ~10 km s⁻¹ [O III]–WIM offset may be partly instrumental (Sect. 4.1, item 18 of the table), which would only make this worse.

### 4.3 Size and height (checklist 3)
Sizes use b = 19.44°, z⊙ = 20.8 pc (Bennett & Bovy 2019) and the geometry v1.0 trace (`data/derived/silk_trace.json`).

| distance | 1° | 10° (paper) | 15° (press) | 17° (NSNS arc) | 16.9° spine | 0.28° strand FWHM | z above midplane |
|---|---|---|---|---|---|---|---|
| 160 pc (paper) | 2.8 pc | 28.0 pc | 42.1 pc | 47.8 pc | 47.3 pc | 0.77 pc | 74 pc |
| 288 pc (PSR B0656+14) | 5.0 | 50.4 | 75.8 | 86.1 | 85.1 | 1.39 | 117 |
| 300 pc | 5.2 | 52.5 | 79.0 | 89.7 | 88.7 | 1.45 | 121 |
| 1.6 kpc (`kd` kpc) | 27.9 | 280 | 421 | 478 | 473 | 7.7 | 553 |

The paper's 30 pc "remnant radius" corresponds to 10.6° at 160 pc and is not derived from any fit. The geometry
teammate's small-circle fit to the full spine gives a 21.3° radius (59.5 pc at 160 pc; 595 pc at 1.6 kpc).

### 4.4 Energetics (checklist 4)
- Dopita (1979), n_[S II] ≈ 45 n0 (V_s/100)²: n_e = 64, n0 = 1 → V_s = 119.3 km s⁻¹ (paper 119). n_e = 150 → 182.6 km s⁻¹.
- Cioffi+88 / Boumis+22: E = 6.8×10⁴³ n0^1.16 V_s^1.35 R_s^3.16 = 2.01×10⁵¹ erg (paper 2.0). R_PDS = 14.0 E51^(2/7) n0^(−3/7) = 17.1 pc
  (paper 18 ± 2.8). Age with the correct "+1/3": 7.96×10⁴ yr (R_PDS = 17.1) or 6.75×10⁴ yr (R_PDS = 18). The printed "+1" gives
  8.7×10⁴ / 7.5×10⁴ yr. As a cross-check, (2/7) R/V = 7.0×10⁴ yr.
- Errors: treating V_s and n0 as independent (paper-like) gives E = 1.3–2.7×10⁵¹, which reproduces ±0.7. Making V_s self-consistent with
  n0 and n_e = 64–150 gives V_s = 127–189 km s⁻¹ and E = 2.2–3.4×10⁵¹ erg.
- Alternative distances (R_s scaled from 30 pc at 160 pc; V_s from Dopita with n_e = 64):

| d | R_s | n0 = 1 | n0 = 0.3 | n0 = 0.1 | n0 = 0.01 |
|---|---|---|---|---|---|
| 160 pc | 30 pc | 2.0×10⁵¹ | 1.1×10⁵¹ | 6.6×10⁵⁰ (R < R_PDS) | 2.2×10⁵⁰ (V_s = 1190!) |
| 288 pc | 54 pc | 1.3×10⁵² | 7.2×10⁵¹ | 4.2×10⁵¹ | 1.4×10⁵¹ |
| 1.6 kpc | 300 pc | **2.9×10⁵⁴** | 1.6×10⁵⁴ | 9.5×10⁵³ | 3.1×10⁵³ |

  At 1.6 kpc the inferred energy is 300–3000 SNe, so a single SNR is unphysical. Lower n0 lowers E but forces V_s ≫ 300 km s⁻¹
  (no longer a radiative PDS shell). At 1.6 kpc the structure would sit at z ≈ 0.55 kpc, in the thick disk.

### 4.5 Line widths (checklist 5)
- Thermal FWHM at 10⁴ K: H 21.4, N 5.74, O 5.37, S 3.79 km s⁻¹. Instrumental FWHM c/R = 40.0 km s⁻¹ (R 7500;
  42.8 / 37.5 for R 7000 / 8000). The observed widths implied by Table 2 are Hα ≈ 51 and [O III] ≈ 59–63 km s⁻¹.
- Eq. 2 (Reynolds+77) is correct: T = 10⁴ (W_H² − W_N²) / (W_H,th²(1 − 1/14)) → coefficient 23.55. Eq. 3 is correct:
  √(14/13) = 1.038, 1/14 = 0.072. Note that V_t is a FWHM, not a dispersion.
- Applied to Table 2: T = 6780 (P1), 6580 (P2), 3220 (P3), 7670 K (P4) with σ_T = 8400, 7300, 10 600, 8000 K and P(T < 0) = 21, 19, 38, 17 %.
  Hα fine structure (two groups 6.4 km s⁻¹ apart, a two-component approximation) adds ~1.2–1.3×10³ K. Neither
  "1–2×10⁴ K on the filament" nor a difference from P4 is supported. Fig. 7a's spread (σ ≈ 6.5–6.8×10³ K) is consistent with noise.
- **[O III] width simulation** (Gaussian + linear fits, 6 km s⁻¹ pixels, ±240 km s⁻¹ window, W_inst = 40 ± 1.5, integrated-flux S/N
  set by the Fisher matrix; the paper's cut FWHM_obs − err > FWHM_inst + err applied):

| true intrinsic FWHM | S/N 6.1: pass fraction | median (16–84 %) of passing fits | S/N 15 | S/N 50 |
|---|---|---|---|---|
| 5.4 (thermal O) | 10 % | 35.8 (31.4–46.2) | 22.8 | 15.6 |
| 20 | 24 % | 36.6 (31.7–46.6) | 25.0 | 20.6 |
| 30 | 48 % | 40.0 (33.0–51.7) | 31.1 | 30.2 |
| 45 | 82 % | 46.9 (37.3–60.3) | 44.7 | 44.8 |

  The paper's 43–49 km s⁻¹ is what passing fits return for any true width ≳ 20–30 km s⁻¹. Even an unresolved line gives
  36 km s⁻¹ in the 10 % of fibres that pass. So "[O III] broader than Hα" is not demonstrated. The Hα/[N II] widths at S/N ~50 are measurable (±2–3 km s⁻¹),
  but they belong to the WIM. The WHAM Hα FWHM_obs of 31–34 km s⁻¹ (instrumental ≈ 12 km s⁻¹ → ~30–32 intrinsic) matches the
  LAMOST Hα 31–32 km s⁻¹.

### 4.6 [S II] density and physical self-consistency (checklist 6)
- The paper's eq. 5 (Proxauf+14) and Sanders+16 ([S II], 10⁴ K) give, for R = 1.03 / 1.10 / 1.14 / 1.19:
  400 / 298 / 247 / 189 and 440 / 328 / 274 / 215 cm⁻³. The T-correction √(T/10⁴) at 7500 K lowers these by 13 %. **n_e = 64 / 100 / 150
  requires R = 1.356 / 1.310 / 1.253.** The paper's 64 cm⁻³ is a median of per-fibre values dominated by fibres near the low-density
  limit (Fig. 7b: >50 % in the first bin). It does not represent Table 2.
- Table 2 intensities give R = 1.08 (P1), 1.13 (P2), 1.10 (P3; table says 1.19, possible only at the edge of rounding) and 1.00 (P4).
  The 16–84 % ranges are 0.87–1.33 / 0.95–1.36 / 0.84–1.45 / 0.75–1.33, i.e. n_e ≈ 60–1400 cm⁻³.
  **The nebula-free P4 has the highest density**, contradicting the text. Fig. 8 (middle) shows background fibres down to
  R ≈ 0.45–0.6, near the high-density limit (≥ 10³–10⁴ cm⁻³), which is impossible for diffuse gas. The per-fibre doublet ratio is
  therefore noise- or systematics-dominated.
- **Emission measure.** WHAM-SS DR1 total Hα: P1 3.17, P2 4.41, P3 2.86, P4 3.40 R (nearest 1° beams 0.27–0.48° away;
  Finkbeiner 2003 map: 3.40, 4.27, 2.74, 3.56 R). This gives EM = 2.75 T4^0.9 I_Hα = 7.9–12.1 pc cm⁻⁶ for the *whole sightline*.
  [O III] at 10–13 R and 10⁴ K needs EM ≥ 7.2–9.4 pc cm⁻⁶ even if all O is O²⁺ (Ω(³P–¹D) = 2.29, branching 0.75, O/H = 4.9×10⁻⁴).
  At 2×10⁴ K it needs ≥ 2.4–3.1 pc cm⁻⁶.
- **Path length.** For EM = 8.3 pc cm⁻⁶ (11.5 R [O III], 10⁴ K): L = EM/n_e² = 2.0×10⁻³ pc (n_e = 64), 3.7×10⁻⁴ (150),
  9.2×10⁻⁵ pc (300). Setting L = strand FWHM gives n_e = 3.3 (160 pc), 2.4 (300 pc) or 1.0 cm⁻³ (1.6 kpc). Gas at n_e = 64 filling
  the paper's 1° width at 160 pc would emit ~4200 R of Hα. A sheet seen edge-on (L > width) makes this worse.
  **n_e ~ 64–150 cm⁻³ is not credible for the [O III] gas.** That gas must have n_e ~ 1–3 cm⁻³, for which [S II] would sit at the
  low-density limit R ≈ 1.45. It is also not credible for the WIM (n_e ~ 0.1–1 cm⁻³).
- Pressure: n_e = 64 (150) at 10⁴ K means P/k = 1.3×10⁶ (3.0×10⁶) K cm⁻³, ~300–1000 × the typical local ISM thermal pressure
  (~10^3.5 K cm⁻³; Jenkins & Tripp 2011). This is self-consistent only inside a strong radiative shock with ram P/k = 2.4×10⁶ (n0 = 1, 119 km s⁻¹),
  which is circular here.

### 4.7 [O III] 5007/4959 (checklist 7)
Measured 1.67/0.73 = 2.29 ± 1.33 (P1), 2.44 ± 1.28 (P2) and 2.53 ± 1.17 (P3), errors from the quoted S/N. These sit 0.4–0.5σ below
the atomic value 2.98 (Storey & Zeippen 2000). The expected λ4959 S/N, given λ5007 S/N 6.1–6.6, is 2.05–2.21, matching the reported 1.8–2.3.
In the simulation, the ratio from independent fits is 1.67–4.94 (16–84 %; median 2.83), 41 % of trials give < 2.5, and λ4959 S/N < 3 in 83 %.
So λ4959 is below the paper's own usability threshold, carries no diagnostic weight, and "strong [O III] λλ4959,5007" (Sect. 4.1) is an overstatement.

### 4.8 Metallicity (checklist 8)
- N2 = log([N II]6583/Hα) = −0.097, −0.153, −0.114, −0.097 (P1–P4). PP04 linear (eq. 6, 8.90 + 0.57 N2): 8.85, 8.81, 8.84, 8.85
  (paper 8.82). The text calls it "third-order", but the PP04 cubic gives 9.19, 9.09, 9.16, 9.19. Marino+13 N2: 8.67–8.70.
  Dopita+16 N2S2Hα: 8.69, 8.64, 8.72, 8.75 (paper 8.70).
- N2 > −0.3 lies **outside the calibrated range** of PP04 (−2.5 < N2 < −0.3) and Marino+13 (−1.6 < N2 < −0.2). Both calibrations are
  built on H II regions photoionized by OB stars. In WIM/DIG and shocks, [N II]/Hα and [N II]/[S II] are driven by T_e and the ionization
  state: the WIM's high [N II]/Hα reflects a higher T_e (Haffner+99; Madsen+06), not a higher N/H. Applied there, the indices return
  "solar-ish" numbers for any Galactic gas.
- Because the P4 background gives the same values, the result describes the WIM along the sightline. "Near-solar, uniform O/H"
  says nothing about the filament and cannot exclude ejecta enrichment. The solar value quoted (8.66, attributed to Allende
  Prieto+2001) should be 8.69 ± 0.05 for that reference; 8.66 is Asplund+2004.

### 4.9 Classification logic vs WIM (checklist 9)
- In WHAM convention (energy units, single lines): [N II]6583/Hα = 0.80 ± 0.09, 0.70 ± 0.09, 0.77 ± 0.10 (P1–P3) vs **0.80 ± 0.10 (P4)**.
  [S II]6716/Hα = 0.47 ± 0.08, 0.46 ± 0.06, 0.42 ± 0.08 vs **0.40 ± 0.09**. Summed [S II]/Hα = 0.90, 0.86, 0.81 vs **0.80**.
  The differences from P4 are ≤ 0.58σ ([S II]) and −1.4…0σ ([N II]).
- Table 2's [N II]/Hα is ([N II]6548 + 6583)/Hα (P1: (0.09 + 0.24)/0.30 = 1.10 vs 1.11). The exception is P4, where 0.95 matches neither
  the sum (1.12) nor 6583 alone (0.80).
- Literature (verified in the PDF text, `data/raw/paper_review_refs/madsen2006.txt`): in the WIM "a typical value of [N II]/Hα is ≈ 0.5,
  but in some cases exceeds 1.0"; "[S II]/Hα is high, ≳ 0.4, in regions of faint Hα emission"; "faintest areas … (l > 150°) have the highest [N II]/Hα (≳ 0.7)";
  WIM [O III]/Hα ≈ 0.06. High-latitude filaments at l = 217–236° (Madsen+06 Table 5) have [N II]/Hα 0.27–0.63, [S II]/Hα 0.23–0.40 and
  [O III]/Hα ≤ 0.004–0.14. Haffner et al. (1999) show [N II]/Hα and [S II]/Hα rising as I_Hα falls, with [S II]/[N II] ~constant, which they read as rising T_e.
- **Implication.** At I_Hα ≈ 3 R the "[S II]/Hα > 0.4 ⇒ SNR" criterion (designed for discrete, bright nebulae) is satisfied by ordinary
  WIM. It classifies the nebula-free background as "shock-excited" too, so the red-line BPT/SMB diagnostics cannot support the SNR verdict.
  The [O III] filament is, however, clearly *not* WIM-like: its [O III]/Hα exceeds 10 (Sect. 4.11) versus ~0.06 for the WIM.

### 4.10 Surface brightness → Rayleigh (checklist 10)
AB: f_ν = 10^(−0.4(m+48.6)) → f_λ × 60 Å; 1 R at 5007 Å = 7.42×10⁻¹⁸ erg s⁻¹ cm⁻² arcsec⁻². P1 23.60 → 12.8 R, P2 23.87 → 9.9 R, P3 23.85 → 10.1 R,
consistent with the paper's "10–13 R". Systematics: ±0.1 mag (Vega/AB ambiguity; Gaia G is Vega, SDSS g ≈ AB) → ±10 %. Effective width
50–70 Å → 8.3–14.9 R. An f/0.95 beam broadens and blue-shifts a 6 nm interference filter, so the effective transmission at 5007 Å is uncertain.
Also unknown: whether Table 3 is background-subtracted (only one value per region is listed despite paired background regions) and how much
star light remains (P1 is visually faintest but photometrically brightest, by 0.27 mag).

### 4.11 New derived constraint: the filament's [O III]/Hα and temperature
- LAMOST relative Hα vs WHAM total Hα: 10.6, 11.9, 11.0, 13.6 R per unit (P1–P4) → 11.8 ± 1.3 R/unit. The red-arm lines trace the WIM seen by WHAM.
- Hα excess over P4: P1 0.59 ± 0.33 R, P2 1.41 ± 0.33 R (P2 lies toward the northern Hα structure), **P3 0.12 ± 0.33 R**
  (WHAM: P3 is 0.54 R *fainter* than P4). 2σ upper limits: 1.25, 2.08, **0.78 R**.
- With the Table 3 [O III] intensities: [O III]5007/Hα (photon units) > 10 (P1), > 4.8 (P2), **> 13 (P3)** at 2σ.
- Collisional-excitation ratio for x(O²⁺) = 1 and solar O/H: 1.7 (8000 K), 3.8 (10⁴), 11.8 (1.5×10⁴), 21 (2×10⁴), 41 (3×10⁴ K).
  **The P3 limit requires T_e ≳ 1.6×10⁴ K** (P1: ≳ 1.4×10⁴ K) even with every O atom in O²⁺.
- Interpretation: ordinary photoionized gas (WIM T_e 6000–9000 K; H II regions ~8000 K) cannot produce [O III]/Hα > 10. Gas hotter than
  ~1.5×10⁴ K with a high O²⁺ fraction is required. That occurs in the cooling zone of a ~100–200 km s⁻¹ shock seen before full
  recombination, or around a very hot ionizing source in very low-density gas. This argues more strongly for shock-like
  excitation of the filament than anything in the paper. It depends on (a) the paper's [O III] calibration (±30–50 %), (b) the Hα
  excess being diluted neither in the LAMOST fibres nor by the WIM, and (c) the [O III] and Hα sampling the same gas. Hβ, a deep Hα image
  of the strand, and [O III]λ4363 (T_e) would test it.

### 4.12 Internal inconsistencies (checklist 11)
1. Table 4 distance units (kpc labelled "100 pc") and everything downstream (160 pc, 28 pc, 30 pc, "100–150 pc" red/blue gap).
2. Extent: 1°×10° (paper) vs 1.5°×15° (press) vs ~17° arc (NSNS) vs 16.9° traced spine (geometry v1.0).
3. "Expansion speed ~140 km s⁻¹, Section 3.4" (Sect. 4.2). No such number appears in Sect. 3.4; the only speed is V_s = 119 km s⁻¹ (Sect. 4.1).
4. Temperature: "about 9000 K" (text), 8727 K (Fig. 7a legend), "≈10,000 K" (adopted). Fig. 7a legend μ = 6471, σ = 6798 K vs caption "σ ≃ 6.47×10³ K".
5. Distance 160 pc vs claimed association with the Monogem Ring (~288–300 pc) via Fesen+2026, without discussion.
6. P1 "appears the faintest" but has the highest surface brightness (Sect. 3.1).
7. λ4959 S/N 1.8–2.3 < the paper's S/N > 3 usability threshold, yet "strong [O III] λλ4959,5007" is invoked.
8. Shock threshold quoted as 0.4 (abstract, Sect. 5), 0.45 (Sect. 1) and 0.45–0.5 (Sect. 4.1).
9. Eq. 6 called "third-order polynomial" but linear.
10. n_e 64–150 vs Table 2 ratios (190–440); background P4 has the highest density, while the text says backgrounds are lower.
11. V_t in Table 4 vs eq. 3 applied to Table 2.
12. Age formula printed with "+1" (Cioffi inversion: "+1/3"). R_PDS 18 pc vs 17.1 pc computed from the paper's E.
13. Table 4 V_LSR − V_helio inconsistent between lines of one region (Sect. 4.1).
14. Fig. 5 (right): dashed [O III] marker at ≈ +9 km s⁻¹, whereas the histogram centroid is ≈ +14 and Table 4 gives +11…+16.
15. Fig. 4: x-axis "RA (degrees)" with hour ticks. Fig. 2 caption mentions a coloured [O III] overlay on a greyscale figure.
16. Solar O/H 8.66 attributed to Allende Prieto+2001 (8.69).
17. Radio: interferometric non-detection presented as supporting an old SNR.
18. Table 2 P3 [S II] ratio 1.19 vs intensities 0.11/0.10 (= 1.10). P4 [N II]/Hα 0.95 vs intensities (1.12 sum / 0.80 single).

---

## 5. Figures (all 200 dpi, `silk_style`)
| file | caption |
|---|---|
| `figures/paper_review_kinematics.png` | (a) V_LSR vs distance toward P1 for five rotation models with ±7 km s⁻¹ cloud motions, [O III] and red-line bands, and 160 pc / 288 pc / 1.6 kpc markers. (b) MC distance PDFs for P1 [O III] (σ_pec = 0, 5, 10). (c) Table 4 velocities by region vs WHAM Hα and rotation at 160 pc, 300 pc and 1.6 kpc |
| `figures/paper_review_density_em.png` | (a) [S II] ratio → n_e (paper eq. 5, Sanders+16) with Table 2 points and the paper's 64–150 cm⁻³ band. (b) Required path length vs n_e for the [O III] and WHAM emission measures, against the strand FWHM at three distances |
| `figures/paper_review_line_ratios.png` | (a) [S II]6716/Hα vs [N II]6583/Hα: LAMOST P1–P4 (P4 nebula-free) with Madsen+06 high-latitude filaments and O-star H II mean. (b) Hα budget (WHAM, LAMOST-scaled, excess over P4) vs [O III] surface brightness, with 2σ [O III]/Hα limits |
| `figures/paper_review_fit_sims.png` | (a) Intrinsic [O III] FWHM returned by fits passing the paper's cut at S/N 6.1 for true widths 5–45 km s⁻¹. (b) Simulated 5007/4959 ratio at S/N 6.1/~2 with P1–P3 values |
| `figures/paper_review_energetics.png` | (a) Cioffi/Boumis explosion energy vs distance for n0 = 1–0.01 (Dopita V_s) and the paper's inputs, with the single-SN band. (b) Spine length, 10° extent, strand FWHM and height above the midplane vs distance |
| `figures/paper_review_sky_context.png` | Finkbeiner (2003) Hα (RA/Dec) with the [O III] spine/envelope (geometry v1.0), Hα patch, LAMOST plate, P1–P4 and the WHAM beams used (labelled with I_Hα) |

Paper renderings for reference: `sources/arxiv/figs/pages/page-02.png` (Fig. 1), `page-05.png` (Fig. 2), `page-07.png` (Fig. 3, Table 2),
`page-09.png` (Fig. 4), `page-11.png` (Fig. 5), `page-12.png` (Fig. 6), `page-13.png` (Table 4), `page-14.png` (Fig. 7), `page-16.png` (Fig. 8).

---

## 6. Methods, sources and data (accessed 2026-11-16)
- **Paper**: arXiv:2610.13329v1, https://arxiv.org/abs/2610.13329 (PDF, text and HTML already in `sources/`). Pages rendered with
  `pdftoppm -r 110` (all pages) and `-r 300` (pp. 7, 11, 14, 16). Images extracted with `pdfimages -png`.
- **kd source** (Wenger et al. 2018, ApJ 856, 52): https://github.com/tvwenger/kd, files `kd/reid14_rotcurve.py`, `kd/reid19_rotcurve.py`,
  `kd/rotcurve_kd.py`, `kd/pdf_kd.py` and `kd/kd_utils.py` fetched from https://raw.githubusercontent.com/tvwenger/kd/master/kd/ as text
  into `data/raw/github_kd/*.txt`. They were read only and re-implemented in `scripts/paper_review_lib.py`; nothing was executed, and the pickle file was not fetched.
- **Madsen, Reynolds & Haffner 2006**, ApJ 652, 401: https://arxiv.org/pdf/astro-ph/0609558 → `data/raw/paper_review_refs/madsen2006_astro-ph0609558.pdf`
  and `madsen2006.txt` (quotes and Table 5 verified there).
- **Haffner, Reynolds & Tufte 1999**, ApJ 523, 223: abstract at https://arxiv.org/abs/astro-ph/9904143 (trend statements only).
- **WHAM-SS DR1** (Haffner et al. 2003, ApJS 149, 405; release tag DR1, version 161116, file date 2017-09-12):
  `data/raw/wham_ss_dr1/wham-ss-DR1-v161116-170912.fits`, downloaded by a teammate (see their note for the URL). The nearest pointings were
  fitted with a single Gaussian over V_LSR −60…+70 km s⁻¹.
- **Finkbeiner 2003** full-sky Hα (ApJS 146, 407; HEALPix N_side 1024, 6′ FWHM, v1.1): `data/raw/lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits`,
  downloaded by a teammate from LAMBDA.
- **Geometry**: `data/derived/silk_trace.json` v1.0 (geometry teammate). It supplies the spine (16.9° arc), strand HWHM 0.138°, envelope, Hα patch,
  LAMOST plate and the curvature fit.
- **Formulae/literature used**: Dopita 1979 (ApJS 40, 455); Cioffi, McKee & Bertschinger 1988 (ApJ 334, 252); Boumis et al. 2022 (MNRAS 512, 1658);
  Proxauf et al. 2014 (A&A 561, A10); Sanders et al. 2016 (ApJ 816, 23); Pettini & Pagel 2004 (MNRAS 348, L59); Marino et al. 2013 (A&A 559, A114);
  Dopita et al. 2016 (Ap&SS 361); Reid et al. 2014 (ApJ 783, 130) and 2019 (ApJ 885, 131); Schönrich, Binney & Dehnen 2010 (MNRAS 403, 1829);
  Bennett & Bovy 2019 (MNRAS 482, 1417); Storey & Zeippen 2000 (MNRAS 312, 813); Osterbrock & Ferland 2006 (Ω(³P,¹D) of O²⁺ = 2.29);
  Haffner et al. 2003 (EM = 2.75 T4^0.9 I_Hα); Jenkins & Tripp 2011 (ApJ 734, 65; ISM pressure).
- **Scripts**: `scripts/paper_review_lib.py` (physics), `scripts/paper_checks.py` (all checks → `data/derived/paper_checks.json`; MC samples
  → `data/derived/paper_review_mc_samples.npz`; ~50 s), `scripts/paper_review_figures.py` (figures). The RNG is seeded (20261116).

## 7. Limitations
- Paper values were transcribed from the PDF tables. Per-fibre MRS-N measurements are not public through any route I could reach,
  so fits could not be redone on real spectra. The line-fit MC uses white noise, a flat continuum, a 6 km s⁻¹ pixel and
  *integrated-flux* S/N. If the paper's S/N is a peak S/N, the real S/N is higher and the width bias smaller (see the S/N 15 column).
- The Reid+2019 linear-model parameters are quoted from memory. The `kd` URC models use the source code values. The `kd` MC uses
  KDE-correlated Reid+19 parameter draws (pickle) that I did not reproduce; independent Gaussian errors were used instead.
- [O III] emissivity uses Ω = 2.29, branching 0.75 and solar O/H. Atomic-data and abundance uncertainties are ~20–30 %, which shifts T_min by
  ~10 %. WIM n_e and T_e are literature values.
- The WHAM beam (1°) is ~3.6× the strand FWHM, so WHAM limits on a strand-confined Hα excess are diluted. The LAMOST-based excess avoids
  this but assumes a common intensity scale for fibres in all four regions (supported to ±11 % by the WHAM comparison).
- The geocoronal-Hα explanation of the LAMOST Hα velocity offset is a hypothesis; it was not tested on spectra.

## 8. Data gaps
- LAMOST MRS-N per-fibre nebular spectra and line catalogue for plate MN074209 (2024-02-01). `dr10.lamost.org` is not reachable from this
  environment. A teammate has `data/raw/lamost_dr12/` stellar MRS catalogues; nebular products were not checked here.
- No Hβ or [O III]λ4363 for this field: no [O III]/Hβ, BPT or direct T_e. No absolute flux calibration in MRS-N.
- Fesen et al. 2026 (arXiv:2608.26666) and the Monogem association were not reviewed here (context teammate).
- No independent check of the SWONS [O III] photometric calibration (Table 3). All [O III]-in-Rayleigh statements inherit it.
