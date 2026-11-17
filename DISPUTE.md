# Comment on "The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament" (Wu et al. 2026, AJ 172, 298; arXiv:2610.13329)


*Draft for discussion, 2026-11-17. Prepared from the Silk Nebula study in this folder (`REPORT.md`). Every number, figure and table can be reproduced from the scripts listed in Section 9.*

**Conventions.**
- "Table n", "Fig. n", "Sect. n" and "Eq. n" refer to Wu et al. (2026). Page numbers refer to arXiv:2610.13329v1. The AJ version has the same tables and footnotes, and Crossref lists no erratum as of 2026-11-17.
- "Figure Cn", "Table Cn" and "Section n" refer to this Comment.
- Velocities are heliocentric (V_helio) unless marked V_LSR. Surface brightness is given in rayleighs (R).
- Every velocity discussed in Sections 2–4 is tabulated side by side, with its reference frame and the Earth-motion and LSR corrections, in the companion [VELOCITY.md](VELOCITY.md) ([VELOCITY.html](VELOCITY.html); 中文 [VELOCITY_zh.md](VELOCITY_zh.md)).
- The Summary is given in English and Chinese (概要). In the HTML version, the EN | 中文 switch at the top right selects the language of the Summary; the rest of the document is in English.

---

<!-- lang:en -->
## Summary

Wu et al. (2026) present LAMOST MRS-N spectroscopy of the Silk Nebula, a 17°-long [O III] filament in Gemini. From it they conclude that the filament is probably an old supernova remnant at ~160 pc, ~28 pc long, with E ≈ 2×10⁵¹ erg and an age of ~6.8×10⁴ yr. The [O III] detection is real and is confirmed by DESI, SDSS and narrow-band imaging. The paper is also right that the red-arm lines trace foreground/background diffuse ionized gas.

The quantitative chain behind the remnant interpretation does not survive inspection. Table C1 lists the disputed claims, and Sections 1–6 give the evidence:
- **Section 1.** The distance is a unit misreading of Table 4. It is shown from the paper's own numbers, not by appeal to software defaults.
- **Section 2.** No kinematic distance is possible in this direction.
- **Section 3.** The [O III] velocity offset is a wavelength-reference artefact.
- **Section 4.** The Hα velocities are red-biased.
- **Section 5.** Table 3 overestimates the surface brightness.
- **Section 6.** The remaining diagnostics either describe the warm ionized medium or are not constrained by the data.

None of this argues against shock excitation of the filament. Independent spectra favour an incompletely cooled radiative shock (`REPORT.md` §5). The point is that the paper's evidence does not establish its distance, size, energy, age or kinematics.

**Table C1.** Disputed claims of Wu et al. (2026), with our findings.

| # | Claim (location) | Finding | Key evidence | Effect | Section |
|---|---|---|---|---|---|
| 1 | Fiducial distance ~160 pc; Table 4 distances "in units of 100pc" (Table 4 footnotes; Sect. 3.3; Abstract; Sects. 4.1, 5) | **Wrong unit.** Table 4 is in kpc, so the paper's own fiducial value is ~1.6 kpc | The printed values reproduce as kpc (1.785 / 1.743–1.803 / 1.154 vs 1.79 / 1.85 / 1.20). At 160 pc, rotation gives only 1.5–1.7 km s⁻¹, not +15.5. Read as ×100 pc, the Monte-Carlo widths imply gradients of 37–47 km s⁻¹ kpc⁻¹ (rotation allows ≤ 10.7). The red-line lower bounds sit at the 0.01-kpc grid floor | All sizes ×10. Energy and age invalid | 1 |
| 2 | A kinematic distance is usable near the anticentre boundary (Sect. 3.3) | **No.** Kinematics carry no distance information here | The rotation gradient is 7.6–9.1 km s⁻¹ kpc⁻¹ toward P1–P4 and 4.3–9.8 along the spine, 8.1° of which lies inside the 160° < ℓ < 200° anticentre band. Peculiar motions of 5–10 km s⁻¹ therefore mean Δd = 0.55–1.32 kpc. Rotation contributes ≤ 3.1 km s⁻¹ within 300 pc. The local clouds LIC and Gem (< 15 pc) would get 0.74 and 2.70 kpc | The distance is unconstrained by the paper's data | 2 |
| 3 | [O III] V_helio +23.8…+26.8 km s⁻¹, 10–16 km s⁻¹ off the red lines; "two-layer shell" (Table 4; Sects. 3.2.3, 4.2) | **Artefact (−10.1 ± 1.0 km s⁻¹).** The blue-arm scale is air and topocentric | In the paper's Fig. 3a, the [N I] 5198/5200 airglow sits at its air rest wavelengths. [O III] measured against that airglow gives +13.6 ± 1.2. DESI, SDSS and Fig. 3a combine to +13.3 ± 1.0, equal to [N II], [S II] and WHAM | The offset, the two-layer interpretation and the [O III] kinematic distances all go | 3 |
| 4 | Hα V_helio +16.3…+19.7 km s⁻¹; Hα used for velocities, maps and distances (Table 4; Fig. 5; Sect. 3.2.3) | **Red-biased by +4.5…+6.7 km s⁻¹** | Hα − [N II] = +5.5 ± 1.1 km s⁻¹ (4.9σ) and Hα − WHAM = +4.8 ± 0.7, while [N II] − WHAM = −0.8 ± 1.0. Geocoronal over-subtraction is plausible: a joint model explains one-half to three-quarters of the offset | The Hα velocities and Hα-based distances are unusable | 4 |
| 5 | 23.60–23.87 mag arcsec⁻², i.e. ~10–13 R, "well above" the Fesen et al. (2024) limits (Table 3; Sect. 4.1) | **Overestimated 4.6–8.6×** | DESI gives 1.62 R; SDSS fibres 1.49 / 1.78 / 1.87 R; SDSS imaging 2.31 R; NSNS 2.0–2.2 R. No multiplicative systematic accounts for it; an additive background error of 7.9–11.3 R would | The brightness argument fails. The GALEX FUV/[O III] ratio contradicts Table 3 | 5 |
| 6 | BPT/SMB and [S II]/Hα shock classification; n_e 64–150 cm⁻³; T_e, V_t, O/H; E, R_PDS, age; [O III] FWHM 43–49 km s⁻¹; radio non-detection (Sects. 3.4–4.1) | **Describe the WIM, or not supported by the data** | The ratios are the same in nebula-free P4. Table 2 gives n_e 225–631 cm⁻³ (169–400 with Eq. 5). T_e has σ_T 7,300–10,600 K. The age formula has "+1" where it should read "+1/3", and R_PDS is 17.1 pc. The widths are not recoverable at S/N 6. Interferometers resolve out a 17′ filament | The SNR classification is not supported by these diagnostics | 6 |
<!-- /lang -->

<!-- lang:zh -->
## 概要

Wu et al. (2026) 发表了丝绸星云（Silk Nebula）的 LAMOST MRS-N 光谱观测。丝绸星云是双子座中一条长 17° 的 [O III] 丝状结构（filament）。作者据此得出结论：该丝状结构很可能是一个年老的超新星遗迹（supernova remnant，SNR），距离 ~160 pc，长 ~28 pc，E ≈ 2×10⁵¹ erg，年龄 ~6.8×10⁴ yr。[O III] 探测是真实的，并已得到 DESI、SDSS 和窄带成像的证实。论文认为红端（red arm）谱线示踪的是前景/背景弥漫电离气体（diffuse ionized gas），这一点也是正确的。

支撑这一遗迹解释的定量论证链经不起检验。表 C1 列出了受质疑的论断，本文第 1–6 节给出证据：
- **本文第 1 节：** 该距离源于对原文表 4 单位的误读。这一点依据论文自身的数值即可证明，而非诉诸软件的默认设置。
- **本文第 2 节：** 在这一方向上无法得到运动学距离（kinematic distance）。
- **本文第 3 节：** [O III] 速度偏移是波长参考假象（wavelength-reference artefact）。
- **本文第 4 节：** Hα 速度存在偏红偏差。
- **本文第 5 节：** 原文表 3 高估了面亮度。
- **本文第 6 节：** 其余诊断要么描述的是温电离介质（warm ionized medium，WIM），要么无法由数据约束。

以上各点均不构成反对该丝状结构受激波激发的论据。独立光谱（`REPORT.md` §5）倾向于支持未完全冷却的辐射激波（incompletely cooled radiative shock）。问题在于，论文的证据不足以确立其距离、尺寸、能量、年龄或运动学。

**表 C1.** Wu et al. (2026) 中受质疑的论断及我们的核查结果。

| # | 论断（出处） | 核查结果 | 关键证据 | 影响 | 本文章节 |
|---|---|---|---|---|---|
| 1 | 基准距离（fiducial distance）~160 pc；原文表 4 的距离“in units of 100pc”（原文表 4 脚注；原文第 3.3 节；原文摘要；原文第 4.1、5 节） | **单位错误。** 原文表 4 的单位是 kpc，因此论文自身的基准值实为 ~1.6 kpc | 按 kpc 解读，印出值可以复现（1.785 / 1.743–1.803 / 1.154，对比 1.79 / 1.85 / 1.20）。在 160 pc 处，银河系自转仅给出 1.5–1.7 km s⁻¹，而非 +15.5。若按 ×100 pc 解读，蒙特卡罗（Monte-Carlo）区间宽度所对应的隐含梯度（implied gradient）为 37–47 km s⁻¹ kpc⁻¹（自转仅允许 ≤ 10.7）。红端谱线的距离下限落在 0.01 kpc 的网格下限（grid floor）上 | 所有尺寸均 ×10。能量与年龄无效 | 1 |
| 2 | 在反银心方向（anticentre）的边界附近，运动学距离可用（原文第 3.3 节） | **不可用。** 在此方向上，运动学不携带距离信息 | 朝 P1–P4 方向的自转速度梯度为 7.6–9.1 km s⁻¹ kpc⁻¹，沿主脊线（spine）为 4.3–9.8，其中有 8.1° 的主脊线位于 160° < ℓ < 200° 的反银心带内。因此，5–10 km s⁻¹ 的本动速度（peculiar motion）意味着 Δd = 0.55–1.32 kpc。在 300 pc 以内，自转的贡献 ≤ 3.1 km s⁻¹。本地星际云 LIC 和 Gem（< 15 pc）会被分别赋予 0.74 和 2.70 kpc | 论文的数据无法约束该距离 | 2 |
| 3 | [O III] V_helio 为 +23.8…+26.8 km s⁻¹，与红端谱线相差 10–16 km s⁻¹；“two-layer shell”（原文表 4；原文第 3.2.3、4.2 节） | **假象（−10.1 ± 1.0 km s⁻¹）。** 蓝端（blue arm）波长标度为空气波长（air wavelength），且为站心（topocentric）参考系 | 在原文图 3a 中，[N I] 5198/5200 气辉（airglow）位于其空气静止波长（rest wavelength）处。以该气辉为参考测量，[O III] 为 +13.6 ± 1.2。DESI、SDSS 与原文图 3a 合并得到 +13.3 ± 1.0，与 [N II]、[S II] 及 WHAM 相同 | 该速度偏移、双层壳层解释以及基于 [O III] 的运动学距离均不再成立 | 3 |
| 4 | Hα V_helio 为 +16.3…+19.7 km s⁻¹；Hα 被用于速度测量、分布图和距离估计（原文表 4；原文图 5；原文第 3.2.3 节） | **存在 +4.5…+6.7 km s⁻¹ 的偏红偏差** | Hα − [N II] = +5.5 ± 1.1 km s⁻¹（4.9σ），Hα − WHAM = +4.8 ± 0.7，而 [N II] − WHAM = −0.8 ± 1.0。地冕 Hα 的过度扣除（geocoronal over-subtraction）是一个合理的解释：联合模型可解释该偏移的一半至四分之三 | Hα 速度以及基于 Hα 的距离均不可用 | 4 |
| 5 | 23.60–23.87 mag arcsec⁻²，即 ~10–13 R，“well above”（远高于）Fesen et al. (2024) 的极限（原文表 3；原文第 4.1 节） | **高估 4.6–8.6×** | DESI 给出 1.62 R；SDSS 光纤（fibre）1.49 / 1.78 / 1.87 R；SDSS 成像 2.31 R；NSNS 2.0–2.2 R。任何乘性系统误差（multiplicative systematic）都无法解释这一差异；7.9–11.3 R 的加性背景误差（additive background error）则可以 | 亮度论据不成立。GALEX FUV/[O III] 比值与原文表 3 相矛盾 | 5 |
| 6 | 基于 BPT/SMB 与 [S II]/Hα 的激波分类；n_e 64–150 cm⁻³；T_e、V_t、O/H；E、R_PDS、年龄；[O III] FWHM 43–49 km s⁻¹；未探测到射电辐射（原文第 3.4–4.1 节） | **描述的是 WIM，或未得到数据支持** | 在无星云的 P4 中，这些谱线比与之相同。由原文表 2 得出 n_e 为 225–631 cm⁻³（若采用原文式 (5) 则为 169–400）。T_e 的不确定度 σ_T 为 7,300–10,600 K。年龄公式中写作“+1”之处应为“+1/3”，且 R_PDS 为 17.1 pc。在 S/N 6 下无法恢复谱线宽度。干涉阵列会将宽 17′ 的丝状结构滤除（resolve out） | 这些诊断并不支持将其归类为 SNR | 6 |
<!-- /lang -->

---

## 1. The 160 pc distance: Table 4 is in kpc, not "100 pc"

**Claim in the paper.**
- Kinematic distances come from `kd` (Wenger et al. 2018), with "the rotation curve ... upgraded by M. Reid et al. (2014)" (Sect. 3.3, p. 10).
- Table 4 (p. 13) gives Monte-Carlo and "traditional" distances "in units of 100pc" (footnotes 1, 2).
- The [O III] rows read 1.63/1.58/0.98 (Monte Carlo) and 1.79/1.85/1.20 (traditional). The text calls these "∼100–190 pc" (p. 10).
- The paper adopts "a fiducial distance of ∼160 pc (the approximate median of the six [O iii]-based kinematic distance estimates; see Table 4)" (Sect. 3.3, p. 11). This value is used again in the Abstract and in Sects. 4.1 and 5.

Our tests use only the paper's Tables 1 and 4, Fig. 5b and text, plus the rotation curve it names; `kd`'s documented unit is corroboration only. Pages refer to arXiv:2610.13329v1. The AJ footnotes are identical, and Crossref lists no erratum (2026-11-17). 

**Why it is wrong.**

1. **Read as kpc, the printed values are reproduced to within 0.005–0.05 kpc.**
   - **Method.** We inverted the printed V_LSR along the Table 1 sightlines with the `kd`-form Reid et al. 2014 curve (Persic et al. 1996 functional form). V_helio-based values use the paper's solar motion (Sect. 3.2.3). `kd` itself was not run: its curve and `pdf_kd` were re-implemented from the source.
   - **[O III] traditional distances** (V_LSR → ours vs Table 4):
     - P1: 15.5 km s⁻¹ → **1.785** vs 1.79;
     - P2: 14.5 km s⁻¹ → 1.743 (**1.803** from V_helio) vs 1.85;
     - P3: 11.2 km s⁻¹ → **1.154** vs 1.20.
   - **[O III] Monte Carlo.** A `kd` emulation gives modes of 1.66/1.67/1.16 kpc. Table 4 prints 1.63/1.58/0.98.
   - **Velocity residuals.** We compared V_LSR with the rotation velocity at the printed distance (χ², kpc vs ×100 pc):
     - [O III]: 0.05 vs **28.9** (3 dof, p = 2×10⁻⁶);
     - Hα: 6.3 (p = 0.18) vs 42.3 (p = 1×10⁻⁸);
     - all 15 rows: 15.9 (p = 0.39) vs 73.0 (p = 1×10⁻⁹);
     - [S II]/[N II] (V_LSR = 0 ± 2 km s⁻¹): 9.5 vs 1.7, which does not discriminate.
   - **The unit result does not depend on the curve.** All seven variants, including `kd`'s default (Reid et al. 2019 URC), give 1.74–2.15 kpc for [O III] P1–P2 and 1.15–1.43 kpc for P3, at least 9× the ×100 pc values (0.12–0.19 kpc).
   - **Only P1 pins the curve.** P1 is matched to 0.005 kpc only by the `kd`-form Reid et al. 2014 curve (others 1.91–2.15). For P2 and P3 that curve is within 0.05–0.11 kpc and its ×cos b variant within 0.02 kpc.

2. **At 160 pc, rotation cannot produce the paper's velocities.**
   - **Rotation curve.** At 160 pc, rotation gives V_LSR = **1.5–1.7 km s⁻¹** toward P1–P3 (0.0–1.7 over all seven variants). The Oort term A sin 2l cos²b × d agrees, with A = 15.3 km s⁻¹ kpc⁻¹ (Bovy 2017).
   - **Distance needed.** Within 300 pc, rotation stays ≤ 3.1 km s⁻¹ (Figure C1a, inset). Reaching 15.5 km s⁻¹ requires **1.79–2.15 kpc**.
   - **The paper's own Fig. 5b** (our Figure C1b) shows this for the whole [O III] sample:
     - 97.5 % of the fibres have V_LSR > +2.3 km s⁻¹, i.e. beyond 0.23 kpc on the P1 curve;
     - 88 % have V_LSR > +6.6 km s⁻¹, i.e. beyond 0.69 kpc.
     - Bin heights are multiples of 1/81, so N is 81 or a multiple. This is our inference; N is not printed, and the fractions do not depend on it.

3. **The Monte-Carlo widths encode the velocity gradient.**
   - **The bound.** A 68 % velocity width W_V maps to W_d ≈ W_V/(dV/dd). Curve resampling, KDE smoothing and extra error terms only widen W_d, so W_V/W_d cannot exceed the largest true dV/dd, 5.8–10.7 km s⁻¹ kpc⁻¹ along P1–P4 at 0.05–2.5 kpc (seven variants plus the Oort term). This assumes the tabulated σ_v (3–4 km s⁻¹ for [O III], Sect. 3.2.3) or larger was used.
   - **Example, P1 [O III].** 8.6/(1.01 + 0.84) = **4.65** km s⁻¹ per unit: allowed as kpc, but **46.5 km s⁻¹ kpc⁻¹** as ×100 pc.
   - **All rows** (Figure C1d, Table C2). As kpc, 0.9–4.7: all 15 allowed. As ×100 pc, 9.5–46.5: 13 of 15 exceed 10.7, [O III] by **3.5–4.4×** and Hα by only 1.3–1.7×. The truncated P1/P2 [N II] rows (9.9, 9.5) do not discriminate.
   - **Emulation check.** Against our `pdf_kd` emulation with the tabulated σ_v, the printed [O III] intervals are 1.5–2.3× wider if kpc (as extra error terms would make them) and 4.4–6.6× narrower if ×100 pc.
   - **Limits.** This is a bound that supports step 1 but is not needed for it. It fails only if the Monte Carlo used velocity errors several times smaller than Table 4 lists.

4. **The Monte-Carlo lower bounds sit at `kd`'s grid floor.**
   - **Printed values.** In all 12 red-line rows, mode − error = **0.00–0.02**.
   - **Read as kpc,** this is `pdf_kd`'s `dist_min = 0.01` kpc. Our emulation puts the lower edge of all eight [S II]/[N II] rows at 0.010–0.017 kpc, the d → 0 edge expected for V ≈ 0 ± 2 km s⁻¹ once negative-velocity samples (no solution) are dropped.
   - **Read as ×100 pc,** the bounds would be 0–2 pc, below that floor. The sightline DIG (Sect. 4.2) would also sit at 6–80 pc, which is unexpected for a layer with a 1.83-kpc scale height (Gaensler et al. 2008). That last point is plausibility only.

5. **"Fiducial 160 pc" is the kpc median multiplied by 100.**
   - **The median.** median(1.63, 1.58, 0.98, 1.79, 1.85, 1.20) = **1.605**.
   - **The range.** "∼100–190 pc" is 0.98–1.85 × 100.
   - **The line-to-line offset.** The [O III] − red-line distance differences in Table 4 are 0.62–1.68 units. Multiplied by 100, they span the "approximately 100 pc" (p. 11) and "100–150 pc" (p. 16) of the text.
   - **The implied gradient.** Sect. 4.2 pairs that "100–150 pc" with "V_LSR differences of 10–16 km s⁻¹". Matching ends gives 100–107 km s⁻¹ kpc⁻¹ in pc (~10× rotation) but 10.0–10.7 in kpc (top of the rotation band).

6. **Downstream propagation.** Within the paper's own method, every physical scale inherits the factor of 10:
   - "∼28 pc" (Abstract; pp. 11, 17) → **281 pc** at 1.605 kpc;
   - "∼3 pc × 28 pc" (Sect. 4.1, p. 15) → 28 × 281 pc;
   - "R ∼ 30 pc" (pp. 11, 16) → ~300 pc;
   - the height above the plane (not quoted) goes from 53 pc to 0.54 kpc;
   - E, R_PDS and age (Sect. 4.1): see Section 6.

**Corroboration from `kd`** (commit 13a8e3f).
- Every docstring and README example gives distances in kpc.
- `pdf_kd` sets `dist_min = 0.01` kpc, the floor in step 4.
- The function default is Reid et al. (2019). The README examples use `reid14_rotcurve`, which is the curve the paper names and the one that matches P1.

![Four panels: rotation V_LSR(d) with Table 4 values under both unit readings; the paper's Fig. 5b [O III] velocity histogram; velocity-space parity of all Table 4 rows; gradients implied by the Monte-Carlo widths.](figures/dispute_distance_fig1.png)

*Figure C1. Table 4 of Wu et al. (2026) compared with the Reid et al. (2014) curve the paper names (`kd` form), on the Table 1 sightlines. All inputs are the paper's printed numbers; no new observations are used.*
- *(a) V_LSR(d) toward P1–P3 (lines), with the range of seven curve variants in grey. [O III] (circles) and Hα (squares) traditional distances are plotted at the printed value read as kpc (filled) and as ×100 pc (open). Dotted lines mark the velocity gap. Inset: 0–300 pc.*
- *(b) [O III] V_LSR histogram digitized from the paper's Fig. 5b (300-dpi page crop; 9 equal bins from −2.0 to 36.7 km s⁻¹; percent axis as in the paper). Orange: velocities rotation can reach within 300 pc. Dotted: 160 pc. Right axis: P1 kinematic distance (`kd`-form Reid et al. 2014).*
- *(c) Rotation velocity at the printed distance against the printed V_LSR ± σ_v, for all 15 rows. Filled = kpc, open = ×100 pc. [S II]/[N II] rows are faded.*
- *(d) Implied gradient W_V/W_d.*
  - *W_d is the printed Monte-Carlo interval (+err plus −err).*
  - *W_V is the 68.3 % width of N(V_LSR, σ_v), truncated at V(0.01 kpc), below which `kd` has no solution. For truncated rows this approximates `kd`'s KDE interval.*
  - *Green: the true gradient, 5.8–10.7 km s⁻¹ kpc⁻¹ (curve variants and Oort term, 0.05–2.5 kpc).*
  - *Red: forbidden, i.e. narrower than σ_v alone allows.*
- *Code: `scripts/dispute_distance_kd.py`. Values: `data/derived/dispute_distance_kd.json` (emulation: `verify_paper_kd.json`).*

**Table C2.** Row-by-row recomputation of Wu et al. (2026) Table 4.
- **Recomputed d:** traditional distance from the printed V_LSR (`kd`-form Reid et al. 2014). Values in parentheses ([O III] and Hα rows) use V_LSR rebuilt from the printed V_helio with the paper's solar motion.
- **z:** (V_LSR − V_rot(d_printed))/σ_v.
- **g_imp:** W_V/W_d in km s⁻¹ kpc⁻¹. The true gradient is 5.8–10.7.
- The two values in the z and g_imp columns are the kpc reading / the ×100 pc reading.

| Region | Line | V_LSR ± σ_v (km/s) | Table 4 MC | Table 4 trad. | Recomputed d (kpc) | z (kpc / ×100 pc) | g_imp (kpc / ×100 pc) | Verdict |
|---|---|---|---|---|---|---|---|---|
| P1 | [O III] | +15.5 ± 4.3 | 1.63 +1.01 −0.84 | 1.79 | 1.785 (1.679) | −0.01 / +3.19 | 4.7 / 46.5 | kpc; ×100 pc rejected |
| P1 | [S II] | −1.0 ± 2.5 | 0.06 +1.17 −0.05 | 0.16 | no solution (V < 0) | −1.04 / −0.46 | 1.7 / 16.7 | not discriminating; MC floor 0.01 |
| P1 | [N II] | −0.9 ± 1.6 | 0.11 +1.12 −0.10 | 0.18 | no solution (V < 0) | −1.68 / −0.68 | 1.0 / 9.9 | not discriminating; MC floor 0.01 |
| P1 | Hα | +4.7 ± 1.2 | 0.47 +1.18 −0.45 | 0.61 | 0.486 (0.432) | −0.95 / +3.41 | 1.5 / 14.7 | kpc; ×100 pc rejected |
| P2 | [O III] | +14.5 ± 4.1 | 1.58 +1.23 −0.83 | 1.85 | 1.743 (1.803) | −0.18 / +3.11 | 4.0 / 39.8 | kpc (via V_helio); ×100 pc rejected |
| P2 | [S II] | +0.5 ± 1.6 | 0.21 +1.06 −0.20 | 0.17 | 0.052 | −0.70 / +0.21 | 1.4 / 14.3 | not discriminating; MC floor 0.01 |
| P2 | [N II] | −0.5 ± 1.6 | 0.08 +1.35 −0.07 | 0.17 | no solution (V < 0) | −1.32 / −0.41 | 0.9 / 9.5 | not discriminating; MC floor 0.01 |
| P2 | Hα | +4.0 ± 1.2 | 0.31 +1.44 −0.31 | 0.54 | 0.430 (0.558) | −0.81 / +2.90 | 1.4 / 13.7 | kpc; ×100 pc disfavoured |
| P3 | [O III] | +11.2 ± 3.3 | 0.98 +1.08 −0.69 | 1.20 | 1.154 (1.156) | −0.12 / +3.01 | 3.7 / 37.3 | kpc; ×100 pc rejected |
| P3 | [S II] | +1.7 ± 1.9 | 0.21 +1.12 −0.21 | 0.23 | 0.161 | −0.37 / +0.77 | 2.1 / 21.4 | not discriminating; MC floor 0.00 |
| P3 | [N II] | +0.4 ± 1.9 | 0.09 +1.02 −0.08 | 0.14 | 0.038 | −0.57 / +0.13 | 1.9 / 18.6 | not discriminating; MC floor 0.01 |
| P3 | Hα | +6.1 ± 1.4 | 0.34 +1.20 −0.33 | 0.58 | 0.600 (0.572) | +0.14 / +3.92 | 1.8 / 18.3 | kpc; ×100 pc rejected |
| P4 | [S II] | +1.1 ± 2.1 | 0.22 +1.61 −0.21 | 0.22 | 0.110 | −0.52 / +0.42 | 1.4 / 14.5 | not discriminating; MC floor 0.01 |
| P4 | [N II] | +1.0 ± 2.3 | 0.28 +0.93 −0.27 | 0.49 | 0.100 | −1.63 / +0.22 | 2.3 / 23.2 | not discriminating; MC floor 0.01 |
| P4 | Hα | +4.5 ± 1.4 | 0.54 +1.18 −0.52 | 0.80 | 0.463 (0.774) | −2.18 / +2.64 | 1.6 / 16.5 | kpc (via V_helio); ×100 pc disfavoured |

Three tabulated [S II]/[N II] velocities are negative and have no outer-Galaxy solution, so their printed "traditional" values cannot come from the tabulated V_LSR. The P2 [O III] and P4 Hα distances are closer to the V_helio-based values than to the printed V_LSR, which suggests that the V_LSR column and the distances were computed separately. Neither point affects the unit test.

**Consequence.**
- **Read correctly**, the paper's kinematics give a fiducial **~1.6 kpc**, not 160 pc. No measurement in the paper supports 160 pc; other evidence does not exclude it, but it does not follow from the paper's data.
- **Within the paper's own method**, every size quoted at 160 pc is 10× too small, and the energetics must be redone (Section 6). The "one of the nearest" claim (Sect. 4.1) loses its basis.
- **We do not claim 1.6 kpc either.**
  - With dV/dd ≤ 10.7 km s⁻¹ kpc⁻¹, a 7 km s⁻¹ non-circular motion shifts the solution by 0.7–1.2 kpc (Section 2).
  - The [O III] V_LSR values are themselves shifted by a wavelength-reference artefact (Section 3).
  - The kinematics do not determine the distance.

**Requested correction.**
1. Table 4, footnotes 1–2: change "in units of 100pc" to "in kpc".
2. Sect. 3.3: change "∼100–190 pc" to "∼1.0–1.9 kpc", and "approximately 100 pc" to "∼0.6–1.7 kpc".
3. Sect. 4.2: change "100–150 pc" to "∼1.0–1.5 kpc".
4. Withdraw the "fiducial distance of ∼160 pc" from the Abstract and Sects. 3.3, 4.1 and 5. Any distance adopted should be non-kinematic and stated as an assumption.
5. Recompute or remove 28 pc, 3 × 28 pc, R_s ∼ 30 pc, E, R_PDS and the age.
6. Publish an erratum to AJ 172, 298.

---

## 2. Why no kinematic distance is possible toward the Silk Nebula

**Claim in the paper.** Sect. 3.3 (p. 10) derives kinematic distances with the `kd` package (Wenger et al. 2018) and the Reid et al. (2014) rotation curve. It excludes objects near the Galactic anticentre (GAC) "owing to the negligible contribution of circular rotation to VLSR along these sightlines" and concedes that the Silk "lies near its boundary, which introduces significant uncertainties into the derived kinematic distance". Even so, it adopts "a fiducial distance of ∼160 pc (the approximate median of the six [O III]-based kinematic distance estimates; see Table 4)" (Sect. 3.3, p. 11). The abstract and Sect. 5 (p. 17) repeat it, and Sects. 3.3 and 4.1 (pp. 11, 15–16) scale size, energy and age from it. For P1 [O III], Table 4 (p. 13) lists RV = 26.8 ± 4.3 km s⁻¹, V_LSR = +15.5 km s⁻¹ (no error printed; we adopt the RV error), a Monte Carlo distance of 1.63 (+1.01/−0.84) and a traditional distance of 1.79. Section 1 deals with the units. Here we show that, in either unit, the velocity carries no useful distance information. The paper flags the uncertainty; our point is that it is total.

**Why it is wrong.** Unless another source is named, every number below comes from `scripts/dispute_kinematic_kd.py` → `data/derived/dispute_kinematic_kd.json`, which re-implements the paper's rotation model from the `kd` source (Persic et al. 1996 curve, Reid et al. 2014 A5 parameters, `kd` projection, IAU standard solar motion). This model reproduces Table 4 (Section 1).

1. **Rotation produces almost no velocity here** (Figure C2a).
   - Across P1–P4 the mean slope dV_LSR/dd over 0–2 kpc is **7.6–9.1 km s⁻¹ kpc⁻¹** (`kd` form, ×cos b projection or flat curve); as d → 0 it is 9.1–10.7.
   - Along the whole [O III] spine it is 4.3–9.8, lowest at the northern end (ℓ = 191.4°).
   - Wenger et al. (2018) define the GAC as ℓ = 180°, where "LSR velocities due to circular motion tend towards zero". Their uncertainty maps mask 160° < ℓ < 200° (captions of their Figs. 5, 8, 19). P1–P4 lie only 2.5–5.5° outside that band. The spine crosses ℓ = 200° at Dec +20.75°: **8.1° of its 16.9° lies inside the band**, 5.9° of it within the paper's own Fig. 1 (Dec ≤ +26.3°; `notes/geometry.md`). Read as a 20° radius about (180°, 0°), the zone misses the spine (24.5–32.8° away), but the gradient is equally small.

2. **Ordinary peculiar velocities equal the whole rotational signal.**
   - Giant molecular clouds within 3 kpc have a one-dimensional cloud-to-cloud dispersion of 7.8 (+0.6/−0.5) km s⁻¹, only ~20 % less without small-scale streaming (Stark & Brand 1989). Anderson et al. (2012) adopt 7 km s⁻¹ streaming motions, which "may be 5 to 10 km s⁻¹". "The vast majority" of masers have peculiar motions ≲ 10 km s⁻¹ (Reid et al. 2019); super-bubbles can accelerate clouds to ≈ 20 km s⁻¹ (Reid et al. 2014).
   - At the slopes of step 1, σ_pec = 5 / 7 / 10 km s⁻¹ means **Δd = 0.55–0.66 / 0.77–0.92 / 1.10–1.32 kpc**.
   - Rotation contributes only 1.4–1.7 km s⁻¹ at 160 pc and 2.6–3.1 km s⁻¹ at 300 pc toward P1–P4 (≤ 3.4 anywhere on the spine). A nearby object's V_LSR is essentially its peculiar velocity.

3. **Gas at a known small distance gets absurd kinematic distances** (Figure C2b).
   - In the Redfield & Linsky (2008) model the P1 sightline crosses the LIC and Gem clouds, both within 15 pc (rotation there: 0.15 km s⁻¹). Their projected V_helio = +19.1 ± 1.1 and +33.9 ± 1.1 km s⁻¹ (calculator output, `data/raw/spectra_ism_lism/dynlism_P1.html`) are V_LSR = +7.0 and +21.8 km s⁻¹.
   - The paper's curve puts them at **0.74 and 2.70 kpc**, 50 and 180 times their maximum distance.
   - The LIC velocity is observed: Ca II K at V_helio +19.1 and +19.4 km s⁻¹ toward 8 Cnc (63 pc, 4.7° from the Silk; `notes/spectra_ism.md`). Gas within 63 pc thus gets a formal distance ≈ 12 times too large. The Gem velocity is a model projection, not detected in the available Ca II/Na I spectra.

4. **The posterior is set by σ_pec, not by the data** (Figure C2c, d; Table C3).
   - For the paper's +15.5 ± 4.3 km s⁻¹, the 16–84 % range widens from 1.28–2.51 kpc (σ_pec = 0) to 1.01–3.31 (σ_pec = 7) and 0.94–3.85 kpc (σ_pec = 10). For the corrected +1.2 ± 1.0 km s⁻¹ (Section 3) it widens from 0.06–0.23 to 0.18–1.26 and 0.26–1.90 kpc. At σ_pec = 7–10 every interval spans a factor of 3.3–7.3.
   - For +1.2 at σ_pec = 7, 44 % of Monte Carlo draws fall below V_rot(0) = 0 and have no solution. The median is 0.12 kpc if they are kept at d = 0 and 0.60 kpc in the posterior; `kd` discards them, biasing it upward.
   - The `kd` Monte Carlo behind Table 4 adds no peculiar-velocity term: Wenger et al. (2018) argue that the Reid et al. (2014) parameter errors "inherently" include streaming. That ensemble argument does not cover the random velocity of one nebula.

5. **The velocity field along the filament does not follow rotation.**
   - At a common distance of 0.3–1.6 kpc, rotation predicts the northern spine (Dec +24°) to be 0.8–4.1 km s⁻¹ *lower* than P1.
   - DESI measures it *higher*: ΔV_LSR ≈ +7…+11 (± ~2.5) km s⁻¹, i.e. +4…+9 km s⁻¹ heliocentric at 2–3.5σ per estimator (`notes/w3_northsouth.md`). SDSS gives a consistent but insignificant +6.7 ± 6.5 km s⁻¹, and the surveys disagree at Dec +17.6…+21°, so this is a ~2–3σ result.
   - Read kinematically, one continuous arc would lie at 0–0.12 kpc in the south (DESI −2.8…−0.2 has no solution; corrected P1 +1.2 gives 0.12) and at 0.86–1.68 kpc in the north.

6. **The paper's own text implies the same.** Sect. 3.3 (p. 10) lists "non-circular motions" as a systematic. Sect. 4.2 (p. 17) proposes that the "∼5–15 km s⁻¹" [O III]/low-ionization offset arises in "different layers of the same expanding shell" (expanding at "∼140 km s⁻¹") or in an inhomogeneously expanding bubble. Section 3 shows the offset is an artefact, but on the paper's own reading either mechanism is a non-circular motion as large as the entire rotational V_LSR at 1.6 kpc toward P1 (14.1 km s⁻¹). Sect. 5 (p. 17) allows "certain systematic uncertainties" yet still adopts the kinematic value.

![Four-panel figure: (a) map of the rotation-velocity gradient over l 180-230, b 0-30 with the Silk spine, P1-P4 and the anticentre exclusion boundaries; (b) rotation-only V_LSR versus distance toward P1 and the northern spine with the paper's, corrected and northern [O III] velocities, and the LIC/Gem local clouds with their formal kinematic distances; (c) and (d) distance posteriors for V_LSR = +15.5 and +1.2 km/s for sigma_pec = 0, 3, 7 and 10 km/s.](figures/dispute_kinematic_fig2.png)

*Figure C2a–d. Kinematic distances toward the Silk Nebula.*
- *Model and code (all panels): the Persic et al. (1996) universal curve with the Reid et al. (2014) A5 parameters as tabulated in `kd`; the `kd` projection; the IAU standard solar motion. Code: `scripts/dispute_kinematic_kd.py`. Numbers: `data/derived/dispute_kinematic_kd.json`.*
- *(a) Mean circular-rotation gradient dV_LSR/dd over 0–2 kpc for ℓ = 180–230°, b = 0–30°, with contours in km s⁻¹ kpc⁻¹. The orange line is the project's [O III] spine (`data/derived/silk_trace.json`), and the white dots are the paper's P1–P4 (Table 1). The dashed line is ℓ = 200°, the edge of the 160° < ℓ < 200° band that Wenger et al. (2018) mask (their Figs. 5, 8 and 19). The dotted line is a 20° radius about (ℓ, b) = (180°, 0°).*
- *(b) Rotation-only V_LSR(d) toward P1 (solid line; grey band ±7 km s⁻¹) and toward the northern spine at Dec +24° (dashed line). The horizontal bands are the paper's P1 [O III] V_LSR (Table 4), the corrected [O III] V_LSR (Section 3; `verify_velocity_summary.json`) and the northern DESI [O III] V_LSR range (`notes/w3_northsouth.md`). The diamonds are the LIC and Gem clouds of the Redfield & Linsky (2008) model, both within 15 pc. Dotted lines lead to their formal (traditional) kinematic distances (open circles).*
- *(c, d) Posterior P(d | V_LSR) toward P1 for σ_pec = 0, 3, 7 and 10 km s⁻¹. The prior is flat on 0–6 kpc and the likelihood Gaussian with σ² = σ_V² + σ_pec². The rotation-curve parameters are marginalised over 800 draws, and the curves are normalised to unit peak. Bars mark 16–84 % intervals and ticks the medians. The orange bar in (c) is the paper's Table 4 Monte Carlo value for P1 [O III], 1.63 (+1.01/−0.84), read as kpc.*

**Table C3.** Kinematic-distance ranges toward P1 (ℓ = 203.8°, b = +19.5°) for each velocity input and peculiar-velocity dispersion σ_pec. "Posterior" uses a flat prior on 0–6 kpc and marginalises over the Reid et al. (2014) parameters. Moving the prior edge to 4 or 10 kpc changes the +15.5 rows at σ_pec = 7 to 1.92 [0.96–2.99] or 2.05 [1.02–3.37], and at σ_pec = 10 to 1.95 [0.84–3.14] or 2.28 [0.97–4.11]. Every other row changes by ≤ 0.05 kpc. "`kd`-style MC" is Wenger et al. (2018) Method C plus a Gaussian σ_pec term; velocities below V_rot(0) are kept at d = 0, whereas `kd` discards them. The lower block lists single-valued (traditional) distances for gas whose location is known independently, and for the two ends of the [O III] spine. All values are computed in `scripts/dispute_kinematic_kd.py`.

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
  - **A physical tie to the Monogem rim.** At Dec < 17.6° the spine lies 12.10 ± 0.12° from the Knies et al. (2024) two-remnant centre, implying 252–318 pc (68 %) only if the association holds.

**Requested correction.**
1. Withdraw the kinematic distance as a result in the abstract, Sect. 3.3, Sects. 4.1–4.2 and Sect. 5. State instead:

   > "Toward (ℓ, b) ≈ (204°, +19°) the circular-rotation gradient is only ≈ 8–9 km s⁻¹ kpc⁻¹ (≈ 4 km s⁻¹ kpc⁻¹ at the northern end, which lies inside the 160° < ℓ < 200° anticentre zone). This is comparable to typical interstellar peculiar velocities of 5–10 km s⁻¹. Kinematic distances are therefore not meaningful for the Silk Nebula, and its distance is not constrained by the present data."

2. In Table 4, either delete both distance columns or label them "kpc; formal values without peculiar motions; not physically meaningful".
3. Give size, radius, energy and age as explicit functions of d. For example, the spine length is ≈ 0.30 d (16.9° of arc).
4. Delete the "100–150 pc" blue/red distance discrepancy (Sect. 4.2).

---

## 3. The [O III] velocity offset is a wavelength-reference artefact

**Claim in the paper.**
- Sect. 3.2.3 (p. 10) says the MRS-N spectra "have been corrected for heliocentric motion using the LAMOST official pipeline", so the emission-line RVs "represent heliocentric velocities". Sect. 2.2 (p. 6) quotes a precision of "< 1.5 km s⁻¹" after the secondary wavelength calibration of Ren et al. (2021).
- The paper's Table 4 (p. 13) gives [O III] RVs for P1/P2/P3 of +26.8 ± 4.3, +26.8 ± 4.1 and +23.8 ± 3.3 km s⁻¹ (V_LSR 15.5/14.5/11.2). For [N II]/[S II] it gives +11.7…+13.3.
- The resulting [O III]–red-line offset is "∼5−15 km s⁻¹" (Sect. 3.2.3), or "∼10–16 km s⁻¹" in V_LSR (Sect. 3.3, p. 11). It underpins three things:
  - a "distance discrepancy of 100–150 pc" between blue- and red-arm distances (Sect. 4.2, p. 16);
  - the "different layers of the same expanding shell" interpretation (Sect. 4.2, p. 17);
  - the [O III]-only kinematic distance (Abstract; Sect. 3.3; Sect. 5).

**Why it is wrong.** The red-arm [N II]/[S II] velocities and the LAMOST stellar heliocentric correction are both sound (steps 5–6). The blue-arm [O III] velocities are not. Numbers come from `scripts/dispute_oiii_values.py` → `data/derived/dispute_oiii_values.json`, or from the project files named. A side-by-side table of all [O III] and red-line velocities, with their frames, is in [VELOCITY.md](VELOCITY.md) §5. 

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
   - Per region (Table C4), the paper minus DESI within 0.5° is +14.0 ± 2.4 km s⁻¹ (5.7σ), and +3.9 ± 2.5 (1.5σ) after correction.
   - The weighted mean of DESI, Fig. 3a and SDSS is **+13.3 ± 1.0 km s⁻¹** (χ² = 2.2 for 2 dof). Other estimator choices give +12.4 ± 0.9, hence a ±1 km s⁻¹ systematic. In V_LSR this is +1.2 ± 1.0, or +1.5 with Schönrich et al. (2010).
   - The velocity varies along the filament. DESI gives +9.7…+12.2 km s⁻¹ at Dec < 17.6° and +17…+23 at Dec 17.6–21°, where SDSS gives +10.0 ± 2.2 (`w3_northsouth.md`). P1 and P3 lie in the southern zone; P2 (Dec 17.74°) is at its edge, with a region-matched DESI value of +15.2 ± 2.9 (n = 9). The 10.1 km s⁻¹ error belongs to the wavelength scale and is unaffected.
   - The corrected [O III] agrees with the paper's own [N II] (P1–P4 mean +12.2 ± 0.9) and with WIM Hα from WHAM (Haffner et al. 2003). WHAM gives +11.5…+13.4 in the nearest beams, and V_LSR +1.3 on the southern filament against +1.2 ± 1.0 for [O III] (`verify_hiwham.md`).
   - [O III] − [N II] becomes **+1.1 ± 1.3 km s⁻¹** (+1.4 ± 1.5 from Fig. 3a alone), against +13.7 ± 2.4 from the paper's P1–P3 values.

![Three panels. (a) The [N I] 5198/5200 airglow doublet from the paper's Fig. 3a spectrum, with markers for the doublet positions on air or vacuum and topocentric or heliocentric scales; the data sit on the air-topocentric markers. (b) [O III] 5007 from Fig. 3a and the Fig. 3b P1-P3 spectra, with air and vacuum rest wavelengths and the plotted-scale positions of lines with true heliocentric velocity +13.3 and +25 km/s; the line peaks at the +13.3 position. (c) Velocity ladder: the paper's published Table 4 [O III], the paper's plotted-scale values, corrected values, DESI, SDSS, the paper's red lines, WHAM H-alpha and the paper's H-alpha, with the adopted +13.3 +- 1.0 km/s band.](figures/dispute_oiii_fig3.png)

*Figure C3a–c. The paper's [O III] velocities are on a topocentric, air-referenced scale.*
- *All panels: spectra re-extracted from the vector polylines of arXiv:2610.13329v1, p. 7 (`scripts/verify_velocity_fig3_vectors.py`). Values from `scripts/dispute_oiii_values.py`; plot from `scripts/dispute_oiii_figure.py`. v_corr = −8.32 km s⁻¹ (astropy, Xinglong, 2024-02-01).*
- *(a) [N I] λλ5198, 5200 in Fig. 3a, continuum-subtracted (white), with a two-Gaussian fit (blue). Vertical lines show where the doublet falls on each candidate scale: air-topocentric (green), air-heliocentric (orange dashed), vacuum-topocentric (cyan dotted) and vacuum-heliocentric, the native LAMOST frame (purple dash-dot). Rest wavelengths are from NIST ASD level energies; air values follow Morton (2000).*
- *(b) [O III] λ5007 in Fig. 3a (white; fit in blue) and the Fig. 3b P1–P3 spectra (thin, each normalised by its own fit). Top axis: velocity read off the plotted scale, as in the paper's Table 4. Green and red lines: where lines with true V_helio = +13.3 and +25 km s⁻¹ fall on this scale.*
- *(c) V_helio (bottom axis) and V_LSR ≈ V_helio − 12.2 (top axis; ±0.4 km s⁻¹ over P1–P4).*
  - *Region rows: the paper's Table 4 (red), Fig. 3b as plotted (orange squares), the paper's Table 4 − 10.1 (green), DESI positive-residual science fibres within 0.5° (blue; `verify_desi_velocity.json`).*
  - *Fig. 3a row: the plotted-scale and airglow-referenced values, with the velocity budget.*
  - *DESI row: on-ridge stack and joint regression (16–84 %). SDSS row: all-strand ON−OFF stack and the Dec < 20° stack of the second re-analysis.*
  - *Lower block: the paper's [N II] (purple) and [S II] (×), the WHAM-SS main Hα component in the nearest beam (`spectra_ism_wham_results.json`), and the paper's Hα (grey; treated in Section 4).*
  - *Green band: the adopted +13.3 ± 1.0 km s⁻¹.*

**Table C4.** [O III] velocities by region (km s⁻¹): the paper's Table 4, the paper's Fig. 3b re-measured, corrected, and compared with DESI, the red lines and WHAM.
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

---

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
   - **Per region.** In the nearest 1° beams of WHAM-SS DR1 (Haffner et al. 2003), the main component lies at V_helio +11.6…+13.4 (V_LSR −0.3…+1.3; FWHM_int 28–31 km s⁻¹; total intensity 2.9–4.4 R). Hα − WHAM = **+4.8 ± 0.7 km s⁻¹**, whereas [N II] − WHAM = −0.8 ± 1.0 (Table C5).
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
   - **Sign test.** With the correct sky intensity or too little subtracted, no δ in −16…+20 km s⁻¹ shifts the centroid by more than +0.4 km s⁻¹. The pipeline placement (δ = 0) gives −0.0…−0.9 km s⁻¹. The literal "6563/6554 × RV(λ6554)" placement gives −0.6…−1.7 km s⁻¹. These ranges cover both normalisations of step 6. A redward shift needs a negative residual on the blue wing (Figure C4a, b).
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

*Figure C4a–d. (a) Forward model for P1 at the adopted joint fit, which fits centroids and widths with a self-consistent sky/nebula ratio (step 6). It shows: the WIM Hα (V_helio +11.8, the Table 4 [N II] velocity, with intrinsic FWHM 28.7 km s⁻¹ from WHAM); the true geocoronal line at v_corr − 2.6 km s⁻¹ (NIST 2s–3p, 2:1); the MRS-N template, 30 % too bright, 2.5 km s⁻¹ from RV(λ6554); the residual; and the Gaussian + quadratic fit, whose centroid is shifted by +3.8 km s⁻¹. All curves are convolved with the 0.32 Å LSF; the ticks mark the centroids. (b) Fitted centroid shift against over-subtraction (model/true − 1) for P1–P4, at the joint-fit template position. The sky/nebula ratio is made self-consistent at each point; below −30 % the iteration does not converge. The red band is the observed Table 4 Hα − [N II] offset (dark: weighted mean ± 1σ; light: P1–P4 range). The green band and line show the joint fit (1σ range, best value). The dotted line shows the centroid-only fit, whose template position (+2.25 km s⁻¹) is close to that of the curves. (c) Fitted intrinsic Hα FWHM for the same models, with the Table 2 Hα widths (weighted mean ± 1σ) and the WHAM widths used as the true WIM width. (d) V_helio of LAMOST Hα, [N II] and [S II] (Table 4); Hα minus the joint-fit model shift; WHAM-SS DR1 Hα from the nearest beam (error: beam-to-beam scatter within 1°, floor 0.5 km s⁻¹); and the [O III] velocity adopted in Sect. 3. Top axis: V_LSR ≈ V_helio − 12.2 km s⁻¹ (mean correction for P1–P4). Model code: `scripts/verify_velocity_halpha_sim.py` (reused) and `scripts/dispute_halpha_analysis.py`; figure: `scripts/dispute_halpha_figure.py`.*

**Table C5.** Per-region Hα velocities and offsets, in km s⁻¹ (V_helio unless stated). LAMOST values are from Wu et al. Table 4; sky/nebula ratios use Table 2. WHAM values are the main-component velocity of the nearest WHAM-SS DR1 beam (`data/derived/verify_hiwham_wham_beams.csv`), converted with the paper's solar motion (V_LSR − V_helio = −11.9…−12.6). The quoted error is the beam-to-beam scatter within 1°; a floor of 0.5 is applied in the differences. Model rows and distances come from `data/derived/dispute_halpha_results.json`. Distances use our re-implementation of the `kd` (Wenger et al. 2018) traditional solution with the Reid et al. (2014) curve (`scripts/verify_paper_kd.py`).

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
- **The Hα distances come from the offset.** At the local gradient of 8.9–9.8 km s⁻¹ kpc⁻¹, the Hα − [N II] offset alone is worth 0.49–0.69 kpc. The printed Hα distances (0.31–0.80, in kpc per Section 1) owe their size to it. With the [N II] or WHAM velocity, the same rotation model gives no solution or 0.01–0.13 kpc, which carries no information here (Section 2).
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

---

## 5. Table 3 overestimates the [O III] surface brightness by a factor of 4.6–8.6

**Claim in the paper.** Table 3 (p. 8) lists the [O III]-band surface brightness of P1, P2 and P3 as 23.60 ± 0.028, 23.87 ± 0.035 and 23.85 ± 0.039 mag arcsec⁻². The footnote says these were "derived by interpolation from Gaia and SDSS photometry". The zero points, 11.23–11.32, have no stated units.

Sect. 3.1 (p. 6) describes the method. The team measured "star-free sub-regions within P1–P3 and adjacent nebula-free background regions" and calibrated them on field stars whose Gaia and SDSS magnitudes were transformed to a common system following Jordi et al. (2010). It notes that P1 "appears the faintest in the narrowband image" yet "yields the highest measured surface brightness", and attributes this to stellar light.

Sect. 4.1 (p. 15) converts the values to "∼10–13 Rayleigh — well above the 0.5–5 Rayleigh detection limit of the Fesen et al. (2024) [O III] images". The survey camera (Sect. 2.1, p. 4) is a 50 mm f/0.95 lens with a 6.0 nm FWHM filter in its aperture and ~50″ pixels, under an SQM ≈ 21.5 mag arcsec⁻² sky; a 180 mm f/2.8 reflector was used for follow-up.

**What is right.** The arithmetic is correct. An AB magnitude (Oke & Gunn 1983) in a 60 Å band that contains only λ5007 gives 12.75 / 9.94 / 10.12 R [calc]. The filament is real, and Fesen et al. (2026) imaged it.

**Why it is wrong.**

1. **Conversion choices cannot rescue the values** (Figure C5a). We varied three choices:
   - band width, 50–70 Å;
   - magnitude system, ±0.1 mag (Blanton & Roweis 2007, Table 1: AB − Vega = −0.08 mag in SDSS g, +0.02 mag in V);
   - λ4959 leakage, from none to full.

   Over these choices Table 3 spans 5.7–16.3 R [calc]. Even a 30 Å band gives 5.0–6.4 R. The paper does not say which instrument produced Table 3. With the 5 nm filter that QHYCCD (2024) reports for the f/2.8 follow-up, the values are 8.3–10.6 R.
2. **Independent measurements give 1.4–2.5 R** (Table C6, Figure C5b). Like Table 3, these are excesses over the local background; the SDSS P-region means are absolute, but the off-strand level is only 0.07 R.
   - DESI DR1 (DESI Collaboration 2026): 1.62 R.
   - SDSS/BOSS sky fibres (DR17, DR19; Abdurro'uf et al. 2022; SDSS Collaboration 2026): 1.49 / 1.78 / 1.87 R at P1 / P2 / P3.
   - SDSS imaging sky levels (DR18; Almeida et al. 2023): 2.31 ± 0.38 R on 17′ scales (3.9σ against coherent nulls).
   - NSNS DR0.2 maps (Ziegenbalg 2025; calibration not re-verified), measured here: 2.04 / 2.06 / 2.18 R, at a 51.5″ pixel scale.

   Over the same regions Table 3 is **4.6–8.6× higher** [calc]. That is 1.7–2.3 mag, or 43–83 times its quoted errors. Without the star-affected P1 the factor is 4.6–6.2×.
3. **Angular scale does not explain the gap.** The Table 3 subregion averages are 2.6–4.3× above even the brightest single measurements [calc]: the brightest SDSS fibre (3.65 R), the brightest DESI node (3.88 R) and the NSNS 1′ peaks (2.95–3.64 R).
4. **The team's own earlier value is fainter.** QHYCCD (2024) reports the SWONS result as 24.3 mag arcsec⁻² for the *brightest* region, calibrated on V-band stars. That is 3.3–6.7 R for a 30–60 Å band, and 0.43–0.70 mag fainter than the Table 3 *region* values [calc].
5. **No identified physical multiplicative factor is large enough** (Figure C5d). As I_true/I_Table3 [calc]:
   - AB vs Vega: 0.91–1.10.
   - Band width: 0.83–1.17.
   - λ4959 leakage: 0.75–1.00.
   - Filter in the lens aperture (Sect. 2.1): 1.01–1.27 at field angles 0–5°, and 0.93–2.8 at 10–15°. These come from a generic 60 Å 3-cavity filter model, not the real SWONS curve.
   - Filter in an f/2.8 beam: 0.97–1.08.
   - An untuned filter in the f/0.95 converging beam gives 2.0–4.0, which would make Table 3 an *under*estimate.

   As a deliberately generous bound we stacked the separate lowering extremes: 50 Å, +0.1 mag and full λ4959. No single filter satisfies these together, and the filter model already contains the width and λ4959 terms. We added ~0.3 mag of P1 starlight. The result, 5.5–5.8 R, is still **2.6–3.7×** above the independent values. Matching them needs a factor of 0.12–0.22.

   A common factor would also keep the regions in order, but they are not in order.
   - Table 3 ranks P1 brightest, by 0.25–0.27 mag.
   - SDSS ranks P1 faintest, by 0.29–0.38 R (errors ≤ 0.09 R).
   - NSNS finds the three regions equal within 0.14 R.
   - The paper's own image shows P1 faintest.

   So at least the P1 excess, 2.8–3.1 R, is additive, as Sect. 3.1 itself suggests.
6. **An additive error of the required size is plausible** (Figure C5c). The excess is 7.9–11.3 R [calc]. Two sources could supply it.
   - **Background or flat-field residual.** The excess is 6–20 % of an in-band sky of 21.1–22.0 mag arcsec⁻², or 9–13 % at the paper's SQM reading of 21.5. SDSS/BOSS sky spectra give 21.1–21.5 AB at 5070 Å; the 22.0 limit allows for a darker narrow-band continuum. A residual this size is large for nearby apertures. However, the paper notes that dome and twilight flats "introduced significant vignetting" in this f/0.95 system (Sect. 2.1).
   - **Mismatched starlight.** Gaia DR3 stars with G > 14 within 0.25° of P1–P3 add up to 3.0–4.4 R-equivalent (Gaia Collaboration 2023; calc). Even with a perfectly star-free background aperture, that is only 27–56 % of the excess, so brighter stars or their halos would also be needed. P1 has 15.5–20.6 R-equivalent more starlight than P2/P3, and its extra 2.8–3.1 R is 13–20 % of that surplus.

   Without the frames and aperture definitions we cannot say which additive term dominates. Nor can we exclude a bookkeeping error: aperture area, pixel scale, exposure normalisation or the unitless zero points. The identified calibration factors, however, cannot account for the excess.

![Four-panel check of Table 3 against independent [O III] data](figures/dispute_brightness_5.png)

*Figure C5. [O III] λ5007 surface brightness of the Silk Nebula: Wu et al. (2026) Table 3 compared with independent data.*
- *(a) Conversion of an [O III]-band AB surface brightness to I(5007) for 30–70 Å effective widths, with the band containing only the line. The white line is 60 Å, the stated SWONS FWHM. Grey shading is ±0.1 mag (AB/Vega) at 60 Å. Red symbols are Table 3: P1 circle, P2 square, P3 triangle. Orange is the QHYCCD-note 24.3 mag arcsec⁻² (bar 30–60 Å, diamond 50 Å). Green is the independent 1.4–2.5 R.*
- *(b) All measurements on one scale, with the angular scale of each at right. Symbol shape gives the region (P1/P2/P3) or a ridge-wide value (diamond; the orange diamond is the QHYCCD value at 50 Å). Filled symbols are region means or medians. Open symbols are peaks or alternative estimators: the star is the brightest SDSS fibre, and the DESI open square is the W-strand node at Dec 18° near P2. The NSNS P1 and P2 symbols overlap. The hatched band is the 0.5–5 R range quoted in Sect. 4.1.*
- *(c) Additive background error, as a percentage of the in-band sky, needed to raise I_true to the Table 3 mean of 10.9 R: (I_T3 − I_true)/I_sky. The box spans the measured I_true and plausible skies, 7–17 % for the mean and 6–20 % region by region (P1–P3 matched to NSNS or SDSS). The dashed line is the paper's SQM reading.*
- *(d) Multiplicative systematics as I_true/I_Table3. Filter cases use a generic 3-cavity, 60 Å FWHM interference-filter model (n\* 1.8–2.05, CWL 5007–5030 Å), not the measured SWONS curve. "At stop" is the filter in the lens aperture (Sect. 2.1), where each field point sees a collimated beam. The green band is the required factor. The stack combines the separate width, magnitude-system and λ4959 lowering extremes with ~0.3 mag of P1 starlight.*

*Data: the paper's Sects. 2.1, 3.1, 4.1 and Table 3; DESI DR1; SDSS DR17/DR19 sky fibres and DR18 imaging-sky levels; NSNS DR0.2; Gaia DR3. Computed by `scripts/dispute_brightness_calc.py` and plotted by `scripts/dispute_brightness_fig.py`.*

**Table C6.** [O III] λ5007 brightness of the southern Silk Nebula by method. The last column is the Table 3 value divided by the measurement, region-matched where possible. Values for P1 / P2 / P3 are separated by "/". Sources: the paper; `data/derived/dispute_brightness_calc.json`; the project notes `verify_desi.md`, `verify_sdss.md`, `verify_precovery.md`, `imaging_optical.md` and `spectra_sdss_lamost.md`.

| Measurement | Method | Scale | I(5007) (R) | Table 3 / value |
|---|---|---|---|---|
| Wu et al. (2026) Table 3 | SWONS narrow-band photometry, Gaia/SDSS zero points; converted here (AB, 60 Å, line only) | star-free subregions, ~50″ px | 12.75 / 9.94 / 10.12 (all variants 5.7–16.3) | 1 |
| Table 3, all lowering extremes stacked | 50 Å, +0.1 mag, full λ4959, P1 −0.3 mag, applied together | — | 5.50 / 5.66 / 5.76 | 2.3 / 1.8 / 1.8 |
| QHYCCD note (SWONS team result, 2024) | 24.3 mag arcsec⁻², V-calibrated | brightest region | 3.3–6.7 (30–60 Å) | 1.5–3.8 |
| LAMOST MRS-N relative × sky-line scale | Table 2 ratios × [N I] 5198 / OH 6553 sky-line scale, APO or UVES (Hanuschik 2003) ratio | fibres | 2.8–3.2 (APO) or 14–16 (UVES) | indeterminate; not used |
| NSNS DR0.2 | calibrated narrow-band survey, on-ridge minus off-ridge (not re-verified) | ≤ 0.25° of ridge, 51.5″ px | 2.04 / 2.06 / 2.18 (1′ peaks 2.95 / 3.64 / 3.37) | 6.3 / 4.8 / 4.6 (peaks 4.3 / 2.7 / 3.0) |
| SDSS imaging sky, g (19 runs, 2002–07) | ridge-template fit to field sky levels; 3.9σ vs coherent nulls | 17′ | 2.31 ± 0.38 (2.39, 2.55 other estimators) | 3.9–5.5 |
| SDSS/BOSS sky fibres | total sky, fibre-acceptance corrected, mean within 0.3° | 2–3″ fibres | 1.49 / 1.78 / 1.87 (±0.05–0.09; ±0.2 zero point); brightest 3.65 | 8.6 / 5.6 / 5.4 |
| SDSS/BOSS south ON−OFF | strand vs 0.45–2° background | r_on 0.15° | 1.39 ± 0.10 | 7.2–9.2 |
| DESI DR1 south ridge | flux-calibrated, GFA fibre acceptance | 1.5″ fibres | 1.62 (±8 % sys; nodes 1.4–2.3; W node 3.88) | 6.1–7.9 |
| Fesen et al. 2024 limits (as quoted) | — | — | 0.5–5 | — |

**Consequence.**
- **The Sect. 4.1 comparison fails.** At 1.4–2.5 R the Silk Nebula lies *inside* the quoted 0.5–5 R range of Fesen et al. (2024), near its logarithmic midpoint of 1.6 R. It is not "well above" that range. Its appearance in Fesen et al. (2026) is consistent with this brightness and does not require a bright filament.
- **Downstream quantities.** The paper's energy uses an adopted n₀ = 1.0 ± 0.3 cm⁻³, not Table 3. Anything scaled from Table 3 would be inflated [calc]: emission measure and [O III] luminosity by 4.6–8.6×, an emission-measure density by 2.2–2.9×, and the Sect. 4.1 energy (E ∝ n₀^1.16) by 6–12× if n₀ were tied to I. The "tenuous" of Sect. 3.1 is, if anything, strengthened.
- **Under the paper's own shock interpretation, Table 3 conflicts with the GALEX FUV.**
  - Regressing GALEX FUV (Martin et al. 2005; Morrissey et al. 2007) on NSNS gives 92 ± 10 CU (photons s⁻¹ cm⁻² sr⁻¹ Å⁻¹) per NSNS rayleigh in the south (`verify_xrayfuv_galex_measure.json`).
  - If Table 3 were right, that would be 15–20 CU per true rayleigh [calc]. This is 2.9–3.8× below (2.6× at +1σ) the lowest non-resonant FUV floor of the tabulated shock models at E(B−V) = 0.03.
  - That floor is 57 CU R⁻¹ for the incomplete shocks of Raymond et al. (1988). For the complete shocks of Raymond (1979), Cox & Raymond (1985) and Dopita & Sutherland (1996) it is 85–281 CU R⁻¹.
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

---

## 6. Further problems: red-line diagnostics, [S II] density, temperature/turbulence/metallicity, energetics, line widths

Section numbers and page numbers (p.) refer to the arXiv v1 PDF of Wu et al. (2026). Unless another project file is named, recomputed values come from `scripts/dispute_other_fig6.py`, which writes `data/derived/dispute_other_values.json`. That script reuses the project's 5-level S⁺ atom (`verify_paper_sii.py`), the energetics formulae (`verify_paper_energetics.py`) and the line-fit Monte Carlo stored in `data/derived/paper_checks.json`.

Several parts of the paper are correct and we reproduce them:
- the coefficients of Eqs. 2 and 3 (23.5 vs 23.55 derived);
- the energetics arithmetic;
- the statistical [O III] velocity errors of 3–4 km s⁻¹;
- the caution in Sect. 4.2 that the red-arm lines are dominated by diffuse ionized gas (DIG).

The problems are in how these numbers are interpreted.

**Claim in the paper.**
1. The red-arm ratios establish shock excitation. The abstract (p. 1) says [S II]/Hα "significantly exceeds the shock-excitation threshold of 0.4"; Sect. 5 (p. 17) says it exceeds "the classical threshold of 0.4". Sect. 4.1 (p. 15) says the ratio "exceeds 0.5 in all four regions, confirming that the ionized gas is shock-excited". It adds that [N II]/Hα = 0.5–1.5 is "a condition satisfied by regions P1, P2, P3, and P4". Fig. 8 places the points in the SNR zone of the Sabin et al. (2013) and Kopsacheili et al. (2020) diagrams.
2. The [S II] doublet gives a median n_e ≈ 64 cm⁻³, rising to ~150 cm⁻³ toward the centre, with "lower electron densities" in the background (Sect. 3.4, Eqs. 4–5, Fig. 7b, pp. 13–14). Sect. 4.1 calls this "the electron density of the [O iii]-emitting gas" (p. 15) and uses it as "the post-shock shell electron density" (p. 16).
3. Line widths give the following (Sect. 3.4, pp. 11–13; Table 4, p. 13):
   - T_e from Eq. 2 (the Reynolds et al. 1977 method): "about 9000 K", "1-2 ×10⁴ K" for most pixels on the [O III] structure against ≈ 7.5×10³ K in P4 (p. 12), and an adopted "approximately 10,000 K" (p. 13);
   - V_t = 29–33 km s⁻¹ from Eq. 3, attributed to "shock-driven gas flows" (p. 13).
4. Near-solar, uniform O/H: 8.82 from N2 (Eq. 6, described as "a third-order polynomial") and 8.70 from N2S2Hα (Eqs. 7–8) (Sect. 3.5, pp. 14–15).
5. The energetics are V_s = 119 ± 5 km s⁻¹, E = (2.0 ± 0.7)×10⁵¹ erg, R_PDS = 18 ± 2.8 pc and t = (6.8 ± 1.8)×10⁴ yr (Sect. 4.1, p. 16).
6. The intrinsic [O III] FWHM is 43.1–48.9 km s⁻¹ (Table 2, p. 7). Sect. 4.1 (p. 15) refers to "strong [O III] λλ4959,5007" emission.
7. "No radio counterpart" in NVSS or TGSS is "consistent with the old age" (Sect. 4.1, p. 16; Sect. 5).

**Why it is wrong.**

1. **The red-line ratios are those of the nebula-free background and of the faint warm ionized medium (WIM).**
   - We recompute the ratios from the Table 2 intensities. On the filament, [S II](6716+6731)/Hα = 0.90 ± 0.11, 0.86 ± 0.09 and 0.81 ± 0.13 (P1–P3); in P4 it is 0.80 ± 0.13. For [N II]6583/Hα the values are 0.80, 0.70 and 0.77 (±0.09–0.10), against 0.80 ± 0.10 in P4.
   - The on-minus-P4 differences are +0.04 to +0.58σ for [S II] and −0.72 to 0.00σ for [N II] (Figure C6a).
   - By the paper's own criteria, the nebula-free P4 is therefore "shock-excited" as well.
   - These are standard values for faint high-latitude WIM. Madsen et al. (2006) find [S II]6716/Hα ≳ 0.4 "in regions of faint Hα emission". That corresponds to ≳ 0.68 for the doublet sum at the low-density limit (derived). They also find [N II]/Hα ≳ 0.7 in the faintest WIM at l > 150°.
   - O-star H II regions average only 0.11 ([S II]) and 0.27 ([N II]) (Madsen et al. 2006). P1–P4 lie at l = 202.5–205.5°, b = +18.0° to +19.6° (Table 1 positions, converted), where the nearest WHAM-SS beams give just 2.9–4.4 R of Hα.
   - Independent data show no Hα filament on top of this WIM. WHAM Hα is the same on and off the strand (paired on−off ΔI = +0.04 ± 0.35 R south, +0.14 ± 0.27 R north; intrinsic FWHM 30.4–31.9 km s⁻¹ on and off). At fibre scale, SDSS sky fibres give an on−off Hα excess of 0.08 ± 0.22 R at Dec < 17.6° (`notes/verify_hiwham.md`).
   - Both ratios rise as I(Hα) falls, which Haffner et al. (1999) and Madsen et al. (2006) attribute to WIM temperatures of 6000–9000 K (for a general review of the WIM see Haffner et al. 2009).
   - The 0.4–0.5 thresholds separate SNRs from H II regions. The faint WIM meets them anyway, so they cannot separate an SNR from the WIM at a few rayleighs.
   - In all three panels of Fig. 8 the grey background points overlap the nebula points (visual inspection of the published figure). The Sabin et al. (2013) diagrams and the Kopsacheili et al. (2020) panel (log [N II]/Hα vs log [S II]/Hα) use the same red-line ratios and inherit the same limitation.
   - The paper's own Sect. 4.2 (pp. 16–17) lists three reasons why the red-arm lines are DIG. We agree with it, but it contradicts the use of these lines as shock evidence in Sect. 4.1 and the Conclusion.
2. **The quoted densities do not follow from Table 2, and densities of this order are unphysical for this gas.**
   - Table 2 gives R = I(6716)/I(6731) = 1.10, 1.14, 1.19 and 1.03 (P1–P4).
   - The paper's own Eq. 5 (Proxauf et al. 2014) gives 298, 247, 189 and 400 cm⁻³ at 10⁴ K, or down to 169 cm⁻³ with a √T scaling to 8000 K.
   - A 5-level S⁺ atom gives 243–631 cm⁻³ at 8000–10⁴ K, or 225–631 cm⁻³ including the Ramsbottom et al. (1996) Ω (Figure C6b). It uses Tayal & Zatsarinny (2010) collision strengths with TZ10, Mendoza & Zeippen (1982) or Rynkun et al. (2019) A-values.
   - Reaching 64 or 150 cm⁻³ needs R = 1.33 or 1.23 with Eq. 5, and R = 1.37 or 1.27–1.28 with the 5-level atom. No region in Table 2 reaches R = 1.23.
   - The 5-level low-density limit is 1.437–1.462; the paper quotes ~1.42 (pp. 14–15).
   - **Uncertainty of R.** Table 2 gives no error on R, and two readings are possible.
     - With the [S II] S/N of 44–65 listed in Table 2, σ_R = 0.026–0.033. The Table 2 ratios then lie 5.6–9.2σ below the R = 1.33 that Eq. 5 needs for 64 cm⁻³.
     - If instead the ±0.02 intensity errors are independent for the two lines, σ_R = 0.20–0.30. Each region's density is then unconstrained (16–84 %: 11–1900 cm⁻³; `verify_paper_sii.json`).
     - Neither reading supports the specific values 64–150 cm⁻³.
   - In the Table 2 central values, P4 has the lowest R and therefore the highest density. That is the opposite of the statement that "background regions show lower electron densities" (p. 14). The difference from P1–P3 is 1.6–3.8σ under the S/N reading and insignificant under the other, so Table 2 does not support the statement.
   - Can noise bias the Fig. 7b median? Our simulation of per-fibre noise under Eq. 5 (`verify_paper_sii.json`) shows it cannot turn R = 1.03–1.19 into a median of 64 cm⁻³.
   - **Emission-measure argument.** WHAM gives EM = 6.4–9.9 pc cm⁻⁶ on these sightlines, using 2.25 pc cm⁻⁶ R⁻¹ at 8000 K (Haffner et al. 2003).
     - At the Table 2 densities, all of that emission measure would have to come from a line-of-sight path of only 6–22 AU (5-level atom) or 12–46 AU (Eq. 5).
     - A 1.37-pc filament (0.28° at 280 pc) is 6×10³–5×10⁴ times longer than that path. At the Table 2 densities it would therefore be correspondingly brighter than the few rayleighs of Hα that are seen on and off the strand alike.
     - Even at 64–150 cm⁻³ the path would be only (0.3–2.4)×10⁻³ pc.
     - Spread over a 1–3 kpc path, the same EM implies ⟨n_e²⟩^½ ≈ 0.05–0.10 cm⁻³ (derived). A path of that length corresponds to the ~1 kpc electron scale height (Haffner et al. 1999) seen at b ≈ 20°.
   - **Direct measurement on this sightline.** Our regressions of the WIM component give 6716/6731 = 1.40 (1.34–1.50) and 1.53 (1.42–1.67) from SDSS sky fibres, and 1.40 (1.33–1.49) from DESI (`w3_oii_summary.json`; `notes/w3_oii.md`). These are at the low-density limit.
   - The Table 2 ratios are therefore most plausibly a measurement systematic, not a gas density. We cannot identify its cause from the published information.
   - **The [O III] gas.** The paper assigns this density to "the [O iii]-emitting gas" (p. 15). We tested that with gas at 64 cm⁻³ filling the 0.28° strand (1.37 pc at 280 pc), with x(O²⁺) = 1 and solar O/H. It would emit 6×10³–2×10⁴ R in λ5007 at 1–2×10⁴ K, and 3.5×10⁴–1.2×10⁵ R at 150 cm⁻³. The measured brightness is ≈ 1.4–2.5 R (Section 5). A brightness of 2–4 R requires n_e ≈ 0.1–3 cm⁻³, for geometries from an edge-on sheet to a filament (`notes/w3_canon.md` §7).
3. **T_e is unconstrained, and neither T_e nor V_t differs between the filament and the background.**
   - Eq. 2 applied to the Table 2 region widths gives T = 6780, 6580, 3220 and 7660 K (P1–P4). The P4 value reproduces the paper's own ≈ 7.5×10³ K.
   - The 1σ errors are σ_T = 7300–10 700 K, and P(T < 0) = 18–38 % (Monte Carlo on the quoted width errors; Figure C6c).
   - In Table 2 the on-filament regions are not hotter than P4. The "1–2×10⁴ K" versus "7.5×10³ K" contrast is therefore not supported by the paper's own widths.
   - The spread in Fig. 7a is no larger than this noise. Its legend gives a Gaussian σ = 6798 K. The caption's "σ ≃ 6.47×10³ K" matches the legend's μ = 6471 K instead.
   - Eq. 3 applied to the same widths gives V_t = 26.7, 25.7, 29.9 and 26.5 km s⁻¹. Table 4 lists 33.1, 32.9, 30.6 and 29.9, which are higher by +0.7 to +7.2 km s⁻¹. The paper does not say how the Table 4 values were obtained (for example, as per-fibre medians).
   - The nebula-free P4 has the same V_t (Table 4: 29.9 against 30.6–33.1). The broadening is therefore a property of the line-of-sight WIM, not of a shock.
4. **The abundances describe the WIM excitation, not the filament.**
   - N2 = log([N II]6583/Hα) = −0.153 to −0.097.
   - This is outside the validity ranges of Pettini & Pagel (2004; −2.5 < N2 < −0.3) and Marino et al. (2013; −1.6 < N2 < −0.2).
   - Eq. 6 is the *linear* Pettini & Pagel (2004) relation. It gives 8.81–8.85 (paper 8.82). PP04's third-order polynomial would give 9.09–9.18.
   - Dopita et al. (2016) gives 8.64–8.75 (paper 8.70). The nebula-free P4 gives 8.85 and 8.75. The paper itself notes that the background "yields relatively consistent results" (p. 15).
   - Both calibrations assume gas photoionized by stars. In the WIM, a high [N II]/Hα traces a higher T_e (Haffner et al. 1999; Madsen et al. 2006), not a higher N/H.
5. **The energetics arithmetic is reproduced, but the printed age formula is wrong and the inputs are invalid.**
   - Reproduced: V_s = 119.3 km s⁻¹ (n_e = 64, n₀ = 1) and E = 2.01×10⁵¹ erg.
   - R_PDS = 17.1 pc, not 18; 18 pc would require E = 2.41×10⁵¹ erg.
   - The printed age formula ends in "+1" where Cioffi et al. (1988) has "+1/3". Ages below use V_s = 119.3 (`verify_paper_energetics.json`); V_s = 119 changes them by ≤ 0.01×10⁴ yr.
     - With "+1/3" the age is 7.95×10⁴ yr.
     - As printed, it is 8.72×10⁴ yr.
     - The closest value to the quoted 6.8×10⁴ yr is 6.74–6.75×10⁴ yr, which requires both "+1/3" and the rounded R_PDS = 18 pc.
   - Every input is invalid or unsupported:
     - R_s = 30 pc follows from the 160 pc unit error (Section 1);
     - n₀ = 1.0 ± 0.3 is adopted, not measured;
     - V_s comes from the WIM [S II] density, and that density is itself wrong (item 2).
   - With the paper's n₀ = 1 and V_s = 119 km s⁻¹, E = 1.2–1.7×10⁵² erg at ~280 pc (R_s = 52.5–59.1 pc) and 2.9–4.2×10⁵⁴ erg at 1.6 kpc.
   - **Balmer bound.** A *complete* radiative shock yields at least 0.117 Hβ photons per swept-up H atom. The Hβ non-detection (≤ 0.057 R from SDSS at 2σ, up to a conservative 0.12 R from an older DESI estimate), with limb brightening of 3–18.6, then bounds n₀ ≤ 0.002–0.034 cm⁻³ and E ≤ 1.4×10⁴⁹–2.9×10⁵⁰ erg at 280 pc. These are model scalings, from `verify_paper_energetics.json`. The bound does not apply to an incomplete shock; our column-truncated shock models give n₀ ≈ 0.012–0.042 cm⁻³ (model only).
6. **The [O III] width is not established, and λ4959 is below the paper's own threshold.**
   - We simulated the paper's procedure at S/N 6.1: Gaussian fits, an instrumental FWHM of 40 ± 1.5 km s⁻¹ (R ≈ 7500) and the cut FWHM_obs − err > FWHM_inst + err (Sect. 3.2.2). Our fits use a linear continuum where the paper used a quadratic one.
   - Fits that pass the cut return medians of 35.8, 36.6, 40.0 and 46.9 km s⁻¹ for true widths of 5.4, 20, 30 and 45 km s⁻¹. The pass fractions are 10, 24, 48 and 82 % (Figure C6d).
   - At an integrated S/N of 6.1, a 43–49 km s⁻¹ result therefore cannot distinguish an unresolved line from a 45 km s⁻¹ one.
   - The paper does not define its S/N. If it is a peak S/N, the integrated S/N is higher and the bias smaller: at S/N 15 the medians are 22.8–44.7 km s⁻¹, and 43–49 km s⁻¹ would then point to an intrinsic width of ≈ 45 km s⁻¹. DESI fibre spectra suggest ≤ 33 km s⁻¹, but that width has not been re-verified (`notes/w3_canon.md` §12, C7).
   - Table 2 gives λ4959 S/N = 1.8–2.3, below the "S/N > 3" usability criterion of Sect. 3.2.1 (p. 8). That is the value expected from λ5007 at S/N 6.1–6.6 and a ratio of 2.98 (Storey & Zeippen 2000), namely 2.05–2.21.
   - The resulting 5007/4959 = 2.29 ± 1.33, 2.44 ± 1.28 and 2.53 ± 1.17 carry no information.
   - λ4959 also lies only 9 Å inside the 4950 Å blue-arm edge (Sect. 2.2).
   - Sect. 3.2.1 itself calls λ5007 "relatively weak" (p. 8), and Sect. 5 describes "moderate [O iii] λλ4959,5007 emissions" (p. 17). "Strong" is therefore an overstatement.
7. **The NVSS/TGSS non-detection carries no information about a filament this wide.**
   - NVSS images are mosaics of VLA snapshots in the compact D/DnC configurations (Condon et al. 1998; Condon 1998). The D-configuration largest angular scale at L band is 970″ in full synthesis, "divided by two" for snapshots (NRAO 2025, VLA OSS), i.e. ≈ 8′.
   - TGSS ADR1's "sensitivity to extended emission will … not go beyond angular scales of 10′–20′" (Intema et al. 2017).
   - The strand FWHM is ≈ 17′ and the arc is 17° long. Both surveys resolve out smooth emission on these scales, so a non-detection is expected whether or not the filament emits.
   - A meaningful test needs single-dish data. Our beam-diluted 408 MHz limit (Haslam et al. 1982 survey) is Σ₁GHz < 1.3×10⁻²² W m⁻² Hz⁻¹ sr⁻¹ (2σ). This is our analysis and has not been independently re-verified.
8. **Internal inconsistencies**
   - From `notes/paper_review.md` §4.12:
     - Sect. 4.2 (p. 17) cites a "total expansion speed (∼140 km s⁻¹, Section 3.4)", but no such number appears in Sect. 3.4; the only speed is V_s = 119 km s⁻¹ (Sect. 4.1).
     - T_e is given as "about 9000 K" (p. 12), ~8.7×10³ K (Fig. 7a caption) and "approximately 10,000 K" (adopted, p. 13).
     - The [S II]/Hα shock threshold is 0.4 (abstract, Sect. 5), ≳ 0.45 (Sect. 1) and 0.45–0.5 (Sect. 4.1).
     - Eq. 6 is called "third-order" but is linear.
     - Table 2's P3 doublet ratio is 1.19, whereas its intensities 0.11/0.10 give 1.10.
     - P4's [N II]/Hα = 0.95 matches neither the line sum (1.12) nor λ6583 alone (0.80). The P1–P3 values follow the sum.
     - The age formula has "+1", and R_PDS = 18 vs 17.1 pc (item 5).
     - The solar O/H of 8.66 is attributed to Allende Prieto et al. (2001), who give 8.69 ± 0.05. The value 8.66 is Asplund et al. (2004)'s; Pettini & Pagel (2004) cite the two together. This is immaterial (0.03 dex).
   - From `notes/geometry.md` §4.1:
     - Fig. 4 plots 173 [O III] fibres, against 279 in the text (p. 6).
     - The Fig. 1 caption's "field of view of 1°×10°" describes the nebula. The image itself spans 3.9°×13.5°.
     - The Table 1 positions of P1 and P4 lie 0.48° and 0.42° from the outlines drawn in Fig. 1b.
     - The background images in Figs. 2, 4 and 5 are offset from the plotted axes by 3.1′ (median), up to 10.4′.

![Figure 6: four panels testing the paper's red-arm diagnostics, [S II] density, Eq. 2 temperature and [O III] width](figures/dispute_other_fig6.png)

*Figure C6. LAMOST quantities are from Wu et al. (2026) Table 2, with P4 the nebula-free reference.*
- *(a) Red-arm ratios recomputed from the Table 2 intensities, with 1σ errors propagated from the quoted ±0.02–0.03.*
  - *Shading: Madsen et al. (2006) Table 5 high-latitude WIM filaments. For the [S II] doublet sum, the λ6716 range is scaled by 1 + 1/1.454 (low-density limit; derived).*
  - *Dashed: Madsen et al.'s faint-WIM lower bounds (≳ 0.4; ≳ 0.7).*
  - *Dotted: their O-star H II means.*
  - *Red: the paper's SNR criteria.*
- *(b) I(6716)/I(6731) against n_e for our 5-level S⁺ atom (four atomic-data sets, 8000–10⁴ K) and for the paper's Eq. 5. Points mark the Table 2 ratios on each curve. Table 2 gives no error on R; the note gives σ_R from its [S II] S/N (44–65) and from its ±0.02 intensities. The red band is the paper's 64–150 cm⁻³. The grey band is the WIM doublet ratio measured on this sightline in SDSS sky fibres (16–84 %).*
- *(c) T_e from Eq. 2 applied to the Table 2 Hα and [N II] FWHM. Error bars are the 16–84 % range of a Monte Carlo on the quoted width errors.*
- *(d) Monte Carlo of the paper's [O III] fitting and FWHM cut (`paper_checks.json`): median and 16–84 % of the recovered intrinsic FWHM among fits that pass the cut, at integrated S/N 6.1, 15 and 50.*

*Script: `scripts/dispute_other_fig6.py`.*

**Table C7.** Claims in Sects. 3.4, 3.5 and 4.1 and Tables 2 and 4, against values recomputed from the paper's own inputs. Sources: `dispute_other_values.json`, `verify_paper_sii.json`, `verify_paper_energetics.json`, `paper_checks.json`.

| Claim (location) | Paper value | Reproduced / implied value | Verdict |
|---|---|---|---|
| [S II] n_e (Sect. 3.4, Fig. 7b, pp. 13–14) | median 64, up to 150 cm⁻³; background lower | Table 2 R = 1.03–1.19 (σ_R 0.03 or 0.2–0.3) → **169–400** (Eq. 5), **225–631** (5-level); P4 highest | Not supported by Table 2; unphysical for WIM (EM, sightline LDL) |
| R for 64 / 150 cm⁻³; low-density limit (p. 14) | LDL ~1.42 | 64 → R 1.33 (Eq. 5) / 1.37 (5-level); 150 → 1.23 / 1.27–1.28; LDL 1.437–1.462 | Not reached by any region |
| T_e, Eq. 2 (Sect. 3.4, p. 12) | ≈ 9000 K; 1–2×10⁴ K on filament vs 7.5×10³ K in P4 | 6780 / 6580 / 3220 / 7660 K; σ_T 7300–10 700 K; P(T<0) 18–38 % | Unconstrained; no on/off difference |
| V_t, Eq. 3 (Table 4, p. 13) | 33.1 / 32.9 / 30.6 / 29.9 km s⁻¹ | 26.7 / 25.7 / 29.9 / 26.5 km s⁻¹ | Not reproduced from Table 2 widths; P4 equal (WIM) |
| 12+log(O/H), N2, Eq. 6 (Sect. 3.5, p. 14) | 8.82, "third-order polynomial" | Linear PP04: 8.81–8.85 (cubic 9.09–9.18); N2 = −0.15…−0.10, outside calibration; P4 8.85 | Arithmetic OK; calibration invalid; WIM |
| 12+log(O/H), N2S2Hα, Eqs. 7–8 | 8.70 | 8.64–8.75; P4 8.75 | Arithmetic OK; describes WIM |
| V_s (Dopita 1979), Sect. 4.1 | 119 ± 5 km s⁻¹ | 119.3 km s⁻¹ (n_e 64, n₀ 1), but n_e is wrong and is WIM | Input invalid |
| E (Boumis et al. 2022), Sect. 4.1 | (2.0 ± 0.7)×10⁵¹ erg | 2.01×10⁵¹ (160 pc); 1.2–1.7×10⁵² (280 pc, n₀ 1); 2.9–4.2×10⁵⁴ (1.6 kpc); complete-shock Balmer bound ≤ 1.4×10⁴⁹–2.9×10⁵⁰ (280 pc) | Arithmetic OK; inputs invalid |
| R_PDS (Cioffi et al. 1988) | 18 ± 2.8 pc | 17.1 pc (18 needs E = 2.41×10⁵¹) | Minor error |
| Age (Cioffi et al. 1988) | (6.8 ± 1.8)×10⁴ yr | 7.95×10⁴ ("+1/3"); 8.72×10⁴ (printed "+1"); 6.74–6.75×10⁴ only with "+1/3" and R_PDS = 18 | Formula misprinted; inputs invalid |
| Intrinsic [O III] FWHM (Table 2) | 46.9 / 48.9 / 43.1 km s⁻¹ | Passing fits return 36–47 km s⁻¹ for true 5–45 km s⁻¹ (integrated S/N 6.1); 23–45 km s⁻¹ at S/N 15 | Not established; S/N undefined |
| [O III] λ4959 (Table 2; Sect. 4.1) | "strong λλ4959,5007"; S/N 1.8–2.3 | Expected S/N 2.05–2.21; 5007/4959 = 2.3–2.5 ± 1.2–1.3 | Below the paper's S/N > 3 cut; uninformative |

**Consequence.**
- In the LAMOST data, Hα, [N II] and [S II] are dominated by the foreground and background WIM: their ratios are the same in P4, and independent Hα data show no filament. The following quantities therefore refer to the WIM, not to the filament:
  - the BPT/SMB classification;
  - n_e;
  - T_e and V_t;
  - O/H.
- The only filament-specific LAMOST observable is [O III] λ5007, at S/N ≈ 6. Its width is not established.
- The shock velocity, explosion energy, PDS radius and age have no observational basis. With the paper's own n₀ = 1, E rises to 10⁵² erg at 280 pc and 10⁵⁴ erg at 1.6 kpc. For a complete shock, the Balmer limit allows at most ~3×10⁵⁰ erg at 280 pc.
- The radio statement does not support an old SNR.
- The paper's own fallback in Sect. 4.2 ("DIG or ionized bubble cannot be excluded") is correct for the red lines. A "possible SNR" interpretation must rest on other evidence, such as [O III]/Hβ or [O III]/Hα, not on the evidence presented in Sects. 3.4–4.1.

**Requested correction.**
1. **Shock criterion.** Delete "confirming that the ionized gas is shock-excited" (Sect. 4.1), and the threshold claim in the abstract and Sect. 5. Replace them with: *"The red-arm ratios are indistinguishable in the nebula-free region P4 and are typical of the faint WIM; they do not constrain the excitation of the [O III] filament."* Make the 0.4/0.45/0.5 threshold consistent throughout.
2. **[S II] density.** Give the uncertainty of each Table 2 doublet ratio. Replace the density paragraph with the Table 2-implied values, or with a statement that the doublet ratios are not physically interpretable. Withdraw the [S II] density as a "post-shock" density and as the density of "the [O iii]-emitting gas" (p. 15). Correct the low-density limit to ≈ 1.45, and correct the background-density statement.
3. **T_e and V_t.** Report T_e per region with its uncertainty and delete "1–2×10⁴ K". Recompute Table 4's V_t from Eq. 3, or state how it was derived.
4. **Abundances.** Call Eq. 6 the linear PP04 relation. State that N2 lies outside the calibrated range and that the abundances characterize the line-of-sight WIM.
5. **Energetics.** Correct the age bracket to "+1/3" and R_PDS to 17.1 pc. Either remove E and the age, or present them as illustrative scalings with n₀ and the distance stated as unknown.
6. **[O III] widths.** Define the S/N in Table 2. Then either show by simulation that the widths survive the FWHM cut at that S/N, or replace them with "not constrained at this S/N". Flag λ4959 as below the S/N threshold, and remove "strong".
7. **Radio.** Replace the radio sentence with *"Interferometric surveys are insensitive to smooth structures larger than ≈ 8′ (NVSS) or 10′–20′ (TGSS); they provide no constraint on smooth emission from a 17′-wide filament."*
8. **Inconsistencies.** Correct the items listed in point 8.

---

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

---

## References

This list covers every work cited in this Comment. The full verification log is in `notes/dispute/references.md`. It uses AAS/AJ order (authors, year, journal, volume, first page or article ID, DOI, arXiv ID) and keeps article titles so that each entry can be checked. Works with more than three authors are given as first author et al.

Each entry was checked on 2026-11-17 against Crossref (by DOI), the arXiv API, DataCite, ADS, or the publisher or data-provider page. Records are in `data/derived/dispute_references_records.json` (`scripts/dispute_references_verify.py`) and `data/derived/dispute_references_manual.json` (`scripts/dispute_references_manual.py`). The citation cross-check is `scripts/dispute_references_crosscheck.py` → `data/derived/dispute_references_crosscheck.json`.

- Abdurro'uf, et al. 2022, The Seventeenth Data Release of the Sloan Digital Sky Surveys: Complete Release of MaNGA, MaStar, and APOGEE-2 Data, ApJS, 259, 35, doi:10.3847/1538-4365/ac4414
- Allende Prieto, C., Lambert, D. L., & Asplund, M. 2001, The "Forbidden" Abundance of Oxygen in the Sun, ApJL, 556, L63, doi:10.1086/322874, arXiv:astro-ph/0106360
- Almeida, A., et al. 2023, The Eighteenth Data Release of the Sloan Digital Sky Surveys: Targeting and First Spectra from SDSS-V, ApJS, 267, 44, doi:10.3847/1538-4365/acda98
- Anderson, L. D., et al. 2012, The Green Bank Telescope H II Region Discovery Survey. III. Kinematic Distances, ApJ, 754, 62, doi:10.1088/0004-637X/754/1/62, arXiv:1205.4228
- Asplund, M., et al. 2004, Line formation in solar granulation. IV. [O I], O I and OH lines and the photospheric O abundance, A&A, 417, 751, doi:10.1051/0004-6361:20034328, arXiv:astro-ph/0312290
- Astropy Collaboration, et al. 2022, The Astropy Project: Sustaining and Growing a Community-oriented Open-source Project and the Latest Major Release (v5.0) of the Core Package, ApJ, 935, 167, doi:10.3847/1538-4357/ac7c74
- Blanton, M. R., & Roweis, S. 2007, K-Corrections and Filter Transformations in the Ultraviolet, Optical, and Near-Infrared, AJ, 133, 734, doi:10.1086/510127, arXiv:astro-ph/0606170
- Boumis, P., et al. 2022, Deep optical study of the mixed-morphology supernova remnant G 132.7+1.3 (HB3), MNRAS, 512, 1658, doi:10.1093/mnras/stac412
- Bovy, J. 2017, Galactic rotation in Gaia DR1, MNRAS, 468, L63, doi:10.1093/mnrasl/slx027, arXiv:1610.07610
- Cioffi, D. F., McKee, C. F., & Bertschinger, E. 1988, Dynamics of radiative supernova remnants, ApJ, 334, 252, doi:10.1086/166834
- Condon, J. J. 1998, Radio Surveys, in IAU Symp. 179, 19, doi:10.1017/S0074180900128177
- Condon, J. J., et al. 1998, The NRAO VLA Sky Survey, AJ, 115, 1693, doi:10.1086/300337
- Cox, D. P., & Raymond, J. C. 1985, Preionization-dependent families of radiative shock waves, ApJ, 298, 651, doi:10.1086/163649
- Cui, X.-Q., et al. 2012, The Large Sky Area Multi-Object Fiber Spectroscopic Telescope (LAMOST), RAA, 12, 1197, doi:10.1088/1674-4527/12/9/003
- DESI Collaboration, et al. 2026, Data Release 1 of the Dark Energy Spectroscopic Instrument, AJ, 171, 285, doi:10.3847/1538-3881/ae4c43, arXiv:2503.14745
- Dopita, M. A. 1979, Optical emission from shocks. V. The Magellanic Cloud supernova remnants, ApJS, 40, 455, doi:10.1086/190594
- Dopita, M. A., & Sutherland, R. S. 1996, Spectral Signatures of Fast Shocks. I. Low-Density Model Grid, ApJS, 102, 161, doi:10.1086/192255
- Dopita, M. A., et al. 2016, Chemical abundances in high-redshift galaxies: a powerful new emission line diagnostic, Ap&SS, 361, 61, doi:10.1007/s10509-016-2657-8, arXiv:1601.01337
- Fesen, R. A., et al. 2024, Deep Optical Emission-line Images of Nine Known and Three New Galactic Supernova Remnants, ApJS, 272, 36, doi:10.3847/1538-4365/ad410a
- Fesen, R. A., et al. 2026, Optical Discovery of New and Candidate Galactic Supernova Remnants Plus Optical Imaging of the Monogem Supernova Remnant, arXiv:2608.26666 (v2, 2026 September 7; no journal version found)
- Gaensler, B. M., et al. 2008, The Vertical Structure of Warm Ionised Gas in the Milky Way, PASA, 25, 184, doi:10.1071/AS08004, arXiv:0808.2550
- Gaia Collaboration, et al. 2023, Gaia Data Release 3: Summary of the content and survey properties, A&A, 674, A1, doi:10.1051/0004-6361/202243940
- Guy, J., et al. 2023, The Spectroscopic Data Processing Pipeline for the Dark Energy Spectroscopic Instrument, AJ, 165, 144, doi:10.3847/1538-3881/acb212
- Haffner, L. M., Reynolds, R. J., & Tufte, S. L. 1999, WHAM Observations of Hα, [S II], and [N II] toward the Orion and Perseus Arms: Probing the Physical Conditions of the Warm Ionized Medium, ApJ, 523, 223, doi:10.1086/307734, arXiv:astro-ph/9904143
- Haffner, L. M., et al. 2003, The Wisconsin Hα Mapper Northern Sky Survey, ApJS, 149, 405, doi:10.1086/378850, arXiv:astro-ph/0309117
- Haffner, L. M., et al. 2009, The warm ionized medium in spiral galaxies, RvMP, 81, 969, doi:10.1103/RevModPhys.81.969, arXiv:0901.0941
- Hanuschik, R. W. 2003, A flux-calibrated, high-resolution atlas of optical sky emission from UVES, A&A, 407, 1157, doi:10.1051/0004-6361:20030885
- Haslam, C. G. T., et al. 1982, A 408 MHz all-sky continuum survey. II. The atlas of contour maps, A&AS, 47, 1 (ADS 1982A&AS...47....1H; no DOI)
- Intema, H. T., et al. 2017, The GMRT 150 MHz all-sky radio survey. First alternative data release TGSS ADR1, A&A, 598, A78, doi:10.1051/0004-6361/201628536
- Jordi, C., et al. 2010, Gaia broad band photometry, A&A, 523, A48, doi:10.1051/0004-6361/201015441
- Katz, D., et al. 2023, Gaia Data Release 3: Properties and validation of the radial velocities, A&A, 674, A5, doi:10.1051/0004-6361/202244220
- Knies, J. R., et al. 2024, A new understanding of the Gemini-Monoceros X-ray enhancement from discoveries with eROSITA, A&A, 688, A90, doi:10.1051/0004-6361/202348834, arXiv:2401.17289
- Kopsacheili, M., Zezas, A., & Leonidaki, I. 2020, A diagnostic tool for the identification of supernova remnants, MNRAS, 491, 889, doi:10.1093/mnras/stz2594
- Kramida, A., Ralchenko, Yu., Reader, J., & NIST ASD Team 2024, NIST Atomic Spectra Database (ver. 5.12), NIST Standard Reference Database 78 (Gaithersburg, MD: NIST), https://physics.nist.gov/asd, doi:10.18434/T4W30F (accessed 2026 November 16–17)
- Kurucz, R. L. 2005, New atlases for solar flux, irradiance, central intensity, and limb intensity, MSAIS, 8, 189 (ADS 2005MSAIS...8..189K; no DOI; atlas file at http://kurucz.harvard.edu/sun/fluxatlas2005/)
- Luo, A.-L., et al. 2015, The first data release (DR1) of the LAMOST regular survey, RAA, 15, 1095, doi:10.1088/1674-4527/15/8/002
- Madsen, G. J., Reynolds, R. J., & Haffner, L. M. 2006, A Multiwavelength Optical Emission Line Survey of Warm Ionized Gas in the Galaxy, ApJ, 652, 401, doi:10.1086/508441, arXiv:astro-ph/0609558
- Marino, R. A., et al. 2013, The O3N2 and N2 abundance indicators revisited: improved calibrations based on CALIFA and T_e-based literature data, A&A, 559, A114, doi:10.1051/0004-6361/201321956, arXiv:1307.5316
- Martin, D. C., et al. 2005, The Galaxy Evolution Explorer: A Space Ultraviolet Survey Mission, ApJL, 619, L1, doi:10.1086/426387
- Mendoza, C., & Zeippen, C. J. 1982, Transition probabilities for forbidden lines in the 3p³ configuration, MNRAS, 198, 127, doi:10.1093/mnras/198.1.127
- Morrissey, P., et al. 2007, The Calibration and Data Products of GALEX, ApJS, 173, 682, doi:10.1086/520512
- Morton, D. C. 2000, Atomic Data for Resonance Absorption Lines. II. Wavelengths Longward of the Lyman Limit for Heavy Elements, ApJS, 130, 403, doi:10.1086/317349
- NRAO 2025, VLA Observational Status Summary 2025B: Resolution (Table 3.1.1, note 4), https://science.nrao.edu/facilities/vla/docs/manuals/oss2025b/performance/resolution (accessed 2026 November 17)
- Oke, J. B., & Gunn, J. E. 1983, Secondary standard stars for absolute spectrophotometry, ApJ, 266, 713, doi:10.1086/160817
- Persic, M., Salucci, P., & Stel, F. 1996, The universal rotation curve of spiral galaxies — I. The dark matter connection, MNRAS, 281, 27, doi:10.1093/mnras/278.1.27 (registered DOI; the "278" in its suffix is correct as registered)
- Pettini, M., & Pagel, B. E. J. 2004, [O III]/[N II] as an abundance indicator at high redshift, MNRAS, 348, L59, doi:10.1111/j.1365-2966.2004.07591.x, arXiv:astro-ph/0401128
- Proxauf, B., Öttl, S., & Kimeswenger, S. 2014, Upgrading electron temperature and electron density diagnostic diagrams of forbidden line emission, A&A, 561, A10, doi:10.1051/0004-6361/201322581
- QHYCCD 2024, Silk Nebula Discovered by Enthusiast Astronomers with QHYCCD Cameras, web note (2024 July 12), https://www.qhyccd.com/silk-nebula-discovered-by-enthusiast-astronomers-with-qhyccd-cameras/ (accessed 2026 November 16; not peer reviewed)
- Ramsbottom, C. A., Bell, K. L., & Stafford, R. P. 1996, Effective Collision Strengths for Electron Impact Excitation of Singly Ionized Sulfur, ADNDT, 63, 57, doi:10.1006/adnd.1996.0009
- Raymond, J. C. 1979, Shock waves in the interstellar medium, ApJS, 39, 1, doi:10.1086/190562
- Raymond, J. C., et al. 1988, Spatial and spectral interpretation of a bright filament in the Cygnus Loop, ApJ, 324, 869, doi:10.1086/165945
- Redfield, S., & Linsky, J. L. 2008, The Structure of the Local Interstellar Medium. IV. Dynamics, Morphology, Physical Properties, and Implications of Cloud-Cloud Interactions, ApJ, 673, 283, doi:10.1086/524002, arXiv:0709.4480
- Reid, M. J., et al. 2014, Trigonometric Parallaxes of High Mass Star Forming Regions: The Structure and Kinematics of the Milky Way, ApJ, 783, 130, doi:10.1088/0004-637X/783/2/130, arXiv:1401.5377
- Reid, M. J., et al. 2019, Trigonometric Parallaxes of High-mass Star-forming Regions: Our View of the Milky Way, ApJ, 885, 131, doi:10.3847/1538-4357/ab4a11, arXiv:1910.03357
- Ren, J.-J., et al. 2021, On the radial velocity calibrations in the LAMOST medium-resolution spectroscopic survey of nebulae, RAA, 21, 051, doi:10.1088/1674-4527/21/3/51, arXiv:2008.06236
- Reynolds, R. J., Roesler, F. L., & Scherb, F. 1977, Faint optical emission lines from the interstellar medium: limits on the temperature and ionization state of the emitting gas, ApJ, 211, 115, doi:10.1086/154908
- Rynkun, P., Gaigalas, G., & Jönsson, P. 2019, Theoretical investigation of energy levels and transition data for S II, Cl III, Ar IV, A&A, 623, A155, doi:10.1051/0004-6361/201834931
- Sabin, L., et al. 2013, New Galactic supernova remnants discovered with IPHAS, MNRAS, 431, 279, doi:10.1093/mnras/stt160
- Schlegel, D. J., Finkbeiner, D. P., & Davis, M. 1998, Maps of Dust Infrared Emission for Use in Estimation of Reddening and Cosmic Microwave Background Radiation Foregrounds, ApJ, 500, 525, doi:10.1086/305772
- Schönrich, R., Binney, J., & Dehnen, W. 2010, Local kinematics and the local standard of rest, MNRAS, 403, 1829, doi:10.1111/j.1365-2966.2010.16253.x
- SDSS Collaboration, et al. 2026, The Nineteenth Data Release of the Sloan Digital Sky Survey, ApJS, 285, 9, doi:10.3847/1538-4365/ae4697, arXiv:2507.07093
- Stark, A. A., & Brand, J. 1989, Kinematics of molecular clouds. II. New data on nearby giant molecular clouds, ApJ, 339, 763, doi:10.1086/167334
- Storey, P. J., & Zeippen, C. J. 2000, Theoretical values for the [O III] 5007/4959 line-intensity ratio and homologous cases, MNRAS, 312, 813, doi:10.1046/j.1365-8711.2000.03184.x
- Tayal, S. S., & Zatsarinny, O. 2010, Breit–Pauli Transition Probabilities and Electron Excitation Collision Strengths for Singly Ionized Sulfur, ApJS, 188, 32, doi:10.1088/0067-0049/188/1/32
- Wenger, T. V. 2018, kd v1.0 (kinematic distance software), Zenodo, doi:10.5281/zenodo.1166001 (code read at https://github.com/tvwenger/kd, commit 13a8e3f, 2024 September 12)
- Wenger, T. V., et al. 2018, Kinematic Distances: A Monte Carlo Method, ApJ, 856, 52, doi:10.3847/1538-4357/aaaec8, arXiv:1802.04203
- Wu, C.-J., et al. 2022, The Data Processing of the LAMOST Medium-resolution Spectral Survey of Galactic Nebulae (LAMOST MRS-N Pipeline), RAA, 22, 075015, doi:10.1088/1674-4527/ac7387, arXiv:2110.04547
- Wu, Y.-F., et al. 2026, The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ, 172, 298, doi:10.3847/1538-3881/aeb11c, arXiv:2610.13329
- Zhang, L.-Y., et al. 2025, LAMOST Medium-resolution Spectroscopic Survey of the Rosette Nebula, AJ, 170, 309, doi:10.3847/1538-3881/ae0f94, arXiv:2510.08950
- Zhang, W., et al. 2021, LAMOST medium-resolution spectral survey of Galactic nebulae (LAMOST-MRS-N): subtraction of geocoronal Hα emission, RAA, 21, 280, doi:10.1088/1674-4527/21/11/280, arXiv:2108.08021
- Ziegenbalg, S. 2025, Northern Sky Narrowband Survey DR0.2: First [O III] and [S II] Results, RNAAS, 9, 227, doi:10.3847/2515-5172/adfec7
