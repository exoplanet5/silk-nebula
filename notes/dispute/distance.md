## 1. The 160 pc distance: Table 4 is in kpc, not "100 pc"

**Claim in the paper.**
- Kinematic distances come from `kd` (Wenger et al. 2018), with "the rotation curve ... upgraded by M. Reid et al. (2014)" (Sect. 3.3, p. 10).
- Table 4 (p. 13) gives Monte-Carlo and "traditional" distances "in units of 100pc" (footnotes 1, 2).
- The [O III] rows read 1.63/1.58/0.98 (Monte Carlo) and 1.79/1.85/1.20 (traditional). The text calls these "∼100–190 pc" (p. 10).
- The paper adopts "a fiducial distance of ∼160 pc (the approximate median of the six [O iii]-based kinematic distance estimates; see Table 4)" (Sect. 3.3, p. 11). This value is used again in the Abstract and in Sects. 4.1 and 5.

Our tests use only the paper's Tables 1 and 4, Fig. 5b and text, plus the rotation curve it names; `kd`'s documented unit is corroboration only. Pages refer to arXiv:2610.13329v1. The AJ footnotes are identical, and Crossref lists no erratum (2026-11-17). "Table 4", "Fig. 5b" and "Sect." refer to the paper; "Figure 1", "Table 2" and "Section" refer to this Comment.

**Why it is wrong.**

1. **Read as kpc, the printed values are reproduced to within 0.005–0.05 kpc.**
   - **Method.** We inverted the printed V_LSR along the Table 1 sightlines with the `kd`-form Reid+14 curve (Persic et al. 1996 functional form). V_helio-based values use the paper's solar motion (Sect. 3.2.3). `kd` itself was not run: its curve and `pdf_kd` were re-implemented from the source.
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
   - **The unit result does not depend on the curve.** All seven variants, including `kd`'s default (Reid+19 URC), give 1.74–2.15 kpc for [O III] P1–P2 and 1.15–1.43 kpc for P3, at least 9× the ×100 pc values (0.12–0.19 kpc).
   - **Only P1 pins the curve.** P1 is matched to 0.005 kpc only by the `kd`-form Reid+14 curve (others 1.91–2.15). For P2 and P3 that curve is within 0.05–0.11 kpc and its ×cos b variant within 0.02 kpc.

2. **At 160 pc, rotation cannot produce the paper's velocities.**
   - **Rotation curve.** At 160 pc, rotation gives V_LSR = **1.5–1.7 km s⁻¹** toward P1–P3 (0.0–1.7 over all seven variants). The Oort term A sin 2l cos²b × d agrees, with A = 15.3 km s⁻¹ kpc⁻¹ (Bovy 2017).
   - **Distance needed.** Within 300 pc, rotation stays ≤ 3.1 km s⁻¹ (Figure 1a, inset). Reaching 15.5 km s⁻¹ requires **1.79–2.15 kpc**.
   - **The paper's own Fig. 5b** (our Figure 1b) shows this for the whole [O III] sample:
     - 97.5 % of the fibres have V_LSR > +2.3 km s⁻¹, i.e. beyond 0.23 kpc on the P1 curve;
     - 88 % have V_LSR > +6.6 km s⁻¹, i.e. beyond 0.69 kpc.
     - Bin heights are multiples of 1/81, so N is 81 or a multiple. This is our inference; N is not printed, and the fractions do not depend on it.

3. **The Monte-Carlo widths encode the velocity gradient.**
   - **The bound.** A 68 % velocity width W_V maps to W_d ≈ W_V/(dV/dd). Curve resampling, KDE smoothing and extra error terms only widen W_d, so W_V/W_d cannot exceed the largest true dV/dd, 5.8–10.7 km s⁻¹ kpc⁻¹ along P1–P4 at 0.05–2.5 kpc (seven variants plus the Oort term). This assumes the tabulated σ_v (3–4 km s⁻¹ for [O III], Sect. 3.2.3) or larger was used.
   - **Example, P1 [O III].** 8.6/(1.01 + 0.84) = **4.65** km s⁻¹ per unit: allowed as kpc, but **46.5 km s⁻¹ kpc⁻¹** as ×100 pc.
   - **All rows** (Figure 1d, Table 2). As kpc, 0.9–4.7: all 15 allowed. As ×100 pc, 9.5–46.5: 13 of 15 exceed 10.7, [O III] by **3.5–4.4×** and Hα by only 1.3–1.7×. The truncated P1/P2 [N II] rows (9.9, 9.5) do not discriminate.
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

*Figure 1. Table 4 of Wu et al. (2026) compared with the Reid et al. (2014) curve the paper names (`kd` form), on the Table 1 sightlines. All inputs are the paper's printed numbers; no new observations are used.*
- *(a) V_LSR(d) toward P1–P3 (lines), with the range of seven curve variants in grey. [O III] (circles) and Hα (squares) traditional distances are plotted at the printed value read as kpc (filled) and as ×100 pc (open). Dotted lines mark the velocity gap. Inset: 0–300 pc.*
- *(b) [O III] V_LSR histogram digitized from the paper's Fig. 5b (300-dpi page crop; 9 equal bins from −2.0 to 36.7 km s⁻¹; percent axis as in the paper). Orange: velocities rotation can reach within 300 pc. Dotted: 160 pc. Right axis: P1 kinematic distance (`kd`-form Reid+14).*
- *(c) Rotation velocity at the printed distance against the printed V_LSR ± σ_v, for all 15 rows. Filled = kpc, open = ×100 pc. [S II]/[N II] rows are faded.*
- *(d) Implied gradient W_V/W_d.*
  - *W_d is the printed Monte-Carlo interval (+err plus −err).*
  - *W_V is the 68.3 % width of N(V_LSR, σ_v), truncated at V(0.01 kpc), below which `kd` has no solution. For truncated rows this approximates `kd`'s KDE interval.*
  - *Green: the true gradient, 5.8–10.7 km s⁻¹ kpc⁻¹ (curve variants and Oort term, 0.05–2.5 kpc).*
  - *Red: forbidden, i.e. narrower than σ_v alone allows.*
- *Code: `scripts/dispute_distance_kd.py`. Values: `data/derived/dispute_distance_kd.json` (emulation: `verify_paper_kd.json`).*

**Table 2.** Row-by-row recomputation of Wu et al. (2026) Table 4.
- **Recomputed d:** traditional distance from the printed V_LSR (`kd`-form Reid+14). Values in parentheses ([O III] and Hα rows) use V_LSR rebuilt from the printed V_helio with the paper's solar motion.
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

**References used in this section.**
- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (2026), The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ 172, 298, doi:10.3847/1538-3881/aeb11c; arXiv:2610.13329 [verified via Crossref and arXiv; no update/erratum records on 2026-11-17]
- Wenger, T. V., Balser, D. S., Anderson, L. D., & Bania, T. M. (2018), Kinematic Distances: A Monte Carlo Method, ApJ 856, 52, doi:10.3847/1538-4357/aaaec8 [verified via Crossref]
- Wenger, T. V. (2018), `kd` v1.0: kinematic distance software, Zenodo, doi:10.5281/zenodo.1166001 (the record the `kd` README asks users to cite); code read at https://github.com/tvwenger/kd, commit 13a8e3f (2024-09-12) [verified via DataCite (Zenodo record) and the GitHub source, read 2026-11-16; project copy `data/derived/verify_paper_raw/kd_*`]
- Reid, M. J., Menten, K. M., Brunthaler, A., Zheng, X. W., et al. (2014), Trigonometric Parallaxes of High Mass Star Forming Regions: The Structure and Kinematics of the Milky Way, ApJ 783, 130, doi:10.1088/0004-637X/783/2/130 [verified via Crossref]
- Reid, M. J., Menten, K. M., Brunthaler, A., Zheng, X. W., et al. (2019), Trigonometric Parallaxes of High-mass Star-forming Regions: Our View of the Milky Way, ApJ 885, 131, doi:10.3847/1538-4357/ab4a11 [verified via Crossref]
- Persic, M., Salucci, P., & Stel, F. (1996), The universal rotation curve of spiral galaxies — I. The dark matter connection, MNRAS 281, 27, doi:10.1093/mnras/278.1.27 [verified via Crossref; the registered DOI suffix does read 278.1.27]
- Bovy, J. (2017), Galactic rotation in Gaia DR1, MNRAS 468, L63, doi:10.1093/mnrasl/slx027; arXiv:1610.07610 [verified via Crossref and arXiv; A = 15.3 ± 0.4 km s⁻¹ kpc⁻¹]
- Gaensler, B. M., Madsen, G. J., Chatterjee, S., & Mao, S. A. (2008), The Vertical Structure of Warm Ionised Gas in the Milky Way, PASA 25, 184, doi:10.1071/AS08004; arXiv:0808.2550 [verified via Crossref and arXiv; scale height 1830 (+120/−250) pc]
