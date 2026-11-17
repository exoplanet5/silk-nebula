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
   - The on-minus-P4 differences are +0.04 to +0.58σ for [S II] and −0.72 to 0.00σ for [N II] (Fig. 6a).
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
   - A 5-level S⁺ atom gives 243–631 cm⁻³ at 8000–10⁴ K, or 225–631 cm⁻³ including the Ramsbottom et al. (1996) Ω (Fig. 6b). It uses Tayal & Zatsarinny (2010) collision strengths with TZ10, Mendoza & Zeippen (1982) or Rynkun et al. (2019) A-values.
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
   - The 1σ errors are σ_T = 7300–10 700 K, and P(T < 0) = 18–38 % (Monte Carlo on the quoted width errors; Fig. 6c).
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
   - Fits that pass the cut return medians of 35.8, 36.6, 40.0 and 46.9 km s⁻¹ for true widths of 5.4, 20, 30 and 45 km s⁻¹. The pass fractions are 10, 24, 48 and 82 % (Fig. 6d).
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

*Figure 6. LAMOST quantities are from Wu et al. (2026) Table 2, with P4 the nebula-free reference.*
- *(a) Red-arm ratios recomputed from the Table 2 intensities, with 1σ errors propagated from the quoted ±0.02–0.03.*
  - *Shading: Madsen et al. (2006) Table 5 high-latitude WIM filaments. For the [S II] doublet sum, the λ6716 range is scaled by 1 + 1/1.454 (low-density limit; derived).*
  - *Dashed: Madsen et al.'s faint-WIM lower bounds (≳ 0.4; ≳ 0.7).*
  - *Dotted: their O-star H II means.*
  - *Red: the paper's SNR criteria.*
- *(b) I(6716)/I(6731) against n_e for our 5-level S⁺ atom (four atomic-data sets, 8000–10⁴ K) and for the paper's Eq. 5. Points mark the Table 2 ratios on each curve. Table 2 gives no error on R; the note gives σ_R from its [S II] S/N (44–65) and from its ±0.02 intensities. The red band is the paper's 64–150 cm⁻³. The grey band is the WIM doublet ratio measured on this sightline in SDSS sky fibres (16–84 %).*
- *(c) T_e from Eq. 2 applied to the Table 2 Hα and [N II] FWHM. Error bars are the 16–84 % range of a Monte Carlo on the quoted width errors.*
- *(d) Monte Carlo of the paper's [O III] fitting and FWHM cut (`paper_checks.json`): median and 16–84 % of the recovered intrinsic FWHM among fits that pass the cut, at integrated S/N 6.1, 15 and 50.*

*Script: `scripts/dispute_other_fig6.py`.*

**Table 7.** Claims in Sects. 3.4, 3.5 and 4.1 and Tables 2 and 4, against values recomputed from the paper's own inputs. Sources: `dispute_other_values.json`, `verify_paper_sii.json`, `verify_paper_energetics.json`, `paper_checks.json`.

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

**References used in this section.**
- Allende Prieto, C., Lambert, D. L., & Asplund, M. (2001), The 'Forbidden' Abundance of Oxygen in the Sun, ApJ 556, L63–L66, doi:10.1086/322874, arXiv:astro-ph/0106360 [verified via Crossref; 8.69 ± 0.05 checked in arXiv abstract]
- Asplund, M., Grevesse, N., Sauval, A. J., Allende Prieto, C., & Kiselman, D. (2004), Line formation in solar granulation. IV. [O I], O I and OH lines and the photospheric O abundance, A&A 417, 751–768, doi:10.1051/0004-6361:20034328, arXiv:astro-ph/0312290 [verified via Crossref; 8.66 ± 0.05 checked in arXiv abstract]
- Boumis, P., Chiotellis, A., Fragkou, V., Akras, S., et al. (2022), Deep optical study of the mixed-morphology supernova remnant G 132.7+1.3 (HB3), MNRAS 512, 1658–1676, doi:10.1093/mnras/stac412 [verified via Crossref]
- Cioffi, D. F., McKee, C. F., & Bertschinger, E. (1988), Dynamics of radiative supernova remnants, ApJ 334, 252, doi:10.1086/166834 [verified via Crossref]
- Condon, J. J. (1998), Radio Surveys, IAU Symposium 179, 19–25, doi:10.1017/S0074180900128177 [verified via Crossref; "217,446 snapshot observations in the compact D and DnC configurations" checked in the publisher abstract]
- Condon, J. J., Cotton, W. D., Greisen, E. W., Yin, Q. F., et al. (1998), The NRAO VLA Sky Survey, AJ 115, 1693–1716, doi:10.1086/300337 [verified via Crossref; snapshot construction checked at https://www.cv.nrao.edu/nvss/]
- Dopita, M. A. (1979), Optical emission from shocks. V – The Magellanic Cloud supernova remnants, ApJS 40, 455, doi:10.1086/190594 [verified via Crossref]
- Dopita, M. A., Kewley, L. J., Sutherland, R. S., & Nicholls, D. C. (2016), Chemical abundances in high-redshift galaxies: a powerful new emission line diagnostic, Ap&SS 361, 61, doi:10.1007/s10509-016-2657-8, arXiv:1601.01337 [verified via Crossref and arXiv]
- Haffner, L. M., Reynolds, R. J., & Tufte, S. L. (1999), WHAM Observations of Hα, [S II], and [N II] toward the Orion and Perseus Arms: Probing the Physical Conditions of the Warm Ionized Medium, ApJ 523, 223–233, doi:10.1086/307734, arXiv:astro-ph/9904143 [verified via Crossref; trends, 6000–9000 K and 1.0 kpc scale height checked in arXiv abstract]
- Haffner, L. M., Reynolds, R. J., Tufte, S. L., Madsen, G. J., et al. (2003), The Wisconsin Hα Mapper Northern Sky Survey, ApJS 149, 405–422, doi:10.1086/378850, arXiv:astro-ph/0309117 [verified via Crossref; 1 R = 2.25 pc cm⁻⁶ at 8000 K checked in arXiv full text]
- Haffner, L. M., Dettmar, R.-J., Beckman, J. E., Wood, K., et al. (2009), The warm ionized medium in spiral galaxies, Rev. Mod. Phys. 81, 969–997, doi:10.1103/RevModPhys.81.969, arXiv:0901.0941 [verified via Crossref; cited as a review only]
- Haslam, C. G. T., Salter, C. J., Stoffel, H., & Wilson, W. E. (1982), A 408 MHz all-sky continuum survey. II. The atlas of contour maps, A&AS 47, 1, bibcode 1982A&AS...47....1H, no DOI [verified via the ADS record and the A&A commentary that reprints the citation (Wielebinski 2009, A&A 500, 245, doi:10.1051/0004-6361/200912161); no DOI found in Crossref]
- Intema, H. T., Jagannathan, P., Mooley, K. P., & Frail, D. A. (2017), The GMRT 150 MHz all-sky radio survey. First alternative data release TGSS ADR1, A&A 598, A78, doi:10.1051/0004-6361/201628536 [verified via Crossref; quoted 10′–20′ sentence checked in A&A full-text HTML]
- Kopsacheili, M., Zezas, A., & Leonidaki, I. (2020), A diagnostic tool for the identification of supernova remnants, MNRAS 491, 889–902, doi:10.1093/mnras/stz2594 [verified via Crossref]
- Madsen, G. J., Reynolds, R. J., & Haffner, L. M. (2006), A Multiwavelength Optical Emission Line Survey of Warm Ionized Gas in the Galaxy, ApJ 652, 401–425, doi:10.1086/508441, arXiv:astro-ph/0609558 [verified via Crossref; quoted values and Table 5 checked in arXiv full text]
- Marino, R. A., Rosales-Ortega, F. F., Sánchez, S. F., Gil de Paz, A., et al. (2013), The O3N2 and N2 abundance indicators revisited: improved calibrations based on CALIFA and T_e-based literature data, A&A 559, A114, doi:10.1051/0004-6361/201321956, arXiv:1307.5316 [verified via Crossref; −1.6 < N2 < −0.2 checked in arXiv full text]
- Mendoza, C., & Zeippen, C. J. (1982), Transition probabilities for forbidden lines in the 3p³ configuration, MNRAS 198, 127–139, doi:10.1093/mnras/198.1.127 [verified via Crossref]
- NRAO (2025), VLA Observational Status Summary 2025B, Performance: Resolution (Table 3.1.1, note 4), https://science.nrao.edu/facilities/vla/docs/manuals/oss2025b/performance/resolution [verified via publisher page]
- Pettini, M., & Pagel, B. E. J. (2004), [O III]/[N II] as an abundance indicator at high redshift, MNRAS 348, L59–L63, doi:10.1111/j.1365-2966.2004.07591.x, arXiv:astro-ph/0401128 [verified via Crossref; Eqs. 1–2, validity range and solar-value citation checked in arXiv full text]
- Proxauf, B., Öttl, S., & Kimeswenger, S. (2014), Upgrading electron temperature and electron density diagnostic diagrams of forbidden line emission, A&A 561, A10, doi:10.1051/0004-6361/201322581 [verified via Crossref]
- Ramsbottom, C. A., Bell, K. L., & Stafford, R. P. (1996), Effective Collision Strengths for Electron Impact Excitation of Singly Ionized Sulfur, At. Data Nucl. Data Tables 63, 57–91, doi:10.1006/adnd.1996.0009 [verified via Crossref]
- Reynolds, R. J., Roesler, F. L., & Scherb, F. (1977), Faint optical emission lines from the interstellar medium – Limits on the temperature and ionization state of the emitting gas, ApJ 211, 115, doi:10.1086/154908 [verified via Crossref]
- Rynkun, P., Gaigalas, G., & Jönsson, P. (2019), Theoretical investigation of energy levels and transition data for S II, Cl III, Ar IV, A&A 623, A155, doi:10.1051/0004-6361/201834931 [verified via Crossref]
- Sabin, L., Parker, Q. A., Contreras, M. E., Olguín, L., et al. (2013), New Galactic supernova remnants discovered with IPHAS, MNRAS 431, 279–291, doi:10.1093/mnras/stt160 [verified via Crossref]
- Storey, P. J., & Zeippen, C. J. (2000), Theoretical values for the [O III] 5007/4959 line-intensity ratio and homologous cases, MNRAS 312, 813–816, doi:10.1046/j.1365-8711.2000.03184.x [verified via Crossref]
- Tayal, S. S., & Zatsarinny, O. (2010), Breit–Pauli transition probabilities and electron excitation collision strengths for singly ionized sulfur, ApJS 188, 32–45, doi:10.1088/0067-0049/188/1/32 [verified via Crossref]
- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (2026), The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament, AJ 172, 298, doi:10.3847/1538-3881/aeb11c, arXiv:2610.13329 [verified via Crossref]
