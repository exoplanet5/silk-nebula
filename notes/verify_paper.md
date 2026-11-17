# Wave-2 verification `paper`: distance units, [S II] density, energetics, [O III] self-consistency, mag→R conversion

*Tag `paper` · 2026-11-16 · adversarial re-check of `notes/paper_review.md` / `data/derived/paper_checks.json` using new code only.*
*Wave-1 scripts were not re-run. Their kd emulation was read only to compare results, and no bug was found in it.*

Reproduce with `~/.venvs/astro313/bin/python -I scripts/verify_paper_<x>.py`, where x is:
- `kd`: (a), ~50 s
- `sii`: (b)
- `gaia_isl`: network query
- `energetics`: (c), (d), (e)

Outputs are in `data/derived/verify_paper_{kd,sii,energetics}.json` and `figures/verify_paper_{kd,sii,oiii_conversion}.png`.

Conventions: **measured** = value printed in Wu+2026; **derived** = computed here; **literature** = external value, with its source given.
Distances are in kpc unless stated otherwise. Velocities are in km s⁻¹.

---

## 1. Results table

| # | Claim tested (wave 1) | Verdict | My value (derived) |
|---|---|---|---|
| a1 | `kd` returns distances in kpc | **CONFIRMED** | `kd` source read, not executed (HEAD 13a8e3f, 2024-09-12; re-fetched 2026-11-16T12:05Z, identical to stored copy). Every docstring gives kpc. The traditional grid starts at 0.001 kpc, `pdf_kd` uses `dist_min = 0.01` kpc, and the KDE is bounded at 0 |
| a2 | Table 4 "Distance … in units of 100 pc" is really kpc; fiducial 160 pc → ~1.6 kpc | **CONFIRMED** | Own Reid+14 URC (`kd` form) gives [O III] P1 V=15.5 → **1.785**, P2 14.5 → 1.743 (V_helio→LSR 14.93 → 1.803), P3 11.2 → **1.154**. Table: 1.79 / 1.85 / 1.20. Unit test on the [O III] rows: χ² = 0.05 (3 dof, p = 0.997) for kpc vs **28.9 (p = 2×10⁻⁶) for ×100 pc**. Hα rows: χ² 6.3 (p = 0.18; 1.9 with V_helio-derived V_LSR) vs 42 (p = 10⁻⁸). The median of the six [O III] values is 1.605 kpc, which is the paper's own definition of "~160 pc" |
| a3 | The AJ version did not change the numbers | **CONFIRMED** | Independent Exa extract of the iopscience HTML (`data/derived/verify_paper_raw/aj_aeb11c_iopscience_exa_20261116.txt`, not the wave-1 Tavily copy). It contains "160 pc" ×5, "28 pc" ×4 and "30 pc" ×2. The Table 4 traditional column is identical, and the footnotes still read "in units of 100 pc". The sentence "approximate median of the six [O iii]-based kinematic distance estimates" is unchanged. The MC column did not render in HTML. Crossref has no `update-to`/`updated-by` and 0 records with `filter=updates:` (no erratum, 2026-11-16) |
| a4 | No kinematic distance is meaningful at l≈204°, b≈19° | **CONFIRMED** | The rotation gradient along P1–P3 is 7.6–9.1 km s⁻¹ kpc⁻¹ (0–2 kpc, `kd` form and ×cos b). Therefore ±7 km s⁻¹ of peculiar motion means ±0.8–0.9 kpc. With the *revised* [O III] V_LSR ≈ 0…+2 (DESI/SDSS), the same code gives **0–0.21 kpc**. Either velocity carries essentially no distance information |
| b1 | Table 2 R = 1.03–1.19 implies n_e 190–440, not 64–150 | **REVISED (stronger)** | With the paper's eq. 5: 400/298/247/189 at 10⁴ K (169–358 at 8000 K, √T-scaled). With my **5-level S⁺ atom** (TZ10 Ω; TZ10, MZ82 or RGJ19 A-values; RBS96 cross-check): **225–631 cm⁻³** over the four data sets at 8000–10⁴ K (TZ10+TZ10: 552/407/338/263 at 10⁴ K; 506/373/311/243 at 8000 K). 64/100/150 cm⁻³ needs R = 1.36–1.39 / 1.32–1.35 / 1.26–1.31 (eq. 5: 1.334 / 1.285 / 1.228) |
| b2 | Low-density limit | **REVISED** | R_LDL = **1.454** (TZ10, 10⁴ K), 1.461 (8000 K), 1.437–1.455 across data sets. Eq. 5 → 1.436. The paper's "~1.42" is slightly low. R_HDL = 0.44 (MZ82/RGJ19 A) to 0.52 (TZ10 A) |
| c1 | Paper energetics arithmetic: E = 2.0×10⁵¹, R_PDS 17.1 (not 18), age formula "+1" vs "+1/3" | **CONFIRMED** | V_s = 119.3. **E = 2.01×10⁵¹**. R_PDS = 17.09 pc (18 pc needs E = 2.41×10⁵¹). Age 7.95×10⁴ yr (+1/3) or 8.72×10⁴ (as printed, +1). 6.75×10⁴ yr only with R_PDS = 18 and +1/3. Paper-like errors (n0 = 1 ± 0.3, V_s = 119 ± 5) give E = (1.33–2.73)×10⁵¹ ✓ but age 7.7–8.3×10⁴ (the ±1.8×10⁴ is not reproduced) |
| c2 | At ~288 pc, E = 1.3×10⁵² (n0 = 1) | **CONFIRMED** | 280 pc: R_s = 52.5 pc (paper's 10.74°) or 59.1 pc (Monogem 12.1°), V_s 119 → **E = 1.2–1.7×10⁵² erg**, age 1.4–1.6×10⁵ yr, R/R_PDS 1.85 |
| c3 | At 1.6 kpc the energy is unphysical | **REVISED** | n0 = 1: **2.9–4.2×10⁵⁴ erg** (unphysical). With V_s fixed at 100–150 and E ∝ n0^1.16, n0 = 0.003 gives 2.7–6.8×10⁵¹ and n0 = 0.001 gives 0.8–1.9×10⁵¹, with R_s = 300–338 pc at z ≈ 0.53 kpc and age 0.6–1.1 Myr. "Unphysical for a single SNR" holds. A low-density superbubble-scale shell is formally possible |
| c4 | With the revised [O III] brightness (2–4 R) | **DERIVED** | Brightness enters only through n0. If I ∝ n0, the 3–5× fainter [O III] lowers n0 3–5× and E 3.6–6.5×. A complete shock bound from Balmer (≥0.117 Hβ photons per swept-up H) with I_Hβ ≤ 0.057 R (SDSS south, 2σ) or 0.12 R (DESI-implied) and limb brightening 3–19 gives **n0 ≤ 0.002–0.034 and E ≤ 1.4×10⁴⁹–2.9×10⁵⁰ erg at 280 pc**. Tying n0 to the [O III] yield instead (ε = 2–10 % of the shock energy flux in λ5007) gives n0 = 0.001–0.23 and E ≈ 10⁴⁹–2×10⁵¹ |
| d1 | n_e ~ 1–3 cm⁻³ for the [O III] gas (not 64–150) | **CONFIRMED / refined** | At 280 pc, I = 2–4 R, T = 1–5×10⁴ K, x(O²⁺) = 0.3–1: EM(n_e n_H) = **0.27–10.3 pc cm⁻⁶**. n_e = **0.47–2.9** (filament, L = w = 1.37 pc), **0.11–0.68** (edge-on shell tangent, L = 25.3 pc), 0.08–0.49 (50-pc sheet). n_e ≈ 0.1–3 is consistent; 64–150 would give ~10⁴ R |
| e1 | 23.6–23.9 mag arcsec⁻² ⇒ 10–13 R (AB, 60 Å) | **CONFIRMED (arithmetic)** | 12.75 / 9.94 / 10.12 R (P1/P2/P3). 50–70 Å: 8.3–14.9. Vega/AB ±0.1 mag: 9.1–14.0. If λ4959 is fully transmitted: 7.4–9.5 R |
| e2 | Plausible systematics inflate 10–13 R by 3–5× | **REFUTED as multiplicative** | Stacking every *multiplicative* systematic I can identify in the worst direction lowers the value only to ~5.5–7.5 R (still 1.4–3.7× above 2–4 R). The f/0.95 passband shift goes the *other* way (true/paper = 2.0–4.0 in an untuned converging beam). A 3–5× error needs an *additive* background/flat error of **7–12 % of a 21.5 mag arcsec⁻² sky** (11–20 % of a 22.0 sky), or a background aperture systematically cleaner of starlight than the nebula aperture |

---

## 2. (a) Distance units: method and details

**Own code** (`scripts/verify_paper_kd.py`):
- **Rotation curve.** The Persic+1996 URC with the Reid+2014 parameters as tabulated in `kd/reid14_rotcurve.py` (a1 = 241, a2 = 0.90, a3 = 1.46, R0 = 8.34 kpc, Θ0 = Θ(R0)), plus a flat curve.
- **Projections.** The `kd` form V = R0 sin l (Θ/R − Θ0/R0) has no cos b. A physically projected ×cos b version gives distances 7 % larger.
- **Solar motion.** The IAU standard solar motion (10.27, 15.32, 7.74), as in `kd_utils.py`.
- **Traditional distances.** 1-pc grid, outer-Galaxy branch, 0.1 km s⁻¹ tolerance, NaN when V < V(d_min).
- **Monte Carlo.** The `pdf_kd` emulation draws 10⁴ samples of a1, a2, a3 and R0 (Gaussian) and of V (Table 4 RV error). The grid floor is 0.01 kpc. It uses a linear-combination boundary-corrected KDE with lower bound 0 and h = σ n^−1/5, and reports the mode and 68.3 % minimum-width interval.

**Expected V_LSR toward P1** (`kd` form / ×cos b):

| distance | V_LSR (`kd` form) | V_LSR (×cos b) |
|---|---|---|
| 160 pc | 1.59 | 1.50 |
| 280 pc | 2.76 | 2.60 |
| 1.0 kpc | 9.26 | 8.73 |
| 1.6 kpc | 14.10 | 13.29 |

| row | V_LSR (table) | V from V_helio (IAU) | Table trad. | mine (`kd` form) | Table MC mode [HPD] | my MC mode [HPD] |
|---|---|---|---|---|---|---|
| P1 [O III] | 15.5 | 14.70 | 1.79 | 1.785 | 1.63 [0.79–2.64] | 1.66 [1.13–2.36] |
| P2 [O III] | 14.5 | 14.93 | 1.85 | 1.743 (1.803 from V_helio) | 1.58 [0.75–2.81] | 1.67 [1.12–2.33] |
| P3 [O III] | 11.2 | 11.22 | 1.20 | 1.154 | 0.98 [0.29–2.06] | 1.16 [0.75–1.53] |
| P1 [S II] | −1.0 | −0.3 | 0.16 | **no solution** | 0.06 [0.01–1.23] | 0.02 [0.01–0.22] (66 % unsolved) |
| P4 Hα | 4.5 | 7.32 | 0.80 | 0.46 (0.77 from V_helio) | 0.54 [0.02–1.72] | 0.43 [0.29–0.61] |

**Independent evidence for kpc.**
- **The MC grid floor.** All 12 red-line Monte Carlo intervals in Table 4 have lower bounds (mode − err) of 0.00–0.02. That is `pdf_kd`'s hard floor of 0.01 kpc. Read as "100 pc" units, the floor would be 1 pc, which `kd` cannot produce.
- **The [O III] rows** fall on the rotation curve only in kpc (Fig. `verify_paper_kd.png`, filled vs open symbols).

**New findings (a).**
- **Red-line traditional distances.** These are 0.14–0.49 for [S II]/[N II] and are *not* reproducible from the tabulated V_LSR. Three of those velocities are negative and have no outer-Galaxy solution, yet Table 4 lists 0.16–0.18. The values correspond to V ≈ 1.5–4.8 km s⁻¹. They must therefore be aggregates of per-fibre solutions in which negative-velocity fibres were dropped by `kd`, which is a positive truncation bias.
- **Per-fibre Hα check.** I used the digitized per-fibre Hα V_LSR from Fig. 5 (geometry tag) inside the drawn P1–P4 regions. d(median V) = 0.41 / 0.39 / 0.62 / 0.60 kpc vs Table 4 Hα 0.61 / 0.54 / 0.58 / 0.80. The aggregation rule therefore cannot be identified.
- **P2 inconsistency.** The P2 [O III] distance (1.85) is consistent with V_LSR recomputed from V_helio (14.93 → 1.80) but not with the printed 14.5 (→ 1.74). The V_LSR column and the distance columns were derived separately (cf. wave-1 item 6: V_LSR − V_helio varies between lines).
- **MC interval widths.** The Table 4 MC intervals are 1.5–2.5× wider than `pdf_kd` gives with the tabulated velocity errors. The red-line upper bounds (+0.9 to +1.6 kpc) need V up to ~15 km s⁻¹. The paper probably used per-fibre velocities or larger errors, but this cannot be determined.

---

## 3. (b) [S II] density with a 5-level atom

`scripts/verify_paper_sii.py` solves the full 5-level statistical equilibrium (4S3/2, 2D3/2, 2D5/2, 2P1/2, 2P3/2).

**Atomic data.** These are plain-text tables from the PyNeb repository (transcriptions of the original papers), fetched 2026-11-16 (`data/derived/verify_paper_raw/pyneb_atomic/fetch_log.txt`):
- **Energies:** NIST.
- **Collision strengths:** Tayal & Zatsarinny 2010, ApJS 188, 32 (TZ10; Ω(4S–2D3/2, 4S–2D5/2, 2D3/2–2D5/2) = 2.56, 3.83, 6.89 at 10⁴ K). Ramsbottom, Bell & Stafford 1996 (RBS96) as a cross-check.
- **A-values:** three sets.
  - TZ10: A(6716) = 2.195×10⁻⁴, A(6731) = 6.32×10⁻⁴ s⁻¹.
  - Mendoza & Zeippen 1982 (MZ82): 2.60×10⁻⁴, 8.82×10⁻⁴.
  - Rynkun+2019 (RGJ19): 2.84×10⁻⁴, 9.43×10⁻⁴.
- The Podobedova+2009 set is not in the PyNeb tree, so it was not used.

Check of the low-density limit by hand: TZ10 Ω plus 2P cascades gives 1.457, and the code gives 1.454.

| R = I6716/I6731 (region) | eq. 5, 10⁴ K | eq. 5 × √0.8 (8000 K) | TZ10+TZ10 10⁴ / 8000 K | TZ10+MZ82 10⁴ / 8000 | TZ10+RGJ19 10⁴ / 8000 | RBS96+MZ82 10⁴ / 8000 |
|---|---|---|---|---|---|---|
| 1.03 (P4) | 400 | 358 | 552 / 506 | 568 / 520 | 631 / 578 | 502 / 461 |
| 1.10 (P1) | 298 | 267 | 407 / 373 | 424 / 389 | 472 / 433 | 372 / 344 |
| 1.14 (P2) | 247 | 221 | 338 / 311 | 355 / 326 | 395 / 363 | 310 / 287 |
| 1.19 (P3) | 189 | 169 | 263 / 243 | 278 / 257 | 310 / 286 | 241 / 225 |

- Eq. 5 (Proxauf+14) runs ~25–30 % *below* the 5-level curves with TZ10 Ω (Fig. `verify_paper_sii.png`). Even the paper's own calibration choice is biased low, yet it still gives 3–6× the quoted 64–150.
- **Can per-fibre noise produce a median of 64?** This is a new simulation using eq. 5. It drops fibres above the low-density limit (eq. 5 has no solution there) and takes the median of the rest. A median of 64 cm⁻³ needs the true R ≈ 1.34 (σ_R = 0.03), or R ≳ 1.42 with σ_R = 0.15; in the second case 47 % of fibres are dropped and the kept mean R is 1.31. With σ_R = 0.30 the median never drops below ~150. **No noise model reconciles R = 1.03–1.19 with a per-fibre median of 64 cm⁻³.** The Table 2 ratios and the Fig. 7b histogram cannot describe the same fibres under the same eq. 5.
- **Table 2 intensities (MC, 10⁴ K, TZ10):**

  | region | R 16–84 % | n_e |
  |---|---|---|
  | P1 | 0.88–1.32 | 110–1040 |
  | P2 | 0.95–1.35 | 84–780 |
  | P3 | 0.85–1.44 | 11–1230 |
  | P4 | 0.76–1.34 | 90–1900 |

  The density is essentially unconstrained per region, and the nebula-free P4 still has the highest central value.

---

## 4. (c) Energetics

**Formulas.**
- Dopita 1979: n_[S II] = 45 n0 (V_s/100)².
- Cioffi+88 / Boumis+22: E = 6.8×10⁴³ n0^1.16 V_s^1.35 R_s^3.16; R_PDS = 14.0 E51^(2/7) n0^(−3/7); t = 9.98×10³ E51^(3/14) n0^(−4/7) [(R/R_PDS)^(10/3) + 1/3].
- As a check, the Cioffi PDS shell speed at the derived age is 116 km s⁻¹, consistent with the V_s that was put in.

**Explosion energy E (erg)** (selected; full grid in `verify_paper_energetics.json → energetics.grid`). R_s scales with d at either the paper's 10.74° (30 pc at 160 pc) or the Monogem 12.1°.

| d | R_s (pc) | n0=1, V_s=119 (paper) | V_s=100, n0=0.03 | V_s=100, n0=0.01 | V_s=150, n0=0.01 | V_s=100, n0=0.003 |
|---|---|---|---|---|---|---|
| 160 pc | 30 | 2.0×10⁵¹ | 4.1×10⁴⁹ | 7.6×10⁴⁸ | 1.3×10⁴⁹ | 1.9×10⁴⁸ |
| 280 pc | 52.5–59.1 | **1.2–1.7×10⁵²** | — | 4.5–6.5×10⁴⁹ | 0.8–1.1×10⁵⁰ | 1.1–1.6×10⁴⁹ |
| 1.6 kpc | 300–338 | **2.9–4.2×10⁵⁴** | 3.9–5.7×10⁵² | 1.1–1.6×10⁵² | 1.9–2.8×10⁵² | 2.7–4.0×10⁵¹ |

**Ages and phase.** Ages are ~(2/7) R/V: 1.1–2.1×10⁵ yr at 280 pc and 0.6–1.1 Myr at 1.6 kpc. R/R_PDS = 0.9–2.4, so the shell is near or beyond the PDS transition everywhere.

**Lowering n0.** Lowering n0 while keeping Dopita's relation with n_e = 64 forces V_s = 377 (n0 = 0.1) or 1193 km s⁻¹ (0.01), which is not a radiative shell. The [S II] density belongs to the WIM (wave 1) and should not set V_s.

**Brightness-tied n0 (derived, d = 280 pc, R_s = 59.1 pc).**
- **Complete radiative shock (Balmer bound).** Each swept-up H recombines at least once, giving ≥0.117 Hβ photons (case B). The face-on Hβ is therefore ≥ 1.17 n0 V₁₀₀ R; Dopita & Sutherland 1996's F_Hβ ≈ 7.44×10⁻⁶ V₁₀₀^2.41 n0 erg cm⁻² s⁻¹ corresponds to ~1.6 recombinations per H (formula quoted from memory). The edge-on limb brightening is L_b = 2√(2R_s/w) = 18.6 for w = 1.37 pc; the range L_b = 3–18.6 was used.
  - With I_Hβ ≤ 0.057 R (SDSS south ON−OFF 2σ) or 0.12 R (DESI ridge 2.4 R / ratio 20): **n0 ≤ 0.0017–0.034** (V_s 100–150).
  - This gives **E ≤ 1.4×10⁴⁹–2.9×10⁵⁰ erg**, compatible with the Monogem E0 ~ (0.1–0.4)×10⁵¹ (Knies+2024 as quoted in `notes/context.md`). No 2×10⁵¹ erg is needed.
- **[O III]-tied.** The face-on λ5007 is 295 ε n0 V₁₀₀³ R. For I = 2–4 R, ε = 0.02–0.1 and L_b 3–18.6: n0 = 0.004–0.23 (V 100) or 0.001–0.067 (V 150), giving E ≈ 10⁴⁹–2.4×10⁵¹. Compared with the old 11.5 R, n0 is lower by 2.9–5.8× and E by 3.4–7.6×.
- The shockmodel agent should pin down ε (incomplete shocks). DESI [O III]/Hβ ≈ 20 already implies that the [O III] is not from a *complete* shock.

---

## 5. (d) Self-consistency at d = 280 pc

**O²⁺ λ5007 emissivity** (own 5-level atom; Aggarwal & Keenan 1999 Ω valid to 2×10⁵ K; Froese Fischer+2004 / Storey & Zeippen 2000 A-values; branching 0.751; n_e = 1). Storey, Sochi & Badnell 2014 Ω agree within 2–3 % up to 2×10⁴ K.

| T_e (K) | ε₅₀₀₇ (cm³ s⁻¹ ph per n_e n_O²⁺) | I₅₀₀₇ per EM, x = 1 (R per pc cm⁻⁶) | [O III]/Hα photon ratio, x = 1 / 0.3 | λ4363/λ5007 |
|---|---|---|---|---|
| 1.0×10⁴ | 8.57×10⁻¹⁰ | 1.30 | 3.6 / 1.1 | 0.0065 |
| 1.5×10⁴ | 1.98×10⁻⁹ | 2.99 | 11.9 / 3.6 | 0.020 |
| 2.0×10⁴ | 2.92×10⁻⁹ | 4.42 | 22.9 / 6.9 | 0.034 |
| 3.0×10⁴ | 4.11×10⁻⁹ | 6.22 | 46 / 14 | 0.058 |
| 5.0×10⁴ | 4.96×10⁻⁹ | 7.50 | 88 / 27 | 0.090 |

Assumptions: solar O/H = 4.90×10⁻⁴ (Asplund+2009), n_e/n_H = 1.15, and Hα α_eff = 1.17×10⁻¹³ T4^−0.9 (Haffner+2003 EM relation).

**Line-of-sight paths.** At 280 pc the 0.28° strand FWHM is w = 1.37 pc and the Monogem-limb radius is R_s = 59.1 pc. Paths used:
1. A filament: L = w = 1.37 pc.
2. An edge-on shell seen at tangency: L = 2√(2R_s w − w²) = 25.3 pc.
3. A flat 50-pc sheet.

**Required emission measure and density:**

| I(5007) | x(O²⁺) | EM(n_e n_H), pc cm⁻⁶ | n_e filament | n_e shell-tangent | n_e 50-pc sheet | Hα from same gas (R) |
|---|---|---|---|---|---|---|
| 2–4 R | 1 | 0.27–3.1 | 0.47–1.61 | 0.11–0.38 | 0.08–0.27 | 0.02–0.97 |
| 2–4 R | 0.3 | 0.89–10.3 | 0.86–2.94 | 0.20–0.68 | 0.14–0.49 | 0.07–3.2 |

(Ranges span T = 1–5×10⁴ K; the high-EM ends are the 10⁴ K cases.)

- **Density range.** n_e ≈ 0.1–3 cm⁻³ **is consistent**: 0.1–0.7 for an edge-on sheet or shell, 0.5–3 for a cylinder. The density scales as d^−½ (×1.32 at 160 pc, ×0.42 at 1.6 kpc).
- **Hα.** SDSS ON−OFF Hα is 0.17 ± 0.03 R and Finkbeiner gives < 0.43 R. These exclude 10⁴ K gas with x(O²⁺) ≲ 0.3, which would give 1.6–3.2 R of Hα.
- **[O III]/Hα.** The measured [O III]/Hα = 4.7 ± 1.0 needs T ≥ 1.08×10⁴ K (x = 1) or ≥ 1.68×10⁴ K (x = 0.3).
- **Pressure.** P/k ≈ 1.96 n_e T gives 5×10³–2.8×10⁴ K cm⁻³ for the shell geometry, comparable to the Monogem hot gas (2 n_e T ≈ (0.8–2.9)×10⁴ for n_e = 3–7×10⁻³ and kT 0.11–0.17 keV, Knies+24 via the context note). A cylinder gives 2×10⁴–1.2×10⁵, i.e. up to ~4–10× over-pressured.

**New finding (d), 4363 vs Balmer tension.** The DESI λ4363/λ5007 = 0.056 (0.022–0.077) gives T = **2.9×10⁴ K** (1.6–4.0×10⁴) with this atom. At 2.9×10⁴ K, collisional [O III] with case-B Balmer from the *same* gas gives [O III]/Hα ≈ 43 x(O²⁺) (photons). The measured ratios are ≈ 4.7 (SDSS) and ≈ 5.3 (DESI [O III]/Hβ = 20 converted to Hα photons). These imply x(O²⁺)·(O/H)/(O/H)⊙ ≈ 0.11–0.12. Otherwise, the Balmer lines come from cooler, recombining gas (a stratified shock), or λ4363 is inflated (Hg I 4358; desi agent).

---

## 6. (e) Surface brightness → Rayleigh

1 R at λ5007 = 7.42×10⁻¹⁸ erg s⁻¹ cm⁻² arcsec⁻². In AB magnitudes for a 60 Å band, **2/3/4/10/13 R = 25.61/25.17/24.86/23.86/23.58 mag arcsec⁻²**.

**Multiplicative systematics:**
- **AB vs Vega.** The magnitude system is not stated. Blanton & Roweis 2007 offsets (g −0.08, V +0.02, quoted from memory) bracket ±0.1 mag, i.e. ±10 %.
- **Effective width.** 50–70 Å gives −17/+17 %.
- **λ4959 leakage.** At most −25 % if both lines are fully transmitted.
- **Interference-filter angle shift.** Model: 3-cavity profile with a 60 Å FWHM, n* = 1.8–2.05, CWL 5007/5015/5030 Å. The stellar zero point calibrates ∫T dλ while the line sees T(5007) + T(4959)/2.98. The resulting true/paper factor is:
  - 0.97–1.08 for an f/2.8 beam or near-normal incidence;
  - **2.0–4.0** for an untuned filter in an f/0.95 converging beam (5007 sits at 15–36 % of peak), so SWONS would *under*-estimate;
  - 0.93–2.8 for a collimated beam at 10–15° field angle.

  It is not stated whether Table 3 used the 50-mm f/0.95 or the 180-mm f/2.8 data, or the real filter curve.
- **Stars.** The paper itself says P1 is contaminated by ~0.3 mag.
- **Worst-case stack.** All of these together give 10–13 R → ~5.5–7.5 R, still 1.4–3.7× above 2–4 R.

**Additive systematics:**
- **Sky fraction.** The claimed signal is 11–15 % of a 21.5 mag arcsec⁻² sky (the paper's SQM). The 2–4 R alternative is 2.3–4.5 %. The difference (6–11 R) is **7–12 % of the sky** (11–20 % for a 22.0 narrow-band continuum).
- **Unresolved starlight (Gaia DR3, measured).** VizieR I/355 cones of 0.25° radius around P1–P4 (`data/derived/verify_paper_raw/gaia_isl/`): ~9 500–10 800 stars deg⁻² with G = 17–21.

  | stars included | SB (G mag arcsec⁻²) | ≈ R-equivalent at λ5007 |
  |---|---|---|
  | all stars | 22.2 (P1) / 22.8 (P2, P3) / 23.4 (P4) | — |
  | G > 14 | 24.7–24.9 | ≈4–4.6 R |
  | G > 15 | 25.2–25.4 | ≈2.5–3 R |
  | G > 16 | 25.9–26.0 | ≈1.5 R |
  | G > 17 | 26.5–26.7 | ≈0.8 R |

  The colour term (~0.3 mag fainter at 5007 AB) is ignored in the R-equivalents. The faint-star light is uniform to ≤0.15 mag between P1–P4, so any nearby background aperture removes it to ≲0.5 R. Unresolved starlight is the same order as the claimed signal only if the "star-free" background aperture was effectively cleaner than the nebula aperture, or if bright-star halos are included. The bright-star term is largest at P1 (total 22.2 vs 22.8 mag).
- **Zero points.** For a 50 mm f/0.95 lens and a 6 nm band, 1 e⁻ s⁻¹ corresponds to m_AB ≈ 14.8 (QE 0.7, T 0.8). The tabulated ZPs of 11.23–11.32 cannot be interpreted without units (gain, pixel area, exposure) and were not used.

**Verdict (e).** The conversion arithmetic is right. A 3–5× inflation is not explained by AB/Vega, bandpass, the f/0.95 geometry or [O III]4959. It requires an additive background or flat-field error of ~7–20 % of the sky, or a mismatched background aperture. Without the SWONS frames and aperture definitions this cannot be decided. The independent NSNS imaging (2–3.6 R) indicates an imaging-calibration problem rather than a fibre problem.

---

## 7. New findings (summary)
1. **The `kd` grid floor in Table 4.** The red-line MC intervals stop at 0.00–0.02, i.e. at `pdf_kd`'s 0.01-kpc floor. This is direct, model-independent proof that Table 4 is in kpc.
2. **The red-line "traditional" distances are not reproducible** from the tabulated V_LSR. Negative velocities have no solution. They are per-fibre aggregates with negative-velocity fibres dropped, which biases them upward. The P2 [O III] distance follows V_helio rather than the printed V_LSR.
3. **Eq. 5 is biased low** by ~25–30 % relative to a 5-level atom with TZ10 Ω. Modern atomic data put Table 2 at 225–631 cm⁻³. No per-fibre noise model gives a median of 64.
4. **The Balmer non-detection caps n0.** For a complete radiative shock at 280 pc, n0 ≲ 0.002–0.03 and E ≲ 3×10⁵⁰ erg, consistent with the Monogem E0 rather than the paper's 2×10⁵¹.
5. **Tension between λ4363 and [O III]/Balmer.** T(4363) ≈ 2.9×10⁴ K vs [O III]/Hα ≈ 5 requires x(O²⁺) ≈ 0.1 in isothermal gas, a stratified shock, or a contaminated λ4363.
6. **The f/0.95 geometry would make SWONS under-estimate λ5007**, not over-estimate it. The 10–13 R vs 2–4 R conflict needs an additive photometric error.

## 8. Limitations
- **`kd` inputs unknown.** `kd` was emulated, not executed. The per-fibre velocities and velocity errors behind the Table 4 MC intervals are unknown, so the interval widths and the red-line traditional values are not reproduced. Only their units are tested.
- **Atomic data.** Taken from PyNeb's tabulations, not re-typed from the original papers. Podobedova+2009 A-values were not available. The Proxauf+14 T-correction was approximated as √T.
- **Energetics are scaling arguments.**
  - The Cioffi/Boumis formulae assume a uniform medium and a single explosion.
  - The Hβ yield bound assumes a *complete* shock.
  - ε (the λ5007 yield) and L_b are free parameters.
  - The DS96 Hβ formula and the Blanton & Roweis offsets are quoted from memory, and the Monogem E0/pressures are taken from a teammate's literature summary.
- **The filter model is generic.** The real SWONS filter curve, the mounting position and the camera used for Table 3 are unknown.
- **The Gaia integrated starlight** uses G-band (Vega) fluxes. The completeness limit (G ≈ 21) and the colour term to λ5007 are ignored.
- **Inputs from other agents.** The DESI/SDSS brightness, Hβ and λ4363 values are wave-1 results, not re-measured here.

## 9. Files and data provenance (all accessed 2026-11-16)
- **Scripts:** `scripts/verify_paper_kd.py`, `scripts/verify_paper_sii.py`, `scripts/verify_paper_energetics.py`, `scripts/verify_paper_gaia_isl.py`.
- **Derived data:** `data/derived/verify_paper_kd.json` (+ `verify_paper_kd_mc_samples.npz`), `verify_paper_sii.json`, `verify_paper_energetics.json`.
- **Figures:** `figures/verify_paper_kd.png`, `figures/verify_paper_sii.png`, `figures/verify_paper_oiii_conversion.png`.
- **Raw data** (`data/derived/verify_paper_raw/`):
  - `kd_*` — `kd` source text from https://raw.githubusercontent.com/tvwenger/kd/master/; re-fetched 12:05Z and byte-identical; commit list in `kd_commits.json`.
  - `aj_aeb11c_iopscience_exa_20261116.txt` — AJ HTML via Exa from https://iopscience.iop.org/article/10.3847/1538-3881/aeb11c.
  - `crossref/` — https://api.crossref.org/works/10.3847/1538-3881/aeb11c and the `filter=updates:` query.
  - `pyneb_atomic/` — data tables and `fetch_log.txt`.
  - `gaia_isl/` — VizieR cones and `query_log.txt`. A NOIRLab Data Lab aggregate query returned HTTP 502 and its log is kept.
