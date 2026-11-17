# verify_altmodel: can anything other than a shock make the Silk [O III]?

Wave-2 modelling agent, tag `altmodel`. Written 2026-11-16 (resumed run).

I wrote all the code (`scripts/verify_altmodel_*.py`): a small CHIANTI-based atomic toolkit, a one-zone photoionisation-equilibrium code with thermal balance, an NEI ionisation solver, and line synthesis with GALEX band conversion. I did not use Cloudy, MAPPINGS or pyneb.

Labels used below:
- **measurement**: a teammate's measured value, used as input;
- **derivation**: my calculation;
- **literature**: a published table or number.

Energy ratios are used unless marked "photons".

---

## 1. Bottom line

**No non-shock mechanism tested here reproduces the Silk.** Each scenario fails at least one robust constraint, and most fail several. The only family that reproduces the optical spectrum is gas cooling from ≳ 10⁵ K at high pressure and stopped before H recombines. That is the cooling zone of an *incomplete radiative shock*.

Three model-independent facts drive these verdicts:

1. **Two temperatures** (derivation, my CHIANTI v7.1.4 inversions).
   - T_e([O III]) = 4.6 × 10⁴ K (3.1–5.8 × 10⁴) from 4363/5007 = 0.083 (0.059–0.096).
   - T([O II]) = 1.25 × 10⁴ K (1.0–1.7 × 10⁴) from 7325/3727 = 0.026 (0.017–0.038). At n_e = 100 it is 1.18 × 10⁴ K.
   - So O⁺⁺ emits in gas about 4× hotter than the O⁺ gas.
   - This rules out every **single-temperature** model: photoionisation equilibrium (PIE), photoionisation with extra heating, a warm CIE layer, and a fossil.
   - It is the signature of a cooling flow.
2. **High pressure** (derivation).
   - [O II] 3729/3726 = 1.30 (measurement: DESI, verify_desi; SDSS 1.29). This gives n_e = 88 / 100 / 125 cm⁻³ at T = 1.0 / 1.25 / 2.0 × 10⁴ K, or 47–208 cm⁻³ over the DESI subset range 1.17–1.37.
   - The CHIANTI v7.1.4 low-density limit is 1.46.
   - So P/k ≈ 2.8 × 10⁶ K cm⁻³ (1.4–5.7 × 10⁶). That is 100–400× the Monogem interior pressure (P/k ≈ 1.5 × 10⁴).
   - Any interface, mixing layer, photoionised skin or fossil in pressure balance with its surroundings predicts the low-density ratio of 1.46.
   - The pressure equals a ram pressure ρ₀v² for v ≈ 127 (n₀ / 1 cm⁻³)^−½ km/s.
   - At this density, the 2 R of [O III] comes from a layer only 2–11 × 10⁻⁴ pc (40–230 AU) deep along the line of sight. The volume filling factor is ~10⁻⁵ of a 25-pc limb path.
3. **[O III]/Hβ ≥ 15 needs T ≥ 1.3 × 10⁴ K** (derivation; ≥ 2.0 × 10⁴ K for the DESI value of 47), at solar O/H with H fully ionised.
   - This is the maximum possible ratio for any mix of O⁺ and O⁺⁺ at the observed [O II]/[O III] = 1.37.
   - Photoionisation equilibrium reaches only 7–12 × 10³ K at the required ionisation state.

**Scenario verdicts**

| # | Scenario | Verdict | Main reasons |
|---|---|---|---|
| 1 | Thermal conduction front (evaporating or static) | **ruled out** | 4363/5007 = 0.11–0.14; C IV/[O III] = 1.6–16 photons (FUV 540–5300 CU per R vs 68–92 observed); 7325/3727 = 0.07; face-on I(5007) = 2 × 10⁻⁴ to 3 × 10⁻³ R per front, so ≥ 30 fronts are needed; low-density [O II] doublet. N(O VI) = 0.9–2.3 × 10¹³ cm⁻² per front agrees with the literature scale, which validates the model. |
| 2 | Turbulent mixing layer (Monogem hot gas / HI wall) | **ruled out** | Optical ratios like isobaric cooling (good), but: energy flux 2 × 10⁻⁴ erg cm⁻² s⁻¹ means a hot-gas inflow of 377 km/s (> c_s,hot); at a plausible 10 km/s, I(5007) ≈ 0.05 R; FUV 254 CU per R; 7325/3727 ×2; X-ray filament 2 × 10³–2 × 10⁴× Monogem; low-density doublet. |
| 3 | Non-equilibrium radiative cooling of Monogem hot gas (GS07 isochoric/isobaric) | **ruled out** | Reproduces [O III]/Hβ 25–48, [O II]/[O III] 1.2–1.7 and 4363/5007 0.08–0.10 only if truncated before H recombines (complete flows give [O III]/Hβ = 3–6). Needs 1.5–2 × 10⁻⁴ erg cm⁻² s⁻¹, about 70 % of it radiated above 5 × 10⁵ K, which predicts an X-ray filament ≥ 10³× the Monogem shell. eRASS1 shows none. FUV 176–254 CU per R. |
| 3w | Warm layer held at a single T ≈ 6 × 10⁴ K (CIE, unspecified heating) | **disfavoured** | Surprisingly good: [O III]/Hβ ≈ 100, [O II]/[O III] 1.2–1.6, 4363/5007 0.096–0.098, FUV 86–92 CU per R. But 7325/3727 = 0.08 (×3 too high), NUV 84–93 (×2), [N II]/Hα ≈ 2. Needs 1.7 × 10⁻⁵ erg cm⁻² s⁻¹ of heating with no identified source, at the thermally unstable peak of the cooling curve. |
| 4a | Photoionisation by WD J0736+1444 (81 kK) | **ruled out** | At the required ionisation: [O III]/Hβ ≤ 5–9, 4363/5007 = 0.006, [N II]/Hα = 0.2. The ionisation parameter falls short by 27× (n_rms, coplanar at 6.4 pc) to 9 × 10³× (n from [O II]); at the Monogem-limb distance by 10⁴–4 × 10⁶×. |
| 4b | Photoionisation by Monogem / Local Bubble EUV | **ruled out** | Boundary flux Φ = 340–6200 cm⁻² s⁻¹, so log U ≤ −6.1, short by 630–2 × 10⁵×. Even at the right U: [O III]/Hβ ≤ 10–14, He II/[O III] = 0.06 (> 0.026), 4363/5007 = 0.011. |
| 4c | Photoionisation by the general ISRF / DIG | **ruled out** | [O III]/Hβ ≤ 1.3 at any U. |
| 4d | Photoionisation plus extra heating (any source) | **ruled out** | Needs 78–350× the photo-heating to reach T = 4 × 10⁴ K (4363/5007 = 0.075), plus the U deficit of 4a/4b. 7325/3727 = 0.076 (×3); NUV 94. |
| 5 | Photoionised shock precursor (alone) | **ruled out** | [O III]/Hβ = 6 (≤ 10.7 at any U), 4363/5007 = 0.007, He II/[O III] 0.05. Pre-shock gas (n ~ 0.1–1) gives the low-density doublet. The parent shock would dominate in any case. |
| 6 | Fossil / recombining structure (old ionised shell, PN halo) | **ruled out** | At the doublet density [O III] fades by e in 21–34 yr (4–7 × 10³ yr at n_H = 0.5). It matches the ratios only in a 9–21-yr window, and only if it started ≥ 4 × 10⁴ K (i.e. shock-heated). 7325/3727 is 0.07 in that window. Photoionisation-heated starts (2 × 10⁴ K) never match. |

**What remains.** The optical ratios, the pressure, the two-temperature structure and the FUV level point the same way: the cooling zone of a radiative shock of ~120–150 km/s into n₀ ~ 1 cm⁻³ gas, truncated before full recombination, with C IV partly suppressed. This is the shockmodel agent's territory.
- My cooling-from-3 × 10⁵ K track is that case.
- It still predicts FUV 236 CU per R against 68–92. Its non-resonant part is 98 CU per R, so C IV must be suppressed ≳ 3–5×.
- 7325/3727 is ×2 too high, so the low-T [O II] zone must be supplied by the photoionised post-shock layer, which a pure cooling track lacks.

---

## 2. Observational inputs adopted

From the final wave-2 notes, not re-measured here except where stated.

| Quantity | Adopted | Source |
|---|---|---|
| I(5007) peak | 2 R (1.4–4) | verify_desi (south median 1.6, W strand 3.9), verify_sdss (1.4–2.5), NSNS 2.6–3.6 |
| [O III]/Hβ | ≥ 15 (robust); 20 (SDSS, wave-1 DESI); 47 (verify_desi) | see caveat 1 in §5 |
| [O II]3727/[O III] | 1.37 (1.1–1.9) | DESI 1.37 (south 1.11); SDSS south 1.84 |
| [O III] 4363/5007 | 0.083 (0.059–0.096) | DESI wave 2 (Hg coefficient ≈ 0); SDSS 0.038 ± 0.035 |
| [O II] 7320+30/[O III] | 0.036 (0.033–0.042), so 7325/3727 = 0.026 | DESI |
| [O II] 3729/3726 | 1.30 (1.24–1.33; subsets 1.17–1.37) | DESI; SDSS 1.29 (1.07–1.59) |
| He II/[O III], [Ne III]/[O III] | < 0.026; 0.099 (0.057–0.109) | DESI |
| [N II]/Hα, [S II]/Hα (south) | 1.08, 0.51 (all: 1.21, 0.75) | DESI wave 2 (differential) |
| GALEX FUV, NUV per R(5007) | 68–92 CU, 36–49 CU (up to 80) | verify_xrayfuv |
| Geometry | G = 2√(2R_c/Δ) = 18.7; Δ = 0.28° × 280 pc = 1.37 pc; R_c = 60 pc; LOS path 25.6 pc | wave-1 geometry, Monogem radius |
| Monogem | kT 0.14 keV; n_e 3.1–7.2 × 10⁻³; R 60 pc; P/k ≈ 1.5 × 10⁴ | Knies+2024 via context.md |
| WD J0736+1444 | 81.3 kK, 150.5 pc, Q(H) = 6.6 × 10⁴⁴ s⁻¹, Q(>35 eV) = 1.1 × 10⁴⁴ s⁻¹, 2.44° (6.4 pc) from the spine | context.md |

---

## 3. Methods

### 3.1 Atomic toolkit (`verify_altmodel_atomic.py`, `verify_altmodel_lines.py`)

**CHIANTI v7.1.4** (literature data; downloaded and logged in `data/raw/verify_altmodel/MANIFEST.jsonl`).
- Multi-level ions use the energies, A-values and Burgess–Tully-scaled Maxwellian collision strengths.
- Lines are summed over each multiplet window: [O III], [O II], [N II], [S II], [Ne III], [S III], [Ar III], C IV, O III], N III], N IV], O IV], Si IV, C III], C II], [Ne IV], Mg II, Si III], N V, O VI, C III λ977, and collisional He II 1640.
- I added a 10⁻²⁰ s⁻¹ leak to ground. It only keeps CHIANTI levels without rates (Mg II) non-singular.

**Rates.**
- Collisional ionisation: Dere (2007, VizieR J/A+A/466/771).
- Recombination: CHIANTI RR/DR fit parameters.
- Photoionisation: Verner & Yakovlev outer-shell fits (CHIANTI `verner_short`).
- Charge exchange: O⁺/O⁺⁺ + H at 10⁻⁹ cm³ s⁻¹. This is approximate; its sensitivity is tested in §4.4.

**Hydrogen and helium.**
- H: case-B recombination (Draine 2011 fits), plus collisional Balmer emission of H⁰ from the CHIANTI H I model (n ≤ 5).
- He II 4686 = 12 × j(Hβ) per He⁺⁺ (±20 %); He II 1640 = 7 × 4686 (±30 %).
- Two-photon continuum: Nussbaumer & Schmutz (1984) shape, normalised to 2.000 photons per decay (checked).

**GALEX bands** (own code; SVO FUV/NUV curves already on disk).
- My W_eff values are 354 Å (C IV), 513 (He II), 632 (O III]), 259 (N IV]), 774 (1400 blend) and 1614 (C III]). They agree with xrayfuv's 357/516/636/261/779 to ≤ 1 %, which is an independent check.
- 1 R of C IV = 225 CU.
- Extinction: CCM89 (A_FUV/E(B−V) = 8.15 at 1528 Å), with E(B−V) = 0.03 applied relative to 5007.

**Sanity checks.**
- 5007 emissivity at 10⁴ K: 3.54 × 10⁻²¹ erg cm³ s⁻¹ per O⁺⁺ per n_e. All-O⁺⁺ solar gas then gives [O III]/Hβ = 13.4. This matches the standard abundance formula (12 + log O⁺⁺/H⁺ = log I/Hβ + 6.200 + 1.251/t …), which gives 13.0.
- 4363/5007 = 0.056 inverts to 2.94 × 10⁴ K. The classical 7.90 exp(32900/T) formula gives 2.98 × 10⁴ K.

### 3.2 Model scripts

- **`verify_altmodel_thermo.py`**: the mechanism-independent "thermometer": the maximum [O III]/Hβ vs T for H⁺ + O⁺/O⁺⁺ gas at fixed [O II]/[O III], plus T inversions.
- **`verify_altmodel_photoion.py`**: one-zone, optically thin PIE with thermal balance.
  - Heating: photo-electrons.
  - Cooling: every CHIANTI level of H, He, C, N, O, Ne, S ions, plus recombination, free-free and H collisional ionisation.
  - Spectra: blackbodies of 35, 81 and 150 kK; a 0.14 keV "hot gas" exp(−E/kT)/E form; and 28 and 50 eV precursor-like forms. All are crude spectral shapes.
  - Grid: log U = −5…−0.5.
- **`verify_altmodel_cooling.py`**.
  - Gnat & Sternberg (2007) published NEI ion fractions (isochoric and isobaric, Z = 1 and 2) and their cooling efficiencies.
  - Integrated line energy per H nucleus: E_l = c_p k ∫ L_l (n_tot/n_H) / ((n_e/n_H) Λ) dT, between 1.6 MK (or 3 × 10⁵ K) and 1–5 × 10⁴ K. Variants add the H-recombination tail.
  - Warm-CIE layer: GS07 CIE fractions, cross-checked against my own CIE balance. O⁺⁺ at 5 × 10⁴ K: 0.31 (own) vs 0.34 (GS07); C³⁺ at 10⁵ K: 0.27 vs 0.29.
  - Static radiative conduction front: F dF/dT = κ₀T^{5/2} n_e n_H Λ, with CIE ions and κ₀ = 6 × 10⁻⁷.
  - Steady classical evaporation front: T^{5/2} ∝ x, with my NEI solver run from a WIM-like start. The classical J for R_cl = 1–60 pc is capped at the saturated flux; q_sat = 1.7 × 10⁻⁴ erg cm⁻² s⁻¹ for Monogem gas.
- **`verify_altmodel_budget.py`**.
  - [O II]-doublet densities.
  - rms densities for I(5007) on the 25.6-pc limb path.
  - The U needed vs the U available for each source, and the Strömgren budget.
  - PIE at fixed T (photoionisation plus heating) on a finer grid with Si.
  - Fossil: isochoric, unheated NEI recombination with cooling, at n_H = 100 and 0.5.
- **`verify_altmodel_summary.py`**: the scorecard, TML energetics, and the X-ray contrast of a 1.6 MK cooling flow.

---

## 4. Results by scenario

### 4.1 Thermometer (`verify_altmodel_thermometer.png`)

T_min (K) for which a fully ionised O⁺/O⁺⁺ zone can reach the observed [O III]/Hβ, at [O II]/[O III] = 1.37 (range for 1.1–1.9):

| O/H | [O III]/Hβ ≥ 15 | ≥ 20 | ≥ 47 |
|---|---|---|---|
| 4.9 × 10⁻⁴ (A09) | 13,050 (12,600–13,850) | 14,300 | 19,800 (19,000–21,300) |
| 3.2 × 10⁻⁴ (ISM gas) | 15,000 | 16,700 | 24,100 |
| 9.8 × 10⁻⁴ (2× solar) | 10,700 | 11,600 | 15,100 |

The maximum [O III]/Hβ is 5.5 at 10⁴ K and 11.1 at 1.2 × 10⁴ K.

Wave-2 [N II] 5755/6583 = 0.14 is at or above the infinite-T limit. That reflects marginal 5755 (0.014 [0.005–0.026]) and differential-red-line systematics, so it gives no T constraint.

### 4.2 Photoionisation equilibrium (`verify_altmodel_photoion.png`)

| Spectrum | Max [O III]/Hβ, at any U | At [O II]/[O III] = 1.37: log U | [O III]/Hβ | T_eq (K) | 4363/5007 |
|---|---|---|---|---|---|
| 35 kK (DIG/ISRF) | 1.3 | −2.47 | 1.2 | 7,160 | 0.002 |
| 81 kK (WD J0736) | 9.1 | −3.34 | 5.2 | 9,790 | 0.006 |
| 150 kK | 13.9 | −3.39 | 9.4 | 11,930 | 0.011 |
| 0.14 keV hot gas | 14.0 | −3.31 | 9.9 | 12,240 | 0.011 |
| 28 / 50 eV precursor | 10.7 / 12.6 | −3.37 / −3.38 | 5.9 / 7.8 | 10,190 / 11,170 | 0.007 / 0.009 |

- The finer grid with Si gives WD 4.8, hot gas 8.6 and precursor 5.9.
- In PIE, He II/[O III] is 0.025 (WD), 0.062 (hot gas) and 0.051 (precursor).
- FUV is 49–53 CU per R, NUV 27–33, and [N II]/Hα 0.19–0.27: all too low.
- Charge-exchange sensitivity (81 kK, log U = −3): k_CX = 0 → 3 × 10⁻⁹ moves [O II]/[O III] from 0.11 to 0.86, but the maximum [O III]/Hβ only from 8.5 to 6.7. CX cannot rescue PIE.
- A stratified, radiation-bounded nebula can only lower the one-zone maximum, so these numbers are upper bounds.

### 4.3 Photon and ionisation-parameter budgets (`verify_altmodel_budget.png`, left)

**Density needed for the brightness.** Required n_H,rms on the 25.6-pc path for I = 1.4–4 R:
- 0.22–0.37 cm⁻³ (PIE-like gas, 1.2 × 10⁴ K, x(O⁺⁺) = 0.45);
- 0.10–0.17 cm⁻³ (4 × 10⁴ K, x = 0.6).

**Recombination budget.** The face-on 5007 photon flux is 1.07 × 10⁵ cm⁻² s⁻¹ for 2 R. Each recombination yields 1.7 (for [O III]/Hβ = 15) to 5.3 (for 47) 5007 photons. So 2.0–6.3 × 10⁴ ionising photons cm⁻² s⁻¹ must be absorbed.

**Sources** (Φ_H in cm⁻² s⁻¹; log U at n_rms / at the [O II] density; shortfall in U):

| Source | Φ_H | log U | Shortfall |
|---|---|---|---|
| WD J0736, coplanar (6.4 pc) | 1.35 × 10⁵ | −4.8 / −7.3 | 27× / 8.6 × 10³× |
| WD J0736, Silk at 282 pc (132 pc away) | 317 | −7.4 / −9.9 | 1.2 × 10⁴× / 3.7 × 10⁶× |
| Monogem EUV (F = n_e n_H Λ R/3; Λ(1.6 MK) = 1.11 × 10⁻²² from GS07; ⟨hν⟩ = 30–100 eV) | 340–6,200 (F = 0.5–3.0 × 10⁻⁷ erg cm⁻² s⁻¹) | −6.1 to −7.4 / ≤ −8.6 | 630–1.1 × 10⁴× / 2 × 10⁵–3.6 × 10⁶× |
| Generic field, generous Φ = 10⁶ | 10⁶ | −3.9 / −6.4 | 3.9× / 1.2 × 10³× |

- The WD's raw photon count could meet the recombination budget only if the WD and the sheet were co-distant and the sheet intercepted all of its output.
- The ionisation-parameter test is the binding one.

### 4.4 Photoionisation plus extra heating (`verify_altmodel_budget.png`, middle)

U is set so that [O II]/[O III] = 1.37, for the 0.14 keV spectrum (the 81 kK spectrum gives similar ratios but needs more heating):

| T (K) | [O III]/Hβ | 4363/5007 | 7325/3727 | FUV / NUV (CU per R) | [N II]/Hα | Heating / photo-heating |
|---|---|---|---|---|---|---|
| 1.5 × 10⁴ | 15 | 0.019 | 0.037 | 56 / 46 | 0.35 | 1.4 (81 kK: 2.9) |
| 2 × 10⁴ | 19 | 0.033 | 0.050 | 92 / 78 | 0.28 | 8 (13) |
| 4 × 10⁴ | 45 | 0.075 | 0.076 | 46 / 94 | 1.7 | 78 (351) |
| 6 × 10⁴ | 94 | 0.098 | 0.088 | 62 / 96 | 1.6 | 920 (7,400) |

- Reaching both [O III]/Hβ and 4363 needs 10²–10³ times the photo-heating.
- 7325/3727 then comes out ×3 too high.

### 4.5 Collisional families (`verify_altmodel_cooling.png`)

**GS07 cooling from 1.6 MK to 10⁴ K, no recombination tail:**

| Track | [O III]/Hβ | [O II]/[O III] | 4363 | 7325/3727 | C IV/[O III] (photons) | FUV / non-res. / NUV (CU per R) | f(5007) | Energy flux for 2 R (erg cm⁻² s⁻¹) |
|---|---|---|---|---|---|---|---|---|
| isochoric Z = 1 | 25.5 | 1.42 | 0.091 | 0.055 | 0.63 | 226 / 93 / 77 | 2.4 × 10⁻³ | 1.7 × 10⁻⁴ |
| isobaric Z = 1 | 25.5 | 1.72 | 0.097 | 0.057 | 0.72 | 254 / 102 / 86 | 2.2 × 10⁻³ | 2.0 × 10⁻⁴ |
| isochoric Z = 2 | 44.9 | 1.22 | 0.082 | 0.052 | 0.49 | 176 / 71 / 65 | 2.9 × 10⁻³ | 1.5 × 10⁻⁴ |
| isobaric Z = 2 | 48.0 | 1.57 | 0.090 | 0.055 | 0.60 | 202 / 77 / 72 | 2.7 × 10⁻³ | 1.6 × 10⁻⁴ |

- **Recombination tail.** With the tail (x(H⁺) = 0.86–0.91 at 10⁴ K in GS07), [O III]/Hβ falls to 3.4–6.4 and [N II]/Hα to 0.15–0.31.
- **Truncation earlier.** Stopping at 2–3 × 10⁴ K gives [O III]/Hβ = 40–90 and [O II]/[O III] = 0.6–1.2.
- **[Ne III]/[O III]** is 0.31–0.41 with GS07's enhanced Ne. Rescaled to Asplund+09 Ne/O it is 0.13–0.17 (observed 0.10).
- **FUV composition** (isobaric Z = 1): C IV 141, O IV] 29, 2γ 26, O III] 23, N IV] 13, Si IV 11, He II 7 CU per R.
- **TML energetics** (Ji, Oh & Masterson 2019, eq. 22: (5/2) P v_hot ≈ ΣL). The required 2.0 × 10⁻⁴ erg cm⁻² s⁻¹ means v_hot = 377 km/s at Monogem pressure. That is supersonic, i.e. a shock. At v_in = 10 km/s the layer gives 0.05 R.
- **X-ray contrast.** A cooling flow from 1.6 MK radiates 38 % of its enthalpy above 10⁶ K. With G = 18.7, its 1–1.6 MK emission would be 1.9 × 10³–2.3 × 10⁴ times the Monogem surface brightness (chords of 22–48 pc, n_e 3.1–7.2 × 10⁻³), using bolometric ratios. verify_xrayfuv finds no X-ray filament: the on-spine bump is ≤ few % of the level.

**Warm CIE layer.** [O II]/[O III] lies in 1–2 only at T = 5.8–6.1 × 10⁴ K. There:
- [O III]/Hβ = 84–101, 4363 = 0.096–0.098, 7325/3727 = 0.081–0.082;
- C IV/[O III] = 0.03–0.04 photons; FUV 86–92 CU per R (Si IV 43–46, O III] 23); NUV 84–93 (C III] 57–60);
- [N II]/Hα = 1.9–2.1, [S II]/Hα = 0.38–0.42;
- 65 % of Hβ is collisional; f(5007) = 2.5 %, so the heating needed is 1.7 × 10⁻⁵ erg cm⁻² s⁻¹ face-on.

**Conduction fronts.**

| | Static, P/k = 1.5 × 10⁴ | Evaporating, J = 1.2 × 10⁻²⁰–3.2 × 10⁻¹⁹ g cm⁻² s⁻¹ |
|---|---|---|
| Face-on I(5007) | 2.8 × 10⁻³ R (needs 38 fronts at G = 18.7) | 1.8 × 10⁻⁴–3.4 × 10⁻³ R |
| Thickness | 86 pc (0.65 pc at P/k = 2 × 10⁶) | — |
| 4363/5007 | 0.111 | 0.116–0.139 |
| [O II]/[O III] | 2.0 | 2.5–5.4 |
| [O III]/Hβ | 25 | 22–52 |
| C IV/[O III] (photons) | 1.6 | 3.6–16 |
| FUV (CU per R) | 541 | 970–5,300 |
| N(O VI) (cm⁻²) | 9.2 × 10¹² | 0.5–2.3 × 10¹³ |
| N(C IV) (cm⁻²) | 1.4 × 10¹² | 2.6–3.6 × 10¹² |

### 4.6 Fossil (`verify_altmodel_budget.png`, right)

Isochoric recombination with no heating:

| Start | n_H = 100 | n_H = 0.5 |
|---|---|---|
| T₀ = 2 × 10⁴ K (O⁺⁺ 0.7): 5007 e-fold time | 21 yr | 4.2 × 10³ yr |
| T₀ = 2 × 10⁴ K: matching epochs | never: 4363 ≤ 0.033 | never |
| T₀ = 4 × 10⁴ K: 5007 e-fold time | 34 yr | 6.9 × 10³ yr |
| T₀ = 4 × 10⁴ K: matching window | 9–21 yr | 1.5–4.2 × 10³ yr |

- "Matching" means [O III]/Hβ ≥ 15, [O II]/[O III] 1–2 and 4363 ≥ 0.05 at once.
- 7325/3727 = 0.068 in the matching window.
- Possible flash or heating events are far too old: the Monogem SN was 1.1 × 10⁵ yr ago (Knies+2024 via context.md), and Geminga's characteristic age is ~3 × 10⁵ yr.
- A low-density fossil cannot give 2 R, because that needs n_rms ≥ 0.2.

---

## 5. Claim → verdict

| Claim (source) | Verdict | My value |
|---|---|---|
| No catalogued hot star can photoionise it; WD J0736 supplies 4.9–67 % of the photons needed for 10 R (context) | **CONFIRMED (stronger reason)** | At the revised 2 R the photon count alone is no longer decisive (coplanar Φ_H = 1.35 × 10⁵ against 2–6 × 10⁴ needed if fully intercepted). But U is short by 27× (n_rms) to 9 × 10³× (n from [O II]), and by ≥ 10⁴× if the Silk is at 282 pc. Its PIE spectrum gives [O III]/Hβ ≤ 9 and 4363 = 0.006. |
| Photoionisation by Monogem hot-gas EUV "remains plausible, not quantified" (context) | **REFUTED** | Φ ≤ 6 × 10³ cm⁻² s⁻¹, log U ≤ −6.1 (needs −3.3); [O III]/Hβ ≤ 14 at best; He II/[O III] 0.06 > 0.026. |
| [O III] 4363/5007 = 0.056 → T_e ~ 3 × 10⁴ K (wave-1 DESI) | **CONFIRMED mapping; REVISED value** | 0.056 → 2.94 × 10⁴ K. Wave-2 0.083 (0.059–0.096) → 4.6 (3.1–5.8) × 10⁴ K. |
| [O II] 3729/3726 = 1.30 → n_e ~ 110 (wave-1 DESI) | **CONFIRMED** | 88 (10⁴ K) / 100 (1.25 × 10⁴) / 125 (2 × 10⁴) cm⁻³; 47–208 over subsets. Low-density limit 1.46 (CHIANTI 7.1.4). P/k = 2.8 (1.4–5.7) × 10⁶. |
| Filament is "Balmer-dark and collisionally excited" (provisional) | **CONFIRMED** | PIE max [O III]/Hβ ≤ 14 at any U for any spectrum tried, ≤ 10 at the observed [O II]/[O III]. T ≥ 1.3 × 10⁴ K is required; T_e(O III) = 4.6 × 10⁴ K. |
| [O III]/Hβ = 47 (> 24 at 97.5 %) (verify_desi, input) | **PLAUSIBLE, but adopt ≥ 15** | See caveat 1 below. |
| FUV floor of a 10–13 R shock from O III] + He II ≈ 210–300 CU (imaging_mw) | **CONFIRMED (order)** | My O III] 1663 alone gives 21–23 CU per R(5007) in 10⁵ K cooling gas. |
| 1 R of C IV ≈ 223 CU (verify_xrayfuv) | **CONFIRMED** | 225 CU (W_eff 354 Å, own code). |
| Conduction fronts too faint, ≤ 3 × 10⁻⁴–1.6 × 10⁻³ R face-on (verify_xrayfuv, from BBF90) | **CONFIRMED, independent model** | 1.8 × 10⁻⁴–3.4 × 10⁻³ R face-on (own NEI evaporation and static fronts). |

**Caveat 1, on DESI [O III]/Hβ = 47.**
- In the same verify_desi file, Hγ/[O III] = 0.031 against Hβ/[O III] = 0.021. That makes Hγ/Hβ = 1.5, against 0.47 for case B.
- The null-line coefficients near Hβ and Hγ are of the same size: null4885 = +0.011 (0.008–0.017), null4330 = +0.020, null4354 = −0.020.
- So Hβ carries a ≈ ±0.01–0.02 systematic. If Hγ is right, [O III]/Hβ ≈ 15. If Hβ is right, the value is 47.
- No conclusion here depends on the choice between 15 and 47.

---

## 6. New findings

1. **Two-temperature structure.** T([O III]) ≈ 4.6 × 10⁴ K against T([O II]) ≈ 1.2 × 10⁴ K, at a pressure of ~3 × 10⁶ K cm⁻³. Together these are a near-unique fingerprint of a radiatively cooling, high-pressure (shocked) layer. They exclude every single-zone mechanism here, including ad hoc heating of photoionised gas.
2. **Cooling flows match the optical spectrum when truncated.** NEI cooling flows from ≥ 10⁵ K reproduce the optical ratios ([O III]/Hβ 25–48, [O II]/[O III] 1.2–1.7, 4363 0.08–0.10) only if H has not yet recombined. A complete flow gives [O III]/Hβ = 3–6. This independently supports the "incomplete shock" reading of xrayfuv.
3. **The FUV level argues for collisional excitation.** Photoionised gas gives 28–53 CU per R, too little. Collisional gas at T ≳ 5 × 10⁴ K gives ≥ 86 CU per R from O III], N IV], Si IV and C IV.
   - Cooling flows give 176–254 CU per R, of which 70–102 is non-resonant.
   - The observed 68–92 therefore needs C IV (and partly Si IV) suppressed ≳ 3–5×, consistent with verify_xrayfuv.
4. **A warm CIE layer at ~6 × 10⁴ K fits much of the spectrum.** It reproduces [O III]/Hβ, [O II]/[O III], 4363 and even the FUV level. Only [O II] 7325/3727 (×3), NUV (×2) and [N II]/Hα (×2) reject it. Better 7325 and NUV measurements are the discriminators against any "heated warm layer" idea.
5. **Hot-gas cooling flows / TMLs cannot be invisible in X-rays.** At this brightness a TML or cooling flow at the Monogem boundary would be a soft-X-ray filament 10³–10⁴× brighter than the Monogem shell. At realistic entrainment speeds it would give only ~0.05 R of [O III].
6. **The emitting layer is extremely thin.** At the [O II] density, 2 R needs only ~40–230 AU of emitting path. The 0.28° width must come from the projection of a thin, rippled or curved sheet, not from a 1.4-pc-thick emitting slab.

---

## 7. Decisive observational tests

Already-available data are marked "done". For the rest, the single most decisive future observation is listed.

| Scenario | Single most decisive test |
|---|---|
| 1 Conduction front | Done: 4363/5007 (predicted 0.11–0.14) plus the [O II] doublet (predicted 1.46) plus GALEX FUV (predicted ≥ 540 CU per R). Future: FUV spectroscopy of C IV/O III] (≥ 1.6 photons). |
| 2 TML | Done: no eRASS1 X-ray filament, and the [O II] doublet. Future: O VI 1032/1038 at ~2 photons per [O III] photon (~4 R, ~3 × 10⁵ LU edge-on), against a typical diffuse 10³–10⁴ LU. |
| 3 NEI cooling of hot gas | Done: eRASS1 soft-X-ray filament predicted at ≥ 10³× Monogem, not seen. |
| 3w Warm heated layer | Deep [O II] 7320/7330 and NUV C III] 1909: predicted 7325/3727 = 0.08 and NUV ≈ 90 CU per R, against 0.026 and ~45 observed. |
| 4a–c Photoionisation | Done: 4363/5007 (predicted ≤ 0.012). Future, independent: a Balmer line at [O III]/Hβ ≤ 10 with Hα ≥ 0.2 R per R of [O III] in high-S/N slit spectra. Also, a WD-driven nebula would show an ionisation gradient with distance from WD J0736 (not tested here). |
| 4d Photoionisation plus heating | [O II] 7325/3727 and the [O II] doublet, as for 3w. |
| 5 Precursor | [O II] doublet: pre-shock gas predicts 1.46. Plus 4363 ≤ 0.01. |
| 6 Fossil | Secular fading: e-folding 20–35 yr at the doublet density, i.e. ~3 %/yr. Compare calibrated multi-epoch narrow-band photometry (SDSS 2002–07 precovery, NSNS 2023–24, a new epoch). A constant brightness over ~20 yr excludes it. |

For the surviving shock picture, one physical test stands out: Hα/Hβ line widths from the [O II] zone (≈ 1.2 × 10⁴ K) against [O III]. The DESI [O III] intrinsic FWHM ≤ 33 km/s, which is consistent with T ≤ 10⁵ K plus modest turbulence.

---

## 8. Limitations

- **One-zone, optically thin PIE.** There is no radiative transfer and no stratification. This only lowers the attainable [O III]/Hβ, so the PIE exclusion is conservative. The hot-gas and precursor spectra are crude exp(−E/kT)/E shapes, not line-dominated plasma spectra.
- **Approximate O–H charge exchange.** The rate is approximate; its sensitivity was tested. CX for other ions (N, S) is neglected, which affects [N II]/Hα and [S II]/Hα in PIE by an unknown factor ≲ 2.
- **NEI cooling uses published GS07 fractions.** Their atomic data, plus my CHIANTI v7.1.4 emissivities: a mixed dataset. Cooling-flow line ratios are uncertain by ~20–30 %.
- **Ne abundance.** GS07 Ne is enhanced (Ne/O = 0.43); I rescaled [Ne III] to A09 Ne/O = 0.17.
- **The recombination tail is bracketed, not modelled.** Gas below 10⁴ K needs photoionisation transfer.
- **The TML is represented by an isobaric cooling flow plus the Ji+2019 energy equation, not by a simulation.** Real TMLs weight intermediate T more, which makes C IV/[O III] larger. Conduction fronts are planar and classical (capped at saturation), with no magnetic suppression, μ = 0.6, and CIE ions in the static case.
- **The X-ray contrast is a bolometric ratio.** Band fractions differ by ≲ 2×, and absorption to 280 pc is τ ≲ 1 at 0.25 keV. Neither changes a 10³–10⁴ factor.
- **Brightness and energy conclusions scale as 1/G**, with G = 18.7 for a curved sheet with R_c = 60 pc. An unknown geometry could change G by ×2–3. That does not rescue any excluded case, since the shortfalls are ≥ 27× (and line ratios do not depend on G).
- **CHIANTI v7.1.4 is not the latest version.** The [O II] doublet low-density limit (1.46) differs by a few % between datasets (1.42–1.50), so n_e is uncertain by ~±30 %.
- **Inputs are taken from teammates.** The [O III] calibration (DESI vs NSNS) affects the FUV per R by up to ×1.35–2. I did not re-measure spectra. The only exceptions are the Hγ/Hβ and null-line checks above, read from verify_desi's file.

---

## 9. Files

**Scripts**
- `scripts/verify_altmodel_atomic.py`, `_lines.py`, `_common.py`: toolkit and inputs.
- `scripts/verify_altmodel_thermo.py`, `_photoion.py`, `_cooling.py`, `_budget.py`, `_summary.py`: models.
- `scripts/verify_altmodel_fetch.py`, `_extract.py`, `_extract2.py`: data.

**Derived**
- `data/derived/verify_altmodel_thermo.json`
- `data/derived/verify_altmodel_photoion.json`
- `data/derived/verify_altmodel_coolion.npz`
- `data/derived/verify_altmodel_linetab.npz`
- `data/derived/verify_altmodel_cooling.json`
- `data/derived/verify_altmodel_budget.json`
- `data/derived/verify_altmodel_scenarios.json`

**Figures**
- `figures/verify_altmodel_thermometer.png`
- `figures/verify_altmodel_photoion.png`
- `figures/verify_altmodel_cooling.png`
- `figures/verify_altmodel_budget.png`
- `figures/verify_altmodel_scorecard.png`

**Downloads** (logged with URL, sha256 and date in `data/raw/verify_altmodel/MANIFEST.jsonl`; all accessed 2026-11-16)
- CHIANTI v7.1.4 database, 18 MB.
- GS07 tables 2–15 and the paper PDF.
- Dere 2007 Table 29 (VizieR CfA mirror).
- arXiv PDFs: Gnat, Sternberg & McKee 2010; Kwak & Shelton 2010; Ji, Oh & Masterson 2019; Fielding+2020; Chen, Fielding & Bryan 2023.
- The Si/Mg CHIANTI directories were extracted from the same tarball, so no new download was needed.
