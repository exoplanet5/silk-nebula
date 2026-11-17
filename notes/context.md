# context: the astrophysical neighbourhood of the Silk Nebula

Topic tag `context` · written 2026-11-16 · all downloads accessed 2026-11-16.

The filament geometry comes from the geometry teammate's `data/derived/silk_trace.json` v1.0 (not preliminary). Its E-strand spine has 66 points, RA 114.28–117.04°, Dec +12.00 to +28.25°, l 191.6–208.1°, b +17.5 to +21.7°. The trace also provides the W strands, the northern Hα patch polygon and the curvature fits. "Measured" means computed here from catalogue or map data. "Literature" means copied from the cited paper.

## Summary

1. **The Silk sits on the east/north rim of the Monogem soft-X-ray emission.** I averaged ROSAT RASS maps along curves parallel to the spine (smoothed to FWHM 1°). In 1/4 keV (R1+R2), the band 1–5° west (inside Monogem) is brighter than the band 1–5° east by these factors:
   - Dec 12–16: 1.14 ± 0.03
   - Dec 16–20: 1.44 ± 0.07
   - Dec 20–24: 1.81 ± 0.06
   - Dec 24–28: 1.02 ± 0.08 (no edge)

   The half-level X-ray edge falls at +0.31°, −0.26° and +0.20° from the spine in the first three bands. In 3/4 keV (R4+R5) the west/east ratios are 1.09–1.26. This is consistent with the teammate's eRASS1 result (`notes/imaging_mw.md`: 0.2–0.6 keV W/E 1.32, edge 1.03 ± 0.20° inside). It also matches Fesen+2026 §6.2.1, which never names the Silk but describes "the nearly vertical 15 degree long group of [O III] filaments" that "lie at the easternmost extent of the Monogem's soft X-ray emission".
2. **Only the southern part of the Silk is concentric with Monogem.**
   - South of the kink (Dec 12–17.5, 23 points) the spine sits at a nearly constant 12.10 ± 0.12° (rms) from the Knies+2024 two-SNR Monogem centre (104.15, +15.30; R = 11.0 ± 1.1°), and at 11.36 ± 0.16° from PSR B0656+14. That is the Monogem outer edge.
   - North of Dec ≈ 17.6° the arc bends away, reaching 16° from those centres (rms 1.0–1.3°).
   - Relative to the single-SNR centre (107.08, +13.50; R = 11.97°), the spine distance rises from 9.1° to 16.2° (rms 2.2°). 56 % of the spine lies inside that circle, for Dec ≤ 21.
   - Caveat: the southern section is short (5.5°) and nearly straight, so a low rms alone does not single out one centre. The informative part is that its mean distance equals the Monogem radius.
3. **Centre of curvature.** It is at (RA, Dec) = (94.5, +16.76), i.e. (l, b) = (194.26°, +0.44°), with R = 21.3° (bootstrap 20.5–22.4°). The single-circle fit is poor (χ²/dof = 3744/63) because the arc is S-shaped. This centre is 12.6° from the Monogem centre, 0.8° from the Origem Loop centre (a projection coincidence: the Origem Arc is at ~1.7 kpc), 3.9° from present-day Geminga and 10.4° from PSR B0656+14.
4. **No pulsar birth locus coincides with the centre of curvature.** I back-traced 8 pulsars with d < 1 kpc and a measured proper motion, in the LSR frame over 0–1 Myr. The Monte Carlo has 3000 draws over distance, PM and an unknown v_r drawn from N(0, 265 km/s).
   - PSR B0656+14 at τ_c = 111 kyr: born at (l, b) = (200.6, +6.9) ± 0.1°. That is 9.05° from the centre, and 0 % of tracks ever enter the centre region.
   - Geminga at τ_c = 342 kyr: (198.5, −11.9), 13.0° (9.2–20.6°) away. Its track does cross the centre region, but only at t = 0.06–0.08 Myr ago, not near τ_c.
   - The E spine is formally concentric with *present-day* Geminga (17.56 ± 0.22°). This is not physical: Geminga moves by that rms in ~4.5 kyr, and its spin-down energy over τ_c (Ėτ ≈ 3.5×10⁴⁷ erg) can inflate only an R ≈ 17 pc bubble, against the ~75 pc needed.
5. **No known hot star can photoionize the [O III] filament.** The best case is WD J0736+1444 (= UCAC4 524-042847): DOZ, Teff 81.3 ± 1.6 kK, d = 150.5 pc, Q(H) = 6.6×10⁴⁴ s⁻¹. It lies 2.44° (6.4 pc projected) from the spine, inside its Strömgren radius of 12.7 pc at n = 0.1 cm⁻³. Even so it supplies at most 4.9 % of the minimum [O III] λ5007 photon rate with a realistic covering factor, and 67 % even if every ionizing photon reached the filament. All other ≥40 kK stars within 15° and 1 kpc fall below 2.0 % (realistic) and 72 % (4π). There is no PN within 5° of the spine, and no proper-motion "wake" aligns with it.
6. **Dust and the Local Bubble.** Foreground dust is tiny: Edenhofer+2024 E(<300 pc) = 0.0050–0.0061 and E(<1.25 kpc) = 0.023–0.048 (ZGR23 units ≈ E(B−V)) along the spine. The Local Bubble wall toward P1–P4 is at 161–171 pc (O'Neill+2024 shell peak, σ 12–18 pc; inner 133–136, outer 196–201 pc). It is very thin (shell A_G = 0.005–0.007 mag).
   - Wu+2026's fiducial 160 pc therefore coincides with the Local Bubble wall.
   - If the Silk is part of the Monogem shell (D = 288 pc, the PSR parallax), the geometry puts it at 250–320 pc (tangent point 277–284 pc), i.e. 1.6–2× farther. Wu+2026's linear size and energy would then scale by ~1.8 and ~(1.8)^2.5–3 (≈4–6×) respectively; that scaling is not derived here.
7. **Northern Hα patch.** The patch where the [O III] arc meets Hα+[SII] (NSNS O200.0a; Fesen+2026's "NE filaments", Hα/[SII] ≈ 1.0 per Ziegenbalg) lies 0.35° from the Planck cold clump PGCC G199.09+19.95. Fesen+2026 also place it at the western tip of a long Hα cloud extending east. No Sh2/LBN/H II region, reflection nebula or SNR is catalogued within 3° of the spine or the patch. All SIMBAD "HII" entries in the box belong to background galaxies.
8. **New catalogue coincidence, not verified.** Two outer-Galaxy HI shell candidates are centred on the Silk:
   - GS202+17+019 (Suad+2014): ellipse centre (202.4, +16.9), semi-axes 5.7 × 4.5°, v₀ = +18.6 km/s, HI minimum at (203.3, +19.1), which is 0.7° from the Silk centre.
   - GS202+19−019 (Ehlerová & Palouš 2005): (202.3, +18.8), 6.5 × 5.5°, v = −18.8 km/s.

   The spine crosses their interiors rather than tracing a rim. GALFA-HI shows no clean expanding shell there (`notes/imaging_mw.md`).

## 1. Monogem Ring

| quantity | value | source (type) |
|---|---|---|
| Discovery as X-ray ring | HEAO-1, centre (l, b) ≈ (200°, +8°), ~25° | Nousek+1981, *via* Knies+2024 §1 and Thorsett+2003 (literature, secondhand) |
| ROSAT model | kT 0.1–0.3 keV; Sedov age 86 kyr; at 300 pc R = 66 pc, n = 5×10⁻³ cm⁻³ | Plucinsky+1996 *via* Thorsett+2003 §1 |
| ROSAT explosion energy | 1.9×10⁵⁰ erg (as quoted by Thorsett+2003); ~0.2×10⁵⁰ erg (as quoted by Knies+2024); 0.2–1.0×10⁵⁰ erg (Fesen+2026) | secondhand; **the quotes disagree**, original not read |
| PSR B0656+14 parallax | 288 +33/−27 pc; PM 44 mas/yr | Brisken+2003 *via* Thorsett+2003 (arXiv astro-ph/0306462) |
| PSR–SNR birth | pulsar moved ~1.25° in 10⁵ yr, from (l, b) ≈ (200.5, +7.2); 9.2° circle centred on the pulsar fits the ring | Thorsett+2003 |
| eROSITA single SNR | centre (107.08°, +13.50°) = (l, b) (202.7, +9.8); R = 11.97° = 63 pc at 300 pc; kT 0.11–0.17 keV; n f^−½ = 3.1–7.2×10⁻³ cm⁻³; age ≈ 1.1–1.2×10⁵ yr; E₀ ≈ (0.1–0.2)×10⁵¹ erg (Sedov), (0.2–0.4)×10⁵¹ (cloudy ISM) | Knies+2024 A&A, arXiv:2401.17289, Table 4 / §7 |
| eROSITA two-SNR | Monogem (104.15, +15.30), R = 11.00° (58 pc); plus G205.6+12.4 (110.33, +12.94), R = 8.82° (47 pc at 300 pc), kT 0.28–0.45 keV, age ≈ 5×10⁴ yr | Knies+2024 Table 4 |
| FUV (FIMS/SPEAR) | half-ring at low latitude along the distorted X-ray edge; strong C IV λλ1548,1550 where the ring meets the Gemini Hα ring; C IV bright only in a small high-latitude part of the ring | Kim+2007 ApJ 665 L139, **secondhand** (Fesen+2026 §6.2.1, Reich+2020 §1, web-search abstract summary); full text not obtained |
| Optical | [O III] filaments along the E, S and NW limbs. The 15° E/N group (= Silk) lies at the X-ray edge. NE Hα+[SII] filaments sit at the north end of the [O III] filaments and are read as a slower shock in denser gas. Earlier filaments: W06 (Weinberger+2006, 07:28:38 +05:09:55, Hα, no [O III]) and R20 (Reich+2020, G203.7+11.5, 4.5° Hα filament toward the centre; B∥ = 26 ± 5 µG, n_e = 1.6 cm⁻³, RM −86 rad m⁻² at 300 pc) | Fesen+2026 arXiv:2608.26666 §6.2; Reich+2020 arXiv:2007.08221; Weinberger+2006 astro-ph/0511546 |

**Where the Silk lies (measured).** Data: `data/derived/context_rosat_spine.csv` and `context_rosat_profiles.csv`; figure `figures/context_rosat_profiles.png`.

| Dec band | R12 W/E (1–5° W vs 1–5° E) | R12 peak offset | R12 half-level edge | R45 W/E |
|---|---|---|---|---|
| 12–16 | 1.144 ± 0.034 | −2.0° | +0.31° | 1.085 ± 0.045 |
| 16–20 | 1.442 ± 0.070 | −2.5° | −0.26° | 1.224 ± 0.045 |
| 20–24 | 1.807 ± 0.060 | −3.5° | +0.20° | 1.257 ± 0.047 |
| 24–28.25 | 1.016 ± 0.076 | — | (+2.8°, no edge) | 1.099 ± 0.047 |

The errors are the scatter along the band divided by √N_indep, with N_indep = band length / 1°. They are not formal and ignore ROSAT systematics. The 1/4-keV map is sensitive to absorption. However, foreground dust within 300 pc is the same on both sides: E(<300) is 0.0047 at 4° E and 0.0054 at 4° W of P1. The teammate also showed the GALFA-HI column difference changes the 0.2–0.6 keV transmission by only ≲ 15 %.

**Concentricity and rim distance (measured).** Script `scripts/context_concentricity.py`, output `data/derived/context_concentricity.json`, figure `figures/context_concentricity.png`.

| centre | mean spine distance | rms | range |
|---|---|---|---|
| Silk curvature centre (E-strand fit) | 21.4° | 0.15° | — |
| Origem Loop centre (SIMBAD 94.13, +16.07) | 21.84° | 0.163° | 21.55–22.29° |
| Geminga now | 17.56° | 0.222° | 17.28–18.37° |
| Monogem 1-SNR (Knies+24) | 11.74° | 2.19° | 9.08–16.23° |
| Monogem 2-SNR W | 13.27° | 1.18° | 11.98–16.03°; **south of the kink 12.10 ± 0.12°** |
| G205.6+12.4 | 9.83° | 3.08° | — |
| PSR B0656+14 now | 13.04° | 1.61° | 11.14–16.51°; **south 11.36 ± 0.16°** |
| PSR B0656+14 at τ_c | 14.32° | 1.44° | 12.59–17.43°; south 12.80 ± 0.16° |
| Geminga at τ_c | 32.7° | 1.1° | — |

If the Silk is on the Monogem shell (sphere of angular radius 11.97°), the line-of-sight intersections give these distances:

| shell distance D | near side (median) | far side (median) | tangent point |
|---|---|---|---|
| 261 pc | 224 pc | 290 pc | 251–258 pc |
| 288 pc | 248 pc | 320 pc | 277–284 pc |
| 321 pc | 276 pc | 357 pc | 308–317 pc |

**Interpretation (derived).** The Silk is the eastern/northern boundary of the X-ray-bright interior.
- South of the Dec ≈ 17.6° kink it follows a Monogem-centred circle at the ring's radius: 12.1° from the two-SNR centre, 11.4° from the pulsar.
- North of the kink it bends north-west, away from the Monogem centre, and flattens to a ~21° radius of curvature.

So the Silk is a Monogem limb only in the south. Natural readings:
- (a) The Monogem (or the hotter G205.6+12.4) shock is flattened against denser gas to the east and north.
- (b) The Silk marks the interface between Monogem's hot gas and the Local Bubble wall or neighbouring warm gas, possibly a photoionized or conductive mixing layer (cf. `notes/imaging_mw.md`).

Either way, the rim distance is 250–320 pc, not 160 pc.

## 2. Pulsars within 30°

Source: ATNF Pulsar Catalogue **v2.9.0** text database (`psrcat.db`, header `#CATALOGUE 2.9.0`). It was downloaded directly from atnf.csiro.au, which *was* reachable on 2026-11-16 despite the brief. Only the text database was parsed; the C code was never compiled. The VizieR B/psr mirror (2017 snapshot) was fetched as a cross-check. 93 pulsars lie within 30° of the Silk centre and 17 of them have d < 1 kpc. Table: `data/derived/context_pulsars_30deg.csv`.

**d < 1 kpc, with proper motion (back-traced).**

| PSR | sep. from Silk | τ_c | Ė (erg/s) | DM | PM (mas/yr) | distance |
|---|---|---|---|---|---|---|
| B0656+14 (J0659+1414) | 11.5° | 111 kyr | 3.8e34 | 13.95 | (44.1, −2.4) | PX 3.47 ± 0.36 mas → 288 pc (Brisken+03) |
| Geminga (J0633+1746) | 17.4° | 342 kyr | 3.2e34 | 2.89 | (142.2, 107.4) ± 1.2 (Faherty+07); ATNF (138, 97) ± 4 | PX 4.0 ± 1.3 mas (Faherty+07) → 251 (188–373) pc; ATNF PX 6.4 ± 1.8 → 157 pc run separately |
| B0823+26 | 13.8° | 4.9 Myr | 4.5e32 | 19.48 | (63.0, −96.8) | 500 pc |
| B0834+06 | 16.0° | 3.0 Myr | 1.3e32 | 12.86 | (2, 51) | DM 191 pc |
| J0742+4110 (MSP) | 24.8° | — | 8.5e33 | 20.81 | (−14.0, −6.0) | DM 523 pc |
| J0613−0200 (MSP) | 29.4° | — | 1.3e34 | 38.78 | (1.8, −10.4) | 780 pc |
| B0609+37 | 29.4° | 79 Myr | 8.9e31 | 27.15 | (5.1, 4.9) | DM 534 pc |
| B0943+10 | 29.8° | 5.0 Myr | 1.0e32 | 15.32 | (−38, −21) | DM 892 pc |

**d < 1 kpc, no proper motion (cannot be traced).** The table gives separation from the curvature centre.

| PSR | sep. from curvature centre | note |
|---|---|---|
| J0611+1436 | **2.7°** | τ_c 1.07 Myr, Ė 8.0e33, DM distance 0.81 kpc; worth a PM measurement |
| J0630+19 | 4.1° | |
| J0609+2130 | 5.1° | recycled |
| J0658+2936, J0653+0443, J0623+0220, J0811+3729, J0803−0942, J0613+3731 | ≥ 14° | |

**Back-trace method.** Script `scripts/context_pulsar_backtrace.py`.
- Present position minus v_LSR·t, for t = 0–1 Myr, in a frame co-moving with the LSR. Solar motion (U, V, W) = (11.1, 12.24, 7.25) km/s (Schönrich+2010). Galactic differential rotation is neglected (≲ 5 km/s for d < 350 pc).
- 3000 MC draws:
  - distance: Gaussian in parallax, or log-normal σ = 0.4 about the YMW16 DM distance;
  - PM: Gaussian within its errors;
  - v_r ~ N(0, 265 km/s) (Hobbs+2005 1-D σ).

**Centre-of-curvature regions.** Data: `data/derived/context_curvature.json`, `context_curvature_regionA.npz`.
- **Region A** (formal) is the Δχ² ≤ 6.18 region of a single circle fitted to the E spine, with σ_int = 0.139° added so that χ²/dof = 1. It covers RA 91.5–94.5°, Dec 16.0–16.5°, R 21.4–24.3°.
- **Region B** (systematic) is the union of the geometry teammate's 16–84 % bootstrap boxes:

| fit | centre RA | centre Dec | R |
|---|---|---|---|
| E_full | 93.4–95.4° | 16.5–17.0° | 21.3° |
| E_north | 90.6–96.3° | 14.8–17.1° | 21.8° |
| E+W | 84.7–90.9° | 14.6–16.2° | 27.3° |

The southern part alone (Dec ≤ 17.5) is nearly straight (R ≈ 46°, 22–71°) and leaves the centre unconstrained.

**Back-trace results.** Data: `data/derived/context_pulsar_coincidence.json`, `context_pulsar_tracks.csv`, `context_pulsar_birth_mc.csv`.

| pulsar | birth locus at τ_c (l, b) | sep. from curvature centre | fraction of tracks entering region A / B | min. track separation (16/50/84 %) |
|---|---|---|---|---|
| B0656+14 | (200.62, +6.89), 16–84 %: l 200.56–200.67, b 6.72–7.01; d 294 (252–338) pc | 9.05° (8.89–9.16°) | 0 / 0 | 3.73/3.95/4.87° at t ≈ 0.8 Myr |
| Geminga (Faherty) | (198.5, −11.9), 16–84 %: b −19.3 to −8.1 | 13.0° (9.2–20.6°) | 0 / 1.0, but only at t = 0.06–0.08 Myr; 0 % within 0.5–1.5 τ_c | 1.6° |
| Geminga (ATNF) | (198.1, −11.5), 16–84 %: b −27.6 to −6.2 | 12.6° (7.3–30.1°) | B within 0.5–1.5 τ_c: 0.5 % | — |
| B0609+37 | — | — | A 0.5 %, B 4 %, track min 17.9° (median) | — |
| others | — | — | 0 % | ≥ 16° |

For comparison, Knies+2024 and Thorsett+2003 put the B0656+14 birthplace near (200.5, +7.2), which agrees with the result here.

**Shell test.** Suppose a sphere centred on the birth site at τ_c had the Silk on its limb. For B0656+14 the spine would then sit at 14.3 ± 1.4° (rms) from the centre, i.e. R ≈ 73 (65–84) pc with the Silk at 287 (252–331) pc. The rms (1.4°) is 10× the circle-fit rms, so the arc is not centred on the B0656+14 explosion site either. That is consistent with the non-concentricity of item 2 above.

**Verdict.** No back-traced pulsar birth locus coincides with the Silk's centre of curvature. The fitted centre lies near the Galactic plane, toward Gem OB1, IC 443 (5.8°) and the Origem Loop (0.8°). I see no physical reason for that alignment, and a 21° circle is a poor model anyway (S-shaped arc).

## 3. SNRs, candidates, loops and shells (l 175–235, b −15…+40)

Table: `data/derived/context_snrs.csv`. It holds 14 Green (2025) SNRs (VizieR VII/297, 310 SNRs) plus 16 literature entries. Selected entries:

| object | (l, b) | size | distance | min. sep. from spine | notes |
|---|---|---|---|---|---|
| Monogem Ring (G203+12) | (202.7, +9.8) | 23.9° | ~300 pc | 9.1° | 56 % of spine inside |
| G205.6+12.4 (Knies+24) | (204.6, +12.4) | 17.6° | 300–500 pc | 5.9° | 44 % of spine inside; hotter, ~50 kyr |
| G190.4+12.5 | (190.4, +12.5) | 6.6° | >1.5 kpc (Knies) | 9.3° | X-ray + [O III] (Fesen) |
| G191.4+11.1 | (191.3, +11.1) | 3.2 × 2.2° | — | 10.6° | new optical candidate (Fesen) |
| G195.9+2.2 | (195.8, +2.1) | ~1° | — | 19.0° | new optical SNR (Fesen) |
| G205.7−1.7 | (205.6, −1.7) | 1.4° | — | 19.3° | optical candidate (Fesen) |
| G206.6+6.1 (= G206.7+5.9) | — | 3.6° | ~440 pc (Gao+22 via Fesen) | 11.4° | |
| G203.1+6.6 | — | 2.5° | — | 11.9° | |
| Monoceros G205.5+0.5 | — | 3.7° | ~1.5 kpc | 17.4° | |
| PKS 0646+06 | — | — | — | 15.2° | |
| IC 443 | — | — | — | 18.9° | |
| S147 | — | — | — | 26.0° | |
| Gemini Hα ring | (191.5, +5.0) | 5.7° | 200–350 pc | 16.7° | wind-blown bubble; Kim+07, Knies+18 via Fesen |
| R20 (G203.7+11.5) | — | 4.5° long | — | 7.1° | |
| W06 | (212.5, +10.6) | — | — | 8.1° | |
| Origem Loop | (194.7, +0.4) | 6.7° | Origem Arc ~1.7 kpc (Gao & Han 2013, arXiv:1212.1523 abstract) | — | non-thermal arc in the N, H II regions in the S |
| λ Ori complex | (195.1, −12.0) | — | ~400 pc | — | south of the map's main area |
| GS202+17+019 | — | 5.7 × 4.5° | 2.5 kpc (rotation-curve kinematic, meaningless toward the anticentre) | 2.9° | Suad+2014; v₀ +18.6, Δv 15.5 km/s, Reff 221 pc |
| GS202+19−019 | — | 6.5 × 5.5° | — | 1.1° | Ehlerová & Palouš 2005; v −18.8, Δv 21.6 km/s, quality 3 |
| GS208+15+003 | — | 3.5 × 4.0° | — | 2.7° | Ehlerová & Palouš 2005; v +2.5 km/s |

Neither Orion–Eridanus nor Barnard's Loop is near the Silk: both lie at b ≲ −10°, more than 30° from the filament. I did not quantify them further. Gem OB1 (the molecular cloud) is at (94.75, +23.0), 6.2° from the curvature centre and ~20° from the Silk. No catalogued SNR is closer to the spine than the Monogem/G205.6 outlines, and none other than those two contains any part of the spine.

**HI shells.** The two shells catalogued at ±19 km/s at nearly the same position might be the approaching and receding sides of one structure. That is speculative and unverified here. The Silk [O III] velocity (V_LSR +11…+16 km/s) lies within the GS202+17+019 velocity range (18.6 ± 7.8 km/s). These catalogues were built for the outer Galaxy, and their distances should not be used.

## 4. Hot ionizing sources (within 15°, d < 1 kpc)

**Inputs.**
- SIMBAD types O/B0–B2/sdO/DO/DAO/PG1159/WR, WD/hot-subdwarf/PN object types (25° cone).
- Gentile Fusillo+2021 Gaia EDR3 WDs (VizieR J/MNRAS/508/3877): P_WD > 0.75, plx > 1 mas, Teff_H > 20 kK, 20° cone.
- Culpan+2022 known hot subdwarfs with spectroscopic Teff/log g (J/A+A/662/A40/knownhsd).
- Hand-checked spectroscopic values for six stars:
  - Reindl+2023 Table A.4 (arXiv:2307.03721), checked in the text;
  - SIMBAD mesFe_h, checked by TAP (`data/raw/context_simbad_verify/mesfeh_litstars.csv`).

**Method.**
- Blackbody Q(H) and Q(>35.1 eV). The radius comes from log g with M = 0.6 (WD) or 0.47 M☉ (sd), or from a literature value.
- Strömgren radius for n = 0.1, 0.3 and 1 cm⁻³, with α_B = 2.59×10⁻¹³ cm³ s⁻¹.
- [O III] photon budget. The required rate is N_req = 10⁶ × I[R] × A, using I(5007) = 10 R (the low end of Wu+2026's 10–13 R) and A = 3 deg² (≈ 0.3° × 10°) at the star's distance. The available rate is 2 × Q(H) × f_cov, where 2 is a generous upper bound on λ5007 photons per photoionization and f_cov is the covering fraction of a 1° × 10° strip at the star's projected separation. An earlier version of the script used 10⁷ per rayleigh; I fixed this to 1 R = 10⁶/(4π) ph cm⁻² s⁻¹ sr⁻¹.

Tables: `data/derived/context_hotstars.csv` (658 stars) and `context_hotstars_top.csv` (45 with Teff ≥ 40 kK, d < 1 kpc).

| star | type | Teff (source) | d (pc) | Q(H) (s⁻¹) | Q(>35 eV) | R_S at n = 0.1 / 1 (pc) | sep. from spine | budget (f_cov / 4π) |
|---|---|---|---|---|---|---|---|---|
| WD J0736+1444 (UCAC4 524-042847) | DOZ | 81.3 ± 1.6 kK (Reindl+23) | 150.5 ± 1.0 | 6.6e44 | 1.1e44 | 12.7 / 2.7 | 2.44° (6.4 pc) | 0.049 / 0.67 |
| Ton 835 | sdO | 51.2 kK (Culpan) | 945 | 1.0e46 | 3.5e44 | 31.7 / 6.8 | 2.39° (39 pc) | 0.020 / 0.26 |
| GALEX J074507.5+144939 | DO: | 80.0 kK (Vincent+24, Gaia XP) | 594 | 4.6e44 | 7.5e43 | 11.3 / 2.4 | 0.28° (2.9 pc) | 0.015 / 0.030 |
| WDJ073957.35+091439.23 | WD | 92.9 kK (GF21 photometric, unreliable) | 427 | 2.1e45 | 4.9e44 | 18.8 / 4.1 | 3.15° | 0.014 / 0.27 |
| [CW83] 0711+22 | sdOp | 41.1 kK (Culpan) | 372 | 4.4e45 | 4.9e43 | 23.9 / 5.2 | 6.75° | 0.011 / 0.72 |
| WDJ074532.90+225914.34 | WD | 57.0 kK (GF21) | 585 | 2.5e44 | 1.4e43 | 9.3 / 2.0 | 0.14° | 0.009 / 0.017 |
| GALEX J074320.3+160752 (= PN G203.8+18.5 in SIMBAD; no PN object type) | DO | 83.2 kK (Bédard+20; 100 kK Kepler+15) | 1233 | 3.8e45 | 6.9e44 | 22.9 / 4.9 | 0.74° | 0.018 / 0.058 |
| HD 60848 | O8:V:pe | 37.0 kK (mesFe_h 2024A&A...681A.107R); R = 7 R☉ assumed | 1913 | 3.9e48 | 2.3e46 | 230 / 50 | 2.26° (75 pc) | 2.0 / 24 |
| NGC 2392 (CSPN) | O6f | 45 kK (Culpan) | 1835 | 2.9e47 | 5.4e45 | 97 / 21 | 4.04° | 0.07 / 2.0 |

**Notes on the table.**
- Only stars at ≳ 1.8 kpc can close the photon budget: HD 60848 and the NGC 2392 nucleus. For that the Silk would have to be at their distance, i.e. 0.3–0.5 kpc long and 0.6 kpc above the plane.
- HD 60848 is a 37 kK O star with Q(>35 eV)/Q(H) = 0.006. Its H II region would be Hα-dominated, not [O III]-dominated with weak Hα. Its SIMBAD RV (+5.5 km/s heliocentric ≈ −6 km/s LSR) also differs from the [O III] gas (V_LSR +11…+16). I therefore reject it as the ionizing source.
- **Gas kinematics (derived).** Flat rotation (R₀ = 8.18 kpc, Θ₀ = 236 km/s) predicts V_LSR along (204°, +19.4°) of:

  | d (kpc) | 0.16 | 0.3 | 0.5 | 1.0 | 1.5 | 2.0 |
  |---|---|---|---|---|---|---|
  | V_LSR (km/s) | +1.5 | +2.8 | +4.5 | +8.7 | +12.5 | +15.9 |

  So the [O III] velocity matches pure rotation only at 1.3–2 kpc. Peculiar and expansion motions of ~10 km/s make this unusable as a distance. The [O III]-minus-low-ion offset (+10 km/s) is better read as local gas motion.
- **PNe.** None lies within 5° of the Silk centre. The nearest, with Gaia distances from SIMBAD parallaxes:

  | PN | sep. from Silk centre | distance |
  |---|---|---|
  | A66 21 (Medusa, DOZ) | 5.3° | 592 pc |
  | NGC 2392 | 6.1° | 1.83 kpc |
  | HDW 7 | 7.2° | 1.43 kpc |
  | A66 24 | 13.5° | 719 pc |
  | K 1-11 | 14.0° | 854 pc |
  | K 2-2 (CSPN 65.9 kK, Culpan) | 14.7° | 831 pc |
- **Wake test** (LSR proper-motion direction and 0–1 Myr LSR back-track vs the spine). No star with Teff ≥ 40 kK has its past track along the spine: the maximum fraction of spine points within 0.5° of a track is 0.076 (GALEX J075633.8+163358, 1.35 kpc). The best cases are cool (20–31 kK) GF21 WDs with fractions 0.15–0.23 and position-angle offsets of 7–27°. Their Q(H) is far too small, so none is a credible wake source.

**Verdict.** No catalogued hot star within 1 kpc can photoionize the Silk to the observed [O III] surface brightness. The DOZ WD J0736+1444 at 150 pc is geometrically close and hot enough, but falls short by ≥ 1.5× even with 4π coverage and ~20× realistically. Shock or interface excitation, or photoionization by the Monogem hot gas's own EUV/soft X-rays, remains the more plausible source of the O⁺⁺. That last option is not quantified here.

## 5. Local Bubble and 3D dust toward (204°, +19°)

Script `scripts/context_dust_analysis.py`, table `data/derived/context_dust_los.csv`, figure `figures/context_dust_profiles.png`.

**Edenhofer+2024** (A&A 685 A82; Zenodo 10658339 `mean_and_std_healpix.fits`, 3.25 GB). Read by HTTP Range requests only (516 radii × 4 nside-256 NEST sub-pixels per sightline; header: NSIDE 256, ORDERING NEST). Units are E of Zhang, Green & Rix 2023 (≈ E(B−V)). Errors are the summed posterior std, ignoring radial correlation.

| sightline | E(<160) | E(<300) | E(<600) | E(<1.25 kpc) | density peak <400 pc |
|---|---|---|---|---|---|
| P1 | 0.0021 | 0.0054 | 0.0105 | 0.0266 ± 0.0011 | 167 pc |
| P2 | 0.0020 | 0.0051 | 0.0109 | 0.0272 ± 0.0010 | 153 pc |
| P3 | 0.0024 | 0.0061 | 0.0100 | 0.0315 ± 0.0015 | 161 pc |
| P4 | 0.0021 | 0.0055 | 0.0095 | 0.0292 ± 0.0013 | 160 pc |
| 4° W of P1 | 0.0021 | 0.0054 | 0.0087 | 0.096 ± 0.008 (cloud at ~700 pc) | — |
| 4° E of P1 | — | 0.0047 | — | 0.026 | — |

Along the spine (Dec 12–28), E(<300) = 0.0040–0.0079, and the nearest density peak is at 150–182 pc. The one exception is Dec 21, where the peak is at 349 pc. Reading A_V ≈ 2.8 E, the foreground to 300 pc is A_V ≈ 0.015 mag.

**O'Neill+2024** Local Bubble shell (ApJ 973 136; Harvard Dataverse doi:10.7910/DVN/INB1RB, nside 256). Toward P1–P4 the peak is at 161–171 pc (σ 12–18 pc), the inner edge at 133–136 pc and the outer edge at 196–201 pc, with shell A_G = 0.0054–0.0067 mag. Along the spine the peak is at 158–167 pc for Dec ≤ 22 and rises to 172–210 pc (σ up to 65 pc) for Dec 23–28.

**Pelgrims+2020** inner surface (A&A 636 A17; Harvard Dataverse doi:10.7910/DVN/RHPVNC, nside 128 RING): 221–244 pc raw and 225–253 pc at l_max = 10 along the spine.

**Lallement+2022.** The third-party HEALPix-7 table of the A₀ = 0.05 mag distance (van Wolde, Zenodo doi:10.5281/zenodo.22129167, derived from 2022A&A...661A.147L) is NULL for all tiles here. It is reported but not used. The STILISM API was unreachable (TLS error).

**Implication.** The Local Bubble wall here is at ~160–170 pc and is very thin (E ≈ 0.002–0.003 mag). Wu+2026's 160 pc puts the Silk *on* the wall. The Monogem interpretation puts it 80–150 pc beyond the wall, inside the low-dust cavity between the wall and ~600 pc. The 3D dust map cannot tell the two apart: there is no dust feature at 250–320 pc to anchor a distance.

## 6. Catalogued nebulae near the trace and the northern Hα patch (SIMBAD TAP)

Box RA 100–135°, Dec −2…+36°. Script `scripts/context_snr_nebulae.py`, table `data/derived/context_nebulae_near_trace.csv` (51 objects within 3° of the spine or the Hα patch, after removing extragalactic H II regions).

- **No** Sh2, LBN, Galactic H II region, reflection nebula, PN or SNR lies within 3° of the spine or the Hα patch. The nearest PN is NGC 2392, 4.0° from the spine. All 'HII' entries in the box (DDO 47, UGC 3995, UGC 4115, CGCG 148-006, Arp 82) are extragalactic.
- **At the Hα patch** (centroid 115.474, +20.536; 1.45 × 0.31°, PA 163°):
  - PGCC G199.09+19.95 (Planck cold clump, Planck 2016 XXVIII), 0.35° from the centroid;
  - PGCC G200.5–201.0+23.3–23.8 clumps, 3.5–3.8° away;
  - cloud WV 216 (Wakker & van Woerden 1991 HVC catalogue), 1.8° away;
  - [GB2001] AG 2 (Génova & Beckman local-ISM absorption cloud, 2001/2003), 3.2° from the patch and 0.68° from the spine.

  A cold clump at the place where the [O III] arc turns into Hα+[SII] fits Fesen's picture: the shock slows where it meets denser gas.
- **HI on the line of sight.** Many GALFA compact HI clouds ([SPG2012], Saul+2012) and HVCs lie on the line of sight. Their velocities are in the names (+52…+149 km/s; e.g. [SPG2012] 116.5+15.0+149, 0.07° from the spine at +162 km/s) and are kinematically unrelated to the Silk (V_LSR −1…+16 km/s).

## Hypotheses for the team

1. **Silk = Monogem's east/north boundary at 250–320 pc.** Supported by the X-ray edge (ROSAT here, eRASS1 by imaging_mw), by Fesen+2026, and by the absence of a local ionizing source.
   - Size: at 288 pc, 10–15° gives 50–75 pc and the 0.28° E-strand FWHM gives ~1.4 pc.
   - Test: the [O III] V_LSR (+11…+16 km/s) against the Monogem expansion geometry at the limb. Expansion is mostly transverse at the limb, so a small radial velocity is expected.
   - Test: the low-ion lines' velocity equal to P4's DIG (Wu+2026).
2. **Silk = the Local Bubble wall interface at 160 pc** (Wu+2026's distance). This would need the [O III] to be energized by something not seen: no SNR interior, no hot star, and the X-ray-bright side faces Monogem, not the observer's own bubble.
3. **The S-shape, with its kink at Dec ≈ 17.6°,** may mean two physically different sections.
   - The southern section sits at the Monogem radius from the Monogem and pulsar centres (12.10 ± 0.12° and 11.36 ± 0.16°). It is a better match to the Monogem rim than to G205.6+12.4 (6.5 ± 0.5°, well inside that 8.8° radius).
   - The northern section bends away toward the north-west. That could be where the shock wraps around denser gas: the cold clump at the Hα patch, the Hα cloud east of it, and more HI just outside (imaging_mw).
   - Worth testing with separate velocity and line-ratio fits of the southern and northern parts (DESI/SDSS teammates).
4. **PSR J0611+1436** (τ_c 1.07 Myr, DM 0.8 kpc, 2.7° from the curvature centre) has no proper motion. A VLBI or timing PM would settle whether it could trace back to the curvature centre. This is low priority given item 3 of the summary.

## Figures

| file | caption |
|---|---|
| `figures/context_map.png` | Galactic map (l 232–178°, b −13…+35°) on ROSAT 1/4 keV. Shows the Silk spine, W strand and Hα patch, P1–P4, the Monogem 1- and 2-SNR outlines and G205.6+12.4 (Knies+24), Green/Fesen SNRs, the Origem Loop, the small-circle fit and its centre, pulsars < 1 kpc with B0656+14 and Geminga LSR back-tracks (0–1 Myr, 100-kyr ticks, MC loci at τ_c), hot stars ≥ 40 kK, PNe, and the R20/W06 filaments. |
| `figures/context_curvature_pulsars.png` | RA/Dec. Small-circle fits to the E, E-north, E-south, W and E+W strands with bootstrap boxes and Δχ² = 2.3/6.2 contours of the E-strand centre. Pulsar back-tracks and MC birth loci, and literature centres (Monogem, G205.6+12.4, Origem, NSNS). |
| `figures/context_concentricity.png` | Left: angular distance of every spine point from the candidate centres vs Dec, with the Monogem radius ± 10 % band. Right: rms of those distances (log scale) for 11 centres, with the circle-fit rms marked. |
| `figures/context_rosat_profiles.png` | ROSAT R1+R2 and R4+R5 profiles parallel to the spine in 4 Dec bands. The X-ray-bright side is west, i.e. the Monogem interior. |
| `figures/context_hotstars.png` | RA/Dec zoom. Hot stars (Teff ≥ 40 kK, d < 1 kpc) coloured by Teff and sized by Q(H), with Strömgren circles at n = 0.1 cm⁻³ and LSR proper-motion arrows (50 kyr). Labelled stars include WD J0736+1444, GALEX J0745+1449, PN G203.8+18.5, HD 60848, Medusa and NGC 2392. |
| `figures/context_dust_profiles.png` | Left: Edenhofer+24 dust density vs distance toward P1–P4 and 4° E/W of P1. Right: Local Bubble wall distance along the spine (O'Neill+24 peak and inner–outer range, Pelgrims+20 inner surface, Edenhofer nearest density peak) with the 160 pc, 288 pc and 250–320 pc Monogem-rim markers. |

All figures are PNG at 200 dpi, made with `scripts/silk_style.py`, and were inspected after the final render.

## Files

**Scripts** (all in `scripts/`; run with `~/.venvs/astro313/bin/python -I`):

| script | purpose |
|---|---|
| `context_common.py` | helpers: small circle, LSR back-track, blackbody Q, Strömgren radius |
| `context_fetch_catalogs.py` | VizieR / SIMBAD fetches |
| `context_psrcat.py` | parse ATNF psrcat |
| `context_pulsar_backtrace.py` | pulsar back-traces |
| `context_snr_nebulae.py` | SNR and nebula tables |
| `context_hotstars.py` | hot-star photon budget |
| `context_xray_profiles.py` | ROSAT profiles |
| `context_dust_fetch.py` | Edenhofer range reads |
| `context_dust_analysis.py` | dust and Local Bubble tables, figure |
| `context_concentricity.py` | concentricity and rim geometry |
| `context_objects.py` | merged object table |
| `context_figures.py` | map, curvature and hot-star figures |
| `context_arxiv_list.py` | list arXiv API results |

**Derived data** (`data/derived/`):
- `context_objects.csv`: 140 rows; Silk references, SNRs/outlines/shells, pulsars, hot stars, PNe, nebulae, named regions.
- `context_pulsars_30deg.csv`, `context_pulsar_tracks.csv`, `context_pulsar_birth_mc.csv`, `context_pulsar_coincidence.json`
- `context_curvature.json`, `context_curvature_regionA.npz`
- `context_concentricity.json`, `context_concentricity.csv`
- `context_snrs.csv`, `context_nebulae_near_trace.csv`
- `context_hotstars.csv`, `context_hotstars_top.csv`
- `context_rosat_profiles.csv`, `context_rosat_spine.csv`
- `context_dust_los.csv`

**Raw data** (`data/raw/`):

| directory | contents |
|---|---|
| `context_catalogs/` | VizieR / SIMBAD extracts, plus `context_fetch_log.json` with exact URLs |
| `context_psr/` | psrcat v2.9.0, VizieR B/psr cone |
| `context_arxiv/` | papers (PDF + pdftotext), arXiv API queries |
| `context_dust/` | Edenhofer range-read profiles and headers, O'Neill, Pelgrims, Lallement-derived table |
| `context_skyview/` | RASS background bands 1, 2, 4, 5, 6, Galactic CAR, 0.2° px |
| `context_simbad_verify/` | mesFe_h check, identifiers, references for nearby shells and clouds |
| `context_vizier_hi/` | Suad+14 and Ehlerová+05 HI shell cones |
| `context_ads/` | ADS page shell; no abstract |

## Sources (all accessed 2026-11-16)

**Catalogues and data**
- ATNF Pulsar Catalogue v2.9.0: https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz (Manchester+2005)
- VizieR B/psr (2017 snapshot): https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=B/psr/psr&-c=116.68+16.42&-c.r=30&-c.u=deg&-out.all&-out.max=2000&-oc.form=d
- VizieR CfA mirror (vizier.cfa.harvard.edu; exact query URLs in `data/raw/context_catalogs/context_fetch_log.json`):
  - VII/297, Green 2025 (2025JApA...46...14G)
  - J/MNRAS/508/3877/maincat, Gentile Fusillo+2021
  - J/A+A/662/A40/hotsd and knownhsd, Culpan+2022
  - J/A+A/635/A193, Geier 2020
  - V/84 main and dist, Acker+1992
  - J/A+A/564/A116, Suad+2014
  - J/A+A/437/101, Ehlerová & Palouš 2005
- SIMBAD TAP (CfA mirror): https://simbad.cfa.harvard.edu/simbad/sim-tap/sync. Tables: basic, flux, ident, mesFe_h, has_ref, ref.
- SkyView RASS Background 1, 2, 4, 5, 6 (Snowden+1997): https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl (Galactic CAR, 0.2° px; 60 × 50° centred (205, +12.5) and 70 × 60° centred (205, +5))
- Edenhofer+2024 3D dust: https://zenodo.org/records/10658339 (`mean_and_std_healpix.fits`, Range requests only)
- O'Neill+2024 Local Bubble shell: doi:10.7910/DVN/INB1RB (`LB_ShellProperties_MeanDraws.fits`)
- Pelgrims+2020: doi:10.7910/DVN/RHPVNC (`map-inner_final.fits`)
- Lallement+2022-derived HEALPix table: doi:10.5281/zenodo.22129167 (van Wolde, v1.0)

**Literature read**
- Wu+2026, arXiv:2610.13329 (`sources/`)
- Fesen+2026, arXiv:2608.26666v1
- Knies+2024, arXiv:2401.17289
- Thorsett+2003, arXiv:astro-ph/0306462
- Reich, Reich & Sun 2020, arXiv:2007.08221
- Weinberger, Temporin & Stecklum 2006, arXiv:astro-ph/0511546
- Reindl+2023, arXiv:2307.03721 (Table A.4)
- Gao & Han 2013, arXiv:1212.1523 (abstract only)
- NSNS gem-01 annotations, via `silk_trace.json` literature_geometry (simg.de/nebulae3/gem-01.html)

Also downloaded but not used for numbers: Yao+2022 (arXiv:2209.14059) and Gómez-Muñoz+2023 (arXiv:2304.01970).

**Secondhand only:** Nousek+1981, Plucinsky+1996, Brisken+2003, Kim+2007 (ApJ 665 L139), Knies+2018, Faherty+2007 (Geminga PM/PX values), Pellizza+ (Geminga birthplace via NSNS).

## Limitations

- **Curvature.** A single small circle is a poor model of the arc (χ²/dof ≫ 1, S-shaped). "Centre of curvature" regions are therefore indicative only; region B's union of bootstrap boxes is the honest systematic range.
- **Back-traces.** They assume LSR-comoving birth gas, no Galactic potential (fine for ≤ 1 Myr) and an unknown v_r. DM distances have 40 % log-normal scatter. 9 of 17 pulsars within 1 kpc have no PM and cannot be traced.
- **Hot stars.**
  - Blackbody fluxes overestimate Q(>35 eV) for H-rich WDs and underestimate it for some He-rich ones; NLTE models were not used.
  - GF21 photometric Teff above ~30 kK is unreliable.
  - SIMBAD-type-only stars use generic Teff/R.
  - The photon budget uses the lower end of the [O III] intensity and area, and a generous η = 2. It is an upper bound on what a star could do.
  - HASH was not queried; PNe come from SIMBAD.
  - The Culpan+22 Gaia candidate list (no Teff) and Geier 2020 were fetched but not used (superseded by the Culpan known-hsd table).
- **ROSAT profiles.** Smoothed to 1°, so a sharp edge narrower than ~0.5° is unresolved. Errors are approximate. The 1/4-keV contrast is partly sensitive to absorption differences, although the foreground dust is symmetric.
- **Dust.** Edenhofer posterior mean and std only (no samples); 4 sub-pixels averaged (~0.46°). The map has no feature at 250–320 pc, so dust cannot fix the Silk's distance.
- **Literature.** Monogem values from HEAO-1, ROSAT and FIMS are secondhand. The Plucinsky+1996 explosion energy is quoted inconsistently: 1.9×10⁵⁰ (Thorsett+03) vs ~0.2×10⁵⁰ erg (Knies+24).
- **HI shells.** The catalogue "distances" for GS202+17+019 etc. are rotation-curve values and meaningless at l ≈ 202°.

## Data gaps

| source | status |
|---|---|
| Kim+2007 FIMS/SPEAR | not on arXiv; ADS abstract page JS-only (405 to WebFetch); KAIST DSpace refused TLS; only secondhand summaries used |
| Plucinsky+1996, Nousek+1981 originals | not read |
| STILISM (Lallement/Vergely G-Tomo) API, stilism.obspm.fr | TLS connection failed; the Lallement-derived Zenodo table is NULL in this region |
| HASH PN database | not attempted |
| eROSITA DR1 maps | analysed by the imaging_mw teammate, not re-analysed here |
| cdsarc.cds.unistra.fr (CDS FTP) | unreachable; the VizieR CfA mirror was used instead |
| atnf.csiro.au | listed as unreachable in the brief but reachable on 2026-11-16 (HTTP 200); used directly |
