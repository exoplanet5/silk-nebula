# w4_distance: distance leads, checkable now and for the future

Wave 4, tag `distance`. Written 2026-11-17.

I read `w3_canon`, `w3_oii`, `w3_northsouth`, the three `w3_judge_*` notes, `verify_distance`, `verify_xrayfuv` and `context` §5. Precedence: w3 > verify > wave 1. All numbers below are my own, from the scripts listed in §6, unless marked *(literature)* or *(canon)*.

## 0. Summary

| # | Question | Answer | Strength |
|---|---|---|---|
| 1a | Distance of PGCC G199.09+19.95 (the cold clump next to the Hα patch) | **Hint that it lies beyond ~400 pc** (about 400–650 pc). Neither 160 nor 282 pc is supported. The PGCC catalogue's own value is 262 ± 121 pc (SDSS DR7 extinction) | ~2–2.5σ; low–medium |
| 1b | Dust at 300–600 pc north and east of the ridge | It is real, but it sits **beyond** the H1 distance: a cloud 3–5° E at **~367 pc** (Dec 16–20), and a rise at 550–650 pc in the north. There is **no dust feature at 250–320 pc** anywhere along the ridge | medium |
| 2 | Kim+2007 FIMS/SPEAR: does the bright high-latitude C IV region "R1" cover the Silk? | **Yes.** R1 lies along the southern Silk (Dec < 17.6). The digitised C IV map peaks +0.28 ± 0.05° (± 0.1° registration) east of the [O III] spine, and its width equals the 2° map beam. Placebo p = 0.017–0.027. R1 table value: C IV = 8.7 ± 0.9 × 10³ LU; O III] and O VI not detected | medium–high (spatial); medium (photometry) |
| 2b | C IV/[O III] implied by R1, compared with GALEX | Beam-deconvolved C IV/[O III]5007 ≈ **0.15–0.4 photons** (central ~0.25). GALEX allows ≤ 0.41 (all FUV is C IV) and gives 0.16 ± 0.05 after removing the R88 non-resonant floor. Unsuppressed shocks give 0.6–2.2, so the implied suppression is **~2.5–9×** (canon ~4–6×). FIMS is ~4–10× above the warm-layer (H6) prediction of 0.03–0.04 | medium |
| 3 | Target list for a high-ion absorption bracket | **52 Gaia hot stars** (38 WDs, 14 A/late-B) within 0.3° of a strand at 120–800 pc. **18** reach S/N 15 at C IV in ≤ 5 COS orbits. A recommended 7–9-star bracket is in §3.3. None has archival UV spectra. A COS G130M QSO sightline exists 0.34° from the W strand | — |
| 4 | Expected proper motion and a two-epoch plan | 160 / 282 / 1000 pc → **132–198 / 75–112 / 21–32 mas/yr** (v_s 100–150). Epoch 1 = SWONS ε-180ED + QHY600, 2024-01, ~1.55″/px. A 10–12 yr repeat detects 282 pc at ≥ 5σ if the per-epoch edge precision reaches ≤ 0.1–0.15″ | plan |

**Effect on the hypotheses**
- FIMS strengthens "shock (any driver)" against H6 and conduction fronts. It does not separate H1 from H2.
- The PGCC hint and the 3D dust weaken the H1-specific "second coincidence" (Hα patch and PGCC on the Monogem radius). If the patch is a shock hitting this clump, the clump sits at ≳ 400 pc, not 252–318 pc.
- Neither result is decisive. The distance is still unmeasured. The COS bracket (§3) and the proper motion (§4) are the routes to a measurement.

## 1. PGCC G199.09+19.95 and the 300–600 pc dust

### 1.1 What the clump is
- **Catalogue row** (Planck 2016 XXVIII, VizieR J/A+A/594/A28 via the CfA mirror): (l, b) = (199.094, +19.947), (RA, Dec) = (115.2847, +20.8334), SNR 10.1.
  - q_Flux = 3, so no size, temperature, N(H₂) or mass is given.
  - **DIST = 0.262 ± 0.121 kpc**, from DIST_OPTION 2 (optical extinction on SDSS DR7), with q_Dist = 2. I did not re-read the paper's flag definitions.
  - It is the only PGCC within 2°.
- **2-D column** (wave-1 SkyView cut-outs on disk):
  - SFD E(B−V): peak **0.097** at (115.308, 20.840), 0.02° from the PGCC position. Mean 0.082 within r < 0.1°, falling to 0.063 at 0.1–0.2° and 0.050 beyond 0.4°. Background 0.048 ± 0.008.
  - **Clump excess ΔE(B−V) ≈ 0.034** (A_V ≈ 0.1). Compact: r ≲ 0.15°.
  - IRIS 100 µm +1.7 MJy/sr; Planck 857 GHz +2.0 MJy/sr.
- **HI** *(canon)*: CNM at V_LSR −0.7 and +6.5 km/s. This gives no distance.

### 1.2 The Edenhofer+24 3D map does not resolve the clump
- **Fetch.** `w4_distance_eden_fetch.py` range-read the nside-256 sub-cube over l 189–211, b 13.5–27.5 at 90–750 pc. That is 37 nside-16 parents and 16.8 MB, from Zenodo 10658339 (`mean_and_std_healpix.fits`, layout as in `context_dust_fetch.py`).
- **On − off** (r < 0.2° against 0.6–1.5°) integrated to 750 pc: **0.0022 ± 0.0026**, against the SFD excess of 0.034. The map recovers ≲ 10 % of the clump.
- The little excess there is accumulates at 500–660 pc (16/50/84 % at 503/644/661 pc). It is not significant.
- Both on and off profiles show the general dust layer peaking at ~650–670 pc in this direction (z ≈ 220 pc).

### 1.3 Stellar extinction against distance ("SFD tomography")
**Method** (`w4_distance_clump_stars.py`).
- For stars near the clump I regress each star's extinction on the SFD E(B−V) at its position, in distance bins.
- Within r ≲ 0.5° the clump dominates the SFD contrast. The slope dA/dE(SFD) should be ≈ 0 for stars in front of the clump and ≈ 2.8–3.2 (A_V) for stars behind it.
- Stars:
  - Gaia DR3 with parallax/error > 5 and RUWE < 1.4 (Data Lab, ±1° box);
  - StarHorse 2021 (Anders+22, VizieR I/354), 0.4° cone, with fidelity > 0.5, FlagOut = 0 and a fractional distance error < 0.25.

| Distance bin (pc) | StarHorse A_V slope (n; n with SFD > 0.065) | Gaia GSP-Phot A₀ slope (n) |
|---|---|---|
| 150–250 | — (7 stars) | −4.5 ± 6.5 (23) |
| 250–350 | — (7 stars) | −3.3 ± 4.2 (39) |
| **350–500** | **−0.9 ± 1.6 (38; 11)** | 3.3 ± 3.6 (73) |
| **500–800** | **3.7 ± 0.9 (128; 28)** | 3.3 ± 2.1 (207) |
| 800–1500 | 4.1 ± 0.6 (427; 97) | 2.9 ± 0.7 (535) |

- **Reading.** Stars at 500–1500 pc carry the full clump extinction (4–7σ). StarHorse stars at 350–500 pc show none: −0.9 ± 1.6 against ~3, about 2.5σ.
- Gaia GSP-Phot is noisier. Its two bins below 350 pc combine to −3.6 ± 3.5, about 1.9σ below "in front".
- A simple on/off median comparison agrees in sign but has low S/N:
  - StarHorse, r < 0.12° against 0.25–0.4°: ΔA_V = +0.05 ± 0.09 (450–600 pc), +0.01 ± 0.09 (600–800), **+0.10 ± 0.04 (800–1100)** and **+0.13 ± 0.05 (1100–1600)**. The on-clump sample is small (4–24 stars per bin), so this shows the clump is in front of the ≳ 800 pc stars but cannot localise it further.
  - Gaia's step fit is unstable with radius (best d_c 200 pc at r < 0.15°, 300 pc at r < 0.25°, ±1-σ windows 200–780 pc). I do not use it.
- **Clump distance: probably ~400–650 pc** (low–medium confidence, ~2–2.5σ). The PGCC catalogue's 262 ± 121 pc is 1.2–1.6σ away.

**Systematics**
- StarHorse posteriors depend on its priors. GSP-Phot A₀ is biased at low extinction.
- Only 4–11 high-SFD stars per bin fall below 500 pc.
- SFD has 6′ resolution, and its 0.1-mag column is spread over a ~0.2° clump.
- StarHorse distances have up to 25 % errors, so the 350–500 bin leaks into its neighbours.
- **Test:** a deeper Gaia XP or Zhang+23 extinction set. The CfA VizieR mirror returned **no rows** for J/MNRAS/524/1855 on 2026-11-17. StarHorse control cones of 1.0° and four offset 0.4° cones also came back empty or with errors, logged.

**Consequence**
- `w3_judge_skeptic` counts the Hα patch and the PGCC on the Monogem-radius circle as a second coincidence for H1.
- If the clump is at ≳ 400 pc, it is not at the H1 distance (252–318 pc). Then either the patch–clump adjacency is a projection, or the patch (and the Silk) is farther away.
- Either way this coincidence should no longer count for H1. It is not strong evidence against H1.

### 1.4 The 300–600 pc dust north and east of the ridge
From `w4_distance_dust.json` → `ridge_sides_300_600`; figure `w4_distance_dust_clump.png`.
- **At 250–320 pc** (the H1 shell distance): E = 0.0008–0.0015 everywhere: on the ridge, 1–5° E and 1–3° W. **There is no feature**, so there is no dense pre-shock layer at the Monogem-limb distance in the 3D map.
- **At 320–400 pc:**
  - A discrete cloud **3–5° east** at Dec 16–20. E(320–400) = 0.006–0.0065 against 0.0016–0.0020 on the ridge and 0.0011–0.0013 to the west. Density peak **367–370 pc**.
  - At 368 pc it is ~89 pc from a Silk at 282 pc, and ~120 pc from the Monogem centre (288 pc). That is twice the 61-pc shell radius, so it is **not in contact with the Monogem shell**.
- **At 400–600 pc:**
  - Dec 20–26: E(400–600) = 0.011–0.018 at 1–5° E, against 0.007–0.008 on the ridge. Density peaks at 557–648 pc.
  - Dec 26–30: the excess shifts to the west side.
  - This NE complex at 400–650 pc is probably what the PGCC clump belongs to (§1.3).
- **Reading.** The dust rising north and east lies 90–370 pc **beyond** a 282-pc Silk. It is not obviously the "denser gas beyond the limb" invoked for H1.
  - It could be the medium of a Silk at ~350–600 pc, which is part of H2's 0.1–0.8 kpc range.
  - 160 pc has nothing either (the O'Neill wall is thin and featureless here; *canon*).

## 2. Kim et al. 2007 (ApJ 665, L139): FIMS/SPEAR C IV and the Silk

### 2.1 Access record
Log: `data/raw/w4_distance_kim07/fetch_log.tsv`.
- **IOP PDF.** The first probe of `https://iopscience.iop.org/article/10.1086/521441/pdf` returned 200/PDF but was discarded (written to /dev/null). The second request hit a Radware bot challenge. I did **not** retry or bypass it.
- **IOP HTML full text** (`.../fulltext/21923.text.html` and `21923.tb1.html`) was read with WebFetch. This gave the abstract, methods, figure captions and **Table 1 verbatim**.
- **Figure 1** (C IV map, 600×611 jpg) and **Fig. 2** were downloaded from the IOP content CDN (`content.cld.iop.org/.../revision1/fg1.jpg`, `fg2.jpg`). The `fg1.html` page itself was bot-blocked.
- **Not available elsewhere:** the paper is not on arXiv (arXiv:0710.1531 is the Antlia paper), and the ADS full-text gateway redirects to the same bot wall.

### 2.2 What the paper reports
All from the HTML text and Table 1.
- **Data.**
  - FIMS L-band 1340–1750 Å (analysed 1360–1680), λ/Δλ ~550, 5′ resolution.
  - Observations 2003-12-13 to 2004-01-19.
  - Analysis area 6h12m–8h, Dec −3 to +27. This includes the Silk to Dec 27.
- **C IV map.**
  - HEALPix Nside 64 (55′ pixels), 120′ FWHM Gaussian smoothing.
  - Doublet fitted with a fixed 2:1 ratio. No extinction correction.
  - Stars with continuum > 20,000 CU masked.
- **Table 1** (10³ LU):

  | Line | Whole region | R1 |
  |---|---|---|
  | C IV | 7.1 ± 0.2 | **8.7 ± 0.9** |
  | C III 977 | 4.9 ± 0.4 | 6.6 ± 2.1 |
  | O VI | 4.0 ± 0.5 | not detected |
  | Si II* | 1.9 ± 0.2 | not detected |
  | He II 1640 | 2.3 ± 0.3 | not detected |
  | O III] 1666 | 2.5 ± 0.3 | not detected (> 2σ threshold) |

  - Other regions: R2 C IV 9.3, O III] 4.9 ± 1.0; R4 C IV 14.3 (peak 2.1 × 10⁴ at (196.2, +4.8)).
  - Si IV is not reported.
- **Interpretation in the paper.** R1 is "the C IV enhancement feature in the high-latitude region … likely due to the interaction of the blast wave with an isolated cloud". The paper says X-rays are fainter in R1's middle and bright around its edges.
  - Distance assumed: 300 pc.
  - **No region coordinates are given in the text.** R1's position below is from the figure.

### 2.3 R1 is the southern Silk
`w4_distance_kim07_digitize.py`, `_placebo.py`, figure `w4_distance_kim07_civ.png`.
- **Registration.** The dotted graticule was traced (42 intersections) and fitted with a cubic (x, y) ↔ (l, b) map. Residual rms 0.10° in l and 0.02° in b. The colour bar was inverted to LU, linear from 0 to 21091.
- **R1 outline:** l ≈ 199.5–210, b ≈ 16.7–21.6.
- **Overlay.** The **south spine (Dec 12–17.6; l 208.1 → 202.5) runs along R1's major axis.** The W strand, the Hα patch and the PGCC sit at R1's western tip. The north arc (l < 201, b ≈ 21–22) runs along the top edge of the FIMS coverage.
- **Cross-spine profile of the digitised map:**

  | Section | Amplitude (LU) | Centre offset | Fitted FWHM / background |
  |---|---|---|---|
  | South (23 spine points) | **5.7 ± 0.3 × 10³** (5.4–5.9 for beam FWHM 1.9–2.4°) | **+0.28 ± 0.05°** (east; ± 0.1° registration) | free FWHM **2.34°** ≈ beam; background 4.9 × 10³ |
  | North | 1.9 ± 0.2 × 10³ | −0.9° (toward W2 and the patch) | ratio N/S 0.33 |

  - South: the source is unresolved, with intrinsic width ≲ 1°.
  - North: the N/S ratio of 0.33 matches the N/S [O III] ratio (~0.3).
- **Placebo.** I placed the rigid south track at 3794 random positions and orientations on the map, using the contrast against ±3° side tracks.
  - Observed contrast 5.3 × 10³ LU.
  - **p = 0.027** (all placements), **p = 0.017** (720 placements at b ≥ 12).
  - This is a posteriori (Kim defined R1 from the same map, and I chose the south section), so read it as ~2σ.
- **Consistency with Table 1.** If the Silk filament (0.28°, ~7° long inside R1) fills ≈ 9 % of the R1 ellipse, then C IV/[O III] = 0.25 predicts an R1 excess of ≈ 3.6 × 10³ LU. Table 1 gives 8.7 ± 0.9 against a map background of ~5 × 10³, an excess of ≈ 3.7 × 10³ LU.

### 2.4 C IV/[O III] compared with GALEX
- **FIMS line integral.** Across the filament, S_CIV = A·√(2π)·σ_beam = **1.28 × 10⁴ LU·deg** (1.19–1.39 × 10⁴).
- **[O III]5007 line integral.**
  - South I = 1.4–2.5 R (DESI 1.62 R; 1 R = 79,577 LU) × 1.064 × FWHM 0.28° (E strand) = 3.3–5.9 × 10⁴ LU·deg.
  - The W strand and the envelope (half-width 0.22°) can add up to ~×1.5.
- **Photon ratio: C IV/[O III] ≈ 0.15–0.4**, central ~0.25. This assumes all of the R1 C IV excess belongs to the Silk.
- **Comparison with GALEX** (`verify_xrayfuv`): FUV 92 ± 10 CU per NSNS-R (south); 1 R of C IV = 223 CU.
  - If all FUV were C IV: ≤ 0.41.
  - After the R88 non-resonant floor of 57 CU/R: 0.16 ± 0.05.
  - FIMS sits at or above the GALEX template value. That is expected if part of the C IV is broader than the [O III] ridge, which a C IV halo outside the front would be. The FIMS beam keeps it; the GALEX template regression does not.
  - The 1σ ranges overlap.
- **Implications**
  - **Suppression factor.** Unsuppressed shocks give C IV/[O III] ≈ (196–541 − 57)/223 = 0.6–2.2, so FIMS implies **×2.5–9** suppression. This is consistent with the canon C9 (~4–6×; `verify_shockmodel` 3.7–6×). It sits at the low end of `verify_xrayfuv`'s 6.5–26×.
  - **H6 warm layer** (C IV/[O III] 0.03–0.04, `w3_judge_skeptic` D3): FIMS is ×4–10 higher. **Disfavoured** if the attribution holds.
  - **Conduction fronts** (1.6–16): FIMS is ×5–60 lower. **Disfavoured.**
  - **FUV/C IV ridge position.** The +0.28 ± 0.05 (± 0.1 systematic) ° eastward offset has the predicted sign (FUV/C IV outside the [O III]; GALEX gives +0.04…+0.09°). But the 2° beam and the 0.1° registration error make it at most a ~2σ consistency, not a test.
  - **Distance:** none. Surface brightness is distance-independent.
- **Caveats**
  - The values are digitised from a colour JPEG (JPEG artefacts; nearest-colour inversion, ±10–15 %).
  - The doublet was fitted assuming a 2:1 ratio. An optically thick 1:1 doublet changes the total by < 10 %.
  - No extinction correction: the foreground E(B−V) to 300 pc ≈ 0.005, so < 5 %.
  - R1 may include unrelated C IV from Monogem's hot-gas boundary. The spatial match (unresolved width, centred within 0.3° along 5.5°) argues that most of it is the Silk.
- **Recommended follow-ups**
  - (a) Ask the KAIST/KASI FIMS team for the R1 spectrum binned along the Silk, to obtain Si IV, O III] and the C IV doublet ratio.
  - (b) The FIMS all-sky line maps (Jo et al. 2019-era releases), if public. Not checked here.

## 3. High-ion absorption bracket: target list

`w4_distance_targets.py` → `data/derived/w4_distance_targets.csv` (87 rows), `_targets_mast.csv`, `_targets_summary.json`; figure `w4_distance_targets_pm.png` (left).

### 3.1 Selection
- Gaia DR3 (Data Lab): BP−RP < 0.10, parallax/error > 5, 1/ϖ = 100–1000 pc, ≤ 0.5° from the E spine, W strand or W2 strand.
  - Data Lab returns NaN parallaxes for `parallax > x`, because PostgreSQL treats NaN as large. These were removed.
- Class from M_G.
- GALEX GR6/7 from **MAST Catalogs** (6″). The CfA VizieR mirror returned no rows for II/312 and II/335.
- Cross-identifications from Gentile Fusillo+21, Culpan+22, Geier+20 and SIMBAD files already on disk.
- MAST UV spectra from the `spectra_ism` CAOM inventory on disk. A fresh astroquery portal query failed through the proxy.
- **COS time:** S/N = 15 per resolution element at C IV (log N(C IV) ≈ 12.7 at 3σ; `w3_judge_phys` T1) with G160M.
  - f_λ is taken from GALEX FUV, assuming a flat f_ν.
  - **A_eff = 1000 cm² (±30 %)** and a 0.0735 Å resel. Both are my approximations; **verify with the COS ETC**.
  - The first orbit gives 2000 s, later orbits 2700 s.
- **Bright-object flag:** a global count rate > 15,000 c/s (FUV ≲ 12). This is a rough screen, not the official COS bright-object protection (BOP) check.

### 3.2 Totals
- **52 primary** targets (≤ 0.3°, 120–800 pc): 38 WDs and 14 "early-type".
  - The early-type stars are mostly A0–A3 (BP−RP −0.03…+0.09, M_G 0.1–1.8). Their 1550 Å continua are line-blanketed and weak, so they are poor C IV backgrounds. Their COS times are underestimated by the flat-f_ν assumption.
- **39** of the 52 have GALEX FUV. **18** reach S/N 15 in ≤ 5 orbits.
- **No primary target has an archival HST, IUE or FUSE UV spectrum.**

### 3.3 Recommended bracket (hot WDs; all on the concave/west side or within the [O III] FWHM)
Offsets are signed from the E spine (+ = east). DA Teff values are from Gentile Fusillo+21 where available.

| Role | Star (Gaia DR3) | d (pc) | Offset | Strand | FUV | COS orbits | Notes |
|---|---|---|---|---|---|---|---|
| Foreground, south | GALEX J074642.6+180022 (668214413526637568) | 154 | −0.11° | E | 16.9 | 3 | DA |
| Foreground / H3 test, south | WDJ074543.00+150756.14 (3164849460692252032) | 234 | −0.18° (0.04° from W) | W | 15.0 | **1** | DA 29.7 kK; best S/N |
| Inside the H1 window, south | GALEX J074507.1+143229 (3164680303403541248) | 294 | −0.26° | E | 17.7 | 5 | DA |
| Just beyond H1 | GALEX J074529.3+184002 (671671827839519744) | 389 | −0.23° (W) | E/W | 17.8 | 5 | DA, junction |
| Background, south | WDJ074604.28+173616.89 (668165420835570304) | 445 | −0.13° | E | 17.2 | 3 | DA 22.4 kK; inside the [O III] FWHM |
| Background, south | WDJ074507.53+144939.84 (3164741777272354816) | 594 | −0.28° | W 0.26° | 15.7 | **1** | **DO, 84 kK**: photospheric C IV/N V. Separate by velocity |
| Background, north | WDJ073925.47+262657.73 (869130097415182848) | 407 | +0.04° | E | 16.8 | 2 | DA 24.7 kK |
| Background, north | WDJ074242.51+242208.78 (867858550874512256) | 588 | −0.07° | E | 17.2 | 3 | DA 28.9 kK |
| Optional, north | HS 0742+2306 (674965174402448000) | 585 | +0.14° | E | (none) | — | DA 57 kK; no GALEX FUV match (coverage gap?); predicted FUV ~15–16 |

- **Total ≈ 23 orbits** for the 8 stars with FUV. A minimal south-only bracket at 234 / 294 / 445 / 594 pc costs ≈ **10 orbits**.
- **Additional 6–8 orbit options:** GALEX J074712.0+150732 (388 pc, +0.15° E side), WDJ074635.97+185654.87 (388 pc), SDSS J074534.68+160709.3 (430 pc, inside the W strand).
- **A-type options** (FUV 12.3–14.4; 1 orbit nominal, contamination risk): Gaia 668163320594868480 (G = 11.1, 711 pc, 0.05° from the spine), 867221143366894208 (396 pc, 0.08°), 867221590043492608 (309 pc, 0.09°).
  - These need STIS E140H if they are too bright for COS. Use them only if no WD is available.
- **Predictions** (`w3_judge_phys` T1, `w3_judge_skeptic` D6):
  - H1: extra C IV/Si IV (log N ~13.5–14.5 edge-on, or 1–3 × 10¹³ per face-on crossing) at V_LSR 0 ± 20 km/s only beyond ~250–320 pc.
  - 160 pc: already present at 234 pc.
  - ≥ 1 kpc: absent in every star listed.
- **Geometry.** Most stars lie inside the limb (west), where a sightline crosses a curved shell twice, at ±√(2Rδ) ≈ ±11 pc about the tangent distance for δ = 0.2° and R = 61 pc. The bracket resolution is therefore ~±10–20 pc.

### 3.4 Existing archival UV spectra near the ridge (≤ 0.5°)
- **2MASS J07452815+1919525 = [VV2006] J074528.2+191952** (an AGN per its VV2006 designation): HST/COS FUV.
  - 1150–1450 Å, 23 exposures, 86 ks. This covers **Si IV 1393/1402 and N V 1238/1242**.
  - Also 1230–2050 Å, likely G140L (low resolution), and 31–2440 Å products.
  - 0.34° from the W strand, 0.62° W of the E spine.
  - It is a d = ∞ reference through the Monogem interior, available now; **not analysed here**.
- **SDSS J074511.56+194926.5:** COS NUV and broad-range products (no z = 0 C IV at medium resolution), 0.43° from W.
- **WD J0747+2106** (99.5 pc): COS G130M, already used (*canon*).
- **Not usable for interstellar high ions** (no hot-star or high-ion coverage):
  - HD 62721: STIS NUV/CCD only.
  - HD 61603: STIS CCD.
  - HD 62345: IUE SWP/LWP; not in the Gaia hot-star selection.

## 4. Proper motion and a two-epoch imaging plan

`w4_distance_pm.py` → `w4_distance_pm.json`, `_pm_table.csv`; figure `w4_distance_targets_pm.png` (right).

### 4.1 Expected motion
μ = v_s / (4.74 d). The [O III] zone moves with the front; on a limb the motion is fully transverse (× cos i otherwise).

| d (pc) | μ for v_s 100 / 120 / 150 km/s (mas/yr) | Shift in 10 yr (″) | Shift in 12 yr (″) |
|---|---|---|---|
| 160 | 132 / 158 / 198 | 1.32–1.98 | 1.58–2.37 |
| 252 / 282 / 318 | 84–126 / **75–112** / 66–100 | 0.75–1.12 (282) | 0.90–1.35 (282) |
| 1000 | 21 / 25 / 32 | 0.21–0.32 | 0.25–0.38 |

- **Reflex term.** Pre-shock gas at rest in the LSR (Schönrich+10 solar motion) gives, at 282 pc, μ ≈ (−6.1, −7.4…−9.7) mas/yr in (α*, δ).
  - The component normal to the filament is −2 to −8.4 mas/yr (westward), i.e. −2.8 to −11 km/s.
  - Its parallel component is 4.6–9.8 mas/yr.
  - It scales as 1/d, like the signal, so it is a fixed ~2–10 % (velocity-equivalent) bias, not a distance-dependent one.
  - The same model gives V_helio = +9.9 to +12.7 km/s for LSR-static gas, consistent with the measured +13.3.
- **Direction.**
  - South: the outward (east) normal has PA 68–115° along the S-shaped spine, close to radial from Monogem (PA 104 → 76° for Dec 12 → 18).
  - North (Dec 24–27): the normal is at PA 80–83°, while radial from Monogem is 38–50° and radial from the curvature centre 58–67°. So the motion's PA there separates the drivers, as `w3_judge_phys` T2 says.
  - Only the component normal to the filament is measurable on a smooth filament.

### 4.2 Precision (statistical only; my Fisher estimates)
**Epoch-1 setup.** SWONS 2024-01: Takahashi ε-180ED (180 mm, 500 mm f/2.8; the focal length is the nominal product spec, not from the paper) + QHY600 (3.76 µm) → **1.551″/px**, 4.1° × 2.8° per field. 36 h of [O III] through a 6 nm filter.
- My model assumes sky V = 21.3 mag/arcsec² and throughput 0.55.
- It gives **S/N ≈ 1.1 per pixel** at 2 R (south) and 0.34 at 0.6 R (north).

| Feature | Per-epoch position σ |
|---|---|
| Sharp edge, width 10″, length 0.5–2° | 0.05–0.10″ |
| Edge, width 60″ | 0.13–0.25″ |
| Edge, width 300″ | 0.28–0.57″ |
| Smooth south Gaussian profile (FWHM 0.28°), 1–5° long | 0.23–0.51″ |
| North edge (60″, 2°) | 0.42″ |

- A 0.6 m f/4 telescope at 0.6″/px with a 3 nm filter for 30 h gives ~3× better precision per unit length.
- **Systematics** will likely set a floor of ~0.1–0.2″: flat-field gradients, star residuals, different filters and continuum subtraction, PSF and seeing. Registration is not the limit: Gaia DR3/DR4 gives thousands of stars per field, ≪ 0.05″.
- **5σ in 12 yr needs μ ≥ 59 / 118 / 295 mas/yr for σ_epoch = 0.1 / 0.2 / 0.5″.**
  - 0.1″ detects d ≲ 340–500 pc (v_s 100–150).
  - 0.2″ detects d ≲ 170–250 pc.
  - 1 kpc needs ≥ 20–30 yr or σ ≤ 0.05″.
- The edge width below the 20″ NSNS sampling is unknown. This is the main uncertainty.

### 4.3 Concrete plan
1. **Epoch 1 (2024.03).** Request the SWONS ε-180ED/QHY600 raw [O III] subs with timestamps, calibrations and filter ID from the SWONS team (Wu, Zhang et al.). The paper uses 36 h of [O III] + 24 h of Hα. Measurements must be made on FITS, not JPEG.
   - Ask NSNS (Ziegenbalg, 2022-12 → 2025-04) for native-scale frames. The public HiPS products are 20″/40″ per pixel, useful only as a coarse check of the 160-pc case over ≥ 20 yr.
   - The NHZ 216 h mosaic and the Fesen+26 data are other possible epoch-1 sets (dates unknown).
2. **Epoch 1b (2026–27, recommended now).** A 0.5–1 m telescope at ≤ 0.6–1″/px with a 3 nm [O III] filter, 20–30 h per 0.5–1° field, on three fields:
   - (A) the south core, Dec 13–16 (P1–P3, NSNS 1′ peaks of 3.0–3.6 R);
   - (B) the north sharp east edge, Dec ~21–26 (the one-sided limb of `w3_judge_obs` §4);
   - (C) the Dec 17–19 junction where the W strand separates (sub-filaments).
   - Also image Hα/[S II] for the patch.
3. **Epoch 2 (2034–2036).** Repeat with the **identical** SWONS setup (≥ 10–12 yr baseline from 2024). Repeat the 1b fields (≥ 8 yr). A third epoch around 2040 covers the 1-kpc case.
4. **Measurement.**
   - Gaia-anchored SIP astrometry with Gaia stars propagated to each epoch.
   - PSF-matched, continuum-subtracted, star-masked images.
   - Cross-correlate profiles in strips every ~5′ perpendicular to the filament.
   - Fit a common normal shift plus segment terms.
   - **Placebo/nulls:** background galaxies and blank-sky strips (zero motion); the component parallel to the filament (should equal only the reflex parallel term, 4.6–9.8 mas/yr at 282 pc).
5. **Decision.**
   - A normal motion > 0 at ≥ 5σ establishes a moving front, i.e. a shock. Combined with v_s from resolved line profiles (predicted FWHM 23–38 km/s; needs R ≳ 20,000 spectroscopy), d = v_s / 4.74μ gives an expansion parallax good to ~±25 %.
   - A motion PA in the north of ≈ 40–50° favours Monogem; ≈ 80° (along the normal) or ≈ 60° favours a western driver at the curvature centre.
   - μ ≈ 0 (≤ 0.2″/decade) excludes H1 and 160 pc (`w3_judge_phys` T2).

## 5. Open issues and limitations
- **PGCC distance (2–2.5σ hint) needs better stellar extinctions.** Use Zhang+23 XP extinctions (the mirror gave no rows; try NOIRLab Data Lab `gaia_dr3_aux` or Zenodo), Green+19 Bayestar (STILISM/Bayestar endpoints were not tried), or LAMOST/DESI spectroscopic stars behind the clump.
- **FIMS:** the ratio is digitised, not a measurement from spectra. The R1 spectrum (Kim+07 Fig. 3) has no Si IV or C IV doublet ratio. The FIMS team or archive could supply the data along the Silk.
- **COS exposure times** come from a simplified A_eff and a flat f_ν. Run the COS ETC with a WD model per star. Check the BOP for the A stars.
- **The archival COS spectrum of the QSO 2MASS J0745+1919** (Si IV/N V at d = ∞) has not been analysed. It is a quick next step.
- The **SWONS optics and pixel scale** are nominal product values. Confirm them from FITS headers.

## 6. Files

**Scripts** (all run with `~/.venvs/astro313/bin/python -I`)

| Script | Purpose |
|---|---|
| `w4_distance_fetch.py` | PGCC, Gaia boxes, StarHorse, hot-star box |
| `w4_distance_eden_fetch.py` | Edenhofer sub-cube |
| `w4_distance_dust.py` | 3D dust, clump on/off, SFD/IRIS/Planck column, ridge sides |
| `w4_distance_clump_stars.py` | Gaia/StarHorse on/off and SFD tomography |
| `w4_distance_kim07_digitize.py` | Kim+07 map registration and profiles |
| `w4_distance_kim07_placebo.py` | Placebo test of the R1 coincidence |
| `w4_distance_kim07_fig.py` | Kim+07 overlay figure |
| `w4_distance_targets.py` | COS target list |
| `w4_distance_pm.py` | Proper-motion expectations and precision |
| `w4_distance_figs.py` | Dust and target/proper-motion figures |

**Derived data** (`data/derived/`)
- `w4_distance_dust.json`, `w4_distance_dust_pixels.npz`
- `w4_distance_clump_stars.json`
- `w4_distance_kim07_digitize.json`, `w4_distance_kim07_profiles.csv`, `w4_distance_kim07_map.npz`, `w4_distance_kim07_placebo.json`
- `w4_distance_targets.csv`, `w4_distance_targets_mast.csv`, `w4_distance_targets_summary.json`, `w4_distance_targets_run.log`
- `w4_distance_pm.json`, `w4_distance_pm_table.csv`

**Raw data** (each folder has a `fetch_log.tsv`)
- `data/raw/w4_distance_kim07/`: Fig. 1 and Fig. 2 jpg, thumbnails, and the bot-challenge page kept as the record of the failed PDF fetch.
- `data/raw/w4_distance_cat/`:
  - PGCC TSV and metadata;
  - five Gaia ±1° boxes;
  - the Gaia hot-star box;
  - StarHorse PGCC cone (2085 stars);
  - empty StarHorse and Zhang responses;
  - 87 MAST GALEX cones (`galex_mast_*.csv`) and the empty VizieR GALEX responses.
- `data/raw/w4_distance_eden/eden_subcube.npz`: 16.8 MB of range reads.

**Figures** (200 dpi, inspected)
- `figures/w4_distance_kim07_civ.png`: the Kim+07 C IV map with the Silk overlaid, and digitised cross-profiles.
- `figures/w4_distance_dust_clump.png`: Edenhofer E(320–400) and E(400–600) maps; clump on/off density profile; SFD-tomography slopes.
- `figures/w4_distance_targets_pm.png`: candidate stars (Dec against distance, FUV-sized, COS feasibility filled); proper motion against distance with detection thresholds and the reflex term.
