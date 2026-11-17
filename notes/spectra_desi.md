# Silk Nebula in DESI DR1 spectra: sky-subtraction residuals (topic tag `spectra_desi`)

Date: 2026-11-16. All velocities are heliocentric unless marked V_LSR. V_LSR = V_helio − 12.1 km/s toward the
nebula (standard solar motion, the same convention as Wu+2026). DESI barycentric velocities equal heliocentric
velocities to within 0.015 km/s. Line ratios are in energy units unless marked "photon".

## Bottom line

1. **DESI detects and traces the Silk Nebula on its own.** I analysed 5,220 DESI DR1 galaxy/QSO coadds from SPARCL;
   4,974 have clean [O III] windows. Their sky-subtraction residuals show [O III]5007 structure along the whole
   trace from Dec +14° to +28°.
   * The structure has the sign pattern expected from per-petal sky subtraction: positive on the strands, negative
     on fibres whose petals contain the strands.
   * A petal forward model of ridges on the geometry trace improves χ² by **Δχ² = 8,903** for 25 parameters
     (6,415 after clipping |χ| > 8).
   * The same model with the trace displaced by ±2–4° in RA gives only Δχ² = 18–69, and ±1° gives 319–707 because
     of partial overlap. This is a > 60σ detection that is specific to the mapped filament.
   * The Koposov/DESI hearsay (discovery note) is therefore technically plausible. What Koposov actually did
     remains unverified.
2. **The [O III] filament is 4–5× fainter than Wu+2026 state.** DESI is flux-calibrated.
   * Ridge peaks (forward model, in-fibre R): south (Dec < 19.5) median **2.4 R**. Examples: E strand 2.52 ± 0.08 R
     at Dec 16 and 2.53 ± 0.12 R at Dec 17; W strand 2.42 ± 0.09 R at Dec 17 and 4.98 ± 0.38 R at Dec 18; E strand
     4.0 ± 0.16 R at Dec 14.
   * North (Dec ≥ 21.5) median **0.86 R**: 0.86–1.19 ± 0.16–0.20 R at Dec 25–27.
   * Model-free strict lower limit (a positive residual = fibre minus petal mean ≤ true intensity):
     3.05 ± 0.14 R in a 0.15° bin at (116.27, +17.28), i.e. ≥ 2.8 R (2σ).
   * Systematic: ±15 % from the fibre-acceptance and fibre-area conversion.
   * These values agree with the NSNS-calibrated peaks from the imaging teammate (2.6–3.6 R south, 0.9–1.2 R north).
     The paper's 23.6–23.9 mag arcsec⁻² ≈ 10–13 R is not supported.
3. **The filament's own spectrum is not DIG and not photoionised gas.** The differential (filament-correlated)
   component, from all fibres with healpix-block bootstrap 68 % ranges, gives:
   * **[O III]5007/Hβ = 20 (15–26)**, and > 11.1 at 97.5 %;
   * **[O II]3727/[O III]5007 = 1.49 (1.33–1.75)**; 1.12 (1.03–1.21) for Dec < 19;
   * **R23 = ([O II] + [O III])/Hβ = 58**, > 45 (84 %) and > 35 (97.5 %);
   * [Ne III]3869/[O III] = 0.10 (0.07–0.13), [S III]9069/[O III] = 0.031 (0.024–0.039),
     [O II]7325/[O III] = 0.037 (0.033–0.044), [Ar III]7136/[O III] = 0.018 (0.012–0.025);
   * He II4686/[O III] = −0.013 (−0.029…+0.013), < 0.034 (97.5 %), so He II/Hβ ≲ 0.7;
   * [O III]4363/5007 = 0.056 (0.022–0.077). This gives T_e([O III]) ≈ 3×10⁴ K, but it is only a ~2σ hint and
     the 95 % interval is unbounded below;
   * [O II]3729/3726 = 1.30 (1.24–1.39), so n_e ≈ 110 (40–160) cm⁻³ and < 220 cm⁻³ (97.5 %) (Sanders+16, 10⁴ K);
   * [O III]5007/4959 = 3.11 (3.00–3.28). The atomic value is ≈ 3.0, which serves as a consistency check.

   R23 ≈ 60 is several times above anything photoionisation equilibrium at T_e ~ 10⁴ K produces (R23 ≲ 10).
   DIG/WIM has [O III]/Hβ ≲ 1. These numbers point to collisionally heated gas: an incompletely cooled
   ("truncated") radiative shock whose recombination zone is weak, or a hot-gas/cool-gas interface.
4. **Low-ionisation lines of the filament itself are weaker than the DIG ratios the paper used.**
   * The differential Hα is biased: Hα/Hβ = 1.9 (1.45–2.5) is unphysical. Positive-residual fibres give
     Hα/[O III] = 0.165 (0.13–0.21); negative-residual fibres give −0.02 (−0.08…+0.04). Hβ, [N II] and [S II]
     are symmetric between the two sets.
   * Using Hβ instead of Hα (Hα ≈ 2.9 Hβ), the filament has **[N II]6583/Hα ≈ 0.5–0.7** and
     **[S II]/Hα ≈ 0.3–0.55**. The positive-residual fibres give 0.53 (0.47–0.64) and 0.40 (0.29–0.60) directly.
   * For comparison, LAMOST gives 0.88–1.11 and 0.83–0.93 on the nebula and 0.95 / 0.78 on background P4.
   * The [S II]/Hα > 0.4 SNR criterion is therefore met only marginally by the filament's own gas. The paper's high
     ratios belong to the smooth DIG/WIM, which DESI removes.
5. **Kinematics: no [O III]–red-line velocity offset in DESI. The [O III] velocity differs between strands.**
   * Flux-weighted [O III] V_helio = **+14.5 km/s** (bootstrap 13–16), i.e. V_LSR ≈ +2.4.
   * E strand: +13.0 (12.5–13.5) overall and +10.5 (10.0–11.0) at Dec 15.5–17.5.
   * W strand: +18.5 (17.5–20.0) overall and +17.5 (16.0–18.5) at Dec 15.5–17.5.
   * By Dec band: +15.5 (13–16) at 12–15.5; +11.0 (11–15) at 15.5–16.5; +13.5 (12–15) at 16.5–17.5; +24.0 (19–26)
     at 17.5–19.5; +19.5 (16–22) at 19.5–22; +15.0 (14–17) at 22–25; +11.0 (8–14) at 25–28.6.
   * The Hα-structure component lies at **+15.0 (14–17) km/s**, so [O III] − Hα = −0.5 ± 2 km/s overall.
   * LAMOST [O III] (Table 4: P1 26.8 ± 4.3, P2 26.8 ± 4.1, P3 23.8 ± 3.3) agrees with DESI on the W strand near P2
     (Δ ≈ −3 ± 5 km/s). It is higher than DESI at P1/E strand by ≈ 13–16 km/s (≈ 3σ) and at P3 by ≈ 8 km/s
     (≈ 2.3σ).
   * The paper's [O III]-based V_LSR = 11–16 km/s, and the kinematic distance built on it (a units slip,
     ~1.6 kpc rather than 160 pc per paper_review), is not reproduced. DESI's V_LSR ≈ 0–+6 km/s is what any local
     gas shows.
   * Frame checks:
     - [O III] velocity vs exposure barycentric correction: slope +0.20 ± 0.22 (n = 129, v_corr −17…+30 km/s).
       A missing correction (slope −1) is excluded at 5σ.
     - Stellar cz: DESI main dark/bright vs Gaia ≈ +2.6/+1.2 km/s, via LAMOST LRS (LRS − Gaia = −4.8 ± 0.03).
       The two main/backup tiles observed on 2022-03-16 are off by −19.1 ± 0.6 km/s. No [O III] fibre comes from
       those tiles.
6. **[O III] line width.** The intrinsic FWHM is ≤ 33 km/s: 32.5 km/s (stat. 28–36; 26–38 for a ±2 % LSF error)
   over all fibres, 39 (34–44) for Dec < 19, and 18 (0–27) for the positive-residual (on-ridge) fibres alone.
   DESI's LSF FWHM is 98 km/s at 5008 Å, so this depends on the SPARCL `wave_sigma` LSF being right to 2–3 %. It is
   consistent with, or narrower than, LAMOST's 43–49 km/s. Thermal O II/O III at 10⁴ K adds only 5 km/s.
7. **The Hα structure in the same fibres is ordinary WIM, plus a low-ionisation northern patch.**
   * Hα-correlated component (all fibres): [N II]/Hα = 0.47 (0.40–0.51), [S II]/Hα = 0.64 (0.55–0.71),
     Hα/Hβ = 2.97 (2.75–3.26), [O III]/Hα = 0.14 (0.09–0.18), [S II]6716/6731 = 1.47 (1.42–1.54).
     The [S II] doublet is at the low-density limit.
   * Without the northern patch: [N II]/Hα = 0.38 and [S II]/Hα = 0.54.
   * Patch region (242 fibres): [N II]/Hα = 0.65 (0.50–0.71), [S II]/Hα = 0.83 (0.56–0.91),
     [O III]/Hα = 0.05 (−0.03…+0.10).
   * Stack on the patch (91 fibres): Hα +1.52 R, [S II] +1.63 R, [O III] +0.21 R in fibre residual.
   * So the [O III] filament and the Hα/[S II] structures are physically different components.
   * The northern arc (Dec ≥ 22) has a lower-ionisation [O III]-correlated spectrum: [O II]/[O III] = 6.2 (5.5–6.8),
     [N II]/[O III] = 0.67, [S II]/[O III] = 0.70, Hα/[O III] = 0.68.

### Verdict and hypothesis

DESI supports the "shock" side of Wu+2026's verdict, but with different evidence. The paper's [S II]/Hα and BPT
arguments use red lines that are DIG-dominated. DESI isolates the filament's own spectrum: very high
[O III]/Hβ ≈ 20 together with strong [O II] (R23 ≈ 60), weak He II, and a hint of T_e ≈ 3×10⁴ K. This combination
is characteristic of an incomplete-cooling radiative shock at ~100–150 km/s, of the kind seen in [O III]-bright
filaments of old SNRs. A conductive or turbulent mixing layer at a hot-gas boundary is the main alternative.
The DESI data argue against two of the paper's options:
* a DIG origin for the [O III], since DIG has [O III]/Hβ ≲ 1;
* a photoionised "ionised bubble", since R23 ≈ 60 is impossible at ~10⁴ K and He II is not detected.

The kinematics are those of a structure moving mostly in the plane of the sky: V_LSR ≈ 0–+6 km/s with a 5–8 km/s
E/W strand difference and narrow lines. This fits the tangential limb of a nearby shell, such as the NE rim of the
Monogem Ring proposed by Fesen+2026, better than a ~1.6 kpc object. It is consistent with, but does not establish,
a distance of ≲ 300 pc; DESI gives no distance.

Tests that would decide:
* deeper Hβ, [O III]4363 and [N II]5755 for T_e;
* He II and [Ne V] for the shock speed;
* UV C IV/O VI from the imaging teammate;
* proper motions of the filament across multi-epoch narrowband imaging.

## Results with sources of each number

* Forward model, detection and null tests: `data/derived/spectra_desi_results.json` (`forward_model`, `detection`,
  `null_tests`, `OIII_v_vs_vcorr`).
* Line ratios, velocities, widths and stacks: `data/derived/spectra_desi_ratios.json`, using the runs `model_all`,
  `nbr_all`, `model_south`, `model_mid`, `model_north`, `model_nopatch`, `model_posP`, `model_negP`, `nbrHa_all`,
  `nbrHa_nopatch` and `nbrHa_patch`, plus `v_by_dec`, `v_by_subset` and `stacks`.
  - Bootstrap arrays: `spectra_desi_ratios_boot.npz`.
  - Spectra: `spectra_desi_regspec_full.npz`.
* Derived limits (R23, n_e, ratio limits, ridge peaks, strict lower limits): `data/derived/spectra_desi_summary.json`.
* Per-fibre table (RA, Dec, distance to strand, [O III] in R, free velocity, model prediction, neighbour mean):
  `data/derived/spectra_desi_fiber_table.csv`. All per-fibre line fits are in `spectra_desi_fits.csv` / `.npz`.
* Velocity-frame checks: `data/derived/spectra_desi_rvcheck.json`, `spectra_desi_vcorr.npz`.

Null and validation numbers:
* Line-free windows: amplitude/error MAD-σ = 0.955–0.971; std 0.98–1.09; fraction > 5σ ≤ 0.14 %.
* [O III]5007: std(amplitude/error) = 1.92, with 1.4 % > +5σ and 1.2 % < −5σ.
* The rms significance of [O III] falls from 2.8 at the trace to ≈ 1.0 at 2.5–3° (validation figure). Subtracting
  the ridge model lowers it near the trace. Some unmodelled sub-structure remains: χ²/dof = 7276/4907.
* The 5577 Å airglow line is not correlated with the filament: regression coefficient −0.005 (−0.12…+0.08) per unit
  [O III].
* The on/off stacks are mirror images in [O III]: top-10 % predicted fibres (n = 498) give +1.089 ± 0.017 R; bottom
  10 % give −0.867 ± 0.016 R.


## Data and access (all accessed 2026-11-16)

| What | Where / how | Version | Local copy |
|---|---|---|---|
| Spectra (flux, ivar, Redrock model; wave_sigma for the pilot) | NOIRLab SPARCL REST API, `https://astrosparcl.datalab.noirlab.edu/api/` (`POST /find/`, `POST /spectras/?include=sparcl_id,flux,ivar,model&format=json&dataset_list=DESI-DR1`). I took the endpoints from the sparclclient source (github.com/astro-datalab/sparclclient, `client.py`) and the server's OpenAPI schema (`/api/docs/?format=openapi`, API v14.2). The client's old `/sparc/` base path returns 404. I parsed JSON only; nothing was unpickled. | DESI-DR1 (iron) healpix coadds | ±40 Å windows only: `data/derived/spectra_desi_windows/*.npz` |
| SPARCL id list in RA 111–122, Dec +6…+28 | SPARCL `/find/` | DESI-DR1 | `data/raw/sparcl_desi_dr1/find_box.csv` |
| Redshift catalogue | Data Lab TAP `https://datalab.noirlab.edu/tap/sync`, `desi_dr1.zpix` | DR1 | `data/raw/datalab_tap/zpix_box.csv` |
| Tile / petal / fibre assignments | TAP `desi_dr1.fiberassign`, `desi_dr1.tile`, `desi_dr1.exposure` | DR1 | `data/raw/datalab_tap/{fiberassign,tiles,exposure}_box.csv` |
| Stellar velocity checks | TAP `desi_dr1.mws`; `gaia_dr3.gaia_source` (radial_velocity) via Data Lab; LAMOST DR12 LRS/MRS catalogues downloaded by the SDSS+LAMOST teammate (`data/raw/lamost_dr12/`, read only) | DR1 / DR3 / DR12 | `data/raw/datalab_tap/mws_box.csv`, `gaia_dr3_rv_box_withrv.csv` |
| Pipeline facts (vacuum, barycentric, 1.52″ mean fibre, sky model) | desispec source on GitHub (`py/desispec/scripts/proc.py`, `scripts/extract.py`, `heliocentric.py`, `sky.py`, `skygradpca.py`, `skymag.py`) | main branch | `data/raw/desispec_src/*.txt` (saved as text, never executed) |

Transfer: about 0.2 MB per spectrum (gzip JSON with three float arrays). 5,220 spectra were downloaded (1.06 GB compressed transfer, logged per batch in `data/derived/spectra_desi_windows/*_log.jsonl`), and only the ±40 Å windows were kept (2,186 pixels per spectrum, 122 MB on disk). The download queue (`spectra_desi_sel_queue.csv`, 21,237 ids; second process `spectra_desi_sel_queue2.csv`) is resumable with `scripts/spectra_desi_fetch.py`. It was stopped at 18:45 JST for the final analysis.

Not reachable or not public, recorded as data gaps:
* `data.desi.lbl.gov` and `legacysurvey.org` (NERSC) are blocked by the proxy. The raw DESI sky-fibre spectra and
  sky models (`sky-*.fits`, `cframe-*.fits`) are not served by SPARCL or Data Lab, so the actual petal sky model
  could not be inspected. Saydjari et al. 2026 (arXiv:2608.28740) used these files.
* `sparclclient` cannot be installed. Its REST calls were reimplemented in `scripts/spectra_desi_sparcl.py`.
* SPARCL throughput was 0.1–1.5 spectra/s through the proxy, with frequent `RemoteDisconnected` errors and
  automatic retries. This rate, not the 15 GB cap, limited the sample.

## Method

1. **Selection** (`scripts/spectra_desi_select.py`). I joined the SPARCL ids to `zpix`. I kept GALAXY/QSO coadds
   with ZWARN = 0, because stars have Balmer/Mg b model mismatch at z ≈ 0. I computed distances to the E/W/W2
   strands of the geometry teammate's trace (`silk_trace.json` v1.0). The selections were:
   * pilot: the Dec 16.1–16.7 strip through P1 and P4, plus 0.2° discs on P2 and P3;
   * strip: all fibres with d < 1° and a random 1–3° subsample;
   * control: fibres more than 5° from every strand.

   I downloaded them in the priority order of `spectra_desi_sel_queue.csv` (`scripts/spectra_desi_fetch.py`).
   Only a part of the queue was obtained (see Limitations).
2. **Residuals** (`scripts/spectra_desi_measure.py`). For each spectrum, r = flux − Redrock model (SPARCL `model`).
   * Fit a local linear continuum, masking ±6 Å around every catalogued line and the object's own redshifted lines.
   * Fit Gaussian lines with σ² = σ_LSF² + (20 km/s)², using σ_LSF from SPARCL `wave_sigma`
     (`spectra_desi_lsf.json`; 0.698 Å = 41.8 km/s at 5008 Å). Doublets are tied ([O III] 1:2.98, [N II] and
     [O I] 1:3.05).
   * Each line is fitted two ways: at fixed velocity (LAMOST: +25 km/s for [O III] and He II, +17 km/s for the
     rest) and at free velocity on a −150…+150 km/s grid, with the error taken from the χ² curvature.
   * Four line-free null windows (4885, 4990, 5030, 6615 Å) calibrate the error scale.
     Their median MAD-σ of amplitude/error is 0.95, so the ivar errors are realistic.
3. **Surface brightness** (`scripts/spectra_desi_analyze.py`). DESI fluxes are calibrated for point sources.
   * In-fibre flux of uniform emission = calibrated flux × MEAN_PSF_TO_FIBER_SPECFLUX (median 0.77).
   * A 1.52″ fibre is 1.815 arcsec². At 5008 Å, 1 R = 7.43×10⁻¹⁸ erg s⁻¹ cm⁻² arcsec⁻², so 10⁻¹⁷ erg s⁻¹ cm⁻²
     in the fibre corresponds to 0.743 R.
4. **Petal forward model** (`scripts/spectra_desi_petalmodel.py`). DESI subtracts a per-exposure, per-petal sky
   built from that petal's sky fibres.
   * Model sky: Gaussian ridges along E, W and W2. Widths come from the trace HWHM, and amplitudes are free at
     1° Dec nodes with linear interpolation between them.
   * Each fibre's predicted residual is T(x_i) − Σ_t w_t ⟨T⟩_(tile t, petal p).
   * ⟨T⟩ is the mean over all fibres of that tile-petal in `desi_dr1.fiberassign`. Sky-fibre positions are not
     public; this mean is a proxy.
   * The weights w_t are proportional to the summed EFFTIME_SPEC of each tile.
   * The fit is a weighted least-squares fit with 4σ clipping.
   * Null test: refit the same model with the whole trace shifted by ±1…4° in RA.
5. **Differential line ratios** (`scripts/spectra_desi_regspec.py`, `scripts/spectra_desi_ratios.py`). The
   continuum-subtracted residual spectra c_i(λ) (in-fibre R Å⁻¹) are combined into a regression spectrum
   S(λ) = Σ w_i P_i c_i / Σ w_i P_i².
   * The instrument P_i is the expected [O III] residual of fibre i (R). It is either the forward-model prediction
     (`model`) or the mean [O III] residual of the other fibres within 0.12° (`nbr`).
   * Integrating a line in S gives (line / [O III]5007) for the emission that is spatially correlated with the
     filament. Smooth DIG and sky cancel, because the same petal-mean subtraction applies to every line.
   * An Hα-structure instrument (neighbour-mean Hα residual) gives the spectrum of whatever Hα structure exists
     within the petals.
   * Lines are fitted at the DESI [O III] velocity.
   * Uncertainties come from a bootstrap over DESI nside-64 healpix blocks (~0.9°, 200–400 resamples). This
     captures per-petal correlated sky-subtraction errors that ivar does not.
6. **Velocity frame checks** (`scripts/spectra_desi_rvcheck.py`; v-vs-v_corr in `spectra_desi_analyze.py`).
   * desispec applies `--barycentric-correction` at extraction (`proc.py` → `extract.py`).
   * Empirically, I regressed per-fibre [O III] velocity on the EFFTIME-weighted barycentric correction of the
     contributing exposures (astropy, KPNO).
   * I compared DESI stellar cz with Gaia DR3 and with LAMOST LRS/MRS for stars in the same box.

7. **Regression with intercept.** For all runs over the full sample, I fit c_i(λ) = a(λ) + S(λ) P_i per pixel.
   The intercept absorbs any additive offset common to all fibres. The sign-split runs (`model_posP`, `model_negP`)
   use S = Σ w P c / Σ w P² without an intercept, because P has a single sign within each subset.

## Figures

| File | Caption |
|---|---|
| `figures/spectra_desi_coverage.png` | DESI DR1 coadd density in the box, with the 5,220 analysed fibres, the trace and P1–P4. |
| `figures/spectra_desi_oiii_map.png` | [O III]5007 residual (R) per fibre, binned 0.2°, ridge-model prediction after petal subtraction, and data − model, Dec 12–28.6. |
| `figures/spectra_desi_oiii_map_zoom.png` | The same around P1–P4 (Dec 13.8–19.2, 0.1° bins), showing the +/− petal sign pattern across the E and W strands. |
| `figures/spectra_desi_profile.png` | [O III] residual vs RA offset from the E spine in seven Dec bands, with the forward-model curve. |
| `figures/spectra_desi_halpha_map.png` | Binned Hα, [N II]6583 and [S II] residuals, with the NSNS Hα patch outline. |
| `figures/spectra_desi_regspec.png` | Regression spectrum of the [O III]-correlated component in 12 line windows, with bootstrap band and fixed-velocity fits (photon units per unit [O III]). |
| `figures/spectra_desi_ratios.png` | Forest plot of differential line ratios (energy, bootstrap 68/95 %) for the [O III] component (all, neighbour instrument, south, north) and the Hα-structure component, with LAMOST P1–P4. |
| `figures/spectra_desi_velocity.png` | [O III] V_helio vs Dec (DESI bands, fibres, LAMOST P1–P3, LAMOST red-line range, WHAM Hα, DESI Hα structure); right: velocity vs barycentric correction, with stellar RV checks. |
| `figures/spectra_desi_ridge_sb.png` | Forward-model ridge peak [O III] (R) vs Dec for the E/W strands, NSNS peak excess and Wu+2026 Table 3. |
| `figures/spectra_desi_stacks.png` | Stacks of the top-10 % predicted, bottom-10 % predicted and Hα-patch fibres around [O III], Hβ, Hα+[N II] and [S II]. |
| `figures/spectra_desi_validation.png` | Significance histograms (null windows vs [O III]), measured vs predicted [O III] and 5577 airglow, and rms significance vs distance from the trace. |

## Scripts (run in this order; all with `~/.venvs/astro313/bin/python -I`)

`spectra_desi_tap.py` (TAP downloads) → `spectra_desi_sparcl.py` (`find`) → `spectra_desi_select.py`
(catalogue and selections) → `spectra_desi_fetch.py` (SPARCL windows) → `spectra_desi_lsf.py` →
`spectra_desi_measure.py` → `spectra_desi_analyze.py` (forward model, petal model in `spectra_desi_petalmodel.py`) →
`spectra_desi_regspec.py` (helpers) / `spectra_desi_ratios.py` → `spectra_desi_summary.py` →
`spectra_desi_rvcheck.py` → figures: `spectra_desi_figs2.py all`, `spectra_desi_figs.py pilotmap` (zoom map),
`spectra_desi_fig_coverage.py` (older coverage version). Line list and grid: `spectra_desi_lines.py`.

## Limitations

* **Sample.** Only 5,220 of the 21,237 queued spectra (and of ~275,000 DR1 coadds in the box) were obtained,
  because SPARCL through the proxy delivered about 1 spectrum/s.
  * Coverage along the trace is patchy (coverage figure). P3 and Dec 12–15 are sparse.
  * No control fibres more than 3° from the trace were downloaded. The far-field null comes from line-free windows,
    the displaced-trace test and the decline of rms significance with distance (to ≈ 1.0 at 2.5–3°), not from a
    > 5° control sample.
  * Running `spectra_desi_fetch.py` again resumes the queue.
* **Differential measurement only.** Emission that is uniform on petal scales (≳ 1°) is invisible, so DESI cannot
  give the total line intensity of smooth components.
  * Ridge peaks depend on the assumed Gaussian cross-section (trace HWHM). The petal mean is approximated by the
    mean over the tile-petal's science fibres, because sky-fibre positions are not public.
  * χ²/dof = 1.48 shows unmodelled structure.
  * Some nodes are degenerate where the E and W strands merge (e.g. E@19 = 0.06 ± 0.23 R).
* **Hα is the least reliable differential line.** Its unphysical Hα/Hβ and its sign asymmetry suggest sky residuals
  (geocoronal Hα and OH 6553/6568/6577 Å) or Hα structure that does not follow [O III]. I therefore quote
  Hβ-based ratios as primary. [O I]6300 is blended with airglow, which is not separable at DESI resolution.
  It is reported but not interpreted.
* **Line widths** depend on the 2–3 % accuracy of the SPARCL `wave_sigma` LSF, measured from 75 pilot spectra.
  Velocity structure between strands broadens the stacked profiles.
* **[O II]3726/3729 deblending** at R ≈ 2,700 is sensitive to LSF and velocity. n_e differs between subsets: the
  on-ridge ratio is 1.45, near the low-density limit, while the off-ridge ratio is 1.12.
* **Flux scale.** Converting point-source calibrated flux to surface brightness uses PSF_TO_FIBER_SPECFLUX and a
  1.52″ mean fibre. The fibre size varies by about ±3 % across the focal plane. I estimate ±15 % systematic.
* **Velocity zero-point.** DESI dark/bright stars agree with Gaia to ≈ +1 to +3 km/s via LAMOST LRS. The −19 km/s
  offset of two backup tiles (one night) is unexplained and is a warning that tile-level frame problems exist in
  DR1, although none is seen in the [O III] v-vs-v_corr test.
* All ratios are **observed**. Reddening toward b = +19° is small: E(B−V) ≈ 0.03–0.05 changes [O II]/[O III] by
  ≲ 5 %.

## Data gaps

* DESI sky-fibre spectra and per-petal sky models (`sky-*.fits`, `cframe-*.fits`) are at NERSC
  (`data.desi.lbl.gov`, blocked) and not in SPARCL. They would allow absolute [O III] maps like those of
  Saydjari et al. 2026, arXiv:2608.28740.
* `legacysurvey.org` is blocked. It was not needed here.
* The remaining 16,000 queued SPARCL spectra, including controls more than 5° from the trace, were not fetched.
* Sergey Koposov's reported DESI detection: no primary source was found (see `notes/discovery.md`). This analysis
  shows the signal is easily detectable in public DR1 data.

