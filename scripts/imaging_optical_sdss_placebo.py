"""Global placebo / injection test for the SDSS sky-level detection.
Refit all core runs (ridge crossing Dec 12-20.5, bright-star fields excluded, fields within WIN deg of the ridge) with
the NSNS [OIII] 'ridge' template (tapered to 0 beyond 1.2 deg) shifted in RA by s/cos(dec), s = -3..+3 deg, and combine
per-run amplitudes (formal-error weights). The same is done for simulated data = 1.0 x unshifted template (times the
expected A) + the actual per-run nuisance structure removed (i.e. pure signal), giving the curve a real signal at
zero shift must produce; the observed curve is compared to it. Usage: python ... [WIN]
Output: data/derived/imaging_optical_sdss_placebo.json
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import numpy as np
from imaging_optical_common import DER, save_json, along_across
from imaging_optical_nsns import ridge_spine
from imaging_optical_sdss_sky import load_fields, template_sampler, bright_star_fields, fit_run, EXPECT

WIN = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
a = load_fields(); sra, sdec = ridge_spine(); bad, _ = bright_star_fields(a)
across, along, beyond = along_across(a['ra'], a['dec'], sra, sdec)
samp = template_sampler('oiii', variant='ridge')
runs = []
for run in np.unique(a['run']):
    s = (a['run'] == run) & (np.abs(across) <= WIN) & ~beyond & (a['quality'] >= 2) & ~bad
    acr = across[s]
    if s.sum() < 40 or not ((acr < -0.7).sum() >= 6 and (acr > 0.7).sum() >= 6 and (np.abs(acr) < 0.3).sum() >= 4):
        continue
    dc = np.median(a['dec'][s][np.abs(acr) < 0.3])
    if 12 <= dc < 20.5:
        runs.append((run, s))
out = {'window_deg': WIN, 'n_runs': len(runs), 'shifts_deg': [], 'A_over_expected': [], 'err_over_expected': [], 'injected_A_over_expected': []}
for sh in np.round(np.arange(-3.0, 3.01, 0.2), 2):
    A, E, Ai = [], [], []
    for run, s in runs:
        sub = np.zeros(s.sum(), dtype=[('field', 'i4'), ('camcol', 'i4'), ('sky', 'f8')])
        sub['field'] = a['field'][s]; sub['camcol'] = a['camcol'][s]; sub['sky'] = a['sky_g'][s]
        T = samp(a['ra'][s] + sh / np.cos(np.radians(a['dec'][s])), a['dec'][s])
        r = fit_run(sub, T)
        T0 = samp(a['ra'][s], a['dec'][s])
        sim = sub.copy(); sim['sky'] = 1.0 + EXPECT['g'] * np.nan_to_num(T0) + 1e-5 * (sub['field'] % 7)
        ri = fit_run(sim, T, clip=99.0)
        if r is not None and np.isfinite(r['A']) and r['T_range'] > 0.2:
            A.append(r['A']); E.append(r['eA_formal']); Ai.append(ri['A'] if ri else np.nan)
    A = np.array(A); E = np.array(E); Ai = np.array(Ai); w = 1 / E ** 2
    m = np.sum(w * A) / np.sum(w); e = 1 / np.sqrt(np.sum(w)); mi = np.nansum(w * Ai) / np.sum(w[np.isfinite(Ai)])
    out['shifts_deg'].append(float(sh)); out['A_over_expected'].append(float(m / EXPECT['g'])); out['err_over_expected'].append(float(e / EXPECT['g']))
    out['injected_A_over_expected'].append(float(mi / EXPECT['g']))
    print(f'{sh:+.1f}  obs {m / EXPECT["g"]:+.2f} +- {e / EXPECT["g"]:.2f}   injected(1.0) {mi / EXPECT["g"]:+.2f}  n={len(A)}')
x = np.array(out['A_over_expected']); s_ = np.array(out['shifts_deg'])
xi = np.array(out['injected_A_over_expected'])
off = np.abs(s_) >= 2.0          # template (taper 1.2 deg) no longer overlaps the true ridge position
out['placebo_mean_absshift_ge2'] = float(np.mean(x[off])); out['placebo_std_absshift_ge2'] = float(np.std(x[off]))
out['zero_shift'] = float(x[s_ == 0][0])
# best scale of the injected curve to the observed curve (shape match), and residual chi2
ee = np.array(out['err_over_expected'])
k = float(np.sum(x * xi / ee ** 2) / np.sum(xi ** 2 / ee ** 2))
out['curve_scale_obs_over_injected'] = k
out['curve_rms_resid'] = float(np.sqrt(np.mean((x - k * xi) ** 2)))
print({k: v for k, v in out.items() if not isinstance(v, list)})
save_json(out, DER / f'imaging_optical_sdss_placebo_w{WIN:.1f}.json')
