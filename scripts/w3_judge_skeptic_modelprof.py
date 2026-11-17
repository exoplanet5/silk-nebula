"""w3_judge_skeptic: what cross-filament ordering does the consensus truncated-shock model actually predict,
measured with the SAME estimators used on the data in w3_judge_skeptic_tests.py?

For each of the 25 'consistent' southern models in verify_shockmodel_summary.json (v_s 100-150 km/s), the
cumulative photon fluxes F_O3 (5007), F_O2 (3727 sum) and F_CIV (1549) in verify_shockmodel_tracks.npz are
projected onto a curved sheet exactly as verify_shockmodel_analysis.profile() does (chord z^2 = 2 R_c (h - s),
R_c = 280 pc * 11.9 deg, d = 280 pc, n0 = the model's width-closure n0), converted to sky offset x = x_front - h
(east positive, the front is the eastern edge), then convolved with a Gaussian kernel that brings the [O III]
FWHM to the observed stacked NSNS FWHM (0.41 deg for Dec 12-15; trace jitter + along-spine structure).
Then: split-normal fit of [O III] (asymmetry sigma_W/sigma_E), template-shift of C IV, of an FUV mix
(C IV scattered + O III] 1664, the latter assumed to follow 5007), and of [O II] relative to [O III].

Output: data/derived/w3_judge_skeptic_modelprof.json
Run: ~/.venvs/astro313/bin/python -I scripts/w3_judge_skeptic_modelprof.py
"""
import sys
import json
import numpy as np

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from verify_xrayfuv_common import DER
from w3_judge_skeptic_tests import fit_split, shift_fit, splitnorm

PC = 3.0857e18
N0REF = 0.1
D_PC = 280.0
RC = D_PC * np.radians(11.9) * PC
FWHM_OBS_STACK = {'dec12_15': 0.406, 'dec12_17p6': 0.528}   # from w3_judge_skeptic_tests.json NSNS fits


def proj(s_ref, Fc_ref, istar, n0):
    s = s_ref[:istar + 1] * N0REF / n0
    Fc = Fc_ref[:istar + 1] * n0 / N0REF
    ds = np.diff(s)
    j = np.diff(Fc) / np.where(ds > 0, ds, np.inf)
    sm = 0.5 * (s[1:] + s[:-1]); sstar = s[-1]
    hs = np.linspace(0, 60 * sstar, 6000)[1:]
    I = np.zeros(hs.size)
    for k, h in enumerate(hs):
        u = np.linspace(np.sqrt(max(h - sstar, 0.0)), np.sqrt(h), 300)
        sv = h - u ** 2
        jj = np.interp(sv, sm, j, left=j[0], right=0.0)
        jj[sv > sstar] = 0.0
        I[k] = 2 * np.sqrt(2 * RC) * np.trapezoid(jj, u)
    return hs, I, sstar


def main():
    summ = json.loads((DER / 'verify_shockmodel_summary.json').read_text())
    T = np.load(DER / 'verify_shockmodel_tracks.npz')
    xg = np.arange(-1.5, 1.5001, 0.005)
    out = []
    for m in summ['consistent']:
        tag = f"v{m['v']}_x{m['xn']}_b{m['b']}_asplund09_trn"
        N = T[tag + '__N']; s_ref = T[tag + '__s']
        istar = int(np.argmin(np.abs(np.log10(np.maximum(N, 1)) - m['logN'])))
        prof = {}
        for key in ('F_O3', 'F_O2', 'F_CIV'):
            hs, I, sstar = proj(s_ref, T[tag + '__' + key], istar, m['n0'])
            hdeg = np.degrees(hs / (D_PC * PC))
            prof[key] = (hdeg, I)
        # put the front at x = 0.15 deg east of the [O III] peak (any constant; shifts are relative)
        rec = dict(tag=tag, v=m['v'], b=m['b'], xn=m['xn'], logN=m['logN'], sstar_deg=float(np.degrees(sstar / (D_PC * PC))))
        hpk = prof['F_O3'][0][np.argmax(prof['F_O3'][1])]
        rec['h_peak_O3_deg'] = float(hpk)
        for seg, fw in FWHM_OBS_STACK.items():
            # intrinsic profiles on the sky grid (x east-positive; front at x = hpk so the [O III] peak sits at 0)
            sky = {}
            for key, (hd, I) in prof.items():
                sky[key] = np.interp(hpk - xg, hd, I, left=0.0, right=0.0)
            # FWHM of intrinsic [O III]
            f0 = fit_split(xg, sky['F_O3'], xmax=1.2)
            fw0 = f0['fwhm'] if f0 else np.nan
            sig_k = np.sqrt(max(fw ** 2 - fw0 ** 2, 1e-6)) / 2.3548 if np.isfinite(fw0) else 0.0
            from scipy.ndimage import gaussian_filter1d
            conv = {k: gaussian_filter1d(v, sig_k / 0.005, mode='constant') for k, v in sky.items()}
            res = {}
            for label_, arr in (('intrinsic', sky), ('convolved', conv)):
                fO3 = fit_split(xg, arr['F_O3'], xmax=1.2)
                tmpl = splitnorm([fO3['A'], fO3['x0'], fO3['sig_w'], fO3['sig_e'], 0, 0], xg)
                d_c4, _ = shift_fit(xg, arr['F_CIV'] / arr['F_CIV'].max(), xg, tmpl / tmpl.max(), xmax=1.2)
                d_o2, _ = shift_fit(xg, arr['F_O2'] / arr['F_O2'].max(), xg, tmpl / tmpl.max(), xmax=1.2)
                # FUV mix: scattered C IV (fraction fC of the FUV) + O III] 1664 following 5007
                mixes = {}
                for fC in (0.3, 0.5, 0.7):
                    mix = fC * arr['F_CIV'] / arr['F_CIV'].sum() + (1 - fC) * arr['F_O3'] / arr['F_O3'].sum()
                    mixes[f'fC{fC}'] = shift_fit(xg, mix / mix.max(), xg, tmpl / tmpl.max(), xmax=1.2)[0]
                res[label_] = dict(O3_fit=fO3, shift_CIV=d_c4, shift_O2=d_o2, shift_FUVmix=mixes)
            res['kernel_sigma_deg'] = float(sig_k); res['intrinsic_O3_fwhm'] = float(fw0)
            rec[seg] = res
        out.append(rec)
        print(tag, json.dumps({k: rec['dec12_15']['convolved'][k] for k in ('shift_CIV', 'shift_O2', 'shift_FUVmix')}),
              'asym', round(rec['dec12_15']['convolved']['O3_fit']['asym_w_over_e'], 2),
              'intrinsic asym', round(rec['dec12_15']['intrinsic']['O3_fit']['asym_w_over_e'], 2), flush=True)

    def rng_(seg, lab, key, sub=None):
        vals = []
        for r in out:
            v = r[seg][lab][key] if sub is None else r[seg][lab][key][sub]
            vals.append(v)
        vals = np.array(vals, float)
        return dict(min=float(vals.min()), median=float(np.median(vals)), max=float(vals.max()))
    summary = {}
    for seg in FWHM_OBS_STACK:
        summary[seg] = {}
        for lab in ('intrinsic', 'convolved'):
            d = {}
            for key in ('shift_CIV', 'shift_O2'):
                d[key] = rng_(seg, lab, key)
            for fC in ('fC0.3', 'fC0.5', 'fC0.7'):
                d['shift_FUVmix_' + fC] = rng_(seg, lab, 'shift_FUVmix', fC)
            asy = np.array([r[seg][lab]['O3_fit']['asym_w_over_e'] for r in out])
            d['O3_asym_w_over_e'] = dict(min=float(asy.min()), median=float(np.median(asy)), max=float(asy.max()))
            summary[seg][lab] = d
    (DER / 'w3_judge_skeptic_modelprof.json').write_text(json.dumps(dict(models=out, summary=summary), indent=1, default=float))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
