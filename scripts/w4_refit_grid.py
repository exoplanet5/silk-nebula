"""w4_refit_grid: re-run the verify_shockmodel 1-D truncated/complete radiative-shock grid (code imported read-only)
on an extended grid, and add the UV lines that wave 2 did not model (full GALEX FUV and NUV line lists + H 2-photon),
computed from the flow's own NEI ion fractions with the verify_altmodel CHIANTI v7.1.4 line toolkit (read-only).

Physics and approximations are exactly those of scripts/verify_shockmodel_flow.py (no radiative transfer; truncated
runs stop at 1e4 K; complete runs hold 8000 K until x(H+) < 0.01).  Tracks are computed at n0 = 0.1 cm^-3 and are
scale-free in n0 at fixed b = B0/sqrt(n0) (s, t ~ 1/n0; fluxes ~ n0).

Grid: v_s = 60..200 km/s (13 values) x pre-shock H0 fraction X = 0 / 0.05 / 0.3 / 1 x b = 1,2,3,5,7,10,14,20,30,45
(uG cm^1.5), solar (Asplund+09) abundances; plus complete runs (b = 3, 10; X = 0, 1).
Outputs: data/derived/w4_refit_tracks.npz (per-run resampled tracks: wave-2 keys + F_uv_<line> cumulative UV fluxes),
         data/derived/w4_refit_grid.json (run metadata), data/derived/w4_refit_grid_log.txt
Run: ~/.venvs/astro313/bin/python -I scripts/w4_refit_grid.py
"""
import sys, json, time, itertools, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_fit as FIT          # read-only import: run(), resample(), GROUPS
import verify_shockmodel_flow as F
import verify_altmodel_lines as AL           # read-only import: synth(), e_line tables (CHIANTI v7.1.4)

ROOT = '/Users/mickey/astronomy/silk-nebula'
N0_REF = FIT.N0_REF
ZOF = {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Ne': 10, 'Mg': 12, 'Si': 14, 'S': 16}
UV_LINES = ['OIII5007', 'CIV1549', 'OIII1663', 'NIII1750', 'NIV1486', 'OIV1400', 'SiIV1400', 'NV1240', 'HeII1640',
            'CIII1909', 'CII2326', 'NeIV2423', 'MgII2800', 'SiIII1892', 'OII2470', 'OII7325', 'OII3726', 'OII3729',
            'Hb', '2g']


def uv_cumulative(sh):
    """Cumulative face-on UV line fluxes (erg cm^-2 s^-1; '2g' in decays cm^-2 s^-1) along the run, from the flow's
    ion fractions and the altmodel CHIANTI emissivities."""
    fr = {ZOF[e]: sh.f[:, sh.net.slices[e]].T for e in ZOF}
    L = AL.synth(sh.T, fr, ne_nH=sh.ne / sh.nH, ne=np.maximum(sh.ne, 1e-3))
    out = {}
    for k in UV_LINES:
        j = np.asarray(L[k]) * sh.nH**2
        seg = 0.5 * (j[:-1] + j[1:]) * np.diff(sh.s)
        out[k] = np.concatenate([[0.0], np.cumsum(seg)])
    return out


if __name__ == '__main__':
    t0 = time.time()
    log = open(ROOT + '/data/derived/w4_refit_grid_log.txt', 'w')
    net = F.Network('asplund09')
    VS = [60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200]
    XN = [0.0, 0.05, 0.3, 1.0]
    BB = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 14.0, 20.0, 30.0, 45.0]
    runs = [dict(v=v, xn=x, b=b, complete=False) for v, x, b in itertools.product(VS, XN, BB)]
    runs += [dict(v=v, xn=x, b=b, complete=True) for v, x, b in itertools.product(VS, (0.0, 1.0), (3.0, 10.0))]
    tracks, meta = {}, []
    for k, rp in enumerate(runs):
        vA = rp['b'] * 1e-6 / np.sqrt(4 * np.pi * net.mass_per_H) / 1e5
        if rp['v'] < 1.5 * vA:
            continue
        try:
            sh, j, C, colC4, colH0, _ = FIT.run(net, rp['v'], rp['xn'], rp['b'], rp['complete'])
        except Exception as e:  # noqa
            print('FAILED', rp, repr(e), file=log, flush=True)
            continue
        tag = f"v{rp['v']}_x{rp['xn']}_b{rp['b']}_{'cmp' if rp['complete'] else 'trn'}"
        tr = FIT.resample(sh, j, C, colC4, colH0)
        uv = uv_cumulative(sh)
        for kk, arr in uv.items():
            tr['F_uv_' + kk] = np.interp(tr['N'], sh.N, arr)
        # C IV- and [O III]-weighted cumulative depth (for FUV-ridge offsets in a flat sheet)
        for kk, arr in tr.items():
            tracks[f'{tag}__{kk}'] = arr
        meta.append(dict(tag=tag, **rp, vA=float(vA), Ts=float(sh.Ts), x1=float(sh.x1), N_end=float(sh.N[-1]),
                         t_end_n0=float(sh.t[-1] * N0_REF), s_end_n0=float(sh.s[-1] * N0_REF), T_end=float(sh.T[-1]),
                         xH_end=float(sh.frac('H', 1)[-1])))
        print(f"{k+1}/{len(runs)} {tag} Ts={sh.Ts:.3g} N_end={sh.N[-1]:.3g} T_end={sh.T[-1]:.3g} "
              f"O3synth/O3flow={uv['OIII5007'][-1]/C['OIII5007'][-1]:.3f} ({time.time()-t0:.0f}s)", file=log, flush=True)
    np.savez_compressed(ROOT + '/data/derived/w4_refit_tracks.npz', **tracks)
    json.dump({'n0_ref': N0_REF, 'grid': dict(VS=VS, XN=XN, BB=BB), 'uv_lines': UV_LINES, 'runs': meta,
               'note': 'tracks from verify_shockmodel_flow (imported read-only); UV lines from verify_altmodel_lines.synth'},
              open(ROOT + '/data/derived/w4_refit_grid.json', 'w'), indent=1)
    print('done', len(meta), 'runs', f'{time.time()-t0:.0f}s', file=log, flush=True)
