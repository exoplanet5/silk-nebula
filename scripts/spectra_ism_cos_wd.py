#!/usr/bin/env python
"""spectra_ism: HST/COS G130M spectrum of WD J074735.98+210635.83 (0.16 deg from the [O III] spine, Dec +21.1).

Product: lfac2c010_x1dsum.fits (HST program 17420), https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/lfac2c010_x1dsum.fits
Distance: Gaia DR3 parallax from VizieR I/355/gaiadr3 (CfA mirror), 5" cone.
Measures apparent-optical-depth (AOD) profiles of interstellar Si II 1260, C II 1334, C II* 1335.7, Si III 1206,
Si IV 1393/1402, N V 1238/1242, S II 1253/1259, O I 1302 in V_helio; column densities in velocity windows.
Outputs: data/derived/spectra_ism_cos_wd.json, figures/spectra_ism_cos_wd.png
Run: ~/.venvs/astro313/bin/python -I scripts/spectra_ism_cos_wd.py
"""
import json
import sys

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter

sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import RAW, DER, FIG, fetch  # noqa: E402
from spectra_ism_vizier_parse import read_asu, fnum  # noqa: E402
from silk_style import new_sky_figure, finish, label, OIII_C, HA_C, FGDIM, WARN, OK  # noqa: E402

C = 299792.458
OUT = RAW / 'spectra_ism_mast'
RA0, DE0 = 116.89992, 21.10995
fx = OUT / 'lfac2c010_x1dsum.fits'
if not fx.exists():
    fetch('https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/lfac2c010_x1dsum.fits', fx, timeout=300,
          min_bytes=10000)
fg = OUT / 'gaia_dr3_wdj0747.tsv'
if not fg.exists():
    fetch('https://vizier.cfa.harvard.edu/viz-bin/asu-tsv', fg,
          params={'-source': 'I/355/gaiadr3', '-c': f'{RA0} {DE0:+}', '-c.rs': '6', '-out': 'Source,RA_ICRS,DE_ICRS,Plx,e_Plx,Gmag,BP-RP,pmRA,pmDE',
                  '-oc.form': 'd'}, timeout=120)
g = list(read_asu(fg).values())[0] if fg.exists() else []
gaia = g[0] if g else {}
plx, eplx = fnum(gaia.get('Plx')), fnum(gaia.get('e_Plx'))

h = fits.open(fx)
hdr0, d = h[0].header, h[1].data
w = np.concatenate([np.asarray(r['WAVELENGTH'], float) for r in d])
f = np.concatenate([np.asarray(r['FLUX'], float) for r in d])
e = np.concatenate([np.asarray(r['ERROR'], float) for r in d])
dq = np.concatenate([np.asarray(r['DQ_WGT'], float) for r in d])
o = np.argsort(w)
w, f, e, dq = w[o], f[o], e[o], dq[o]
good = (dq > 0) & np.isfinite(f) & (e > 0)
w, f, e = w[good], f[good], e[good]
# COS x1d wavelengths are heliocentric (vacuum)
LINES = {'SiIII1206': (1206.500, 1.67), 'NV1238': (1238.821, 0.156), 'NV1242': (1242.804, 0.0780),
         'SII1253': (1253.805, 0.0109), 'SiII1260': (1260.422, 1.18), 'OI1302': (1302.168, 0.0480),
         'CII1334': (1334.532, 0.128), 'CII*1335': (1335.708, 0.115), 'SiIV1393': (1393.760, 0.513), 'SiIV1402': (1402.773, 0.254)}
res = {'target': 'WD J074735.98+210635.83', 'ra': RA0, 'dec': DE0, 'gaia_dr3': gaia, 'd_pc': 1000 / plx if plx > 0 else None,
       'd_pc_err': 1000 / plx ** 2 * eplx if plx > 0 else None, 'cos_header': {k: str(hdr0.get(k)) for k in
       ('TARGNAME', 'PROPOSID', 'OPT_ELEM', 'CENWAVE', 'DATE-OBS', 'EXPTIME')}, 'lines': {}}
WIN = {'local_m20_p20': (-20, 20), 'oiii_p15_p40': (15, 40), 'neg_m100_m25': (-100, -25), 'pos_p40_p100': (40, 100)}
prof = {}
for name, (lam0, fval) in LINES.items():
    v = (w / lam0 - 1) * C
    m = (v > -450) & (v < 450)
    if m.sum() < 50:
        continue
    vv, ff, ee = v[m], f[m], e[m]
    # continuum: running median over 250 km/s with the line core (|v| < 120) excluded, then linear fit across it
    cmask = (np.abs(vv) > 120) & (np.abs(vv) < 420)
    if cmask.sum() < 20:
        continue
    pcoef = np.polyfit(vv[cmask], median_filter(ff, 25)[cmask], 2)
    cont = np.polyval(pcoef, vv)
    nf = ff / cont
    ne = ee / cont
    tau = -np.log(np.clip(nf, 0.03, None))
    dv = np.gradient(vv)
    Na = 3.768e14 * tau / (fval * lam0)       # cm^-2 (km/s)^-1
    out = {'lam0': lam0, 'f': fval, 'snr_cont': float(np.nanmedian(1 / ne[cmask]))}
    for wn, (v1, v2) in WIN.items():
        mm = (vv >= v1) & (vv < v2)
        N = float(np.sum(Na[mm] * dv[mm]))
        eN = float(np.sqrt(np.sum((3.768e14 / (fval * lam0) * ne[mm] / np.clip(nf[mm], 0.03, None) * dv[mm]) ** 2)))
        ew = float(np.sum((1 - nf[mm]) * dv[mm]) * lam0 / C * 1000)
        out[wn] = {'logN': round(float(np.log10(N)), 2) if N > 0 else None, 'N': N, 'eN': eN, 'EW_mA': round(ew, 1),
                   'sat_frac': float(np.mean(nf[mm] < 0.1))}
    cen = (vv > -60) & (vv < 80)
    out['v_centroid_tau'] = float(np.sum(tau[cen] * vv[cen]) / np.sum(tau[cen])) if np.sum(tau[cen]) > 0 else None
    res['lines'][name] = out
    prof[name] = (vv, nf, ne)
json.dump(res, open(DER / 'spectra_ism_cos_wd.json', 'w'), indent=1, default=float)
print(json.dumps({k: res[k] for k in ('d_pc', 'd_pc_err', 'cos_header')}, default=float))
for k, v in res['lines'].items():
    print(k, 'S/N', round(v['snr_cont'], 1), 'vcen', None if v['v_centroid_tau'] is None else round(v['v_centroid_tau'], 1),
          {wn: (v[wn]['logN'], v[wn]['EW_mA']) for wn in WIN})

# Figure
show = [k for k in ('SiII1260', 'CII1334', 'CII*1335', 'SiIII1206', 'SiIV1393', 'SiIV1402', 'NV1238') if k in prof]  # O I 1302 omitted: geocoronal airglow
fig, axs = new_sky_figure(11, 1.25 * len(show) + 1.4, 'Silk Nebula · COS UV absorption toward WD J0747+2106',
                          f"0.16° from the [O III] spine; Gaia DR3 d = {res["d_pc"]:.1f} ± {res["d_pc_err"]:.1f} pc; G130M"
                          if res['d_pc'] else 'G130M', nrows=len(show), ncols=1, sharex=True)
for ax, k in zip(axs, show):
    vv, nf, ne = prof[k]
    ax.step(vv, nf, where='mid', color='#e8eaed', lw=0.8)
    ax.fill_between(vv, nf - ne, nf + ne, step='mid', color='#9aa4ae', alpha=0.25, lw=0)
    ax.axhline(1, color=FGDIM, lw=0.5)
    ax.axvspan(23.8, 26.8, color=OIII_C, alpha=0.35, lw=0)
    ax.axvline(12.5, color=HA_C, lw=0.8, ls='-')
    ax.axvline(19.2, color=OK, lw=0.8, ls=':')
    ax.axvline(33.9, color=WARN, lw=0.8, ls=':')
    ax.set_ylim(-0.1, 1.4)
    ax.set_xlim(-200, 200)
    label(ax, -195, 1.12, k.replace('II*', 'II*'), size=8)
axs[-1].set_xlabel('V_helio (km/s)')
axs[0].plot([], [], color=OIII_C, lw=6, alpha=0.5, label='[O III] V (LAMOST)')
axs[0].plot([], [], color=HA_C, lw=0.8, label='red lines / WHAM Hα')
axs[0].plot([], [], color=OK, lw=0.8, ls=':', label='LIC')
axs[0].plot([], [], color=WARN, lw=0.8, ls=':', label='Gem cloud')
axs[0].legend(loc='lower right', ncol=4, fontsize=7)
finish(fig, FIG / 'spectra_ism_cos_wd.png')
