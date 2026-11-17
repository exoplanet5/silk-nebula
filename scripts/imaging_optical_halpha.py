"""Narrowband / Halpha precovery tests for the Silk Nebula.

Data (all in data/raw/):
  lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits   Finkbeiner 2003 composite, HEALPix Nside 1024 NESTED, Galactic, R
  wham_ss_dr1/wham-ss-DR1-v161116-170912.fits        WHAM-SS DR1 kinematic table (1-deg beams, 12 km/s res.)
  skyview/SHASSA_{H,C,CC,Sm}/                         SHASSA (Gaustad+2001) via SkyView, dR, Dec < +15
  nsns_hips/{oiii,halpha,sii,sn_*}/Norder3            NSNS DR0.2 HiPS (Ziegenbalg), R
Outputs:
  data/derived/imaging_optical_wham_beams.csv, imaging_optical_halpha_measurements.json
  figures/imaging_optical_halpha_wide.png, imaging_optical_shassa_south.png, imaging_optical_wham_spectra.png
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import gaussian_filter, median_filter
from imaging_optical_common import *
from imaging_optical_nsns import ridge_spine
import silk_style as ss
import matplotlib.pyplot as plt

RIDGE = json.loads((DER / 'imaging_optical_nsns_ridge.json').read_text())


def ridge_overlay(ax, tr, lw=0.9, color=ss.OIII_C, label='[OIII] ridge'):
    sra, sdec = ridge_spine()
    ax.plot(sra, sdec, '-', color=color, lw=lw, alpha=0.85, transform=tr, label=label)


def teammate_trace_overlay(ax, tr):
    t = trace_info()
    if t is None:
        return None
    sra, sdec, sig, ver, prelim = spine()
    o = np.argsort(sdec)
    ax.plot(sra[o], sdec[o], ls='--', color=ss.SUN, lw=0.8, alpha=0.8, transform=tr,
            label=f'geometry trace v{ver}' + (' (prelim.)' if prelim else ''))
    return ver, prelim


def common_marks(ax, tr, p4=True, patch=True, decmax=99):
    for k, (ra, dec) in sg.REGIONS.items():
        if (k == 'P4' and not p4) or dec > decmax:
            continue
        c = ss.WARN if k == 'P4' else '#ffffff'
        ax.plot([ra], [dec], 'o', ms=5, mfc='none', mec=c, mew=1.0, transform=tr)
        ss.label(ax, ra + 0.35, dec, k, color=c, size=7.5, transform=tr, va='center', ha='right')


# ------------------------------------------------------------------ WHAM
def wham_beams():
    t = fits.getdata(RAW / 'wham_ss_dr1/wham-ss-DR1-v161116-170912.fits', 1)
    l = t['GAL-LON'].astype(float); b = t['GAL-LAT'].astype(float)
    c = SkyCoord(l * u.deg, b * u.deg, frame='galactic').icrs
    ra, dec = c.ra.deg, c.dec.deg
    m = (ra > 104) & (ra < 129) & (dec > 1) & (dec < 33)
    vel = t['VELOCITY'][m].astype(float); dat = t['DATA'][m].astype(float); var = t['VARIANCE'][m].astype(float)
    out = dict(l=l[m], b=b[m], ra=ra[m], dec=dec[m], inten=t['INTEN'][m].astype(float),
               err=t['ERROR'][m].astype(float), star=t['STAR'][m].astype(bool), vel=vel, data=dat, var=var)
    wins = {'v_m80_p80': (-80, 80), 'v_m10_p30': (-10, 30), 'v_m40_m10': (-40, -10), 'v_p30_p60': (30, 60),
            'v_p5_p25': (5, 25), 'v_m10_p10': (-10, 10)}
    for k, (v0, v1) in wins.items():
        I = np.zeros(m.sum()); E = np.zeros(m.sum())
        for i in range(m.sum()):
            v = vel[i]; d = dat[i]; s2 = var[i]
            ok = np.isfinite(v) & np.isfinite(d) & (v >= v0) & (v <= v1)
            if ok.sum() < 2:
                I[i] = np.nan; E[i] = np.nan; continue
            dv = np.median(np.diff(v[np.isfinite(v)]))
            I[i] = np.nansum(d[ok]) * dv; E[i] = np.sqrt(np.nansum(s2[ok])) * dv
        out[k] = I; out[k + '_err'] = E
    return out


def wham_analysis(W):
    sra, sdec = ridge_spine()
    across, along, beyond = along_across(W['ra'], W['dec'], sra, sdec)
    W['across'] = across; W['beyond'] = beyond
    res = {}
    for key in ('v_m80_p80', 'v_m10_p30', 'v_p5_p25', 'v_m10_p10', 'v_m40_m10', 'v_p30_p60'):
        rows = []
        for d0 in np.arange(11, 28, 2.0):
            sel = (np.abs(W['dec'] - (d0 + 1)) <= 1.0) & ~W['star'] & ~beyond
            on = sel & (np.abs(across) <= 0.5)
            off = sel & (np.abs(across) >= 1.0) & (np.abs(across) <= 3.0)
            if on.sum() == 0 or off.sum() < 3:
                continue
            I_on = np.nanmean(W[key][on]); I_off = np.nanmedian(W[key][off])
            sc = 1.4826 * np.nanmedian(np.abs(W[key][off] - I_off))
            rows.append(dict(dec_bin=[float(d0), float(d0 + 2)], n_on=int(on.sum()), n_off=int(off.sum()),
                             on_mean_R=float(I_on), off_median_R=float(I_off), excess_R=float(I_on - I_off),
                             off_scatter_R=float(sc)))
        res[key] = rows
    # north Halpha patch: beams whose centre lies within 0.6 deg of (115.8, 20.6)
    dpat = np.hypot((W['ra'] - 115.50) * np.cos(np.radians(20.55)), W['dec'] - 20.55)
    pa = dpat < 0.6
    ring = (dpat > 1.2) & (dpat < 3.0) & (np.abs(across) > 0.7)
    res['Ha_patch_north'] = {k: dict(n=int(pa.sum()), on_mean_R=float(np.nanmean(W[k][pa])) if pa.any() else None,
                                     ring_median_R=float(np.nanmedian(W[k][ring])),
                                     ring_scatter_R=float(1.4826 * np.nanmedian(np.abs(W[k][ring] - np.nanmedian(W[k][ring])))))
                             for k in ('v_m80_p80', 'v_m10_p30', 'v_p5_p25')}
    # mean spectra: on-ridge Dec 13-19.5 vs matched off-ridge
    sel = (W['dec'] > 13) & (W['dec'] < 19.5) & ~W['star'] & ~W['beyond']
    on = sel & (np.abs(across) <= 0.5); off = sel & (np.abs(across) >= 1.0) & (np.abs(across) <= 3.0)
    vgrid = np.arange(-60, 60.1, 2.0)
    def mean_spec(mask):
        S = [np.interp(vgrid, W['vel'][i][np.isfinite(W['vel'][i])], W['data'][i][np.isfinite(W['vel'][i])],
                       left=np.nan, right=np.nan) for i in np.where(mask)[0]]
        S = np.array(S)
        return np.nanmean(S, 0), np.nanstd(S, 0) / np.sqrt(np.sum(np.isfinite(S), 0)), len(S)
    son, eon, non = mean_spec(on); soff, eoff, noff = mean_spec(off)
    res['spectra'] = dict(v=vgrid.tolist(), on=son.tolist(), on_err=eon.tolist(), n_on=non,
                          off=soff.tolist(), off_err=eoff.tolist(), n_off=noff)
    # patch spectrum
    sp, ep, npch = mean_spec(pa & ~W['star']); res['spectra']['patch'] = sp.tolist(); res['spectra']['n_patch'] = npch
    return res


def save_wham_csv(W):
    keys = ['ra', 'dec', 'l', 'b', 'across', 'inten', 'err', 'v_m10_p30', 'v_m10_p30_err', 'v_p5_p25', 'v_m40_m10', 'v_p30_p60']
    with open(DER / 'imaging_optical_wham_beams.csv', 'w') as f:
        f.write('# WHAM-SS DR1 beams near the Silk Nebula; intensities in R integrated over LSR windows (km/s); '
                'across = signed distance (deg, + east) from the NSNS [OIII] ridge\n')
        f.write(','.join(keys + ['star']) + '\n')
        for i in range(len(W['ra'])):
            f.write(','.join(f'{W[k][i]:.4f}' for k in keys) + f',{int(W["star"][i])}\n')


# ------------------------------------------------------------------ SHASSA
def shassa_analysis():
    out = {}
    hd = fits.open(RAW / 'skyview/SHASSA_Sm/SHASSA_Sm_wide.fits')[0]
    w = WCS(hd.header); sm = hd.data.astype(float) / 10.0     # dR -> R
    cc = fits.getdata(RAW / 'skyview/SHASSA_CC/SHASSA_CC_wide.fits').astype(float) / 10.0
    yy, xx = np.mgrid[0:sm.shape[0], 0:sm.shape[1]]
    ra, dec = w.wcs_pix2world(xx, yy, 0)
    sra, sdec = ridge_spine()
    across, along, beyond = along_across(ra, dec, sra, sdec)
    # NSNS [OIII] on the same grid for comparison
    o3 = hips_mosaic('oiii', 3, w, sm.shape)
    edges = np.arange(-2.0, 2.0001, 3 / 60); cen = 0.5 * (edges[1:] + edges[:-1])
    for nm, (d0, d1) in {'S12.5_14': (12.5, 14.0), 'S14_15': (14.0, 15.0), 'P3_14.4_15.0': (14.42, 15.02)}.items():
        sel = (dec >= d0) & (dec < d1) & ~beyond & np.isfinite(sm) & (sm != 0)
        prof = np.array([np.nanmedian(sm[sel & (across >= edges[k]) & (across < edges[k + 1])]) for k in range(len(cen))])
        prof_o3 = np.array([np.nanmedian(o3[sel & (across >= edges[k]) & (across < edges[k + 1])]) for k in range(len(cen))])
        on = sel & (np.abs(across) <= 0.25); off = sel & (np.abs(across) >= 0.8) & (np.abs(across) <= 2.0)
        offE = sel & (across >= 0.8) & (across <= 2.0); offW = sel & (across <= -0.8) & (across >= -2.0)
        exc = np.nanmedian(sm[on]) - np.nanmedian(sm[off])
        sys_ = 0.5 * abs(np.nanmedian(sm[offE]) - np.nanmedian(sm[offW]))
        # scatter of 0.25-deg-wide strips in the off region as an empirical error for a 0.5-deg band
        strips = [np.nanmedian(sm[sel & (across >= a) & (across < a + 0.5)]) for a in np.arange(-2.0, -0.6, 0.25)] + \
                 [np.nanmedian(sm[sel & (across >= a) & (across < a + 0.5)]) for a in np.arange(0.6, 1.6, 0.25)]
        out[nm] = dict(dec_range=[d0, d1], Ha_on_median_R=float(np.nanmedian(sm[on])), Ha_off_median_R=float(np.nanmedian(sm[off])),
                       Ha_excess_R=float(exc), sys_EW_R=float(sys_), strip_scatter_R=float(np.nanstd(strips)),
                       OIII_nsns_excess_R=float(np.nanmedian(o3[on]) - np.nanmedian(o3[off])),
                       across_deg=cen.tolist(), Ha_profile_R=prof.tolist(), OIII_nsns_profile_R=prof_o3.tolist(),
                       coverage_fraction_on=float(np.mean(np.isfinite(sm[(dec >= d0) & (dec < d1) & (np.abs(across) <= 0.25)]))))
    return out, (w, sm, cc, o3, ra, dec)


# ------------------------------------------------------------------ figures
def fig_wide(W, maps_wide, w):
    fink, nsns_ha, nsns_o3 = maps_wide
    fig, axs = ss.new_sky_figure(16, 6.5, 'Silk Nebula · Hα and [OIII] wide field',
                                 'epochs: SHASSA 1997–2000 · WHAM-SS 1997–98 · NSNS DR0.2 2018–2025 ([OIII] 2022-12 on)',
                                 1, 4, subplot_kw={'projection': w})
    panels = [('Finkbeiner Hα [R]', fink, (0.5, 6.0), 'magma'),
              ('WHAM-SS Hα −10…+30 km/s [R]', None, (0.5, 4.0), 'magma'),
              ('NSNS Hα [R]', gaussian_filter(np.nan_to_num(nsns_ha), 1.0), (0.8, 6.0), 'magma'),
              ('NSNS [OIII] [R]', gaussian_filter(np.nan_to_num(nsns_o3), 1.0), (-0.4, 3.0), 'bone')]
    for ax, (ttl, img, (v0, v1), cm) in zip(axs, panels):
        tr = ax.get_transform('world')
        if img is not None:
            im = ax.imshow(img, origin='lower', vmin=v0, vmax=v1, cmap=cm, interpolation='nearest')
        else:
            ax.set_facecolor('#05070a')
            sel = ~W['star']
            im = ax.scatter(W['ra'][sel], W['dec'][sel], c=W['v_m10_p30'][sel], s=95, marker='h', vmin=v0, vmax=v1,
                            cmap=cm, transform=tr, edgecolors='none')
            ax.set_xlim(-0.5, img_shape[1] - 0.5); ax.set_ylim(-0.5, img_shape[0] - 0.5); ax.set_aspect('equal')
        if '[OIII]' in ttl:
            overlay_v1(ax, tr, spine=False, wstrand=False, patch=True)
        else:
            overlay_v1(ax, tr)
        common_marks(ax, tr)
        cb = plt.colorbar(im, ax=ax, orientation='horizontal', fraction=0.04, pad=0.07)
        cb.ax.tick_params(labelsize=7, colors=ss.FGDIM); cb.outline.set_edgecolor('#2a323c')
        ax.set_title(ttl, color=ss.FG, fontsize=9)
        ax.coords[0].set_axislabel('RA'); ax.coords[1].set_axislabel('Dec')
        ax.coords[0].set_ticks(spacing=3 * u.deg); ax.coords[0].set_major_formatter('d')
        ax.coords[1].set_major_formatter('d'); ax.coords[1].set_ticks(spacing=4 * u.deg)
        ax.coords.grid(color='#2a323c', alpha=0.6, ls=':', lw=0.5)
        for c in ax.coords:
            c.set_ticklabel(size=7, color=ss.FGDIM)
    h, l = axs[0].get_legend_handles_labels()
    axs[0].legend(h, l, loc='lower left', fontsize=6.5, ncol=1, framealpha=0.7)
    # annotate source boundaries on the Finkbeiner panel
    ss.label(axs[0], 110.3, 7.0, 'SHASSA', size=7, color=ss.FGDIM, transform=axs[0].get_transform('world'))
    ss.label(axs[0], 117.0, 26.0, 'WHAM 1°', size=7, color=ss.FGDIM, transform=axs[0].get_transform('world'))
    ss.label(axs[0], 123.6, 16.5, 'VTSS', size=7, color=ss.FGDIM, transform=axs[0].get_transform('world'))
    ss.finish(fig, FIG / 'imaging_optical_halpha_wide.png')


def fig_shassa(sh, shassa_meas):
    w, sm, cc, o3, ra, dec = sh
    fig = plt.figure(figsize=(13, 6.2))
    fig.text(0.012, 0.985, 'Silk Nebula · SHASSA Hα, south end', color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
    fig.text(0.012, 0.945, 'SHASSA 1997–2000, continuum-subtracted, 4′ smoothed · NSNS [OIII] contours 1, 2 R',
             color=ss.FGDIM, fontsize=9.5, va='top')
    ax = fig.add_subplot(1, 2, 1, projection=w)
    tr = ax.get_transform('world')
    img = np.where(sm == 0, np.nan, sm)
    im = ax.imshow(img, origin='lower', vmin=0.5, vmax=5.0, cmap='magma')
    o3s = gaussian_filter(np.nan_to_num(o3), 4.0)
    ax.contour(o3s, levels=[1.0, 2.0], colors=[ss.OIII_C, '#b3e5fc'], linewidths=[0.8, 1.0], alpha=0.9)
    common_marks(ax, tr, patch=False, decmax=15.6)
    ax.set_xlim(150, 830); ax.set_ylim(60, 825)
    ax.coords[0].set_axislabel('RA'); ax.coords[1].set_axislabel('Dec')
    ax.coords[0].set_major_formatter('d'); ax.coords[1].set_major_formatter('d')
    ax.coords[0].set_ticks(spacing=2 * u.deg); ax.coords[1].set_ticks(spacing=2 * u.deg)
    ax.coords.grid(color='#2a323c', alpha=0.6, ls=':', lw=0.5)
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label('Hα [R]', color=ss.FGDIM)
    cb.ax.tick_params(labelsize=7, colors=ss.FGDIM)
    ax2 = fig.add_subplot(1, 2, 2)
    for nm, col in [('S12.5_14', ss.HA_C), ('S14_15', ss.WARN)]:
        r = shassa_meas[nm]
        x = np.array(r['across_deg'])
        ax2.plot(x, np.array(r['Ha_profile_R']), '-', color=col, lw=1.4, label=f"Hα Dec {r['dec_range'][0]}–{r['dec_range'][1]}")
        ax2.plot(x, np.array(r['OIII_nsns_profile_R']) + 1.0, '--', color=col, lw=1.0, alpha=0.8,
                 label=f"[OIII] NSNS +1 R")
    ax2.axvspan(-0.25, 0.25, color=ss.OIII_C, alpha=0.08)
    ax2.set_xlabel('distance from [OIII] ridge [deg]  (+ east)'); ax2.set_ylabel('intensity [R]')
    ax2.grid(alpha=0.3); ax2.legend(fontsize=7.5, ncol=2, loc='upper left')
    ax2.set_xlim(-2, 2)
    fig.tight_layout(rect=(0, 0.02, 1, 0.9))
    ss.finish(fig, FIG / 'imaging_optical_shassa_south.png', tight=False)


def fig_wham_spectra(res):
    s = res['spectra']
    fig, ax = ss.new_sky_figure(8.5, 5.0, 'Silk Nebula · WHAM-SS Hα spectra', 'mean of 1° beams, Dec 13–19.5, on vs off ridge')
    fig.texts[1].set_y(0.938)
    v = np.array(s['v'])
    ax.fill_between(v, np.array(s['on']) - np.array(s['on_err']), np.array(s['on']) + np.array(s['on_err']), color=ss.HA_C, alpha=0.25)
    ax.plot(v, s['on'], color=ss.HA_C, lw=1.5, label=f"on ridge |Δ|≤0.5° (n={s['n_on']})")
    ax.plot(v, s['off'], color=ss.FGDIM, lw=1.2, label=f"off ridge 1–3° (n={s['n_off']})")
    ax.plot(v, s['patch'], color=ss.SII_C, lw=1.0, ls='--', label=f"north Hα patch (n={s['n_patch']})")
    ax.plot(v, np.array(s['on']) - np.array(s['off']), color='#ffffff', lw=1.0, label='on − off')
    ax.axvspan(-1.0, 6.1, color=ss.HA_C, alpha=0.08); ax.axvspan(11.2, 15.5, color=ss.OIII_C, alpha=0.12)
    ss.label(ax, 2.5, -0.008, 'LAMOST red', size=7, color=ss.HA_C, ha='center')
    ss.label(ax, 13.3, -0.008, '[OIII]', size=7, color=ss.OIII_C, ha='center')
    ax.set_ylim(-0.012, None)
    ax.axhline(0, color='#2a323c', lw=0.8)
    ax.set_xlabel('V_LSR [km/s]'); ax.set_ylabel('Hα [R / (km/s)]'); ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc='upper left')
    ss.finish(fig, FIG / 'imaging_optical_wham_spectra.png')


if __name__ == '__main__':
    W = wham_beams()
    wres = wham_analysis(W)
    save_wham_csv(W)
    w, shp = wide_grid(2.0)
    img_shape = shp
    fm = hp.read_map(str(RAW / 'lambda_finkbeiner_halpha/Halpha_fwhm06_1024.fits'))   # converted to RING
    fink = healpix_sample(fm, w, shp, 'G', nest=False)
    nsns_ha = hips_mosaic('halpha', 3, w, shp); nsns_o3 = hips_mosaic('oiii', 3, w, shp)
    # Finkbeiner on/off ridge
    ra, dec = pixel_radec(w, shp)
    sra, sdec = ridge_spine()
    across, along, beyond = along_across(ra, dec, sra, sdec)
    fres = {}
    for d0 in np.arange(11, 28, 2.0):
        sel = (np.abs(dec - (d0 + 1)) <= 1) & ~beyond
        on = sel & (np.abs(across) <= 0.25); off = sel & (np.abs(across) >= 0.8) & (np.abs(across) <= 2.0)
        offE = sel & (across >= 0.8) & (across <= 2.0); offW = sel & (across <= -0.8) & (across >= -2.0)
        fres[f'{d0:.0f}-{d0 + 2:.0f}'] = dict(
            fink_excess_R=float(np.nanmedian(fink[on]) - np.nanmedian(fink[off])), fink_off_R=float(np.nanmedian(fink[off])),
            fink_EW_asym_R=float(np.nanmedian(fink[offE]) - np.nanmedian(fink[offW])),
            nsns_ha_excess_R=float(np.nanmedian(nsns_ha[on]) - np.nanmedian(nsns_ha[off])),
            nsns_o3_excess_R=float(np.nanmedian(nsns_o3[on]) - np.nanmedian(nsns_o3[off])))
    shres, sh = shassa_analysis()
    out = dict(wham=wres, finkbeiner_vs_ridge=fres, shassa=shres,
               notes='on = |across|<=0.25 deg (maps) or <=0.5 deg (WHAM beam centres); off = 0.8-2 deg (maps) or 1-3 deg (WHAM).')
    save_json(out, DER / 'imaging_optical_halpha_measurements.json')
    for k, v in fres.items():
        print('fink', k, {kk: round(vv, 2) for kk, vv in v.items()})
    for k in ('v_m80_p80', 'v_m10_p30', 'v_p5_p25'):
        for r in wres[k]:
            print('wham', k, r['dec_bin'], r['n_on'], round(r['excess_R'], 2), '+-', round(r['off_scatter_R'], 2))
    print('patch', wres['Ha_patch_north'])
    for k, v in shres.items():
        print('shassa', k, {kk: (round(vv, 2) if isinstance(vv, float) else None) for kk, vv in v.items() if not isinstance(vv, list)})
    fig_wide(W, (fink, nsns_ha, nsns_o3), w)
    fig_shassa(sh, shres)
    fig_wham_spectra(wres)
