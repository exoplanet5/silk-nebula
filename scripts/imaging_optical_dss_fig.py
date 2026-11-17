"""Figure + summary of the DSS photographic precovery search (reads imaging_optical_dss.py outputs).

Images: per 2x2 deg tile the star-masked 1' block medians are turned into a density contrast (D - D_bg)/D_sky, with
D_bg a 2-D quadratic fitted to |across| >= 0.7 deg from the NSNS ridge, then smoothed (Gaussian sigma 2') and
mosaicked onto one TAN grid. Under a linear density-intensity response near the sky level this contrast equals the
fractional sky-intensity excess (systematic factor ~2, see notes).
Limits: template amplitudes from imaging_optical_dss_profiles.json, rescaled so that the measured plate sky matches the
literature sky (B_J 22.3, R_F 21.0, POSS-I O 22.4, E 20.9 mag/arcsec^2), expressed as 'true line R per NSNS-[OIII] R'.
Outputs: figures/imaging_optical_dss_strips.png, data/derived/imaging_optical_dss_summary.json
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import astropy.units as u
from reproject import reproject_interp
from scipy.ndimage import gaussian_filter
from imaging_optical_common import DER, FIG, save_json, along_across, tan_wcs
from imaging_optical_nsns import ridge_spine
import silk_style as ss
import silk_geometry as sg

BANDS = [('DSS2_Blue', 'DSS2 J  IIIaJ', '1987–94'), ('DSS2_Red', 'DSS2 F  IIIaF', '1987–96'),
         ('DSS1_Blue', 'DSS1 O  103aO', '1950–55'), ('DSS1_Red', 'DSS1 E  103aE', '1950–55')]
SKY_LIT = {'DSS2_Blue': 22.3, 'DSS2_Red': 21.0, 'DSS1_Blue': 22.4, 'DSS1_Red': 20.9}


def m1R(band):
    e = json.loads((DER / 'imaging_optical_expected_signal.json').read_text())['cases']
    if band == 'DSS2_Blue':
        return e['OIII_12R_paper']['POSSII_J']['mag_arcsec2'] + 2.5 * np.log10(12), '[OIII]'
    if band == 'DSS1_Blue':
        return e['OIII_12R_paper']['POSSI_O']['mag_arcsec2'] + 2.5 * np.log10(12), '[OIII]4959 only'
    if band == 'DSS2_Red':
        return e['Halpha_10R']['POSSII_F']['mag_arcsec2'] + 2.5, 'Halpha'
    return e['Halpha_10R']['POSSI_E']['mag_arcsec2'] + 2.5, 'Halpha'


def load_rows(band):
    """Binned tiles grouped by row: list of (row, [(img, wcs, skyD, k), ...])."""
    z = np.load(DER / f'imaging_optical_dss_binned_{band}.npz')
    keys = sorted({k.rsplit('_', 1)[0] for k in z.files})
    rows = {}
    for k in keys:
        r = k.split('_')[0]
        rows.setdefault(r, []).append((z[k + '_img'], WCS(str(z[k + '_wcs'])), float(z[k + '_sky'][0]), float(z[k + '_sky'][1])))
    return rows


def row_fields(tiles, sra, sdec):
    """Concatenate the E+W tiles of a row with ONE zero point (the tiles are cut from the same plate scans, so their
    densities are continuous across the seam, which lies near the ridge). Returns ra, dec, across, along, contrast
    (D - D_bg)/D_sky with D_bg a quadratic in (across, along) fitted to |across| >= 0.7 deg, and the flux scale."""
    ra_l, de_l, v_l, idx = [], [], [], []
    for j, (img, tw, skyD, k) in enumerate(tiles):
        yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
        ra, dec = tw.wcs_pix2world(xx, yy, 0)
        ra_l.append(ra.ravel()); de_l.append(dec.ravel()); v_l.append(img.ravel()); idx.append(np.full(img.size, j))
    ra = np.concatenate(ra_l); dec = np.concatenate(de_l); v = np.concatenate(v_l); idx = np.concatenate(idx)
    acr, alo, bey = along_across(ra, dec, sra, sdec)
    ok = np.isfinite(v)
    skyD = float(np.nanmedian(v)); kmean = float(np.nanmean([t[3] for t in tiles]))
    bgm = ok & (np.abs(acr) >= 0.7) & (np.abs(acr) <= 2.2)
    a0 = np.nanmedian(alo[ok])
    X = lambda x, y: np.c_[np.ones(x.size), x, y, x * x, x * y, y * y]
    c, *_ = np.linalg.lstsq(X(acr[bgm], alo[bgm] - a0), v[bgm], rcond=None)
    con = (v - X(acr, alo - a0) @ c) / skyD
    return dict(ra=ra, dec=dec, across=acr, along=alo, contrast=con, ok=ok, idx=idx, skyD=skyD, k=kmean,
                sky_mag=float(-2.5 * np.log10(skyD / (kmean * 36.0))))


def contrast_mosaic(band, w, shp, sra, sdec):
    acc = np.zeros(shp); cnt = np.zeros(shp)
    for r, tiles in load_rows(band).items():
        R = row_fields(tiles, sra, sdec)
        n0 = 0
        for j, (img, tw, skyD, k) in enumerate(tiles):
            con = R['contrast'][R['idx'] == j].reshape(img.shape)
            ok = np.isfinite(con)
            f = np.where(ok, con, 0.0); wgt = gaussian_filter(ok.astype(float), 2.0)
            sm = gaussian_filter(f, 2.0) / np.maximum(wgt, 1e-3)
            sm[wgt < 0.3] = np.nan
            arr, fp = reproject_interp((sm, tw), w, shape_out=shp, order='bilinear')
            g = (fp > 0) & np.isfinite(arr)
            acc[g] += arr[g]; cnt[g] += 1
    return np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)


def row_profile2(R, o3s, nw):
    """Cross-filament profile (2' bins, % of sky) of one row + NSNS [OIII] template fit with placebo errors."""
    x = R['across']; v = 100 * R['contrast']; ok = R['ok'] & np.isfinite(v) & (np.abs(x) <= 2.2)
    px, py = nw.wcs_world2pix(R['ra'], R['dec'], 0)
    ix = np.clip(np.round(px).astype(int), 0, o3s.shape[1] - 1); iy = np.clip(np.round(py).astype(int), 0, o3s.shape[0] - 1)
    t = o3s[iy, ix]
    edges = np.arange(-2.2, 2.2001, 2 / 60); cen = 0.5 * (edges[1:] + edges[:-1])
    prof = np.full(len(cen), np.nan); pt = np.full(len(cen), np.nan)
    for i in range(len(cen)):
        s = ok & (x >= edges[i]) & (x < edges[i + 1])
        if s.sum() >= 5:
            prof[i] = np.median(v[s]); pt[i] = np.median(t[s])
    g = np.isfinite(prof) & np.isfinite(pt)
    bgm = g & (np.abs(cen) >= 0.7)
    t3 = pt - np.polyval(np.polyfit(cen[bgm], pt[bgm], 2), cen)
    def fitA(tt):
        X = np.c_[tt[g], np.ones(g.sum()), cen[g], cen[g] ** 2]
        cf, *_ = np.linalg.lstsq(X, prof[g], rcond=None)
        return cf
    cf = fitA(t3)
    t3f = np.where(np.isfinite(t3), t3, 0.0)
    pl = [fitA(np.interp(cen - sh, cen, t3f, left=0, right=0))[0] for sh in np.arange(-1.6, 1.61, 0.1) if abs(sh) >= 0.8]
    bg = cf[1] + cf[2] * cen + cf[3] * cen ** 2
    return dict(across_deg=cen.tolist(), profile_pct=(prof - bg).tolist(), template_R=t3.tolist(),
                A_pct_per_R=float(cf[0]), A_err_placebo=float(np.std(pl)), n_placebo=len(pl),
                bin_scatter_pct=float(1.4826 * np.nanmedian(np.abs((prof - bg - cf[0] * t3)[g & (np.abs(cen) > 0.7)]))),
                sky_mag=R['sky_mag'])


def summary(sra, sdec):
    z = np.load(DER / 'imaging_optical_nsns_maps.npz')
    nw = WCS(str(z['wcs_header'])); o3s = gaussian_filter(np.nan_to_num(z['oiii']), 1.0)
    plates = json.loads((DER / 'imaging_optical_dss_profiles.json').read_text())
    out = {}
    for band, *_ in BANDS:
        mm, line = m1R(band)
        pct_1R = 100 * 10 ** (-0.4 * (mm - SKY_LIT[band]))     # % of the literature sky per R of the line
        rows = {}; A = []; E = []
        for r, tiles in load_rows(band).items():
            R = row_fields(tiles, sra, sdec)
            P = row_profile2(R, o3s, nw)
            P['A_trueR_per_NSNS_R'] = P['A_pct_per_R'] / pct_1R
            P['err'] = P['A_err_placebo'] / pct_1R
            P['plates'] = sorted({p for t in plates[band]['rows'].get(r, {}).get('tiles', []) for p in t['plates']})
            rows[r] = P; A.append(P['A_trueR_per_NSNS_R']); E.append(P['err'])
        A = np.array(A); E = np.array(E); keys = list(rows)
        def comb(sel):
            w = 1 / E[sel] ** 2; m = np.sum(w * A[sel]) / np.sum(w); e = 1 / np.sqrt(np.sum(w))
            chi = np.sum(w * (A[sel] - m) ** 2) / max(sel.sum() - 1, 1)
            return dict(A=float(m), err=float(e * max(1, np.sqrt(chi))), chi2_dof=float(chi), n_rows=int(sel.sum()))
        core = np.array([k in ('d13', 'd15', 'd17', 'd19') for k in keys])
        out[band] = dict(plate=plates[band]['plate'], line_for_conversion=line, mag_arcsec2_of_1R=mm,
                         pct_of_sky_per_R=pct_1R, rows=rows,
                         combined_all=comb(np.ones(len(A), bool)), combined_d13_d19=comb(core))
        c = out[band]['combined_d13_d19']
        out[band]['upper3sigma_trueR_per_NSNS_R_d13_d19'] = float(max(c['A'], 0) + 3 * c['err'])
    out['_notes'] = ('Row-level processing (E+W tiles with one zero point; quadratic background in (across, along) '
                     'fitted to 0.7-2.2 deg from the ridge); contrast (D-D_bg)/D_sky taken as fractional sky intensity '
                     '(linear density response assumed; systematic factor ~2). A in true line R per R of the NSNS [OIII] '
                     'template (line_for_conversion); red plates are a control ([OIII] out of band). Supersedes the '
                     'per-tile zero points of imaging_optical_dss_profiles.json, which created a step at the tile seam '
                     'near the ridge. SDSS sky-level analysis gives 1.35 +- 0.19 (imaging_optical_sdss_sky.json).')
    save_json(out, DER / 'imaging_optical_dss_summary.json')
    return out


def main():
    sra, sdec = ridge_spine()
    S = summary(sra, sdec)
    for b, v in S.items():
        if not b.startswith('_'):
            print(b, v['combined_d13_d19'], v['upper3sigma_trueR_per_NSNS_R_d13_d19'])
    w, shp = tan_wcs(116.0, 18.8, 5.0, 18.4, 1.0)
    z = np.load(DER / 'imaging_optical_nsns_maps.npz')
    nw = WCS(str(z['wcs_header']))
    o3, _ = reproject_interp((gaussian_filter(np.nan_to_num(z['oiii']), 1.5), nw), w, shape_out=shp)
    fig = plt.figure(figsize=(15.5, 9.6))
    fig.text(0.012, 0.985, 'Silk Nebula · DSS photographic plates', color='#cfd6dd', fontsize=12, fontweight='bold', va='top')
    fig.text(0.012, 0.955, 'star-masked 1′ medians, quadratic plate background removed, 2′ smoothing · contrast (D−D_bg)/D_sky',
             color=ss.FGDIM, fontsize=9.5, va='top')
    W0, dx = 0.03, 0.128
    panels = [('NSNS [OIII] 2022–25', o3, (-0.3, 3.0), 'bone', 'R')]
    for band, name, ep in BANDS:
        panels.append((f'{name} {ep}', 100 * contrast_mosaic(band, w, shp, sra, sdec), (-1.5, 1.5), 'RdBu_r', '% of sky'))
    for i, (ttl, img, (v0, v1), cm, unit) in enumerate(panels):
        ax = fig.add_axes([W0 + i * (dx + 0.012), 0.08, dx, 0.82], projection=w)
        tr = ax.get_transform('world')
        im = ax.imshow(img, origin='lower', vmin=v0, vmax=v1, cmap=cm, interpolation='nearest')
        if i > 0:
            ax.contour(gaussian_filter(np.nan_to_num(o3), 3), levels=[1.0], colors=[ss.SUN], linewidths=0.6, alpha=0.8)
        for k, (ra0, de0) in sg.REGIONS.items():
            c = ss.WARN if k == 'P4' else ('#ffffff' if i == 0 else '#111111')
            ax.plot([ra0], [de0], 'o', ms=4, mfc='none', mec=c, mew=0.9, transform=tr)
            if i == 0:
                ss.label(ax, ra0 + 0.75, de0, k, color=c, size=7, transform=tr, va='center', ha='right')
        ax.set_title(ttl, color=ss.FG, fontsize=8.5)
        ax.coords[0].set_major_formatter('d'); ax.coords[1].set_major_formatter('d')
        ax.coords[0].set_ticks(spacing=2 * u.deg); ax.coords[1].set_ticks(spacing=2 * u.deg)
        ax.coords[0].set_axislabel('RA')
        ax.coords[1].set_axislabel('Dec' if i == 0 else '')
        if i > 0:
            ax.coords[1].set_ticklabel_visible(False)
        for c in ax.coords:
            c.set_ticklabel(size=7, color=ss.FGDIM)
        cax = fig.add_axes([W0 + i * (dx + 0.012) + 0.01, 0.035, dx - 0.02, 0.012])
        cb = plt.colorbar(im, cax=cax, orientation='horizontal'); cb.set_label(unit, color=ss.FGDIM, fontsize=7)
        cb.ax.tick_params(labelsize=6.5, colors=ss.FGDIM); cb.outline.set_edgecolor('#2a323c')
    # profile panel
    ax = fig.add_axes([0.70 + 0.07, 0.56, 0.21, 0.34])
    for band, col, lab in [('DSS2_Blue', ss.OIII_C, 'DSS2 J'), ('DSS2_Red', ss.HA_C, 'DSS2 F')]:
        P = [np.array(S[band]['rows'][r]['profile_pct']) for r in ('d13', 'd15', 'd17', 'd19')]
        x = np.array(S[band]['rows']['d15']['across_deg'])
        ax.plot(x, np.nanmean(P, 0), '-', color=col, lw=1.2, label=lab)
    t = np.nanmean([np.array(S['DSS2_Blue']['rows'][r]['template_R']) for r in ('d13', 'd15', 'd17', 'd19')], 0)
    f1 = S['DSS2_Blue']['pct_of_sky_per_R']
    ax.plot(x, 1.35 * f1 * t, '--', color='#ffffff', lw=1.0, label='SDSS-scaled [OIII]')
    ax.plot(x, 2.7 * f1 * t, ':', color=ss.DANGER, lw=1.2, label='paper 10–13 R (min.)')
    ax.axhline(0, color='#2a323c', lw=0.8)
    ax.set_xlim(-1.9, 1.9); ax.set_ylim(-2.5, 2.5); ax.set_xlabel('distance from ridge [deg]', fontsize=8); ax.set_ylabel('% of sky', fontsize=8)
    ax.set_title('mean profile, Dec 12–20', color=ss.FG, fontsize=8.5)
    ax.grid(alpha=0.3); ax.legend(fontsize=6.5, loc='lower left')
    ax.tick_params(labelsize=7)
    # limits text panel
    ax2 = fig.add_axes([0.70 + 0.07, 0.08, 0.21, 0.38]); ax2.axis('off')
    lines = ['3σ limits, Dec 12–20', '(× NSNS [OIII] intensity)', '']
    for band, name, ep in BANDS:
        v = S[band]; c = v['combined_d13_d19']
        lines.append(f"{name.split()[0]} {name.split()[1]}: {c['A']:+.1f}±{c['err']:.1f}  <{v['upper3sigma_trueR_per_NSNS_R_d13_d19']:.0f}")
    lines += ['', 'SDSS sky levels: 1.35±0.19', 'paper 10–13 R: 2.7–6.2']
    for j, s in enumerate(lines):
        ax2.text(0.0, 1 - j * 0.075, s, color=ss.FG if j else '#cfd6dd', fontsize=8, transform=ax2.transAxes, va='top')
    ss.credit(fig)
    fig.savefig(FIG / 'imaging_optical_dss_strips.png', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    main()
