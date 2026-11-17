"""Is the Silk [O III] spine concentric with any candidate explosion/driving centre?  Where does it sit on the Monogem rim?

Inputs : data/derived/silk_trace.json (geometry teammate, v1.0: E-strand spine, 66 points, Dec 12-28.25)
         data/derived/context_pulsar_coincidence.json (MC birth loci at tau_c), literature centres typed below
Outputs: data/derived/context_concentricity.json  per-centre mean / rms / range of the spine's angular distance,
                                                   Monogem-rim geometry, Geminga energetics check
         data/derived/context_concentricity.csv   spine angular distance from each centre vs Dec
         figures/context_concentricity.png
Method : angular distance of every spine point to each centre (spherical); 'concentric' would mean rms of the
         distances comparable to the trace's intrinsic wiggle (single-circle fit rms 0.15 deg).
         Monogem rim: for a sphere of radius R at distance D (Knies+24: D = 300 pc, R = 63 pc = 11.97 deg), the
         line of sight at angular offset theta from the centre meets the shell surface at
         d = D cos(theta) -/+ sqrt(R^2 - D^2 sin^2(theta)); this is the distance the Silk would have IF it is part of the
         Monogem shell (near/far intersection; tangent point d = D cos(theta)).
Usage  : ~/.venvs/astro313/bin/python -I scripts/context_concentricity.py
"""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from context_common import *
from silk_style import *

D = ROOT / 'data/derived'
tr = load_trace()
sp = tr['spine']
SRA = np.array([q['ra'] for q in sp]); SDE = np.array([q['dec'] for q in sp])
cv = tr['curvature']
pc = json.loads((D / 'context_pulsar_coincidence.json').read_text())['pulsars']
bB = pc['B0656+14']['birth_at_tau']; bG = pc['Geminga']['birth_at_tau']

CENTRES = {   # name: (ra, dec, source)
    'Silk curvature centre (E strand)': (cv['center_ra'], cv['center_dec'], 'silk_trace.json v1.0'),
    'Monogem 1-SNR centre (Knies+24)': (107.08, 13.50, 'Knies+2024 Table 4'),
    'Monogem 2-SNR W centre (Knies+24)': (104.15, 15.30, 'Knies+2024 Table 4'),
    'G205.6+12.4 centre (Knies+24)': (110.33, 12.94, 'Knies+2024 Table 4'),
    'PSR B0656+14 now': (104.9506, 14.2393, 'ATNF v2.9.0'),
    'PSR B0656+14 at tau_c (MC median)': (bB['ra_p16_50_84'][1], bB['dec_p16_50_84'][1], 'context_pulsar_backtrace.py'),
    'Geminga now': (98.4756, 17.7703, 'ATNF v2.9.0'),
    'Geminga at tau_c (MC median)': (bG['ra_p16_50_84'][1], bG['dec_p16_50_84'][1], 'context_pulsar_backtrace.py'),
    'Geminga birthplace (Pellizza+ via NSNS)': (84.7088, 7.6303, 'simg.de/nebulae3/gem-01.html'),
    'Origem Loop centre': (94.1342, 16.0724, 'SIMBAD'),
    'GS202+17+019 HI supershell (Suad+14)': (113.6317, 16.7505, 'VizieR J/A+A/564/A116'),
}


def main():
    rows, summ = [], {}
    for nm, (ra, dec, src) in CENTRES.items():
        dd = angsep(SRA, SDE, ra, dec)
        summ[nm] = dict(ra=round(float(ra), 3), dec=round(float(dec), 3), source=src,
                        mean_deg=round(float(dd.mean()), 2), rms_deg=round(float(dd.std()), 3),
                        min_deg=round(float(dd.min()), 2), max_deg=round(float(dd.max()), 2),
                        rms_south_dec_le17p5=round(float(dd[SDE <= 17.5].std()), 3),
                        rms_north_dec_ge17p75=round(float(dd[SDE >= 17.75].std()), 3))
        for de, x in zip(SDE, dd):
            rows.append(dict(centre=nm, dec=de, sep_deg=x))
    pd.DataFrame(rows).to_csv(D / 'context_concentricity.csv', index=False, float_format='%.4f')

    # Monogem rim geometry (Knies+24 single SNR: D=300 pc, R=63 pc; also D = 288 (+33,-27) pc parallax scaling)
    th = angsep(SRA, SDE, 107.08, 13.50)
    rim = {}
    for Dm in (261.0, 288.0, 300.0, 321.0):
        R = Dm * np.tan(np.radians(11.97))          # keep the observed angular radius fixed
        disc = R ** 2 - (Dm * np.sin(np.radians(th))) ** 2
        near = np.where(disc >= 0, Dm * np.cos(np.radians(th)) - np.sqrt(np.clip(disc, 0, None)), np.nan)
        far = np.where(disc >= 0, Dm * np.cos(np.radians(th)) + np.sqrt(np.clip(disc, 0, None)), np.nan)
        tang = Dm * np.cos(np.radians(th))
        rim[f'D{Dm:.0f}'] = dict(R_pc=round(float(R), 1), frac_spine_inside_disc=round(float((disc >= 0).mean()), 2),
                                 tangent_dist_pc_range=[round(float(tang.min()), 0), round(float(tang.max()), 0)],
                                 near_side_pc_median=round(float(np.nanmedian(near)), 0) if np.isfinite(near).any() else None,
                                 far_side_pc_median=round(float(np.nanmedian(far)), 0) if np.isfinite(far).any() else None,
                                 near_side_of_shell_on_centre_los_pc=round(float(Dm - R), 0))
    inside = th < 11.97
    rim_summary = dict(spine_sep_from_monogem_centre_deg=dict(min=round(float(th.min()), 2), max=round(float(th.max()), 2),
                                                              median=round(float(np.median(th)), 2)),
                       dec_range_inside_11p97deg=[float(SDE[inside].min()), float(SDE[inside].max())] if inside.any() else None,
                       frac_inside=round(float(inside.mean()), 2), by_distance=rim)

    # Geminga: could a shell centred on Geminga's *present* position be its own?  (spine rms about Geminga-now)
    g = summ['Geminga now']
    mu = np.hypot(142.2, 107.4)                       # mas/yr (Faherty+07)
    t_rms = g['rms_deg'] * 3.6e6 / mu                 # yr for Geminga to move by the rms
    Edot = 3.25e34; tau = 3.42e5 * 3.156e7
    rho = 1.4 * 1.6726e-24 * 0.01                     # n_H = 0.01 cm^-3 hot/warm-ionised local medium (generous)
    # energy-driven wind bubble (Weaver+77): R = 0.76 (L t^3 / rho)^(1/5)
    R_bub = 0.76 * (Edot * tau ** 3 / rho) ** 0.2 / 3.0857e18
    R_need = 250.0 * np.sin(np.radians(g['mean_deg']))
    gem = dict(spine_mean_sep_deg=g['mean_deg'], spine_rms_deg=g['rms_deg'], pm_total_masyr=round(mu, 1),
               time_to_move_by_rms_yr=round(float(t_rms), 0),
               shell_radius_needed_pc_at_250pc=round(float(R_need), 1),
               spindown_energy_over_tau_erg=float(Edot * tau),
               wind_bubble_radius_pc_n0p01_over_tau=round(float(R_bub), 1),
               verdict='not physical: the present Geminga position moves out of concentricity in ~'
                       f'{t_rms/1e3:.0f} kyr, and its spin-down power (Edot*tau ~ {Edot*tau:.1e} erg) inflates only an'
                       f' R ~ {R_bub:.0f} pc bubble, vs the ~{R_need:.0f} pc needed')
    out = dict(spine_source=f"silk_trace.json v{tr.get('version')} E strand, {len(SRA)} points, Dec {SDE.min()}-{SDE.max()}",
               single_circle_fit_rms_deg=cv['fits']['E_full']['rms_resid_deg'], centres=summ,
               monogem_rim=rim_summary, geminga_present_concentricity=gem)
    (D / 'context_concentricity.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))

    # figure
    fig, axes = new_sky_figure(11, 5.6, 'Silk spine distance from candidate centres',
                               'angular separation of every E-strand spine point · Monogem radius 11.97° (Knies+24)', ncols=2,
                               gridspec_kw=dict(width_ratios=[1.35, 1]))
    ax = axes[0]
    show = [('Monogem 1-SNR centre (Knies+24)', '#e8eaed', '-'), ('Monogem 2-SNR W centre (Knies+24)', '#e8eaed', ':'),
            ('G205.6+12.4 centre (Knies+24)', OK, '-'), ('PSR B0656+14 now', SUN, '-'),
            ('PSR B0656+14 at tau_c (MC median)', SUN, '--'), ('Geminga now', '#ff7043', '-'),
            ('Silk curvature centre (E strand)', OIII_C, '-'), ('GS202+17+019 HI supershell (Suad+14)', '#c792ea', '-')]
    for nm, c, ls in show:
        dd = angsep(SRA, SDE, *CENTRES[nm][:2])
        ax.plot(SDE, dd, ls, color=c, lw=1.6, label=nm.replace(' (Knies+24)', '').replace(' (MC median)', '')
                .replace(' (E strand)', '').replace(' (Suad+14)', ''))
    ax.axhspan(11.97 * 0.9, 11.97 * 1.1, color='#e8eaed', alpha=0.08, lw=0)
    ax.axhline(11.97, color='#e8eaed', lw=0.8, ls='--', alpha=0.6)
    label(ax, 12.2, 10.3, 'Monogem radius ±10%', color=FGDIM, size=7.5, ha='left')
    ax.axvline(17.6, color=FGDIM, lw=0.8, ls=':')
    label(ax, 17.75, 2.0, 'kink', color=FGDIM, size=7.5)
    ax.set_xlabel('Dec along the spine (deg)'); ax.set_ylabel('angular separation (deg)')
    ax.set_ylim(0, 31); ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    ax = axes[1]
    names = list(summ)
    y = np.arange(len(names))
    ax.barh(y, [summ[n]['rms_deg'] for n in names], color=[OIII_C if 'curvature' in n else '#5c6b7a' for n in names])
    ax.axvline(cv['fits']['E_full']['rms_resid_deg'], color=OIII_C, lw=1, ls='--')
    for k, n in enumerate(names):
        label(ax, summ[n]['rms_deg'] + 0.05, k, f"{summ[n]['mean_deg']:.1f}°", color=FGDIM, size=7, va='center')
    ax.set_yticks(y); ax.set_yticklabels([n.replace(' (Knies+24)', '').replace(' (MC median)', '').replace(' via NSNS', '')
                                          for n in names], fontsize=7)
    ax.invert_yaxis(); ax.set_xscale('log'); ax.set_xlim(0.1, 12)
    ax.set_xlabel('rms of spine distance (deg; label = mean)'); ax.grid(alpha=0.25, axis='x')
    finish(fig, ROOT / 'figures/context_concentricity.png')


if __name__ == '__main__':
    main()
