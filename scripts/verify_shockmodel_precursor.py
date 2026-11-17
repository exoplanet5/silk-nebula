"""verify_shockmodel_precursor: order-of-magnitude precursor budget for the consistent truncated shocks.
Energy radiated per H nucleus while T > 5e4 K (from the 1-D model's energy equation), half assumed to go upstream,
a fraction f_EUV = 0.5-0.8 of that in H-ionizing photons of mean energy 20-40 eV -> upstream ionizing photons per
incoming H (phi).  Precursor length, flow, recombination and cooling times at n0 = 0.025 cm^-3.
Output: data/derived/verify_shockmodel_precursor.json"""
import sys, json, numpy as np
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
import verify_shockmodel_flow as F
import verify_shockmodel_fit as V
EV = 1.602e-12
net = F.Network('asplund09')
out = {}
for v in (100, 120, 140):
    sh = F.Shock(net, v, 0.1, B0_uG=10 * np.sqrt(0.1), pre=V.preshock(0.05), complete=False).run()
    E = np.array([sh.Efun(x) for x in sh.x])
    dE = (sh.E0 - E) * net.mass_per_H / EV                    # eV radiated per H nucleus so far
    i5 = int(np.argmax(sh.T < 5e4)); i2 = int(np.argmax(sh.T < 2e4))
    tot = 0.5 * net.mass_per_H * (v * 1e5)**2 / EV
    phi = [0.5 * f * dE[i5] / hv for f, hv in ((0.5, 40.0), (0.8, 20.0))]
    n0 = 0.025
    sig = (6.3e-18 * (13.6 / 30.0)**3, 6.3e-18)               # cross-section at 30 eV / at threshold
    lens = {f'X{X}': [1 / (n0 * X * s) / 3.0857e18 for s in sig] for X in (0.05, 0.3)}
    t_flow = {k: [L * 3.0857e18 / (v * 1e5) / 3.156e7 for L in Ls] for k, Ls in lens.items()}
    t_rec = 1 / (2.6e-13 * n0) / 3.156e7
    t_cool = 3.45 * 1.38e-16 * 1e4 / (n0 * 1e-24) / 3.156e7   # Lambda(1e4 K, photoionised) ~ 1e-24 erg cm3/s
    out[v] = dict(eV_per_H_total=float(tot), eV_radiated_T_gt_5e4=float(dE[i5]), eV_radiated_T_gt_2e4=float(dE[i2]),
                  phi_up_range=[float(p) for p in phi], precursor_length_pc=lens, t_flow_yr=t_flow,
                  t_rec_yr=float(t_rec), t_cool_yr=float(t_cool))
    print(v, json.dumps(out[v]))
json.dump(out, open('/Users/mickey/astronomy/silk-nebula/data/derived/verify_shockmodel_precursor.json', 'w'), indent=1)
