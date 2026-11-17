#!/usr/bin/env python
"""spectra_ism: Gaia DR3 GSP-Spec DIB 862 nm in the box RA 104-128, Dec +4..+34.

ESA Gaia TAP (gea.esac.esa.int) failed with an SSL error through the proxy on 2026-11-16, so the NOIRLab Data Lab
TAP (https://datalab.noirlab.edu/tap/sync) table gaia_dr3.astrophysical_parameters is tried. Output:
data/raw/spectra_ism_gaia/gaia_dr3_dib_box_datalab.csv
"""
import sys
sys.path.insert(0, '/Users/mickey/astronomy/silk-nebula/scripts')
from spectra_ism_common import fetch, RAW  # noqa: E402

q = """SELECT a.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.phot_g_mean_mag,
a.dibew_gspspec, a.dibewuncertainty_gspspec, a.dib_gspspec_lambda, a.dib_gspspec_lambda_uncertainty,
a.dibp0_gspspec, a.dibp2_gspspec, a.dibqf_gspspec, a.distance_gspphot, a.teff_gspspec, a.ag_gspphot
FROM gaia_dr3.astrophysical_parameters AS a JOIN gaia_dr3.gaia_source AS g ON a.source_id = g.source_id
WHERE a.dibew_gspspec IS NOT NULL AND g.ra BETWEEN 104 AND 128 AND g.dec BETWEEN 4 AND 34"""
dest = RAW / 'spectra_ism_gaia' / 'gaia_dr3_dib_box_datalab.csv'
ok, n, st = fetch('https://datalab.noirlab.edu/tap/sync', dest,
                  params={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': q}, timeout=600, retries=1)
print(ok, n, st)
if ok:
    t = dest.read_text(errors='replace')
    print(t[:500])
    print('lines', len(t.splitlines()))
