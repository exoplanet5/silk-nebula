"""verify_altmodel: fetch atomic data and published model tables (small downloads, logged).

Downloads (all accessed 2026-11-16):
  * CHIANTI v7.1.4 database tarball (Dere+1997; Landi+2013) -- 17 MB, energy levels (.elvlc),
    radiative data (.wgfa), Maxwellian-averaged collision strengths (.splups), ionization
    (.diparams), recombination (.rrparams, .drparams) for every ion.
  * Gnat & Sternberg 2007 (ApJS 168, 213) on-line tables 2-15: CIE and time-dependent
    (isochoric / isobaric) ion fractions and cooling efficiencies, from the authors' site.
Downloaded files are untrusted data: they are only parsed as text / tar members, never executed.
Run:  ~/.venvs/astro313/bin/python -I scripts/verify_altmodel_fetch.py
"""
import hashlib
import json
import os
import sys
import time
import urllib.request

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(ROOT, "data", "raw", "verify_altmodel")
os.makedirs(RAW, exist_ok=True)
LOG = os.path.join(RAW, "MANIFEST.jsonl")

GS = "http://wise-obs.tau.ac.il/~orlyg/cooling/"
URLS = [("https://download.chiantidatabase.org/CHIANTI_v7.1.4_database.tar.gz",
         "chianti/CHIANTI_v7.1.4_database.tar.gz")]
URLS += [(GS + "CIEion/tab2.txt", "gs07/tab2.txt")]
URLS += [(GS + f"TDion/tab{i}.txt", f"gs07/tab{i}.txt") for i in range(3, 13)]
URLS += [(GS + "CIEcool/tab13.txt", "gs07/tab13.txt"),
         (GS + "TDcool/tab14.txt", "gs07/tab14.txt"),
         (GS + "TDcool/tab15.txt", "gs07/tab15.txt")]
# Dere (2007, A&A 466, 771) Table 29: spline fits to scaled ionization rate coefficients (VizieR CfA mirror)
URLS += [("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2bA/466/771/table29&-out.all&-out.max=1000",
          "dere2007/table29.tsv")]
# literature (arXiv PDFs; read as text with pdftotext)
for aid, tag in [("1002.1309", "gnat_sternberg_mckee2010"), ("1006.3811", "kwak_shelton2010"),
                 ("1809.09101", "ji_oh_masterson2019"), ("2003.08390", "fielding2020"),
                 ("2211.01395", "chen_fielding_bryan2023")]:
    URLS.append((f"https://arxiv.org/pdf/{aid}", f"lit/{tag}_{aid}.pdf"))
URLS.append((GS + "GS07.pdf", "lit/gnat_sternberg2007_GS07.pdf"))


def fetch(url, rel):
    out = os.path.join(RAW, rel)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print("have", rel)
        return
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib (research script)"})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=300) as r:
        data = r.read()
    with open(out, "wb") as f:
        f.write(data)
    rec = dict(url=url, file=os.path.relpath(out, ROOT), bytes=len(data),
               sha256=hashlib.sha256(data).hexdigest(),
               accessed=time.strftime("%Y-%m-%dT%H:%M:%S%z"), seconds=round(time.time() - t0, 1))
    with open(LOG, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(rec)


if __name__ == "__main__":
    which = sys.argv[1:]
    for url, rel in URLS:
        if which and not any(w in rel for w in which):
            continue
        try:
            fetch(url, rel)
        except Exception as e:  # keep going, record failure
            print("FAILED", url, e)
            with open(LOG, "a") as f:
                f.write(json.dumps(dict(url=url, error=str(e),
                                        accessed=time.strftime("%Y-%m-%dT%H:%M:%S%z"))) + "\n")
