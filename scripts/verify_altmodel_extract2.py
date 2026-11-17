"""verify_altmodel: extract further CHIANTI v7.1.4 ion directories (Si, Mg) with the safe 'data' filter.
Run: ~/.venvs/astro313/bin/python -I scripts/verify_altmodel_extract2.py"""
import os, re, tarfile
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TGZ = os.path.join(ROOT, "data/raw/verify_altmodel/chianti/CHIANTI_v7.1.4_database.tar.gz")
OUT = os.path.join(ROOT, "data/raw/verify_altmodel/chianti/dbase")
pat = re.compile(r"^(si|mg)/(\1)_(\d+)/[^/]+$")
with tarfile.open(TGZ, "r:gz") as tf:
    mem = [m for m in tf.getmembers() if m.isfile() and pat.match(m.name)]
    tf.extractall(OUT, members=mem, filter="data")
print(len(mem), "files extracted to", OUT)
