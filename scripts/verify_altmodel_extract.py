"""verify_altmodel: extract the needed CHIANTI v7.1.4 ion directories (safe 'data' filter)."""
import os, re, tarfile
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TGZ = os.path.join(ROOT, "data/raw/verify_altmodel/chianti/CHIANTI_v7.1.4_database.tar.gz")
OUT = os.path.join(ROOT, "data/raw/verify_altmodel/chianti/dbase")
os.makedirs(OUT, exist_ok=True)
ELS = {"h": 1, "he": 2, "c": 6, "n": 7, "o": 8, "ne": 10, "s": 16, "ar": 18}
pat = re.compile(r"^(%s)/(\1)_(\d+)/[^/]+$" % "|".join(ELS))
with tarfile.open(TGZ, "r:gz") as tf:
    mem = []
    for m in tf.getmembers():
        if not m.isfile():
            continue
        if pat.match(m.name) or m.name.startswith("ioneq/chianti") or m.name.startswith("abundance/sun_photospheric_2009_asplund") or m.name in ("VERSION", "README_CHIANTI"):
            mem.append(m)
    tf.extractall(OUT, members=mem, filter="data")
print(len(mem), "files extracted to", OUT)
