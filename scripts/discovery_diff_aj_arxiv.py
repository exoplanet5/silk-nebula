"""Sentence-level diff between the AJ published text (iopscience capture) and arXiv v1 body text.

Usage: python -I discovery_diff_aj_arxiv.py <aj_text> <arxiv_md> <out_txt>
Normalises whitespace/markup, splits into sentences, and lists sentences present in one version
but without a close match (difflib ratio < 0.85) in the other. Output is for human review only.
"""
import re, sys, difflib

def sentences(t):
    t = re.sub(r'Image \d+: \$[^$]*\$', '', t)
    t = re.sub(r'\$[^$]*\$', '', t)
    t = re.sub(r'[_*]', '', t)
    t = re.sub(r'\s+', ' ', t)
    s = re.split(r'(?<=[.!?])\s+(?=[A-Z(\[])', t)
    return [x.strip() for x in s if len(x.strip()) > 40]

aj_raw = open(sys.argv[1]).read().splitlines()
# keep the article body: from the abstract heading to the acknowledgments paragraph
start = next(i for i, l in enumerate(aj_raw) if l.startswith('We report [O iii]'))
end = next(i for i, l in enumerate(aj_raw) if 'not the sole discoverer' in l) + 1
aj = sentences('\n'.join(aj_raw[start:end]))
ax_raw = open(sys.argv[2]).read()
ax = sentences(ax_raw[ax_raw.find('###### Abstract'):])
norm = lambda x: re.sub(r'[^a-z0-9]', '', x.lower())
axn = [norm(x) for x in ax]; ajn = [norm(x) for x in aj]
out = []
def best(s, pool):
    return max((difflib.SequenceMatcher(None, s, p).quick_ratio() for p in pool), default=0)
out.append('=== Sentences in AJ (published) without close match in arXiv v1 ===')
for s, n in zip(aj, ajn):
    r = best(n, axn)
    if r < 0.85:
        out.append(f'[{r:.2f}] {s}')
out.append('\n=== Sentences in arXiv v1 without close match in AJ ===')
for s, n in zip(ax, axn):
    r = best(n, ajn)
    if r < 0.85:
        out.append(f'[{r:.2f}] {s}')
open(sys.argv[3], 'w').write('\n'.join(out))
print('\n'.join(out))
