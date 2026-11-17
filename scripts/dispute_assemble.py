"""Assemble DISPUTE.md from the section drafts in notes/dispute/ plus the frame text in this script.

- Relabels this Comment's own figures/tables with a 'C' prefix (Figure C1..C6, Table C1..C7) so they cannot be
  confused with the paper's Figures 1-8 and Tables 1-4. Paper references ('Table 4', 'Fig. 3a', 'Sect. 3.3') stay.
- Strips each section's 'References used in this section' list (one consolidated, verified list is appended).
- Expands informal shorthand (Reid+14, Reid+19, Wu+26).
- Bilingual Summary: notes/dispute/frame/summary_en.md + summary_zh.md (the rest is English only).
usage: ~/.venvs/astro313/bin/python scripts/dispute_assemble.py
"""
from pathlib import Path

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
D = ROOT / 'notes/dispute'
FRAME = ROOT / 'notes/dispute/frame'


def section(name, replacements):
    t = (D / f'{name}.md').read_text()
    cut = t.find('**References used in this section.**')
    if cut > 0:
        t = t[:cut].rstrip() + '\n'
    for old, new in replacements:
        if old not in t:
            raise SystemExit(f'[{name}] pattern not found: {old!r}')
        t = t.replace(old, new)
    for old, new in [('Reid+14', 'Reid et al. 2014'), ('Reid+19', 'Reid et al. 2019'), ('Wu+26', 'Wu et al. (2026)')]:
        t = t.replace(old, new)
    return t.strip() + '\n'


SECTIONS = [
    ('distance', [
        ('"Table 4", "Fig. 5b" and "Sect." refer to the paper; "Figure 1", "Table 2" and "Section" refer to this Comment.', ''),
        ('(Figure 1a, inset)', '(Figure C1a, inset)'),
        ('(our Figure 1b)', '(our Figure C1b)'),
        ('(Figure 1d, Table 2)', '(Figure C1d, Table C2)'),
        ('*Figure 1. Table 4', '*Figure C1. Table 4'),
        ('**Table 2.** Row-by-row', '**Table C2.** Row-by-row'),
    ]),
    ('kinematic', [
        ('(Fig. 2a)', '(Figure C2a)'),
        ('(Fig. 2b)', '(Figure C2b)'),
        ('(Fig. 2c, d; Table 3)', '(Figure C2c, d; Table C3)'),
        ('*Figure 2a–d.', '*Figure C2a–d.'),
        ('**Table 3.**', '**Table C3.**'),
    ]),
    ('oiii', [
        ('"The paper\'s Table 4" is Wu et al.\'s table; "Table 4" alone is ours, below.', ''),
        ('Per region (Table 4)', 'Per region (Table C4)'),
        ('*Figure 3a–c.', '*Figure C3a–c.'),
        ('**Table 4.**', '**Table C4.**'),
    ]),
    ('halpha', [
        ('(Table 5)', '(Table C5)'),
        ('(Fig. 4a, b)', '(Figure C4a, b)'),
        ('*Figure 4a–d.', '*Figure C4a–d.'),
        ('**Table 5.**', '**Table C5.**'),
        ('in kpc per Sect. 1)', 'in kpc per Section 1)'),
        ('no information here (Sect. 2)', 'no information here (Section 2)'),
    ]),
    ('brightness', [
        ('## 5. Table 3 overestimates the [O III] surface brightness by a factor of ~5–9',
         '## 5. Table 3 overestimates the [O III] surface brightness by a factor of 4.6–8.6'),
        ('(Fig. 5a)', '(Figure C5a)'),
        ('(Table 6, Fig. 5b)', '(Table C6, Figure C5b)'),
        ('(Fig. 5d)', '(Figure C5d)'),
        ('(Fig. 5c)', '(Figure C5c)'),
        ('*Figure 5.', '*Figure C5.'),
        ('**Table 6.**', '**Table C6.**'),
    ]),
    ('other', [
        ('(Fig. 6a)', '(Figure C6a)'),
        ('(Fig. 6b)', '(Figure C6b)'),
        ('Fig. 6c)', 'Figure C6c)'),
        ('(Fig. 6d)', '(Figure C6d)'),
        ('*Figure 6.', '*Figure C6.'),
        ('**Table 7.**', '**Table C7.**'),
    ]),
]

parts = [(FRAME / 'head.md').read_text().strip() + '\n']
# Bilingual Summary only: the English and Chinese summaries sit between language markers. Markdown viewers show both;
# scripts/make_report_html.py turns the markers into blocks toggled by the EN | 中文 switch.
parts.append('<!-- lang:en -->\n' + (FRAME / 'summary_en.md').read_text().strip() + '\n<!-- /lang -->\n\n'
             + '<!-- lang:zh -->\n' + (FRAME / 'summary_zh.md').read_text().strip() + '\n<!-- /lang -->\n')
for name, repl in SECTIONS:
    parts.append(section(name, repl))
parts.append((FRAME / 'tail.md').read_text().strip() + '\n')
refs = (D / 'references.md').read_text()
refs = refs[:refs.find('### Verification log')].rstrip() + '\n'
refs = refs.replace('This list covers every work cited in Sections 1–6 of this Comment.',
                    'This list covers every work cited in this Comment. The full verification log is in `notes/dispute/references.md`.')
parts.append(refs)
out = ROOT / 'DISPUTE.md'
out.write_text('\n---\n\n'.join(parts))
print('wrote', out, len(out.read_text().split()), 'words')
