"""Render REPORT.md (+ REPORT_zh.md if present) into one self-contained, bilingual REPORT.html (mistune, CommonMark lists).

- Both languages are rendered into separate <article> elements, and a fixed EN | 中文 switch sits at the top right.
- Figures are embedded once (base64, in the English article). The Chinese article's <img data-fig> elements
  borrow those sources via JS, so the page does not double in size. English works without JS.
- Initial language: URL hash (#en / #zh) > saved choice (localStorage) > browser language > English.
- Switching languages keeps the reader at the same section (headings correspond one-to-one).

usage: ~/.venvs/astro313/bin/python scripts/make_report_html.py            # REPORT.md (+ REPORT_zh.md) -> REPORT.html
       ~/.venvs/astro313/bin/python scripts/make_report_html.py DISPUTE    # DISPUTE.md (+ DISPUTE_zh.md) -> DISPUTE.html
"""
import base64
import re
import sys
from pathlib import Path

import html as htmllib

import mistune
from mistune.renderers.html import HTMLRenderer
from markdown.extensions.toc import slugify_unicode   # only the slug helper is used

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
DOC = sys.argv[1] if len(sys.argv) > 1 else 'REPORT'
TITLES = {'REPORT': ('Silk Nebula Study', '丝绸星云研究报告'),
          'DISPUTE': ('Comment on Wu et al. 2026 (Silk Nebula)', '对 Wu 等（2026）丝绸星云论文的评论'),
          'VELOCITY': ('Silk Nebula velocities', '丝绸星云速度对比')}
_t = TITLES.get(DOC, (DOC, DOC))
LANGS = [  # (code, html lang, markdown file, page title, toc label)
    ('en', 'en', f'{DOC}.md', _t[0], 'Contents'),
    ('zh', 'zh-Hans', f'{DOC}_zh.md', _t[1], '目录'),
]

CSS = """
:root{--bg:#fcfcfb;--fg:#0b0b0b;--muted:#52514e;--rule:#e4e3df;--accent:#2a78d6;--code:#f3f2ef;--fig:#101418;--chip:#ffffff}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#1a1a19;--fg:#f0efec;--muted:#c3c2b7;--rule:#383835;--accent:#3987e5;--code:#262624;--fig:#101418;--chip:#232321}}
:root[data-theme="dark"]{--bg:#1a1a19;--fg:#f0efec;--muted:#c3c2b7;--rule:#383835;--accent:#3987e5;--code:#262624;--fig:#101418;--chip:#232321}
html{background:var(--bg)}
body{background:var(--bg);color:var(--fg);font:15px/1.6 -apple-system,BlinkMacSystemFont,"Helvetica Neue",Arial,sans-serif;
     max-width:1120px;margin:0 auto;padding:24px 16px 64px}
article:lang(zh),.lang-block:lang(zh){font-family:-apple-system,BlinkMacSystemFont,"PingFang SC","Hiragino Sans GB","Noto Sans CJK SC","Microsoft YaHei",sans-serif;line-height:1.75}
[hidden]{display:none!important}
h1{font-size:1.65em;line-height:1.25;margin:.2em 0 .6em;padding-right:120px}
h2{font-size:1.3em;margin-top:2.2em;padding-top:.4em;border-top:1px solid var(--rule)}
h3{font-size:1.08em;margin-top:1.6em}
a{color:var(--accent)}
em{color:var(--muted)}
code{font-family:Menlo,Consolas,monospace;font-size:.86em;background:var(--code);padding:.08em .3em;border-radius:4px}
pre{background:var(--code);padding:12px 14px;border-radius:8px;overflow-x:auto;font-size:.85em;line-height:1.45}
pre code{background:none;padding:0}
table{border-collapse:collapse;margin:1em 0;font-size:.88em;display:block;overflow-x:auto;max-width:100%}
th,td{border-bottom:1px solid var(--rule);padding:5px 10px;text-align:left;vertical-align:top}
th{color:var(--muted);font-weight:600}
img{max-width:100%;height:auto;display:block;margin:1.2em auto .3em;background:var(--fig);border-radius:6px}
hr{border:none;border-top:1px solid var(--rule);margin:2em 0}
li{margin:.25em 0}
hr + h2{border-top:none;margin-top:1.2em}
.toc{font-size:.9em;border:1px solid var(--rule);border-radius:8px;padding:.4em 1.2em;margin:1em 0}
.lang-switch{position:fixed;top:12px;right:16px;z-index:20;display:flex;border:1px solid var(--rule);border-radius:999px;
             background:var(--chip);box-shadow:0 1px 4px rgba(0,0,0,.12);overflow:hidden}
.lang-switch button{font:600 13px/1 -apple-system,BlinkMacSystemFont,"PingFang SC","Helvetica Neue",sans-serif;color:var(--muted);
             background:transparent;border:0;padding:8px 14px;cursor:pointer;min-width:44px;min-height:32px}
.lang-switch button[aria-pressed="true"]{background:var(--accent);color:#fff}
.lang-switch button:focus-visible{outline:2px solid var(--accent);outline-offset:-2px}
@media print{.lang-switch{display:none}}
"""

JS = """
(function(){
  var langs=['en','zh'], titles=%(titles)s, htmlLang=%(htmllang)s;
  var arts={}; langs.forEach(function(l){arts[l]=document.getElementById('doc-'+l);});
  if(!arts.zh){return;}
  // Chinese figures borrow the embedded sources from the English article (figures stored once).
  var src={}; arts.en.querySelectorAll('img[data-fig]').forEach(function(i){src[i.getAttribute('data-fig')]=i.getAttribute('src');});
  arts.zh.querySelectorAll('img[data-fig]').forEach(function(i){var s=src[i.getAttribute('data-fig')]; if(s){i.setAttribute('src',s);}});
  var buttons=document.querySelectorAll('.lang-switch button'), current=null;
  function heads(l){return Array.prototype.slice.call(arts[l].querySelectorAll('h1,h2,h3'));}
  function visibleIndex(l){var hs=heads(l), idx=0; for(var k=0;k<hs.length;k++){if(hs[k].getBoundingClientRect().top<=80){idx=k;}} return idx;}
  function setLang(l, keepPlace){
    if(langs.indexOf(l)<0){l='en';}
    if(l===current){return;}
    var idx=(keepPlace&&current)?visibleIndex(current):-1;
    langs.forEach(function(x){arts[x].hidden=(x!==l);});
    document.documentElement.lang=htmlLang[l]; document.title=titles[l];
    buttons.forEach(function(b){b.setAttribute('aria-pressed', b.getAttribute('data-lang')===l?'true':'false');});
    try{localStorage.setItem('silk-'+%(doc)s+'-lang',l);}catch(e){}
    if(history.replaceState){history.replaceState(null,'','#'+l);}
    current=l;
    if(idx>0){var hs=heads(l); if(hs[idx]){hs[idx].scrollIntoView();}} else if(idx===0){window.scrollTo(0,0);}
  }
  buttons.forEach(function(b){b.addEventListener('click',function(){setLang(b.getAttribute('data-lang'),true);});});
  document.addEventListener('click',function(e){var a=e.target.closest&&e.target.closest('a[data-switch]'); if(a){e.preventDefault(); setLang(a.getAttribute('data-switch'),false);}});
  var h=(location.hash||'').replace('#','').replace('lang=',''), stored=null;
  try{stored=localStorage.getItem('silk-'+%(doc)s+'-lang');}catch(e){}
  var nav=(navigator.language||'').toLowerCase().indexOf('zh')===0?'zh':'en';
  setLang(langs.indexOf(h)>=0?h:(stored||nav), false);
})();
"""

PARTIAL_JS = """
(function(){
  // Partial bilingual page: only elements marked data-lang (e.g. the Summary and its contents entry) switch language.
  var blocks=document.querySelectorAll('[data-lang]:not(button)'), buttons=document.querySelectorAll('.lang-switch button');
  function setLang(l){
    if(l!=='zh'){l='en';}
    blocks.forEach(function(b){b.hidden=(b.getAttribute('data-lang')!==l);});
    buttons.forEach(function(b){b.setAttribute('aria-pressed', b.getAttribute('data-lang')===l?'true':'false');});
    try{localStorage.setItem('silk-'+%(doc)s+'-lang',l);}catch(e){}
    if(history.replaceState){history.replaceState(null,'','#'+l);}
  }
  buttons.forEach(function(b){b.addEventListener('click',function(){setLang(b.getAttribute('data-lang'));});});
  var h=(location.hash||'').replace('#','').replace('lang=',''), stored=null;
  try{stored=localStorage.getItem('silk-'+%(doc)s+'-lang');}catch(e){}
  setLang(h==='zh'||h==='en'?h:(stored||'en'));
})();
"""


def b64(path):
    p = ROOT / path
    if not p.exists():
        print('missing figure:', p)
        return None
    return 'data:image/png;base64,' + base64.b64encode(p.read_bytes()).decode()


class Renderer(HTMLRenderer):
    """mistune (CommonMark) renderer: language-prefixed heading ids, and h2/h3 collected for the contents box."""

    def __init__(self, code):
        super().__init__(escape=True)
        self.code, self.heads, self.used, self.cur_lang = code, [], set(), None

    def heading(self, text, level, **attrs):
        plain = htmllib.unescape(re.sub(r'<[^>]+>', '', text))
        hid = f'{self.cur_lang or self.code}-' + (slugify_unicode(plain, '-') or 'section')
        base, k = hid, 2
        while hid in self.used:
            hid, k = f'{base}-{k}', k + 1
        self.used.add(hid)
        if level in (2, 3):
            self.heads.append((level, hid, text, self.cur_lang))
        return f'<h{level} id="{hid}">{text}</h{level}>\n'


def toc_html(heads, label):
    out, open_sub = [f'<nav class="toc"><strong>{label}</strong><ul>'], False
    for level, hid, text, blk in heads:
        dl = f' data-lang="{blk}"' + ('' if blk == 'en' else ' hidden') if blk else ''
        if level == 2:
            if open_sub:
                out.append('</ul></li>'); open_sub = False
            elif len(out) > 1:
                out.append('</li>')
            out.append(f'<li{dl}><a href="#{hid}">{text}</a>')
        else:
            if not open_sub:
                out.append('<ul>'); open_sub = True
            out.append(f'<li{dl}><a href="#{hid}">{text}</a></li>')
    out.append('</ul></li>' if open_sub else '</li>')
    out.append('</ul></nav>')
    return ''.join(out)


def render(code, md_file, toc_label, embed):
    r = Renderer(code)
    md = mistune.create_markdown(renderer=r, plugins=['table'])
    pieces = re.split(r'<!-- lang:(\w+) -->\n?(.*?)<!-- /lang -->', (ROOT / md_file).read_text(), flags=re.S)
    out = []
    for k in range(0, len(pieces), 3):
        r.cur_lang = None
        out.append(md(pieces[k]))
        if k + 2 < len(pieces):
            lang, text = pieces[k + 1], pieces[k + 2]
            r.cur_lang = lang
            hl = 'zh-Hans' if lang == 'zh' else lang
            hide = '' if lang == 'en' else ' hidden'
            out.append(f'<div class="lang-block" data-lang="{lang}" lang="{hl}"{hide}>{md(text)}</div>')
    r.cur_lang = None
    body = ''.join(out)
    global PARTIAL
    PARTIAL = PARTIAL or len(pieces) > 1
    body = body.replace('<hr />', toc_html(r.heads, toc_label) + '<hr />', 1)

    def img(m):
        path = m.group(1)
        if embed:
            data = b64(path)
            return f'data-fig="{path}" src="{data}"' if data else m.group(0)
        return f'data-fig="{path}"'
    body = re.sub(r'src="([^"]+\.png)"', img, body)
    # cross-language links become in-page switches
    body = re.sub(rf'href="{DOC}_zh\.md"', 'href="#zh" data-switch="zh"', body)
    body = re.sub(rf'href="{DOC}\.md"', 'href="#en" data-switch="en"', body)
    body = re.sub(r'href="([^"#:]+)\.md"', r'href="\1.html"', body)
    return body


PARTIAL = False
present = [L for L in LANGS if (ROOT / L[2]).exists()]
articles = []
for i, (code, hl, md_file, title, toc_label) in enumerate(present):
    body = render(code, md_file, toc_label, embed=(code == 'en'))
    hidden = '' if i == 0 else ' hidden'
    articles.append(f'<article id="doc-{code}" lang="{hl}"{hidden}>{body}</article>')

switch = ''
script = ''
if len(present) > 1:
    switch = ('<div class="lang-switch" role="group" aria-label="Language / 语言">'
              '<button type="button" data-lang="en" aria-pressed="true" lang="en">EN</button>'
              '<button type="button" data-lang="zh" aria-pressed="false" lang="zh-Hans">中文</button></div>')
    titles = '{' + ','.join(f'"{c}":"{t}"' for c, _, _, t, _ in present) + '}'
    hlang = '{' + ','.join(f'"{c}":"{h}"' for c, h, _, _, _ in present) + '}'
    script = '<script>' + JS % {'titles': titles, 'htmllang': hlang, 'doc': f'"{DOC.lower()}"'} + '</script>'
elif PARTIAL:
    switch = ('<div class="lang-switch" role="group" aria-label="Summary language / 概要语言">'
              '<button type="button" data-lang="en" aria-pressed="true" lang="en">EN</button>'
              '<button type="button" data-lang="zh" aria-pressed="false" lang="zh-Hans">中文</button></div>')
    script = '<script>' + PARTIAL_JS % {'doc': f'"{DOC.lower()}"'} + '</script>'

html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{present[0][3]}</title><style>{CSS}</style></head>
<body>{switch}{''.join(articles)}{script}</body></html>"""
out = ROOT / f'{DOC}.html'
out.write_text(html)
print('wrote', out, f'{len(html) / 1e6:.1f} MB', 'languages:', [c for c, *_ in present], '(bilingual blocks)' if PARTIAL else '')
