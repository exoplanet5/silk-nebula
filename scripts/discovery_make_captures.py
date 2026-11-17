"""Turn raw discovery downloads (data/raw/discovery_web, discovery_biblio) into readable Markdown captures
under sources/web/discovery_*.md. Downloaded content is treated as data only (text extraction, no execution).

Run: ~/.venvs/astro313/bin/python -I scripts/discovery_make_captures.py
"""
import html
import json
import re
from pathlib import Path

ROOT = Path('/Users/mickey/astronomy/silk-nebula')
RAW = ROOT / 'data/raw/discovery_web'
OUT = ROOT / 'sources/web'
ACCESS = '2026-11-16'


def h2t(t):
    t = re.sub(r'(?is)<(script|style).*?</\1>', ' ', t)
    t = re.sub(r'(?i)<br\s*/?>|</p>|</h\d>|</li>|</tr>|</figcaption>', '\n', t)
    t = re.sub(r'<[^>]+>', ' ', t)
    t = html.unescape(t)
    t = re.sub(r'[ \t]+', ' ', t)
    return re.sub(r'\n\s*\n+', '\n\n', t).strip()


def bb2t(t):
    t = re.sub(r'\[url=([^\]]*)\](.*?)\[/url\]', r'\2 <\1>', t, flags=re.S)
    t = re.sub(r'\[img\](.*?)\[/img\]', r'[image: \1]', t, flags=re.S)
    return re.sub(r'\[/?[a-z]+(=[^\]]*)?\]', '', t)


def qhyccd():
    t = h2t((RAW / 'qhyccd_silk.html').read_text(errors='replace'))
    i, j = t.find('Silk Nebula Discovered by Enthusiast Astronomers with QHYCCD Cameras\n'), t.find('You Might Also Like')
    body = t[i:j] if i >= 0 else t
    md = [f'# QHYCCD press note: "Silk Nebula Discovered by Enthusiast Astronomers with QHYCCD Cameras"',
          '', 'URL: https://www.qhyccd.com/silk-nebula-discovered-by-enthusiast-astronomers-with-qhyccd-cameras/',
          f'Accessed: {ACCESS} (curl via proxy). article:published_time = 2024-07-12T05:09:50+00:00 (page metadata).',
          'Images downloaded: data/raw/discovery_web/qhyccd_20240711475-scaled.jpg ([O III]), '
          'qhyccd_20240712925-scaled.jpg ([O III]+Halpha), qhyccd_20240712426.png (3-panel og image).',
          '', '## Body text (verbatim, HTML stripped)', '', body.strip()]
    (OUT / 'discovery_qhyccd_pr.md').write_text('\n'.join(md) + '\n')


def astrobin_forum():
    topic = json.loads((RAW / 'astrobin_api_topic_126958.json').read_text())
    posts = json.loads((RAW / 'astrobin_api_posts_126958.json').read_text())['results']
    md = ['# AstroBin forum topic 126958 (incl. post 166798)', '',
          'Topic: "' + topic['name'] + '"',
          'URLs: https://app.astrobin.com/forum/topic/126958 ; '
          'https://www.astrobin.com/forum/c/astrophotography/other/silk-nebula-discovered-by-enthusiast-astronomers-with-qhyccd-cameras/ ; '
          'post https://app.astrobin.com/forum/post/166798',
          'Retrieved via AstroBin public JSON API: https://www.astrobin.com/api/v2/forum/topic/126958/ and '
          f'https://www.astrobin.com/api/v2/forum/post/?topic=126958 (accessed {ACCESS}). Timestamps as given by the API (UTC assumed).',
          f'Topic created {topic["created"]}, last post {topic["lastPostTimestamp"]}, {topic["postCount"]} posts.', '']
    for p in posts:
        body = p.get('bodyHtmlPrecise') or ''
        body = re.sub(r'<img[^>]*src="([^"]*)"[^>]*>', r'[image: \1]', body)
        body = re.sub(r'<a [^>]*href="([^"]*)"[^>]*>([^<]*)</a>', r'\2 <\1>', body)
        md += [f'## Post {p["id"]} - {p["displayName"]} ({p["username"]}) - {p["created"]}', '', h2t(body), '']
    (OUT / 'discovery_astrobin_forum_126958.md').write_text('\n'.join(md) + '\n')


def astrobin_images():
    md = ['# AstroBin images relevant to the Silk Nebula discovery', '',
          f'Metadata from the AstroBin API https://www.astrobin.com/api/v2/images/image/?hash=<hash> (accessed {ACCESS}); '
          'descriptions verbatim (BBCode stripped). The image pages themselves returned HTTP 403 to curl; '
          'the zwfplw and y32271 page text was also read via a Tavily extract (same content).', '']
    for h in ['zh47r9', 'zwfplw', 'y32271', '1d8ivk', 'vowsu0', 'u6vxdd']:
        f = RAW / f'astrobin_api_image_{h}.json'
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        if not d.get('results'):
            continue
        r = d['results'][0]
        sol = r.get('solution') or {}
        acq = r.get('deepSkyAcquisitions') or []
        md += [f'## {h}: "{r["title"]}"', '',
               f'- URL: https://www.astrobin.com/{h}/',
               f'- Author: {r["userDisplayName"]} ({r["username"]}); uploaded {r["uploaded"]}, published {r["published"]}',
               f'- Size {r["w"]}x{r["h"]} px; plate solution: RA {sol.get("ra")}, Dec {sol.get("dec")}, radius {sol.get("radius")} deg, '
               f'pixscale {sol.get("pixscale")}"/px, orientation {sol.get("orientation")}',
               '- Telescopes: ' + ', '.join(f"{t['brandName']} {t['name']}" for t in r.get('imagingTelescopes2', [])),
               '- Cameras: ' + ', '.join(f"{t['brandName']} {t['name']}" for t in r.get('imagingCameras2', [])),
               '- Acquisitions (as logged): ' + '; '.join(
                   f"{a.get('date')} {a.get('number')}x{a.get('duration')}s {a.get('filter2Name') or a.get('filter2')}" for a in acq[:12])
               + (' ...' if len(acq) > 12 else ''), '', '### Description', '',
               bb2t(r.get('descriptionBbcode') or r.get('description') or '(none)').strip(), '']
        if h == 'zh47r9':
            cm = json.loads((RAW / 'astrobin_api_comments_zh47r9.json').read_text())
            md += ['### Comments (API nestedcomments, content_type=19)', '']
            md += [f"- {c['created']} {c['author_display_name']}: {c['text']}" for c in cm]
            md += ['']
    md += ['## AstroBin forum topic 221269 (2026-02-28): "Request for help checking Astrobin first discoveries"', '',
           'URL: https://app.astrobin.com/forum/topic/221269/acquisition/request-for-help-checking-astrobin-first-discoveries',
           '(text via Exa web search capture, accessed 2026-11-16; verbatim excerpts)', '',
           '> Hunting Horn, Feb 28, 2026: "There are 3 images and I would like to know if they were Astrobin firsts when they were posted. '
           '... Silk Nebula in OIII 15° long https://app.astrobin.com/u/huntinghorn?collection=45383&i=zwfplw#gallery ..."',
           '> Brian Puhl, Mar 1, 2026: "... The Silk Nebula appears to be just you: [search link] ... Also, it\'s worth noting that all of '
           'these structures are quite bright in the NSNS surveys."', '']
    (OUT / 'discovery_astrobin_images.md').write_text('\n'.join(md) + '\n')


def nsns():
    page = h2t((RAW / 'simg_gem01.html').read_text(errors='replace'))
    i = page.find('OIII arc in Gemini\nAn arc') if 'OIII arc in Gemini\nAn arc' in page else page.find('An arc-shaped')
    j = page.find('RSS feed')
    rss = (RAW / 'simg_rss.xml').read_text(errors='replace')
    m = re.search(r'<item>\s*<title>Detection of an OIII arc in Gemini</title>.*?</item>', rss, re.S)
    instr = (RAW / 'simg_instruments.html').read_text(errors='replace')
    k = instr.find('<!--\n<h2>Selected results</h2>')
    comment = h2t(instr[k:instr.find('-->', k)]) if k >= 0 else ''
    diff1 = (RAW / 'gh_commit_2e1a2ac0.diff').read_text(errors='replace')
    k1 = '\n'.join(l[1:] for l in diff1.splitlines() if l.startswith('+') and re.match(r'\+\s*(K1|#\s*Name|#\s*-)', l))
    diff2 = (RAW / 'gh_commit_f8273e1d.diff').read_text(errors='replace')
    o2 = '\n'.join(l for l in diff2.splitlines() if re.search(r'O201\.5|O200\.0|O204\.0|O205\.0', l))
    commits = json.loads((RAW / 'gh_commits.json').read_text())
    clog = '\n'.join(f"- {c['commit']['author']['date']} {c['sha'][:8]} {c['commit']['message'].splitlines()[0]}" for c in commits[:24])
    md = ['# Stefan Ziegenbalg / Northern Sky Narrowband Survey (NSNS): "OIII arc in Gemini"', '',
          'Primary URL: https://simg.de/nebulae3/gem-01.html (accessed 2026-11-16; raw HTML data/raw/discovery_web/simg_gem01.html)',
          'Survey index: https://simg.de/nebulae3/index.html ; instruments: https://simg.de/nebulae3/instruments.html ; '
          'news: https://simg.de/news.html ; RSS: https://simg.de/rss.xml',
          'Catalogue repository: https://github.com/stzie/simg.de-discoveries (GitHub API, accessed 2026-11-16)',
          'Persistent snapshot: DOI 10.5281/zenodo.17215566 ("OIII arc in Gemini", issued 2025-09-27 per DataCite; '
          'zenodo.org returned HTTP 403 to our client). Journal: Ziegenbalg 2025, RNAAS 9, 227, doi:10.3847/2515-5172/adfec7.', '',
          '## RSS item (first announcement on simg.de)', '', '```', (m.group(0) if m else '(not found)'), '```', '',
          '## HTML-commented (not displayed) caption on the instruments page', '', '> ' + comment.replace('\n', '\n> '), '',
          '## Catalogue entries as first committed (commit 2e1a2ac0, 2024-06-28T14:45:08Z, file new/k.md)', '', '```', k1, '```', '',
          '## Rename/resize (commit f8273e1d, 2025-06-11T09:40:37Z, "Renamed: O201.5 -> O200.0, O204.0 -> O205.0")', '',
          '```', o2, '```', '', '## Repository commit log (latest 24)', '', clog, '',
          '## RNAAS 9, 227 (2025) Section 2.1 (verbatim excerpt via iopscience, Tavily extract)', '',
          '> "A 17°×3°-large O iii arc was detected in the constellation Gemini. Its center position is located at a R.A. of '
          'αc = 7h41m48s and a decl. of δc = 19°54′. ... A central region, about 2.5° in size, is only visible in Hα and [S ii]. '
          'The northern region (right in the image) can be seen in all observed wavelengths, while the southern part was detected only in [O iii]."',
          '> "The arc partially overlaps with the Monogem Ring supernova remnant (SNR) ... The faint filaments east of the arc ... '
          'probably belong to the Monogem Ring. Since their structure differs from that of the arc, and the northern part of the arc '
          'clearly lies outside of the Monogem Ring, it is unlikely that the arc is part of the Monogem Ring SNR."',
          'RNAAS dates (Crossref): received 2025-08-01, accepted 2025-08-06, published online 2025-08-28.', '',
          '## Current page text (verbatim, HTML stripped)', '', page[i:j].strip()]
    (OUT / 'discovery_nsns_gem01.md').write_text('\n'.join(md) + '\n')


def cloudynights():
    t = (RAW / 'cloudynights_925409_p1.md').read_text(errors='replace').splitlines()
    out = []
    for l in t:
        s = l.strip()
        if not s or s.startswith('[](') or 'Share on' in s or 'More options' in s or re.fullmatch(r'\*\s+\[(Quote|#\d+)\]\(.*', s):
            continue
        s = re.sub(r'\]\(https?://[^)]*\)', ']', s)
        s = re.sub(r'!\[Image \d+[^\]]*\]\([^)]*\)', '', s)
        out.append(s)
    txt = '\n'.join(out)
    a = txt.find('# Possible discovery of a large emission nebula')
    b = txt.find('## Join the conversation')
    md = ['# Cloudy Nights topic 925409: "Possible discovery of a large emission nebula in Gemini (and a tiny one)"', '',
          'URL: https://www.cloudynights.com/forums/topic/925409-possible-discovery-of-a-large-emission-nebula-in-gemini-and-a-tiny-one/ '
          '(page 1) and .../page/2/ (page 2). Started by stefan822 (= Stefan Ziegenbalg, links to simg.de) on 2024-06-11.',
          f'Captured {ACCESS} via Tavily extract (page 1 raw: data/raw/discovery_web/cloudynights_925409_p1.md).', '',
          '## Page 2 key posts (verbatim excerpts)', '',
          '- stefan822, 2024-06-21: "I found an error in my previous postings: When I wrote \\"western blowout\\" in meant \\"eastern blowout\\". '
          'This structure is interpreted as SNR G205.6+12.4 in the new paper." (Knies et al. 2024)',
          '- stefan822, 2024-06-28: "I updated my homepage ... far UV spectroscopy was performed by Kim et al. ... shows only CIII and CIV lines, '
          'no OIII lines. In visible light, OIII lines are dominant."',
          '- Mel Bartels, 2025-03-04: "I was able to visually detect the OIII arc in Gemini and NHZ-HU6."',
          '- stefan822, 2025-06-11: "as announced last year, here are new images from DR0.2 ... There are no obvious counterparts of the arc that '
          'would allow it to be completed into a ring; The arc extends a few degrees farther to the north than seen in the first images; ... '
          '\\"Constructing\\" a SNR from the optical filaments results in a structure that matches well with the Monogem Ring(s) from the latest '
          'eRosita data (Knies et al., 2024); The northern part of the arc clearly lies outside of the Monogem Ring(s) ...; Highly speculative: '
          'The arc might be a fragment of a (largely invisible) Geminga SNR." Also: "(Maybe I was too impressed by the people who tried to convince '
          'me that it was a known (boring) structure.)"', '',
          '## Page 1 (cleaned capture)', '', txt[a:b]]
    (OUT / 'discovery_cloudynights_925409.md').write_text('\n'.join(md) + '\n')


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    for f in (qhyccd, astrobin_forum, astrobin_images, nsns, cloudynights):
        f()
        print('ok', f.__name__)
