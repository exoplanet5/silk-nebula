"""List entries of an arXiv API Atom response (data/raw/context_arxiv/*.xml). Usage: python -I context_arxiv_list.py file.xml [grep]"""
import sys, xml.etree.ElementTree as ET
ns = {'a': 'http://www.w3.org/2005/Atom'}
r = ET.parse(sys.argv[1]).getroot()
pat = sys.argv[2].lower() if len(sys.argv) > 2 else None
for e in r.findall('a:entry', ns):
    t = ' '.join(e.find('a:title', ns).text.split())
    s = ' '.join(e.find('a:summary', ns).text.split())
    if pat and pat not in (t + s).lower():
        continue
    print(e.find('a:id', ns).text.strip(), '|', e.find('a:published', ns).text[:10], '|', t[:130], '|',
          e.find('a:author/a:name', ns).text)
    if '-v' in sys.argv:
        print('   ', s[:1500])
