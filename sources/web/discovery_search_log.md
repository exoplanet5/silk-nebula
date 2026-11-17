# Discovery topic: search and access log (2026-11-16)

## Literature / bibliographic searches

| Query / endpoint | Result |
|---|---|
| arXiv API `all:"Silk Nebula"` (https://export.arxiv.org/api/query?search_query=all:%22Silk%20Nebula%22) | 1 hit: 2610.13329v1 only (saved data/raw/discovery_biblio/arxiv_q_silk.xml) |
| arXiv API `abs:Monogem` sorted by date (30 newest) | Newest relevant: 2608.26666 (Fesen+); others are pulsar-halo/TeV papers; nothing else on the optical arc (arxiv_q_monogem.xml) |
| OpenAlex full-text `"Silk Nebula"` | 1 work: AJ 10.3847/1538-3881/aeb11c (2026-11-04), cited_by_count 0 |
| OpenAlex `cites:W7221525878` | 0 citing works |
| OpenAlex full-text `"OIII arc in Gemini"` | 0 works |
| Crossref works/10.3847/1538-3881/aeb11c and /10.3847/2515-5172/adfec7 | AJ and RNAAS metadata incl. received/accepted dates and reference list |
| DataCite 10.5281/zenodo.17215566 | "OIII arc in Gemini", S. Ziegenbalg, issued 2025-09-27, derived from 10.5281/zenodo.16780036 |
| Semantic Scholar Graph API (arXiv:2610.13329, arXiv:2608.26666) | HTTP 429 (rate limited) twice: **not obtained** |
| NASA ADS | **Not queried**: no ADS API token available in this environment. ADS citation and full-text searches are a data gap; OpenAlex/Crossref were used instead |

## Web searches (Tavily, Exa, WebSearch)

- "Silk Nebula" Gemini OIII (English): finds only the arXiv/AJ paper and aggregators (alphaxiv, papers.cool, arxivtldr, pubdb, themoonlight, wispaper.ai 2026-11-02), the QHYCCD note, the AstroBin pages, the NSNS page and the Cloudy Nights thread. **No Sky & Telescope, Universe Today, phys.org or space.com article found** as of 2026-11-16.
- 丝绸星云 LAMOST 双子座: Xinmin/Shangguan, Sina and ifeng (2026-11-16), Sun Guoyou blog (2024-05-21), arXivDaily. **No NAOC (nao.cas.cn) press release found** by these searches.
- Koposov DESI [O III] Gemini / sky fibres; "Sergey Koposov DESI sky spectra map of [O III]": **no statement, paper, talk or post by Koposov or Belokurov about the Silk/Gemini arc found.** The only source is S. Jackson's AstroBin description (zwfplw, 2025-05-02), which is secondhand ("Prof Vasily Belokurov told me about how Prof Sergey Koposov found it by ingenious analysis of DESI data").
- Mel Bartels visual observation: only his one-line Cloudy Nights reply (2025-03-04). The linked observing report was not found.
- Stargazers Lounge topic 420847 "Help me identify some nebulae" (cited by Fesen+ footnote 27): read via Exa (Cloudflare blocked curl). It is about alexbb's April 2024 Monoceros [O III] mosaic (Sigma 135 mm, Rosette to Hu 6), with replies suggesting the Monogem Ring. The mosaic does not reach the Silk (AstroBin vowsu0 plate solution: centre 101.1, +11.4; radius 12.6°, so the Silk is >15° from the centre). **Not a Silk precovery.**

## Access problems (recorded as data gaps)

- zenodo.org: HTTP 403 ("unusual traffic"). Zenodo record metadata was taken from DataCite instead; the snapshot files were not downloaded.
- AstroBin image pages: HTTP 403 to curl. Read via the AstroBin public JSON API (`/api/v2/images/image/?hash=`, `/api/v2/forum/...`, `/api/v2/nestedcomments/...`) and Tavily extracts.
- Stargazers Lounge: Cloudflare challenge for curl; read with Exa.
- vizier.cds.unistra.fr: HTTP 404 "Page unavailable" for an ASU query. The vizier.cfa.harvard.edu mirror worked (HIP 37797: RA 116.25248502, Dec +15.89503528; data/raw/discovery_vizier/).
- China National Astronomy (《中国国家天文》) 2024 No. 4: not available online. Not seen.
- The NSNS full-resolution JavaScript-viewer images and the Hunting Horn FITS mosaic (Google Drive link in AstroBin zwfplw) were not downloaded. They belong to the imaging teammates and the files are large.
