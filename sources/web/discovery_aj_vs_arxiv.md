# Wu et al. 2026: AJ published version vs arXiv v1

## Bibliographic record

- Wu, Y.-F., Jia, S.-Y., Wu, C.-J., Zhang, C., et al. (17 authors), "The Silk Nebula: A Newly Discovered Giant [O III] Emission Filament", **AJ 172(6), 298**, doi:10.3847/1538-3881/aeb11c, CC-BY 4.0.
- Dates from Crossref (`assertion` field) and iopscience: **received 2026-07-30, accepted 2026-10-05, published online 2026-11-04**; print issue 2026-12. Crossref lists 63 references and an `is-referenced-by-count` of 0. Metadata saved as data/raw/discovery_biblio/crossref_aeb11c.json.
- arXiv:2610.13329v1 was submitted 2026-10-09 ("20 pages, 8 figures, accepted for publication in AJ"); the PDF header reads "Draft version October 12, 2026", AASTeX 7.0.1. There is only one arXiv version.
- OpenAlex W7221525878: publication date 2026-11-04, cited_by_count 0. An OpenAlex `cites:` query returned 0 citing works (accessed 2026-11-16).
- AJ full text was read through a Tavily extract of https://iopscience.iop.org/article/10.3847/1538-3881/aeb11c (accessed 2026-11-16; saved as data/raw/discovery_biblio/aj_aeb11c_iopscience_tavily.txt). The figures were not re-downloaded.

## Differences found (AJ vs arXiv v1)

A sentence-level diff (scripts/discovery_diff_aj_arxiv.py, output data/derived/discovery_aj_vs_arxiv_diff.txt) is dominated by changes in citation style (numbered → author-year) and copy-editing. Besides those, I found these content differences:

1. **New acknowledgment of independent discoverers (AJ only)** — verbatim:
   > "The SWONS team was not the sole discoverer of the Silk Nebula. It thanks the team of Fesen, Marcel Drechsler, and Tim Schaeffer, as well as the team of Stefan Ziegenbalg, for their contributions to the discovery of the Silk Nebula; each of these teams independently discovered the Silk Nebula as well."

   This paragraph is absent from arXiv v1 (sources/arxiv/2610.13329v1.txt, ACKNOWLEDGEMENTS).
2. **Reference list**:
   - Ziegenbalg, S. 2025, RNAAS 9, 227 (NSNS DR0.2; first journal description of the 17° "O III arc in Gemini") is in the AJ list (Crossref key ajaeb11cbib63) and absent from v1. I did not locate its in-text citation in the captured HTML.
   - Lumbreras-Calle et al. 2025, A&A 704, A224 is added. In the AJ text it is cited next to Fesen et al. 2020, 2024 for "[O III] traces shock structures more extended than Hα".
   - The garbled v1 entry "Slane, P. O., Badenes, C., Yamaguchi, H., et al. 2013" becomes H. Yamaguchi et al. 2014, ApJ 780, 136.
   - Bond et al. 2025 changes from arXiv:2506.11306 to AJ 170, 137.
3. **Unchanged substantive claims**: "The Silk Nebula was discovered and named by the SWONS team" (Conclusion); "In 2023 November, a candidate north–south filamentary [O III] structure was detected" (Sect. 2.1); 1°×10° extent; 0.5°×1° Hα structure at the north end; P1–P4 coordinates; Tables 2–4 values; "So the Silk Nebula has also been detected in Fesen's project (R. A. Fesen et al. 2026), and it is considered to be part of the Monogem Ring"; the verdict (possible SNR; DIG or ionized bubble not excluded); E = (2.0±0.7)×10^51 erg; age (6.8±1.8)×10^4 yr. In the captured AJ Table 4 the Monte Carlo distance column (superscript-subscript ranges) did not render. The traditional distances (1.79, 1.85, ... ×100 pc) are identical to v1.

## Note on chronology

AJ received the manuscript 2026-07-30. That is 4 weeks before Fesen et al. appeared on arXiv (2026-08-27), yet both the v1 and the AJ text cite "Fesen et al. 2026, arXiv:2608.26666". That identifier did not exist before late August 2026, so the Fesen citation (with the Monogem sentence) was added or updated during revision, before acceptance on 2026-10-05. The acknowledgment of independent discoverers was added after the arXiv v1 posting on 2026-10-09 (proof stage). Whether the submitted version mentioned Fesen's work in another form cannot be checked from public material.
