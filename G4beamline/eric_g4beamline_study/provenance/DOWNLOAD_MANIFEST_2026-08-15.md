# Prebys G4Beamline study download manifest

Downloaded on 2026-08-15 from the main study page and its
`g4beamline_study/` directory at `https://prebys.physics.ucdavis.edu/misc/AAAreadme/`.
The initial crawl was bounded to the main page, its two linked instruction
pages, and every artifact listed in that directory index. On 2026-08-15, the
general G4Beamline notes and the linked converter were added after a separate
relevance review. Archives remain unextracted.

| Relative path | Bytes | SHA-256 |
|---|---:|---|
| `AAA_g4beamline_study.html` | 13,096 | `7c40b61b9a1fbf6a9fb521ae393d16374e177aec4918b3966a50fd7d68fcb203` |
| `AAAREADME_g4beamline_grid.html` | 2,234 | `a3d0dc422e756b3efe2371dc6e342ea7bcaae7b96fc1bc547b57bb2e70544b59` |
| `AAAREADME_g4beamline.html` | 4,259 | `e0ef87910e21953e6161629a2e5f2cf8dc1170626b488f992dba044e8a32369b` |
| `AAAREADME_mad2g4bl.html` | 10,020 | `4bff42ef208639a7e1e472286087d3f2b66ea5a715d2f31dad608496798e2b6b` |
| `g4beamline_study/index.html` | 2,599 | `c4591b2956f99e25bc6251ffe90220f6da8243481a39c777787de3c267759b6a` |
| `g4beamline_study/ExtrBeamTransToCMAG-D.dat` | 49,465,532 | `cf7c836f923dd814f87ca7a6645edff3a02c0c9c05bcf021186002ed9e6b8149` |
| `g4beamline_study/exttracks_rot_1800000.root` | 45,192,824 | `7df4c14fcb651c9259293b33ca9e5ae645d9ad8a6b61a13095af64c36f055edb` |
| `g4beamline_study/exttracks_rot_e30_e40_1800000.root` | 13,577,449 | `583286959bb1cdd94df638a62c52c10412e8e3f6a852d0cc9bb4a2fbf130cf74` |
| `g4beamline_study/g4beamline_study_20170103.tgz` | 79,158 | `7d6c5e6ee0c731be57e7f3b23ba6c164dcd0c5f3c4adf468a37a97e24f784350` |
| `g4beamline_study/g4beamline_study_export.tgz` | 95,474,202 | `7e5e4a30e5a85688a44154ab70fc5ded54f2ddfd8c53b2e62b3af3b09a8a7137` |
| `g4beamline_study/mu2e_downstream.root` | 17,985,327 | `7e5cefcb2fad3f2eb615dbfa4e9d47a804a717ee541f6d2eedef96cd748ffc0a` |
| `g4beamline_study/mu2e_upstream.root` | 18,010,516 | `3853c89068ffd225a358b52309d87f812cf8d8034ecddd8fade34ff75f848968` |
| `g4beamline_study/mu2e_upstream_gaussian.root` | 16,404,038 | `50404b83ff69620be6bfdb06ab0ba8f4d4e403bddb38714461a1efd5f73c4196` |
| `mad2g4bl.py` | 11,362 | `5b457d97606ec86f7788773e5414b2082dcb45f5e69df9adf3fd2df13ae8055a` |

Both `.tgz` files pass `gzip -t` and complete `tar -tzf` listings. The 2017
archive has 41 members. The 2026 export has 560 members; GNU tar emits
non-fatal warnings for macOS/libarchive extended-header keywords.

## Current extraction state

The table above records the files as originally downloaded on 2026-08-15. As
of the 2026-08-17 folder comparison, the two `.tgz` files are no longer present
at those paths. Their extracted payloads are present as
`g4beamline_study/g4beamline_study_20170103/g4beamline_study/` and
`g4beamline_study/g4beamline_study_export/g4beamline_study_export/`. No claim is
made here about when or by whom the extraction and archive removal occurred.
