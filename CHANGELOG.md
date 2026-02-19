# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-02-19

### Changed
- Updated SBML model metadata: id, metaid, name, RDF annotations to v2.0
- Renamed production model to `CreolePotato.xml` (version-agnostic, standard-GEM convention)
- Modernized CI/CD: updated GitHub Actions to latest versions, added PR validation workflow
- Replaced deprecated `actions/upload-release-asset@v1` with `gh release upload`

### Added
- `CHANGELOG.md` (this file)
- `CITATION.cff` for standardized citation metadata
- `LICENSE.md` (CC-BY-4.0)
- `version.txt` for machine-readable version tracking
- `scripts/curation/update_model_metadata.py` for reproducible metadata updates
- `scripts/string/` - STRING DB protein interaction analysis scripts (7 scripts)
- `.github/workflows/pr-validation.yml` for automated PR checks
- FBA smoke test step in release workflow

### Removed
- `models/curated/creole_v1.0.xml` (duplicate of `creole_v1.0_final.xml`)
- `models/curated/creole_with_gprs.xml` (pre-versioned artifact)
- `models/curated/creole_with_gprs_annot.xml` (pre-versioned artifact)

### Fixed
- `.gitignore`: removed blanket `*.md` rule that blocked documentation files
- SBML metadata version mismatch (file was named v1.9 but metadata said v1.3)

## v1.x Series (Deprecated)

> All v1.x releases point to the same commit and are superseded by v2.0.0.
> They are preserved for historical reference.

### [v1.9] - 2025-09-13
- Corrected exchange reaction directionality
- Simplified gene product labels (removed KEGG/GO from labels)
- Curated species and reactions for consistency

### [v1.8] - 2025-09-13
- Added eggNOG-mapper MIRIAM annotations (162K functional annotations)
- EC, KEGG, GO, UniProt cross-references from eggNOG

### [v1.7] - 2025-09-13
- Enriched with MIRIAM/Identifiers.org standard annotations
- Improved interoperability with external tools

### [v1.6] - 2025-09-13
- Gene annotations from external databases
- Gene-to-PGSC identifier mapping

### [v1.5] - 2025-09-13
- Expanded gene-protein-reaction (GPR) associations
- Improved GPR coverage

### [v1.4] - 2025-09-13
- Additional mass balance corrections
- Fixed stoichiometric inconsistencies

### [v1.2] - 2025-09-13
- Gene and metabolite annotation improvements

### [v1.1] - 2025-09-13
- First pass mass balance on reactions

### [v1.0] - 2025-09-13
- Initial curated model from KBase draft reconstruction
- Base model for *Solanum tuberosum* Group Phureja

[2.0.0]: https://github.com/Ampiria/GEM_Creole/releases/tag/v2.0.0
[v1.9]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.9
[v1.8]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.8
[v1.7]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.7
[v1.6]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.6
[v1.5]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.5
[v1.4]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.4
[v1.2]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.2
[v1.1]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.1
[v1.0]: https://github.com/Ampiria/GEM_Creole/releases/tag/v1.0
