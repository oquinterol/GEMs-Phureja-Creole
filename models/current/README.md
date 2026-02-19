# Current Production Model: CreolePotato.xml

This directory contains the current production version of the genome-scale metabolic model, used for simulations and analysis.

## Version: 2.0

This version is a publication-ready release that consolidates all curation work from the v1.x series.

### Model Summary

| Property | Value |
|----------|-------|
| Model ID | `CreolePotato_v2_0` |
| Organism | *Solanum tuberosum* Group Phureja (Cultivar 'Criolla Colombia') |
| Reactions | 1,063 |
| Metabolites | 1,078 |
| Genes | 2,048 |
| SBML Level | 3 Version 1 + FBC v2 |
| FBA Status | Optimal (objective: 361.32) |

### Changes from v1.9

- Updated SBML metadata (id, metaid, name, RDF annotations) to v2.0
- Renamed to version-agnostic filename `CreolePotato.xml` (standard-GEM convention)
- Version tracked in `version.txt` and SBML metadata, not filename

### History

The v1.x series is preserved in `models/curated/` and as GitHub releases (deprecated). Key milestones:

- **v1.0**: Initial curated model from KBase draft
- **v1.1-v1.4**: Mass balance corrections
- **v1.5**: GPR expansion
- **v1.6-v1.7**: Gene annotations and MIRIAM enrichment
- **v1.8**: eggNOG-mapper MIRIAM annotations
- **v1.9**: Exchange reaction directionality fixes, consistency curation
- **v2.0**: Metadata standardization, publication release
