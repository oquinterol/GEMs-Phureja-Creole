# Curated Model History

This directory contains previous curated versions of the metabolic model. Each file represents a specific state in the development and refinement process. The current production model is in `models/current/`.

## Versions

- **creole_v1.0_final.xml**: First stable base model version.

- **creole_v1.1_balanced.xml**: First pass mass balance on reactions.

- **creole_v1.2_annotated.xml**: Gene and metabolite annotation improvements.

- **creole_v1.4_balanced.xml**: Additional mass balance iteration, fixing stoichiometric inconsistencies.

- **creole_v1.5_gprs_improved.xml**: Expansion of gene-protein-reaction (GPR) associations.

- **creole_v1.6_gene_annotated.xml**: Gene annotations from external databases.

- **creole_v1.6_mapped_genes.xml**: Gene-to-PGSC identifier mapping variant.

- **creole_v1.7_miriam_annotated.xml**: MIRIAM/Identifiers.org standard annotations.

- **creole_v1.8_eggnog_miriam.xml**: eggNOG-mapper MIRIAM annotations (162K functional annotations).

## Removed Files

The following files were removed as duplicates in v2.0.0:

- `creole_v1.0.xml` — exact duplicate of `creole_v1.0_final.xml`
- `creole_with_gprs.xml` — pre-versioned artifact, superseded by v1.5
- `creole_with_gprs_annot.xml` — pre-versioned artifact, superseded by v1.2
