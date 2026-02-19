#!/usr/bin/env python3
"""
Enrich FASTA headers with annotation information for better STRING mapping.

STRING recognizes headers with gene names, descriptions, and EC numbers.
This script creates improved FASTA files with enriched headers that include:
- Gene ID
- Preferred gene name (from eggNOG)
- Protein description
- EC numbers
- KEGG KO

Format: >geneID [preferred_name] description [EC:X.X.X.X] [KO:KXXXXX]

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
from pathlib import Path
from typing import Dict, Tuple


def parse_emapper_annotations(emapper_file: Path) -> Dict[str, Dict]:
    """Parse eggNOG-mapper annotations."""
    annotations = {}

    with open(emapper_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 21:
                continue

            gene_id = fields[0]
            annotations[gene_id] = {
                'preferred_name': fields[8] if fields[8] != '-' else '',
                'description': fields[7] if fields[7] != '-' else '',
                'EC': fields[10] if fields[10] != '-' else '',
                'KEGG_ko': fields[11] if fields[11] != '-' else '',
            }

    return annotations


def parse_fasta(fasta_file: Path) -> Dict[str, Tuple[str, str]]:
    """Parse FASTA file."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = (current_header, ''.join(current_seq))

                current_header = line[1:]
                current_id = current_header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            sequences[current_id] = (current_header, ''.join(current_seq))

    return sequences


def create_enriched_header(
    gene_id: str,
    annotations: Dict[str, Dict],
    organism: str = "Solanum tuberosum"
) -> str:
    """
    Create enriched FASTA header for STRING.

    Format examples that STRING recognizes:
    >gene_id gene_name description [EC:X.X.X.X]
    >gene_id description OS=Organism
    """
    if gene_id not in annotations:
        return gene_id

    annot = annotations[gene_id]

    # Build header components
    header_parts = [gene_id]

    # Add gene name if available
    if annot['preferred_name']:
        header_parts.append(f"GN={annot['preferred_name']}")

    # Add description
    if annot['description']:
        # Clean description (remove extra spaces)
        desc = ' '.join(annot['description'].split())
        header_parts.append(desc)

    # Add organism
    header_parts.append(f"OS={organism}")

    # Add EC number if available
    if annot['EC']:
        # EC numbers can be comma-separated in eggNOG
        ec_numbers = annot['EC'].split(',')
        for ec in ec_numbers[:3]:  # Limit to first 3 EC numbers
            header_parts.append(f"EC={ec.strip()}")

    # Add KEGG KO if available
    if annot['KEGG_ko']:
        # Extract KO number (format: ko:K12345)
        ko_ids = [ko.replace('ko:', 'KO:') for ko in annot['KEGG_ko'].split(',')]
        for ko in ko_ids[:2]:  # Limit to first 2 KOs
            header_parts.append(ko.strip())

    return ' '.join(header_parts)


def enrich_fasta_file(
    input_fasta: Path,
    output_fasta: Path,
    annotations: Dict[str, Dict],
    organism: str = "Solanum tuberosum"
):
    """Create enriched FASTA file with better headers."""

    sequences = parse_fasta(input_fasta)

    with open(output_fasta, 'w') as out:
        for gene_id, (old_header, seq) in sequences.items():
            # Create enriched header
            new_header = create_enriched_header(gene_id, annotations, organism)

            # Write enriched FASTA
            out.write(f">{new_header}\n")

            # Write sequence in lines of 70 characters
            for i in range(0, len(seq), 70):
                out.write(f"{seq[i:i+70]}\n")

    print(f"Created enriched FASTA: {output_fasta.name}")


def main():
    parser = argparse.ArgumentParser(
        description="Enrich FASTA headers with annotations for better STRING mapping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Enrich a single file
  python enrich_fasta_headers.py \\
    --input data/string/web_example/creole_web_example_batch01.faa \\
    --annotations data/annotations/emapper_creole.emapper.annotations \\
    --out data/string/web_enriched/creole_enriched_batch01.faa

  # Enrich all files in a directory
  for f in data/string/web_example/*.faa; do
    python enrich_fasta_headers.py \\
      --input "$f" \\
      --annotations data/annotations/emapper_creole.emapper.annotations \\
      --out data/string/web_enriched/$(basename "$f")
  done
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input FASTA file'
    )
    parser.add_argument(
        '--annotations',
        type=Path,
        required=True,
        help='eggNOG-mapper annotations file'
    )
    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Output enriched FASTA file'
    )
    parser.add_argument(
        '--organism',
        type=str,
        default='Solanum tuberosum',
        help='Organism name (default: Solanum tuberosum)'
    )

    args = parser.parse_args()

    # Create output directory if needed
    args.out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading annotations from {args.annotations}...")
    annotations = parse_emapper_annotations(args.annotations)
    print(f"  Loaded {len(annotations):,} annotations")

    print(f"\nEnriching FASTA file: {args.input.name}")
    enrich_fasta_file(args.input, args.out, annotations, args.organism)

    print("\n✅ Done!")


if __name__ == '__main__':
    main()
