#!/usr/bin/env python3
"""
Map generic gene IDs (g*.t*) to PGSC/UniProt IDs for STRING database compatibility

This script extracts ortholog information from eggNOG-mapper annotations
and creates a mapping table compatible with STRING database requirements.

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
import re
from pathlib import Path
from typing import Dict, Tuple, Optional
import sys


def parse_emapper_annotations(emapper_file: Path) -> Dict[str, Tuple[str, str]]:
    """
    Parse eggNOG-mapper annotations to extract PGSC/UniProt IDs

    Args:
        emapper_file: Path to .emapper.annotations file

    Returns:
        Dictionary mapping gene_id -> (ortholog_id, description)
    """
    mappings = {}
    pgsc_count = 0
    uniprot_count = 0
    no_mapping_count = 0

    print(f"📖 Reading eggNOG annotations from: {emapper_file}")

    with open(emapper_file, 'r') as f:
        for line in f:
            # Skip comments and header
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            query_id = fields[0]  # g*.t* ID
            seed_ortholog = fields[1]  # Ortholog ID
            description = fields[7] if len(fields) > 7 else "No description"

            # Skip if no ortholog
            if seed_ortholog == '-' or not seed_ortholog:
                no_mapping_count += 1
                continue

            # Extract taxonomy ID and gene ID from seed_ortholog
            # Format: 4113.PGSC0003DMT400058594 or 4081.Solyc01g020440.2.1
            match = re.match(r'(\d+)\.(.*)', seed_ortholog)
            if not match:
                no_mapping_count += 1
                continue

            tax_id = match.group(1)
            gene_id = match.group(2)

            # Prioritize S. tuberosum (4113) PGSC IDs
            if tax_id == '4113':
                if 'PGSC' in gene_id:
                    pgsc_count += 1
                    mappings[query_id] = (f"4113.{gene_id}", description)
                else:
                    # Other S. tuberosum IDs (UniProt, etc.)
                    uniprot_count += 1
                    mappings[query_id] = (f"4113.{gene_id}", description)
            elif tax_id == '4081':
                # S. lycopersicum - only if no S. tuberosum mapping exists
                if query_id not in mappings:
                    mappings[query_id] = (f"4081.{gene_id}", description)
            else:
                # Other species - lowest priority
                if query_id not in mappings:
                    mappings[query_id] = (f"{tax_id}.{gene_id}", description)

    print(f"✅ Parsed {len(mappings)} mappings:")
    print(f"   - PGSC IDs (4113): {pgsc_count}")
    print(f"   - UniProt/Other (4113): {uniprot_count}")
    print(f"   - No mapping: {no_mapping_count}")

    return mappings


def write_mapping_table(mappings: Dict[str, Tuple[str, str]], output_file: Path, strip_taxon: bool = True):
    """Write mapping table in TSV format"""

    print(f"💾 Writing mapping table to: {output_file}")
    print(f"   Strip taxon prefix: {strip_taxon}")

    with open(output_file, 'w') as f:
        f.write("original_id\tstring_id\tdescription\n")

        for orig_id, (string_id, description) in sorted(mappings.items()):
            # Optionally strip taxon prefix
            if strip_taxon and '.' in string_id:
                clean_id = string_id.split('.', 1)[1]
            else:
                clean_id = string_id

            f.write(f"{orig_id}\t{clean_id}\t{description}\n")

    print(f"✅ Wrote {len(mappings)} mappings")


def create_string_proteome(
    original_fasta: Path,
    mappings: Dict[str, Tuple[str, str]],
    output_fasta: Path,
    strip_taxon: bool = True
):
    """
    Create a new proteome FASTA with STRING-compatible IDs

    Args:
        original_fasta: Original proteome FASTA file
        mappings: ID mappings from parse_emapper_annotations
        output_fasta: Output FASTA with new IDs
        strip_taxon: If True, remove taxon prefix (e.g., '4113.') from IDs
    """
    print(f"🧬 Creating STRING-compatible proteome FASTA")
    print(f"   Input: {original_fasta}")
    print(f"   Output: {output_fasta}")
    print(f"   Strip taxon prefix: {strip_taxon}")

    mapped_count = 0
    unmapped_count = 0

    with open(original_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
        current_id = None

        for line in fin:
            if line.startswith('>'):
                # Extract ID from header
                header = line[1:].strip()
                gene_id = header.split()[0]

                # Check if we have a mapping
                if gene_id in mappings:
                    string_id, description = mappings[gene_id]

                    # Optionally strip taxon prefix for cleaner IDs
                    if strip_taxon and '.' in string_id:
                        # Remove '4113.' prefix, keep only the actual gene ID
                        clean_id = string_id.split('.', 1)[1]
                    else:
                        clean_id = string_id

                    # Write new header with STRING-compatible ID
                    fout.write(f">{clean_id} {description}\n")
                    current_id = clean_id
                    mapped_count += 1
                else:
                    # Keep original ID if no mapping
                    fout.write(line)
                    current_id = None
                    unmapped_count += 1
            else:
                # Write sequence line as-is
                fout.write(line)

    print(f"✅ Created STRING proteome:")
    print(f"   - Mapped sequences: {mapped_count}")
    print(f"   - Unmapped sequences: {unmapped_count}")
    print(f"   - Total: {mapped_count + unmapped_count}")


def generate_statistics(mappings: Dict[str, Tuple[str, str]]) -> Dict[str, int]:
    """Generate mapping statistics by taxonomy ID"""

    stats = {}

    for _, (string_id, _) in mappings.items():
        tax_id = string_id.split('.')[0]
        stats[tax_id] = stats.get(tax_id, 0) + 1

    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Map gene IDs to PGSC/UniProt IDs for STRING compatibility',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate mapping table only
  %(prog)s --emapper data/annotations/emapper_creole.emapper.annotations \\
           --out data/mappings/gene_to_string_ids.tsv

  # Generate mapping table + STRING proteome
  %(prog)s --emapper data/annotations/emapper_creole.emapper.annotations \\
           --proteome data/Proteome/Predicted_Proteins.faa \\
           --out data/mappings/gene_to_string_ids.tsv \\
           --out-proteome data/Proteome/Predicted_Proteins_STRING.faa
        """
    )

    parser.add_argument(
        '--emapper',
        type=Path,
        required=True,
        help='Input eggNOG-mapper annotations file'
    )

    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Output mapping table (TSV format)'
    )

    parser.add_argument(
        '--proteome',
        type=Path,
        help='Original proteome FASTA file (optional)'
    )

    parser.add_argument(
        '--out-proteome',
        type=Path,
        help='Output STRING-compatible proteome FASTA (requires --proteome)'
    )

    parser.add_argument(
        '--keep-taxon',
        action='store_true',
        help='Keep taxon prefix in FASTA IDs (e.g., 4113.PGSC...). Default is to strip it.'
    )

    parser.add_argument(
        '--stats',
        type=Path,
        help='Output statistics file (optional)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not args.emapper.exists():
        print(f"❌ Error: eggNOG file not found: {args.emapper}", file=sys.stderr)
        sys.exit(1)

    if args.out_proteome and not args.proteome:
        print("❌ Error: --out-proteome requires --proteome", file=sys.stderr)
        sys.exit(1)

    if args.proteome and not args.proteome.exists():
        print(f"❌ Error: Proteome file not found: {args.proteome}", file=sys.stderr)
        sys.exit(1)

    # Parse eggNOG annotations
    mappings = parse_emapper_annotations(args.emapper)

    if not mappings:
        print("❌ Error: No mappings found in eggNOG file", file=sys.stderr)
        sys.exit(1)

    # Write mapping table
    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Use same strip_taxon setting for both table and FASTA
    strip_taxon = not args.keep_taxon
    write_mapping_table(mappings, args.out, strip_taxon)

    # Create STRING proteome if requested
    if args.proteome and args.out_proteome:
        args.out_proteome.parent.mkdir(parents=True, exist_ok=True)
        # By default strip taxon (STRING adds it automatically)
        # Use --keep-taxon to preserve it
        strip_taxon = not args.keep_taxon
        create_string_proteome(args.proteome, mappings, args.out_proteome, strip_taxon)

    # Generate statistics
    stats = generate_statistics(mappings)

    print("\n📊 Mapping Statistics by Taxonomy ID:")
    for tax_id, count in sorted(stats.items(), key=lambda x: x[1], reverse=True):
        tax_names = {
            '4113': 'Solanum tuberosum (potato)',
            '4081': 'Solanum lycopersicum (tomato)',
            '33090': 'Viridiplantae (green plants)',
        }
        tax_name = tax_names.get(tax_id, f"Tax ID {tax_id}")
        print(f"   {tax_id}: {count:6d} genes  ({tax_name})")

    if args.stats:
        args.stats.parent.mkdir(parents=True, exist_ok=True)
        with open(args.stats, 'w') as f:
            f.write("taxonomy_id\ttaxonomy_name\tgene_count\n")
            for tax_id, count in sorted(stats.items(), key=lambda x: x[1], reverse=True):
                tax_names = {
                    '4113': 'Solanum tuberosum',
                    '4081': 'Solanum lycopersicum',
                    '33090': 'Viridiplantae',
                }
                tax_name = tax_names.get(tax_id, f"Tax_{tax_id}")
                f.write(f"{tax_id}\t{tax_name}\t{count}\n")
        print(f"\n📈 Statistics saved to: {args.stats}")

    print("\n✅ Done!")


if __name__ == '__main__':
    main()
