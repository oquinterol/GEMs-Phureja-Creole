#!/usr/bin/env python3
"""
Select one representative isoform per gene from complete proteome.

This script selects the longest isoform for each gene to comply with
STRING requirement: "one protein sequence per gene/locus"

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
from pathlib import Path
from collections import defaultdict
import re


def extract_gene_id(transcript_id: str) -> str:
    """
    Extract gene ID from transcript ID.
    E.g., "g1.t1" -> "g1"
    """
    match = re.match(r'(g\d+)\.t\d+', transcript_id)
    if match:
        return match.group(1)
    return transcript_id


def parse_fasta(fasta_file: Path):
    """Parse FASTA file and return sequences grouped by gene."""

    gene_isoforms = defaultdict(list)

    current_id = None
    current_header = None
    current_seq = []

    print(f"Leyendo archivo: {fasta_file.name}")

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    gene_id = extract_gene_id(current_id)
                    seq_str = ''.join(current_seq)
                    gene_isoforms[gene_id].append({
                        'transcript_id': current_id,
                        'header': current_header,
                        'sequence': seq_str,
                        'length': len(seq_str)
                    })

                # Start new sequence
                current_header = line[1:]
                current_id = current_header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id:
            gene_id = extract_gene_id(current_id)
            seq_str = ''.join(current_seq)
            gene_isoforms[gene_id].append({
                'transcript_id': current_id,
                'header': current_header,
                'sequence': seq_str,
                'length': len(seq_str)
            })

    return gene_isoforms


def select_longest_isoforms(gene_isoforms: dict) -> list:
    """Select the longest isoform for each gene."""

    representatives = []

    for gene_id, isoforms in gene_isoforms.items():
        # Sort by length (descending)
        isoforms_sorted = sorted(isoforms, key=lambda x: x['length'], reverse=True)

        # Select longest
        longest = isoforms_sorted[0]
        representatives.append(longest)

    return representatives


def write_fasta(representatives: list, output_file: Path):
    """Write representative isoforms to FASTA file."""

    with open(output_file, 'w') as out:
        for rep in representatives:
            # Write header
            out.write(f">{rep['header']}\n")

            # Write sequence in lines of 70 characters
            seq = rep['sequence']
            for i in range(0, len(seq), 70):
                out.write(f"{seq[i:i+70]}\n")

    print(f"\n✅ Archivo creado: {output_file}")
    print(f"   Total de genes: {len(representatives):,}")


def main():
    parser = argparse.ArgumentParser(
        description="Select one representative isoform per gene (longest)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python select_one_isoform_per_gene.py \\
    --input data/Proteome/Predicted_Proteins_ENRICHED_FULL.faa \\
    --out data/string/proteome_1_per_gene_ENRICHED.faa
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input FASTA file (can have multiple isoforms per gene)'
    )
    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Output FASTA file (one isoform per gene)'
    )

    args = parser.parse_args()

    print("="*70)
    print("Selección de 1 Isoforma por Gen - Para STRING")
    print("="*70)

    # Parse FASTA
    print("\n1. Parseando archivo FASTA...")
    gene_isoforms = parse_fasta(args.input)

    total_isoforms = sum(len(isoforms) for isoforms in gene_isoforms.values())
    print(f"   Total de secuencias: {total_isoforms:,}")
    print(f"   Genes únicos: {len(gene_isoforms):,}")

    # Count genes with multiple isoforms
    multi_isoform = sum(1 for isoforms in gene_isoforms.values() if len(isoforms) > 1)
    print(f"   Genes con múltiples isoformas: {multi_isoform:,}")

    # Select longest isoforms
    print("\n2. Seleccionando isoforma más larga por gen...")
    representatives = select_longest_isoforms(gene_isoforms)

    print(f"   Seleccionadas: {len(representatives):,} proteínas")

    # Create output directory if needed
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # Write output
    print("\n3. Escribiendo archivo de salida...")
    write_fasta(representatives, args.out)

    print("\n" + "="*70)
    print("✅ Proceso completado!")
    print(f"   De {total_isoforms:,} secuencias → {len(representatives):,} genes únicos")
    print(f"   Reducción: {total_isoforms - len(representatives):,} isoformas alternativas eliminadas")
    print("="*70)


if __name__ == '__main__':
    main()
