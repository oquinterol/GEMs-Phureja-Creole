#!/usr/bin/env python3
"""
Filter proteome for STRING analysis based on eggNOG-mapper annotations.

This script:
1. Parses eggNOG-mapper annotations to identify metabolically relevant proteins
2. Selects one representative isoform per gene (longest sequence)
3. Creates filtered FASTA files in batches suitable for STRING upload (max 2000 proteins)

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


def parse_emapper_annotations(emapper_file: Path) -> Dict[str, Dict]:
    """
    Parse eggNOG-mapper annotations file.

    Returns:
        Dict mapping gene_id -> annotation data
    """
    annotations = {}

    with open(emapper_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 21:
                continue

            gene_id = fields[0]

            # Extract relevant annotation fields
            annotations[gene_id] = {
                'description': fields[7],
                'preferred_name': fields[8],
                'GOs': fields[9],
                'EC': fields[10],
                'KEGG_ko': fields[11],
                'KEGG_Pathway': fields[12],
                'KEGG_Reaction': fields[14],
                'BiGG_Reaction': fields[19],
                'score': float(fields[3]) if fields[3] != '-' else 0.0
            }

    return annotations


def parse_fasta(fasta_file: Path) -> Dict[str, Tuple[str, str]]:
    """
    Parse FASTA file.

    Returns:
        Dict mapping sequence_id -> (header, sequence)
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = (current_header, ''.join(current_seq))

                # Start new sequence
                current_header = line[1:]  # Remove '>'
                # Extract gene ID (e.g., "g1.t1 gene=g1" -> "g1.t1")
                current_id = current_header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_id:
            sequences[current_id] = (current_header, ''.join(current_seq))

    return sequences


def extract_gene_id(transcript_id: str) -> str:
    """
    Extract gene ID from transcript ID.
    E.g., "g1.t1" -> "g1"
    """
    match = re.match(r'(g\d+)\.t\d+', transcript_id)
    if match:
        return match.group(1)
    return transcript_id


def filter_by_metabolic_relevance(
    annotations: Dict[str, Dict],
    min_score: float = 50.0,
    require_ec: bool = False,
    require_kegg: bool = False
) -> Set[str]:
    """
    Filter transcripts by metabolic relevance criteria.

    Args:
        annotations: Transcript annotations from eggNOG-mapper (keyed by transcript ID)
        min_score: Minimum annotation score
        require_ec: Require EC number
        require_kegg: Require KEGG annotation (pathway or reaction)

    Returns:
        Set of transcript IDs passing filters
    """
    filtered = set()

    for transcript_id, annot in annotations.items():
        # Check score threshold
        if annot['score'] < min_score:
            continue

        # Check EC requirement
        if require_ec and (not annot['EC'] or annot['EC'] == '-'):
            continue

        # Check KEGG requirement
        if require_kegg:
            has_kegg = (
                (annot['KEGG_Pathway'] and annot['KEGG_Pathway'] != '-') or
                (annot['KEGG_Reaction'] and annot['KEGG_Reaction'] != '-') or
                (annot['KEGG_ko'] and annot['KEGG_ko'] != '-')
            )
            if not has_kegg:
                continue

        filtered.add(transcript_id)

    return filtered


def select_representative_isoforms(
    sequences: Dict[str, Tuple[str, str]],
    filtered_transcripts: Set[str]
) -> Dict[str, Tuple[str, str, str]]:
    """
    Select the longest isoform for each gene.

    Args:
        sequences: Dict of transcript_id -> (header, sequence)
        filtered_transcripts: Set of transcript IDs that passed filtering

    Returns:
        Dict mapping gene_id -> (transcript_id, header, sequence)
    """
    gene_isoforms = defaultdict(list)

    # Group isoforms by gene (only for filtered transcripts)
    for transcript_id, (header, seq) in sequences.items():
        # Only consider transcripts that passed filters
        if transcript_id not in filtered_transcripts:
            continue

        gene_id = extract_gene_id(transcript_id)
        gene_isoforms[gene_id].append((transcript_id, header, seq))

    # Select longest isoform per gene
    representatives = {}
    for gene_id, isoforms in gene_isoforms.items():
        # Sort by sequence length (descending)
        isoforms_sorted = sorted(isoforms, key=lambda x: len(x[2]), reverse=True)
        transcript_id, header, seq = isoforms_sorted[0]
        representatives[gene_id] = (transcript_id, header, seq)

    return representatives


def write_fasta_batches(
    sequences: Dict[str, Tuple[str, str, str]],
    output_dir: Path,
    batch_size: int = 2000,
    prefix: str = "creole_metabolic"
):
    """
    Write sequences to FASTA files in batches.

    Args:
        sequences: Dict of gene_id -> (transcript_id, header, sequence)
        output_dir: Output directory for batch files
        batch_size: Maximum sequences per file
        prefix: Prefix for output filenames
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Sort genes alphabetically for consistency
    sorted_genes = sorted(sequences.keys())

    batch_num = 1
    batch_count = 0
    batch_file = None

    for gene_id in sorted_genes:
        transcript_id, header, seq = sequences[gene_id]

        # Open new batch file if needed
        if batch_count == 0:
            if batch_file:
                batch_file.close()

            batch_filename = output_dir / f"{prefix}_batch{batch_num:02d}.faa"
            batch_file = open(batch_filename, 'w')
            print(f"Creating {batch_filename.name}...")

        # Write sequence
        batch_file.write(f">{header}\n")
        # Write sequence in lines of 70 characters
        for i in range(0, len(seq), 70):
            batch_file.write(f"{seq[i:i+70]}\n")

        batch_count += 1

        # Check if batch is full
        if batch_count >= batch_size:
            batch_file.close()
            print(f"  Wrote {batch_count} sequences to batch {batch_num}")
            batch_num += 1
            batch_count = 0
            batch_file = None

    # Close last batch file
    if batch_file:
        batch_file.close()
        print(f"  Wrote {batch_count} sequences to batch {batch_num}")

    print(f"\nTotal batches created: {batch_num}")


def write_summary_report(
    output_file: Path,
    total_sequences: int,
    total_genes: int,
    filtered_genes: int,
    representative_count: int,
    annotations: Dict[str, Dict],
    representatives: Dict[str, Tuple[str, str, str]]
):
    """Write summary report of filtering process."""

    if representative_count == 0:
        with open(output_file, 'w') as f:
            f.write("# Proteome Filtering Summary for STRING Analysis\n\n")
            f.write(f"Total protein sequences in original file: {total_sequences:,}\n")
            f.write(f"Total unique genes: {total_genes:,}\n")
            f.write(f"Transcripts passing metabolic filters: {filtered_genes:,}\n")
            f.write(f"Representative isoforms selected: {representative_count:,}\n\n")
            f.write("WARNING: No proteins passed the filtering criteria.\n")
            f.write("Consider relaxing the filtering parameters.\n")
        print(f"\nSummary report written to: {output_file}")
        return

    # Count annotations in filtered set (use transcript_id from representatives)
    ec_count = 0
    kegg_pathway_count = 0
    kegg_reaction_count = 0
    go_count = 0

    for gene_id, (transcript_id, header, seq) in representatives.items():
        if transcript_id in annotations:
            annot = annotations[transcript_id]
            if annot['EC'] and annot['EC'] != '-':
                ec_count += 1
            if annot['KEGG_Pathway'] and annot['KEGG_Pathway'] != '-':
                kegg_pathway_count += 1
            if annot['KEGG_Reaction'] and annot['KEGG_Reaction'] != '-':
                kegg_reaction_count += 1
            if annot['GOs'] and annot['GOs'] != '-':
                go_count += 1

    with open(output_file, 'w') as f:
        f.write("# Proteome Filtering Summary for STRING Analysis\n\n")
        f.write(f"Total protein sequences in original file: {total_sequences:,}\n")
        f.write(f"Total unique genes: {total_genes:,}\n")
        f.write(f"Transcripts passing metabolic filters: {filtered_genes:,}\n")
        f.write(f"Representative isoforms selected: {representative_count:,}\n\n")

        f.write("## Annotation Coverage in Filtered Set\n\n")
        f.write(f"Genes with EC numbers: {ec_count:,} ({100*ec_count/representative_count:.1f}%)\n")
        f.write(f"Genes with KEGG pathways: {kegg_pathway_count:,} ({100*kegg_pathway_count/representative_count:.1f}%)\n")
        f.write(f"Genes with KEGG reactions: {kegg_reaction_count:,} ({100*kegg_reaction_count/representative_count:.1f}%)\n")
        f.write(f"Genes with GO terms: {go_count:,} ({100*go_count/representative_count:.1f}%)\n\n")

        f.write("## Filtering Strategy\n\n")
        f.write("- Selected longest isoform per gene\n")
        f.write("- Filtered for metabolic relevance (EC numbers or KEGG annotations)\n")
        f.write("- Minimum annotation score: 50.0\n")
        f.write("- Batched into files of max 2000 sequences for STRING upload\n")

    print(f"\nSummary report written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter proteome for STRING analysis based on metabolic relevance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Filter for proteins with EC numbers or KEGG annotations
  python filter_proteome_for_string.py \\
    --proteome data/Proteome/Predicted_Proteins.faa \\
    --annotations data/annotations/emapper_creole.emapper.annotations \\
    --out data/string/batches

  # Strict filter: require EC numbers
  python filter_proteome_for_string.py \\
    --proteome data/Proteome/Predicted_Proteins.faa \\
    --annotations data/annotations/emapper_creole.emapper.annotations \\
    --out data/string/batches \\
    --require-ec
        """
    )

    parser.add_argument(
        '--proteome',
        type=Path,
        required=True,
        help='Input proteome FASTA file'
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
        help='Output directory for filtered FASTA batches'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=2000,
        help='Maximum sequences per batch file (default: 2000)'
    )
    parser.add_argument(
        '--min-score',
        type=float,
        default=50.0,
        help='Minimum annotation score (default: 50.0)'
    )
    parser.add_argument(
        '--require-ec',
        action='store_true',
        help='Require EC number annotation'
    )
    parser.add_argument(
        '--require-kegg',
        action='store_true',
        help='Require KEGG annotation (pathway/reaction/KO)'
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default='creole_metabolic',
        help='Prefix for output batch files (default: creole_metabolic)'
    )

    args = parser.parse_args()

    print("="*70)
    print("STRING Proteome Filtering Pipeline")
    print("="*70)

    # Parse annotations
    print("\n1. Parsing eggNOG-mapper annotations...")
    annotations = parse_emapper_annotations(args.annotations)
    print(f"   Loaded annotations for {len(annotations):,} genes")

    # Parse FASTA
    print("\n2. Parsing proteome FASTA file...")
    sequences = parse_fasta(args.proteome)
    print(f"   Loaded {len(sequences):,} protein sequences")

    # Count unique genes
    unique_genes = set(extract_gene_id(tid) for tid in sequences.keys())
    print(f"   Found {len(unique_genes):,} unique genes")

    # Filter by metabolic relevance
    print("\n3. Filtering by metabolic relevance...")
    filtered_transcripts = filter_by_metabolic_relevance(
        annotations,
        min_score=args.min_score,
        require_ec=args.require_ec,
        require_kegg=args.require_kegg
    )
    print(f"   {len(filtered_transcripts):,} transcripts passed filters")

    # Select representative isoforms
    print("\n4. Selecting representative isoforms (longest per gene)...")
    representatives = select_representative_isoforms(sequences, filtered_transcripts)
    print(f"   Selected {len(representatives):,} representative proteins")

    # Write batches
    print(f"\n5. Writing FASTA batches (max {args.batch_size} sequences/file)...")
    write_fasta_batches(
        representatives,
        args.out,
        batch_size=args.batch_size,
        prefix=args.prefix
    )

    # Write summary report
    summary_file = args.out / f"{args.prefix}_summary.txt"
    write_summary_report(
        summary_file,
        total_sequences=len(sequences),
        total_genes=len(unique_genes),
        filtered_genes=len(filtered_transcripts),
        representative_count=len(representatives),
        annotations=annotations,
        representatives=representatives
    )

    print("\n" + "="*70)
    print("Filtering complete!")
    print(f"Output files written to: {args.out}")
    print("="*70)


if __name__ == '__main__':
    main()
