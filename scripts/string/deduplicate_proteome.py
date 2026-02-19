#!/usr/bin/env python3
"""
Remove duplicate protein sequences from FASTA, keeping longest isoform per ID.

For STRING database submission, we need ONE protein per gene/locus.
This script keeps the longest sequence when multiple sequences share the same ID.

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
from pathlib import Path
from typing import Dict, Tuple
import sys


def parse_fasta_with_duplicates(fasta_file: Path) -> Dict[str, Tuple[str, str]]:
    """
    Parse FASTA and keep longest sequence per ID.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        Dict mapping gene_id -> (full_header, sequence)
        Only keeps the longest sequence for each ID
    """
    sequences = {}
    current_id = None
    current_header = None
    current_seq = []

    duplicates_found = 0
    total_sequences = 0

    print(f"📖 Reading FASTA file: {fasta_file}")

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    total_sequences += 1
                    seq = ''.join(current_seq)

                    # Check if ID already exists
                    if current_id in sequences:
                        duplicates_found += 1
                        # Keep longer sequence
                        existing_seq = sequences[current_id][1]
                        if len(seq) > len(existing_seq):
                            sequences[current_id] = (current_header, seq)
                    else:
                        sequences[current_id] = (current_header, seq)

                # Start new sequence
                current_header = line[1:].strip()
                # Extract ID (first word before space or pipe)
                current_id = current_header.split()[0].split('|')[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

        # Save last sequence
        if current_id is not None:
            total_sequences += 1
            seq = ''.join(current_seq)

            if current_id in sequences:
                duplicates_found += 1
                existing_seq = sequences[current_id][1]
                if len(seq) > len(existing_seq):
                    sequences[current_id] = (current_header, seq)
            else:
                sequences[current_id] = (current_header, seq)

    print(f"✅ Processed {total_sequences:,} sequences")
    print(f"   - Unique IDs: {len(sequences):,}")
    print(f"   - Duplicates found: {duplicates_found:,} ({duplicates_found*100/total_sequences:.1f}%)")

    return sequences


def write_deduplicated_fasta(sequences: Dict[str, Tuple[str, str]], output_file: Path):
    """Write deduplicated sequences to FASTA"""

    print(f"\n💾 Writing deduplicated FASTA to: {output_file}")

    with open(output_file, 'w') as f:
        for gene_id in sorted(sequences.keys()):
            header, seq = sequences[gene_id]
            f.write(f">{header}\n")

            # Write sequence in 80-character lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    print(f"✅ Wrote {len(sequences):,} unique sequences")


def main():
    parser = argparse.ArgumentParser(
        description='Deduplicate FASTA file for STRING submission (one protein per gene)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Why deduplicate?
  STRING expects ONE protein sequence per gene/locus. If you submit multiple
  isoforms or copies per gene, orthology predictions will be affected and
  confidence scores may be incorrect.

Deduplication strategy:
  When multiple sequences share the same ID, keep the LONGEST sequence.
  This is the standard approach for proteome databases.

Example:
  python deduplicate_proteome.py \\
    --input data/Proteome/Predicted_Proteins_STRING.faa \\
    --out data/Proteome/Predicted_Proteins_STRING_unique.faa
        """
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input FASTA file (may contain duplicates)'
    )

    parser.add_argument(
        '--out',
        type=Path,
        required=True,
        help='Output deduplicated FASTA file'
    )

    parser.add_argument(
        '--stats',
        type=Path,
        help='Output statistics file (optional)'
    )

    args = parser.parse_args()

    # Validate input
    if not args.input.exists():
        print(f"❌ Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print("="*70)
    print("STRING Proteome Deduplication")
    print("="*70)
    print()

    # Parse and deduplicate
    sequences = parse_fasta_with_duplicates(args.input)

    if not sequences:
        print("❌ Error: No sequences found", file=sys.stderr)
        sys.exit(1)

    # Write output
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_deduplicated_fasta(sequences, args.out)

    # Write statistics
    if args.stats:
        args.stats.parent.mkdir(parents=True, exist_ok=True)
        with open(args.stats, 'w') as f:
            f.write("metric\tvalue\n")
            f.write(f"total_sequences_input\t{sum(1 for _ in open(args.input) if _.startswith('>'))}\n")
            f.write(f"unique_sequences_output\t{len(sequences)}\n")
            f.write(f"duplicates_removed\t{sum(1 for _ in open(args.input) if _.startswith('>')) - len(sequences)}\n")
        print(f"\n📊 Statistics saved to: {args.stats}")

    print("\n" + "="*70)
    print("✅ Deduplication completed!")
    print("="*70)
    print(f"\n📌 Summary:")
    print(f"   Input:  {sum(1 for _ in open(args.input) if _.startswith('>')):,} sequences")
    print(f"   Output: {len(sequences):,} unique sequences")
    print(f"   Ready for STRING submission!")


if __name__ == '__main__':
    main()
