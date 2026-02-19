#!/usr/bin/env python3
"""
Split FASTA file into batches of specified size.

Author: GEM_Creole Team
Date: 2025-10-19
"""

import argparse
from pathlib import Path


def split_fasta_into_batches(input_file: Path, output_dir: Path, batch_size: int, prefix: str):
    """Split FASTA file into batches."""

    output_dir.mkdir(parents=True, exist_ok=True)

    batch_num = 1
    seq_count = 0
    current_batch = []

    print(f"Leyendo archivo: {input_file}")

    with open(input_file, 'r') as f:
        current_seq = []

        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_seq:
                    current_batch.extend(current_seq)
                    seq_count += 1

                    # Check if batch is full
                    if seq_count >= batch_size:
                        # Write batch
                        output_file = output_dir / f"{prefix}_batch{batch_num:02d}.faa"
                        with open(output_file, 'w') as out:
                            out.writelines(current_batch)

                        print(f"✅ Creado: {output_file.name} ({seq_count:,} secuencias)")

                        # Reset for next batch
                        batch_num += 1
                        seq_count = 0
                        current_batch = []

                # Start new sequence
                current_seq = [line]
            else:
                current_seq.append(line)

        # Save last sequence
        if current_seq:
            current_batch.extend(current_seq)
            seq_count += 1

    # Write remaining sequences
    if current_batch:
        output_file = output_dir / f"{prefix}_batch{batch_num:02d}.faa"
        with open(output_file, 'w') as out:
            out.writelines(current_batch)

        print(f"✅ Creado: {output_file.name} ({seq_count:,} secuencias)")

    print(f"\n✅ Total de batches creados: {batch_num}")


def main():
    parser = argparse.ArgumentParser(
        description="Split FASTA into batches",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python split_fasta_batches.py \\
    --input data/string/proteome_1_per_gene_ENRICHED.faa \\
    --out data/string/batches_100k \\
    --batch-size 100000 \\
    --prefix proteome_100k
        """
    )

    parser.add_argument('--input', type=Path, required=True, help='Input FASTA file')
    parser.add_argument('--out', type=Path, required=True, help='Output directory')
    parser.add_argument('--batch-size', type=int, default=100000, help='Sequences per batch')
    parser.add_argument('--prefix', type=str, default='batch', help='Output file prefix')

    args = parser.parse_args()

    print("="*70)
    print("División de FASTA en Batches para STRING")
    print("="*70)
    print(f"\nArchivo de entrada: {args.input}")
    print(f"Tamaño de batch: {args.batch_size:,} secuencias")
    print(f"Directorio de salida: {args.out}")
    print()

    split_fasta_into_batches(args.input, args.out, args.batch_size, args.prefix)

    print("\n" + "="*70)
    print("✅ Proceso completado!")
    print("="*70)


if __name__ == '__main__':
    main()
