#!/usr/bin/env python3
"""Download genome .fna.gz files and extract 150bp subsequences with metadata.

Usage:
    python3 make_down_list_and_extract.py --metadata metadata.tsv --sample-size 10
"""

import os
import pandas as pd
from Bio import SeqIO
import gzip
import urllib.request
import argparse
import random


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True,
                        help="Path to the metadata TSV file (uncompressed)")
    parser.add_argument("--sample-size", type=int, default=10,
                        help="Number of genomes to randomly sample and extract")
    args = parser.parse_args()

    # Step 1: Read metadata
    df = pd.read_csv(args.metadata, sep='\t')

    required_cols = ['assembly_accession', 'asm_name']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f'Missing required column: {col}')

    # Step 2: Sample N random genomes
    sample_df = df.sample(n=min(args.sample_size, len(df)), random_state=42)

    # Step 3: Process each genome
    all_data = []
    for _, row in sample_df.iterrows():
        asb, asm, unique_org_name = row['assembly_accession'], row['asm_name'], row['unique_name']
        url = get_fna_path(asb, asm)
        fna_filename = f"{asb}_{asm.replace(' ', '_')}_genomic.fna.gz"

        print(f"Downloading: {url}")
        try:
            urllib.request.urlretrieve(url, fna_filename)
            df = extract_subsequences_from_fna(fna_filename, unique_org_name)
            all_data.append(df)
        except Exception as e:
            print(f"Failed to download or process {url}: {e}")

    # Step 4: Save combined results
    if all_data:
        final_df = pd.concat(all_data, ignore_index=True)
        final_df.to_parquet("subsequences_dataset.parquet", index=False)
        print("✅ Saved to subsequences_dataset.parquet")


def get_fna_path(asb, asm):
    asm = asm.replace(' ', '_').replace('#', '_')
    srv = 'https://ftp.ncbi.nlm.nih.gov/genomes/all'
    ext = 'genomic.fna.gz'
    return ('%s/%s/%s/%s/%s/%s_%s/%s_%s_%s'
            % (srv, asb[:3], asb[4:7], asb[7:10], asb[10:13],
               asb, asm, asb, asm, ext))


def extract_subsequences_from_fna(fna_path, organism_name, label=None, window_size=150, stride=150):
    records = []
    with gzip.open(fna_path, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            seq = str(record.seq)
            for start in range(0, len(seq) - window_size + 1, stride):
                sub_seq = seq[start:start + window_size]
                records.append({
                    "organism": organism_name,
                    "sequence": sub_seq,
                    "location": (start, start + window_size),
                    "label": 1,
                    "contig": record.id
                })
    return pd.DataFrame(records)


if __name__ == '__main__':
    main()
