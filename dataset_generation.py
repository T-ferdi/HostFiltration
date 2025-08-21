#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO
import gzip
import urllib.request
import argparse
import random
from tqdm import tqdm


def tqdm_hook(t):
    last_b = [0]
    def update_to(block_num=1, block_size=1, total_size=None):
        if total_size is not None:
            t.total = total_size
        downloaded = block_num * block_size
        t.update(downloaded - last_b[0])
        last_b[0] = downloaded
    return update_to


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True,
                        help="Path to the metadata TSV file (uncompressed)")
    parser.add_argument("--sample-size", type=int, default=10,
                        help="Number of microbial genomes to sample and extract from")
    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep='\t')

    required_cols = ['assembly_accession', 'asm_name']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f'Missing required column: {col}')

    sample_df = df.sample(n=min(args.sample_size, len(df)), random_state=42)

    all_microbial_data = []

    # Extract microbial sequences
    for _, row in sample_df.iterrows():
        asb, asm, unique_org_name = row['assembly_accession'], row['asm_name'], row['unique_name']
        url = get_fna_path(asb, asm)
        fna_filename = f"{asb}_{asm.replace(' ', '_')}_genomic.fna.gz"

        print(f"Downloading microbial genome: {url}")
        try:
            with tqdm(unit='B', unit_scale=True, desc=fna_filename, leave=True) as t:
                urllib.request.urlretrieve(url, fna_filename, reporthook=tqdm_hook(t))
            df = extract_subsequences_from_fna(fna_filename, unique_org_name, label=1)
            all_microbial_data.append(df)
        except Exception as e:
            print(f"Failed to download or process {url}: {e}")

    if not all_microbial_data:
        print("❌ No microbial data extracted.")
        return

    microbial_df = pd.concat(all_microbial_data, ignore_index=True)
    num_microbial_seqs = len(microbial_df)
    print(f"✅ Extracted {num_microbial_seqs} microbial sequences.")

    # Download and extract matched number of human sequences
    human_url = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
                 "GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz")
    human_fna = "human_GRCh38.fna.gz"
    print(f"Downloading human genome: {human_url}")
    try:
        with tqdm(unit='B', unit_scale=True, desc=human_fna, leave=True) as t:
            urllib.request.urlretrieve(human_url, human_fna, reporthook=tqdm_hook(t))
        human_df = extract_random_subsequences_from_fna(
            human_fna, "Homo sapiens", label=0, window_size=150, n=num_microbial_seqs)
        print(f"✅ Extracted {len(human_df)} human sequences.")
    except Exception as e:
        print(f"Failed to download or process human genome: {e}")
        return

    # Combine and save
    final_df = pd.concat([microbial_df, human_df], ignore_index=True)
    final_df.to_csv("subsequences_dataset.csv", index=False)
    print("✅ Dataset saved to subsequences_dataset.csv")


def get_fna_path(asb, asm):
    asm = asm.replace(' ', '_').replace('#', '_')
    srv = 'https://ftp.ncbi.nlm.nih.gov/genomes/all'
    ext = 'genomic.fna.gz'
    return ('%s/%s/%s/%s/%s/%s_%s/%s_%s_%s'
            % (srv, asb[:3], asb[4:7], asb[7:10], asb[10:13],
               asb, asm, asb, asm, ext))


def extract_subsequences_from_fna(fna_path, organism_name, label, window_size=150, stride=150):
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
                    "label": label,
                    "contig": record.id,
                    "source_file": os.path.basename(fna_path)
                })
    return pd.DataFrame(records)


def extract_random_subsequences_from_fna(fna_path, organism_name, label, window_size=150, n=10):
    records = []
    all_seqs = []

    with gzip.open(fna_path, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            if len(record.seq) >= window_size:
                all_seqs.append((record.id, str(record.seq)))

    while len(records) < n and all_seqs:
        contig_id, seq = random.choice(all_seqs)
        if len(seq) < window_size:
            continue
        start = random.randint(0, len(seq) - window_size)
        sub_seq = seq[start:start + window_size]
        records.append({
            "organism": organism_name,
            "sequence": sub_seq,
            "location": (start, start + window_size),
            "label": label,
            "contig": contig_id,
            "source_file": os.path.basename(fna_path)
        })

    return pd.DataFrame(records)


if __name__ == '__main__':
    main()
