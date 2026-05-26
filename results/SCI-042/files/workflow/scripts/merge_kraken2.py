#!/usr/bin/env python3
"""Merge Kraken2/Bracken output files into a unified abundance table."""
import pandas as pd
import sys
import os

def merge_bracken(input_files, output_file):
    dfs = []
    for f in input_files:
        sample = os.path.basename(f).replace('.bracken', '')
        df = pd.read_csv(f, sep='\t')
        df = df[['name', 'fraction_total_reads']].rename(
            columns={'name': 'Taxon', 'fraction_total_reads': sample}
        )
        dfs.append(df)
    
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on='Taxon', how='outer')
    merged = merged.fillna(0)
    merged.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    merge_bracken(snakemake.input, snakemake.output[0])
