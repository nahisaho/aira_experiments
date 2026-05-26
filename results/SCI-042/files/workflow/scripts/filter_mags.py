#!/usr/bin/env python3
"""Filter MAGs by quality and create unified summary."""
import pandas as pd
import os

def filter_mags(checkm_files, gtdbtk_files, output_file, params):
    all_data = []
    for cf, gf in zip(checkm_files, gtdbtk_files):
        checkm = pd.read_csv(cf, sep='\t')
        try:
            gtdbtk = pd.read_csv(gf, sep='\t')
            merged = checkm.merge(gtdbtk, left_on='Name', right_on='user_genome', how='left')
        except Exception:
            merged = checkm.copy()
        all_data.append(merged)
    
    df = pd.concat(all_data, ignore_index=True)
    # Quality categories
    df['Quality'] = 'Low'
    medium = (df['Completeness'] >= params['min_completeness']) & \
             (df['Contamination'] < params['max_contamination'])
    high = (df['Completeness'] >= params['hq_completeness']) & \
           (df['Contamination'] < params['hq_contamination'])
    df.loc[medium, 'Quality'] = 'Medium'
    df.loc[high, 'Quality'] = 'High'
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    filter_mags(snakemake.input.checkm, snakemake.input.gtdbtk,
                snakemake.output[0], dict(snakemake.params))
