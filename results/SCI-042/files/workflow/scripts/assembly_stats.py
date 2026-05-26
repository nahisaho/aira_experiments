#!/usr/bin/env python3
"""Calculate assembly statistics from FASTA contigs."""
import sys

def calc_stats(fasta_file, output_file):
    lengths = []
    current_len = 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_len > 0:
                    lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_len > 0:
        lengths.append(current_len)
    
    lengths.sort(reverse=True)
    total = sum(lengths)
    cumsum = 0
    n50 = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total / 2:
            n50 = l
            break
    
    with open(output_file, 'w') as out:
        out.write('Metric\tValue\n')
        out.write(f'Total_contigs\t{len(lengths)}\n')
        out.write(f'Total_length\t{total}\n')
        out.write(f'N50\t{n50}\n')
        out.write(f'Largest\t{lengths[0] if lengths else 0}\n')
        out.write(f'Mean_length\t{total/len(lengths) if lengths else 0:.0f}\n')

if __name__ == '__main__':
    calc_stats(snakemake.input[0], snakemake.output[0])
