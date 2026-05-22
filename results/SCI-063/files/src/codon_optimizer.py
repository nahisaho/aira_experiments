#!/usr/bin/env python3
"""
Module 2: Codon Optimization & Genome Stability
"""
import numpy as np
import pandas as pd
import json, os

class CodonOptimizer:
    CODON_TABLE = {
        'F':['TTT','TTC'],'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
        'I':['ATT','ATC','ATA'],'M':['ATG'],'V':['GTT','GTC','GTA','GTG'],
        'S':['TCT','TCC','TCA','TCG','AGT','AGC'],'P':['CCT','CCC','CCA','CCG'],
        'T':['ACT','ACC','ACA','ACG'],'A':['GCT','GCC','GCA','GCG'],
        'Y':['TAT','TAC'],'*':['TAA','TAG'],'H':['CAT','CAC'],'Q':['CAA','CAG'],
        'N':['AAT','AAC'],'K':['AAA','AAG'],'D':['GAT','GAC'],'E':['GAA','GAG'],
        'C':['TGT','TGC'],'W':['TGG','TGA'],'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
        'G':['GGT','GGC','GGA','GGG'],
    }
    PREFERRED = {
        'F':'TTT','L':'TTA','I':'ATT','M':'ATG','V':'GTT','S':'TCT','P':'CCT',
        'T':'ACT','A':'GCT','Y':'TAT','*':'TAA','H':'CAT','Q':'CAA','N':'AAT',
        'K':'AAA','D':'GAT','E':'GAA','C':'TGT','W':'TGA','R':'AGA','G':'GGT',
    }

    def __init__(self, seed=42):
        self.rng = np.random.RandomState(seed)
        self.codon_to_aa = {}
        for aa, codons in self.CODON_TABLE.items():
            for c in codons: self.codon_to_aa[c] = aa

    def generate_synthetic_genome(self, n_genes=380, avg_len=300):
        aas = [a for a in self.CODON_TABLE if a != '*']
        freqs = {'M':0.02,'W':0.01,'C':0.01,'H':0.02,'Y':0.03,'F':0.04,'Q':0.04,'N':0.04,
                 'D':0.05,'E':0.06,'K':0.06,'R':0.05,'S':0.07,'T':0.05,'A':0.08,
                 'G':0.07,'P':0.05,'V':0.07,'I':0.06,'L':0.10}
        aa_list = list(freqs.keys())
        aa_prob = np.array([freqs[a] for a in aa_list]); aa_prob /= aa_prob.sum()
        genes = {}
        for i in range(n_genes):
            length = max(50, int(self.rng.lognormal(np.log(avg_len), 0.5)))
            protein = 'M' + ''.join(self.rng.choice(aa_list, size=length-1, p=aa_prob))
            dna = ''.join(self.rng.choice(self.CODON_TABLE[aa]) for aa in protein)
            dna += self.rng.choice(self.CODON_TABLE['*'])
            genes[f"MG_{i+1:04d}"] = {'protein': protein, 'original_dna': dna, 'length_aa': len(protein), 'length_nt': len(dna)}
        self.genes = genes
        return genes

    def calculate_cai(self, seq):
        codons = [seq[i:i+3] for i in range(0, len(seq)-3, 3)]
        if not codons: return 0
        scores = []
        for c in codons:
            aa = self.codon_to_aa.get(c)
            if aa and aa != '*':
                scores.append(1.0 if c == self.PREFERRED[aa] else 1.0/len(self.CODON_TABLE[aa]))
        return np.exp(np.mean(np.log(np.array(scores)+1e-10))) if scores else 0

    def find_repeats(self, seq, min_len=12):
        repeats = []
        for kl in range(min_len, min(50, len(seq)//2)):
            seen = {}
            for i in range(len(seq)-kl+1):
                km = seq[i:i+kl]
                seen.setdefault(km, []).append(i)
            for km, pos in seen.items():
                if len(pos) > 1: repeats.append({'sequence':km,'length':kl,'count':len(pos),'positions':pos})
        return repeats

    def optimize_gene(self, gene_id, data):
        protein = data['protein']
        orig = data['original_dna']
        codons = [self.PREFERRED[aa] for aa in protein] + [self.PREFERRED['*']]
        prot_stop = protein + '*'
        for _ in range(5):
            seq = ''.join(codons)
            reps = self.find_repeats(seq, 12)
            if not reps: break
            for r in reps:
                ci = r['positions'][-1] // 3
                if ci < len(prot_stop):
                    alts = [c for c in self.CODON_TABLE.get(prot_stop[ci],[]) if c != codons[ci]]
                    if alts: codons[ci] = self.rng.choice(alts)
        final = ''.join(codons)
        gc = lambda s: (s.count('G')+s.count('C'))/max(len(s),1)
        return {'gene_id':gene_id,'optimized_dna':final,'length_nt':len(final),
                'gc_original':gc(orig),'gc_optimized':gc(final),
                'cai_original':self.calculate_cai(orig),'cai_optimized':self.calculate_cai(final),
                'repeats_original':len(self.find_repeats(orig,12)),'repeats_optimized':len(self.find_repeats(final,12))}

    def optimize_genome(self):
        self.optimization_results = pd.DataFrame([self.optimize_gene(gid, gd) for gid, gd in self.genes.items()])
        return self.optimization_results

    def genome_stability_analysis(self):
        df = self.optimization_results
        full = ''.join(df['optimized_dna'].values)
        genome_reps = self.find_repeats(full, 20)
        return {
            'total_genes':len(df),'total_genome_length_bp':len(full),
            'avg_cai_before':float(df['cai_original'].mean()),'avg_cai_after':float(df['cai_optimized'].mean()),
            'cai_improvement':float(df['cai_optimized'].mean()-df['cai_original'].mean()),
            'avg_gc_before':float(df['gc_original'].mean()),'avg_gc_after':float(df['gc_optimized'].mean()),
            'total_gene_repeats_before':int(df['repeats_original'].sum()),
            'total_gene_repeats_after':int(df['repeats_optimized'].sum()),
            'genome_level_repeats_20bp':len(genome_reps),
            'repeat_reduction_pct':float((1-df['repeats_optimized'].sum()/max(df['repeats_original'].sum(),1))*100),
        }

    def generate_plots(self, output_dir='figures'):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        os.makedirs(output_dir, exist_ok=True)
        df = self.optimization_results
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        axes[0].scatter(df['cai_original'], df['cai_optimized'], alpha=0.4, s=10, color='#1565C0')
        axes[0].plot([0,1],[0,1],'k--',alpha=0.3)
        axes[0].set_xlabel('Original CAI'); axes[0].set_ylabel('Optimized CAI')
        axes[0].set_title('Codon Adaptation Index Improvement')
        axes[1].hist(df['gc_optimized'], bins=30, alpha=0.7, color='#2E7D32', edgecolor='white')
        axes[1].axvline(0.317, color='red', linestyle='--', linewidth=2, label='Target (31.7%)')
        axes[1].set_xlabel('GC Content'); axes[1].set_title('GC Distribution (Optimized)'); axes[1].legend()
        axes[2].bar(['Before','After'],[df['repeats_original'].sum(),df['repeats_optimized'].sum()],
                     color=['#E53935','#43A047'], edgecolor='white')
        axes[2].set_ylabel('Repeat Count'); axes[2].set_title('Repetitive Sequences (≥12bp)')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig5_codon_optimization.png', dpi=300, bbox_inches='tight')
        plt.close()
        return [f'{output_dir}/fig5_codon_optimization.png']

def run_module2():
    print("\n"+"="*60); print("MODULE 2: Codon Optimization & Genome Stability"); print("="*60)
    o = CodonOptimizer(42)
    o.generate_synthetic_genome(380, 300)
    print(f"  Generated {len(o.genes)} genes")
    o.optimize_genome()
    o.optimization_results.to_csv('results/codon_optimization_results.csv', index=False)
    s = o.genome_stability_analysis()
    with open('results/module2_stability.json','w') as f: json.dump(s, f, indent=2)
    print(f"  Genome: {s['total_genome_length_bp']:,} bp")
    print(f"  CAI: {s['avg_cai_before']:.4f} → {s['avg_cai_after']:.4f}")
    print(f"  Repeats reduced by {s['repeat_reduction_pct']:.1f}%")
    plots = o.generate_plots()
    for p in plots: print(f"  Saved: {p}")
    return o, s

if __name__ == '__main__': run_module2()
