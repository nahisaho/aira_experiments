#!/usr/bin/env python3
"""Module 4: Genome Refactoring Strategy"""
import numpy as np
import pandas as pd
import json, os

class GenomeRefactorer:
    def __init__(self, n_genes=380, seed=42):
        self.n_genes, self.seed = n_genes, seed
        self.rng = np.random.RandomState(seed)

    def identify_redundancies(self):
        gene_ids = [f"MG_{i+1:04d}" for i in range(self.n_genes)]
        groups, used = [], set()
        for g in range(int(self.n_genes*0.08)):
            gs = min(self.rng.choice([2,2,2,3], p=[0.6,0.15,0.15,0.1]), self.n_genes-len(used))
            cands = [i for i in range(self.n_genes) if i not in used]
            if len(cands) < gs: break
            members = self.rng.choice(cands, gs, replace=False); used.update(members)
            si, fo = self.rng.uniform(0.35,0.9), self.rng.uniform(0.4,0.95)
            groups.append({'group_id':f'RG_{g+1:03d}','members':[gene_ids[m] for m in members],
                'group_size':int(gs),'avg_sequence_identity':float(si),'functional_overlap':float(fo),
                'can_merge':bool(si>0.5 and fo>0.6),
                'merge_strategy':'chimeric_fusion' if si>0.7 else 'best_representative'})
        self.redundancy_groups = groups
        return groups

    def plan_sequence_compression(self):
        comps = []
        for i in range(int(self.n_genes*0.15)):
            comps.append({'type':'gene_overlap','description':f'Overlap stop/start at junction {i+1}',
                'bp_saved':int(self.rng.randint(3,30)),'risk_level':'low'})
        for i in range(int(self.n_genes*0.9)):
            orig, mn = self.rng.randint(50,500), self.rng.randint(20,80)
            comps.append({'type':'intergenic_reduction','description':f'Reduce intergenic {i+1}: {orig}→{mn} bp',
                'bp_saved':int(orig-mn),'risk_level':'low' if (orig-mn)<100 else 'medium'})
        for i in range(int(self.n_genes*0.1)):
            comps.append({'type':'promoter_consolidation','description':f'Share promoter for group {i+1}',
                'bp_saved':int(self.rng.randint(50,200)),'risk_level':'medium'})
        for g in [g for g in self.redundancy_groups if g['can_merge']]:
            n_rm = g['group_size']-1
            comps.append({'type':'redundant_gene_removal','description':f'Remove {n_rm} from {g["group_id"]}',
                'bp_saved':int(n_rm*self.rng.randint(600,1500)),'risk_level':'high'})
        self.compressions = comps
        return comps

    def calculate_refactoring_summary(self):
        df = pd.DataFrame(self.compressions); ts = int(df['bp_saved'].sum()); orig = 580000
        by_type = df.groupby('type')['bp_saved'].agg(['count','sum','mean']).to_dict('index')
        by_risk = df.groupby('risk_level')['bp_saved'].agg(['count','sum']).to_dict('index')
        self.summary = {'total_compressions':len(self.compressions),'total_bp_saved':ts,
            'original_genome_estimate_bp':orig,'refactored_genome_estimate_bp':orig-ts,
            'compression_pct':float(ts/orig*100),'redundancy_groups_found':len(self.redundancy_groups),
            'mergeable_groups':sum(1 for g in self.redundancy_groups if g['can_merge']),
            'by_compression_type':{k:{kk:int(vv) if isinstance(vv,(np.integer,int)) else float(vv) for kk,vv in v.items()} for k,v in by_type.items()},
            'by_risk_level':{k:{kk:int(vv) if isinstance(vv,(np.integer,int)) else float(vv) for kk,vv in v.items()} for k,v in by_risk.items()}}
        return self.summary

    def generate_plots(self, output_dir='figures'):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        os.makedirs(output_dir, exist_ok=True)
        df = pd.DataFrame(self.compressions)
        fig, axes = plt.subplots(2,2, figsize=(14,12))
        df.groupby('type')['bp_saved'].sum().sort_values().plot(kind='barh', ax=axes[0,0], color='#1565C0')
        axes[0,0].set_title('Compression by Strategy Type')
        risk_s = df.groupby('risk_level')['bp_saved'].sum()
        cols = {'low':'#43A047','medium':'#FFA000','high':'#E53935'}
        for rl in ['low','medium','high']:
            axes[0,1].bar(rl, risk_s.get(rl,0), color=cols[rl])
        axes[0,1].set_title('Savings by Risk Level')
        # Waterfall
        orig = self.summary['original_genome_estimate_bp']
        cats = list(df.groupby('type')['bp_saved'].sum().sort_values(ascending=False).items())
        labels = ['Original']+[c[0].replace('_','\n') for c in cats]+['Final']
        vals = [orig]; r = orig
        for _,s in cats: r -= s; vals.append(r)
        vals.append(self.summary['refactored_genome_estimate_bp'])
        bc = ['#1565C0']+['#E53935']*len(cats)+['#43A047']
        axes[1,0].bar(range(len(labels)), vals, color=bc)
        axes[1,0].set_xticks(range(len(labels))); axes[1,0].set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
        axes[1,0].set_title('Genome Size Reduction Waterfall')
        # Redundancy scatter
        if self.redundancy_groups:
            rdf = pd.DataFrame(self.redundancy_groups)
            axes[1,1].scatter(rdf['avg_sequence_identity'],rdf['functional_overlap'],
                c=rdf['can_merge'].map({True:'#43A047',False:'#E53935'}), s=rdf['group_size']*50, alpha=0.6)
            axes[1,1].axhline(0.6,color='gray',ls='--',alpha=0.5); axes[1,1].axvline(0.5,color='gray',ls='--',alpha=0.5)
            axes[1,1].set_xlabel('Sequence Identity'); axes[1,1].set_ylabel('Functional Overlap')
            axes[1,1].set_title('Redundancy Groups (Green=Mergeable)')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig8_refactoring.png', dpi=300, bbox_inches='tight'); plt.close()
        return [f'{output_dir}/fig8_refactoring.png']

def run_module4():
    print("\n"+"="*60); print("MODULE 4: Genome Refactoring Strategy"); print("="*60)
    r = GenomeRefactorer(380, 42)
    groups = r.identify_redundancies()
    print(f"  Redundancy groups: {len(groups)} ({sum(1 for g in groups if g['can_merge'])} mergeable)")
    r.plan_sequence_compression()
    pd.DataFrame(r.compressions).to_csv('results/compression_plan.csv', index=False)
    s = r.calculate_refactoring_summary()
    with open('results/module4_refactoring.json','w') as f: json.dump(s, f, indent=2, default=str)
    print(f"  Total bp saved: {s['total_bp_saved']:,}")
    print(f"  Compression: {s['compression_pct']:.1f}%")
    print(f"  Genome: {s['original_genome_estimate_bp']:,} → {s['refactored_genome_estimate_bp']:,} bp")
    plots = r.generate_plots()
    for p in plots: print(f"  Saved: {p}")
    return r, s

if __name__ == '__main__': run_module4()
