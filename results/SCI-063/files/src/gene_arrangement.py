#!/usr/bin/env python3
"""Module 3: Gene Arrangement Optimization"""
import numpy as np
import pandas as pd
import networkx as nx
import json, os

class GeneArrangementOptimizer:
    CATEGORIES = ['translation','transcription','replication','cell_division','membrane_transport',
                  'energy_metabolism','nucleotide_metabolism','lipid_metabolism','protein_folding',
                  'dna_repair','cell_envelope','regulation','unknown']

    def __init__(self, n_genes=380, genome_size=530000, seed=42):
        self.n_genes, self.genome_size, self.seed = n_genes, genome_size, seed
        self.rng = np.random.RandomState(seed)

    def generate_gene_annotations(self):
        genes, cat_w = [], [.20,.06,.05,.04,.08,.08,.06,.04,.05,.04,.06,.04,.20]
        operon_id, i = 0, 0
        while i < self.n_genes:
            op_size = min(self.rng.geometric(0.4), 8, self.n_genes - i)
            operon_id += 1
            cat = self.rng.choice(self.CATEGORIES, p=cat_w)
            for j in range(op_size):
                gl = max(150, int(self.rng.lognormal(np.log(900), 0.5)))
                es = self.rng.beta(3,1) if cat in ['translation','transcription','replication'] else self.rng.beta(2,2)
                genes.append({'gene_id':f'MG_{i+1:04d}','gene_length_bp':gl,
                    'category': cat if j==0 or self.rng.random()<0.7 else self.rng.choice(self.CATEGORIES, p=cat_w),
                    'operon_id':f'OP_{operon_id:03d}','position_in_operon':j,'essentiality_score':es,
                    'expression_level':self.rng.lognormal(2.0,1.0),'strand':1,'start_pos':0})
                i += 1
        self.genes_df = pd.DataFrame(genes)
        return self.genes_df

    def build_interaction_network(self):
        G = nx.Graph()
        for _, g in self.genes_df.iterrows(): G.add_node(g['gene_id'])
        for op in self.genes_df['operon_id'].unique():
            og = self.genes_df[self.genes_df['operon_id']==op]['gene_id'].tolist()
            for i in range(len(og)):
                for j in range(i+1,len(og)): G.add_edge(og[i],og[j],weight=0.9)
        for cat in self.CATEGORIES:
            cg = self.genes_df[self.genes_df['category']==cat]['gene_id'].tolist()
            for _ in range(min(len(cg)*2, len(cg)*(len(cg)-1)//2)):
                if len(cg)>=2:
                    g1,g2 = self.rng.choice(cg,2,replace=False)
                    if not G.has_edge(g1,g2): G.add_edge(g1,g2,weight=0.4)
        self.interaction_graph = G
        return G

    def optimize_strand_assignment(self):
        genes = self.genes_df.copy()
        cum = 0
        positions = []
        for _, g in genes.iterrows():
            positions.append(cum); cum += g['gene_length_bp'] + self.rng.randint(50,300)
        scale = self.genome_size / cum
        genes['start_pos'] = [int(p*scale) for p in positions]
        ter = self.genome_size // 2
        for idx, row in genes.iterrows():
            leading = 1 if row['start_pos'] <= ter else -1
            genes.at[idx,'strand'] = leading if self.rng.random() < 0.5+0.4*row['essentiality_score'] else -leading
        for op in genes['operon_id'].unique():
            mask = genes['operon_id']==op
            genes.loc[mask,'strand'] = genes.loc[mask].iloc[0]['strand']
        self.genes_df = genes
        return genes

    def calculate_arrangement_metrics(self):
        g = self.genes_df; ter = self.genome_size//2
        leading = sum(1 for _,r in g.iterrows() if r['strand']==(1 if r['start_pos']<=ter else -1))
        ess = g[g['essentiality_score']>0.7]
        ess_lead = sum(1 for _,r in ess.iterrows() if r['strand']==(1 if r['start_pos']<=ter else -1))
        ops = g.groupby('operon_id'); intact = sum(1 for _,o in ops if o['strand'].nunique()==1)
        return {
            'leading_strand_bias':float(leading/len(g)),
            'essential_leading_bias':float(ess_lead/max(len(ess),1)),
            'operon_integrity':float(intact/len(ops)),
            'intact_operons':intact,'total_operons':len(ops),
            'coding_density':float(g['gene_length_bp'].sum()/self.genome_size),
            'n_genes':len(g),
        }

    def generate_plots(self, output_dir='figures'):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        os.makedirs(output_dir, exist_ok=True)
        genes = self.genes_df
        # Circular genome map
        fig, ax = plt.subplots(figsize=(10,10), subplot_kw={'projection':'polar'})
        cats = genes['category'].unique()
        cmap = plt.cm.get_cmap('tab20', len(cats))
        cc = {c:cmap(i) for i,c in enumerate(cats)}
        for _,g in genes.iterrows():
            angle = 2*np.pi*g['start_pos']/self.genome_size
            width = 2*np.pi*g['gene_length_bp']/self.genome_size
            r = 1.0 if g['strand']==1 else 0.85
            ax.bar(angle, 0.12, width=width, bottom=r, color=cc[g['category']], alpha=0.7, edgecolor='white', linewidth=0.2)
        ax.annotate('oriC',xy=(0,1.2),fontsize=12,fontweight='bold',ha='center',color='red')
        ax.annotate('ter',xy=(np.pi,1.2),fontsize=12,fontweight='bold',ha='center',color='blue')
        ax.set_ylim(0.5,1.3); ax.set_yticklabels([])
        ax.set_title('Minimal Genome Map\n(Outer: +strand, Inner: -strand)', fontsize=14, pad=20)
        patches = [mpatches.Patch(color=cc[c],label=c.replace('_',' ')) for c in sorted(cats) if c!='unknown']
        ax.legend(handles=patches, loc='center', fontsize=7, ncol=2)
        plt.savefig(f'{output_dir}/fig6_genome_map.png', dpi=300, bbox_inches='tight'); plt.close()

        # Strand bias
        fig, axes = plt.subplots(1,2, figsize=(14,6))
        sd = genes.groupby(['category','strand']).size().unstack(fill_value=0)
        sd.plot(kind='barh', stacked=True, ax=axes[0], color=['#E53935','#1E88E5'])
        axes[0].set_title('Strand Distribution by Category')
        axes[0].legend(['Lagging (-)','Leading (+)'])
        ter = self.genome_size//2
        genes['on_leading'] = genes.apply(lambda r: r['strand']==(1 if r['start_pos']<=ter else -1), axis=1)
        eb = pd.cut(genes['essentiality_score'], bins=5)
        genes.groupby(eb, observed=True)['on_leading'].mean().plot(kind='bar', ax=axes[1], color='#1565C0')
        axes[1].set_title('Leading Strand Bias vs Essentiality')
        axes[1].axhline(0.5, color='red', linestyle='--', alpha=0.5)
        axes[1].tick_params(axis='x', rotation=45)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig7_strand_bias.png', dpi=300, bbox_inches='tight'); plt.close()
        return [f'{output_dir}/fig6_genome_map.png', f'{output_dir}/fig7_strand_bias.png']

def run_module3():
    print("\n"+"="*60); print("MODULE 3: Gene Arrangement Optimization"); print("="*60)
    o = GeneArrangementOptimizer(380, 530000, 42)
    o.generate_gene_annotations()
    print(f"  Genes: {len(o.genes_df)}, Operons: {o.genes_df['operon_id'].nunique()}")
    G = o.build_interaction_network()
    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    o.optimize_strand_assignment()
    o.genes_df.to_csv('results/gene_arrangement.csv', index=False)
    m = o.calculate_arrangement_metrics()
    with open('results/module3_arrangement.json','w') as f: json.dump(m, f, indent=2)
    print(f"  Leading strand bias: {m['leading_strand_bias']:.3f}")
    print(f"  Essential leading bias: {m['essential_leading_bias']:.3f}")
    print(f"  Operon integrity: {m['operon_integrity']:.3f}")
    plots = o.generate_plots()
    for p in plots: print(f"  Saved: {p}")
    return o, m

if __name__ == '__main__': run_module3()
