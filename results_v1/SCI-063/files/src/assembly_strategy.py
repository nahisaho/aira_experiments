#!/usr/bin/env python3
"""Module 5: Hierarchical Gibson Assembly Strategy"""
import numpy as np
import pandas as pd
import json, os

class GibsonAssemblyDesigner:
    GIBSON_OVERLAP = 40
    MAX_GIBSON_FRAGMENTS = 8

    def __init__(self, genome_size=530000, seed=42):
        self.genome_size, self.seed = genome_size, seed
        self.rng = np.random.RandomState(seed)

    def design_fragment_hierarchy(self):
        pos, l0 = 0, []
        i = 0
        while pos < self.genome_size:
            fs = min(self.rng.randint(6000,9500), self.genome_size - pos + self.GIBSON_OVERLAP)
            l0.append({'fragment_id':f'L0_{i+1:03d}','level':0,'start':pos,
                'end':min(pos+fs,self.genome_size),'size_bp':fs,
                'gc_content':float(self.rng.normal(0.32,0.03)),
                'synthesis_feasibility':float(self.rng.beta(8,2))})
            pos += fs - self.GIBSON_OVERLAP; i += 1
        self.level0 = l0
        fpl1 = self.MAX_GIBSON_FRAGMENTS - 1
        n1 = int(np.ceil(len(l0)/fpl1))
        l1 = []
        for j in range(n1):
            si, ei = j*fpl1, min((j+1)*fpl1, len(l0))
            children = [f['fragment_id'] for f in l0[si:ei]]
            ts = sum(f['size_bp'] for f in l0[si:ei]) - self.GIBSON_OVERLAP*(len(children)-1)
            l1.append({'fragment_id':f'L1_{j+1:03d}','level':1,'n_children':len(children),
                'children':children,'assembled_size_bp':ts,'assembly_method':'Gibson Assembly',
                'expected_efficiency':float(self.rng.beta(6,2))})
        self.level1 = l1
        fpl2 = max(2, int(np.ceil(len(l1)/3)))
        l2 = []
        for k in range(3):
            si, ei = k*fpl2, min((k+1)*fpl2, len(l1))
            if si >= len(l1): break
            children = [f['fragment_id'] for f in l1[si:ei]]
            ts = sum(f['assembled_size_bp'] for f in l1[si:ei])
            l2.append({'fragment_id':f'L2_{k+1:03d}','level':2,'n_children':len(children),
                'children':children,'assembled_size_bp':ts,'assembly_method':'Yeast TAR Cloning',
                'expected_efficiency':float(self.rng.beta(5,3))})
        self.level2 = l2
        ta = sum(f['assembled_size_bp'] for f in l2)
        self.level3 = {'fragment_id':'L3_001','level':3,'n_children':len(l2),
            'children':[f['fragment_id'] for f in l2],'assembled_size_bp':ta,
            'assembly_method':'Genome Transplantation'}
        return {'level0':l0,'level1':l1,'level2':l2,'level3':self.level3}

    def estimate_assembly_costs(self):
        ts_bp = sum(f['size_bp'] for f in self.level0)
        costs = {
            'synthesis':{'total_bp':ts_bp,'cost_per_bp':0.09,'total_cost':ts_bp*0.09,'n_fragments':len(self.level0),'estimated_time_weeks':4},
            'level1_gibson':{'n_reactions':len(self.level1),'cost_per_reaction':50,'total_cost':len(self.level1)*50,'estimated_time_weeks':2},
            'level2_tar_cloning':{'n_reactions':len(self.level2),'cost_per_reaction':200,'total_cost':len(self.level2)*200,'estimated_time_weeks':3},
            'level3_transplantation':{'n_attempts':5,'cost_per_attempt':5000,'total_cost':25000,'estimated_time_weeks':4},
        }
        costs['total'] = {'total_cost_usd':sum(v['total_cost'] for v in costs.values() if 'total_cost' in v),
                          'total_time_weeks':sum(v.get('estimated_time_weeks',0) for v in costs.values())}
        self.costs = costs
        return costs

    def quality_control_checkpoints(self):
        self.checkpoints = [
            {'stage':'Post-synthesis','method':'Sanger sequencing','acceptance':'100% identity'},
            {'stage':'Post-Level1','method':'Restriction digest + gel','acceptance':'Correct size'},
            {'stage':'Post-Level1-Seq','method':'Nanopore sequencing','acceptance':'<1 err/10kb'},
            {'stage':'Post-Level2','method':'PFGE + Southern blot','acceptance':'Correct size/pattern'},
            {'stage':'Post-Transplant','method':'WGS (Illumina+Nanopore)','acceptance':'<5 SNVs'},
            {'stage':'Viability','method':'Growth curve + omics','acceptance':'Doubling <120 min'},
        ]
        return self.checkpoints

    def generate_plots(self, output_dir='figures'):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        os.makedirs(output_dir, exist_ok=True)
        # Assembly hierarchy
        fig, ax = plt.subplots(figsize=(16,10))
        lc = {0:'#90CAF9',1:'#42A5F5',2:'#1565C0',3:'#0D47A1'}
        w0 = 14/max(len(self.level0),1)
        for i in range(len(self.level0)):
            ax.add_patch(plt.Rectangle((1+i*w0,0), w0*0.9, 0.8, facecolor=lc[0], edgecolor='white', lw=0.5, alpha=0.7))
        w1 = 14/max(len(self.level1),1)
        for i,f in enumerate(self.level1):
            ax.add_patch(plt.Rectangle((1+i*w1,3), w1*0.9, 1.0, facecolor=lc[1], edgecolor='white'))
            ax.text(1+i*w1+w1*0.45, 3.5, f['fragment_id'], ha='center', va='center', fontsize=7, color='white', fontweight='bold')
        w2 = 14/max(len(self.level2),1)
        for i,f in enumerate(self.level2):
            ax.add_patch(plt.Rectangle((1+i*w2,6), w2*0.9, 1.2, facecolor=lc[2], edgecolor='white'))
            ax.text(1+i*w2+w2*0.45, 6.6, f['fragment_id'], ha='center', va='center', fontsize=10, color='white', fontweight='bold')
        ax.add_patch(plt.Rectangle((3,9), 10, 1.5, facecolor=lc[3], edgecolor='black', lw=2))
        ax.text(8, 9.75, f"Complete Minimal Genome\n{self.genome_size/1000:.0f} kb", ha='center', va='center', fontsize=14, color='white', fontweight='bold')
        ax.text(0.3,0.4,f'Level 0\n({len(self.level0)} frags)',ha='center',va='center',fontsize=9,fontweight='bold')
        ax.text(0.3,3.5,f'Level 1\nGibson\n({len(self.level1)})',ha='center',va='center',fontsize=9,fontweight='bold')
        ax.text(0.3,6.6,f'Level 2\nTAR\n({len(self.level2)})',ha='center',va='center',fontsize=9,fontweight='bold')
        ax.text(0.3,9.75,'Level 3\nTransplant',ha='center',va='center',fontsize=9,fontweight='bold')
        for yf,yt in [(0.8,3),(4,6),(7.2,9)]:
            ax.annotate('',xy=(8,yt),xytext=(8,yf),arrowprops=dict(arrowstyle='->',lw=2,color='gray'))
        ax.set_xlim(-0.5,16); ax.set_ylim(-0.5,11.5); ax.axis('off')
        ax.set_title('Hierarchical Genome Assembly Strategy', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig9_assembly_hierarchy.png', dpi=300, bbox_inches='tight'); plt.close()

        # Cost breakdown
        fig, axes = plt.subplots(1,2, figsize=(14,6))
        ci = {'DNA Synthesis':self.costs['synthesis']['total_cost'],'Gibson Assembly':self.costs['level1_gibson']['total_cost'],
              'TAR Cloning':self.costs['level2_tar_cloning']['total_cost'],'Transplantation':self.costs['level3_transplantation']['total_cost']}
        cs = ['#90CAF9','#42A5F5','#1565C0','#0D47A1']
        axes[0].pie(ci.values(), labels=ci.keys(), autopct='%1.1f%%', colors=cs, startangle=90)
        axes[0].set_title(f'Cost Breakdown (Total: ${self.costs["total"]["total_cost_usd"]:,.0f})')
        stages = ['Synthesis','Level 1\nGibson','Level 2\nTAR','Level 3\nTransplant','QC']
        starts, durs = [0,4,6,9,13], [4,2,3,4,4]
        for i,(st,s,d) in enumerate(zip(stages,starts,durs)):
            axes[1].barh(i,d,left=s,color=cs[min(i,len(cs)-1)],edgecolor='white',height=0.6)
            axes[1].text(s+d/2,i,f'{d}w',ha='center',va='center',fontsize=10,color='white',fontweight='bold')
        axes[1].set_yticks(range(len(stages))); axes[1].set_yticklabels(stages)
        axes[1].set_xlabel('Time (weeks)'); axes[1].set_title('Assembly Timeline'); axes[1].invert_yaxis()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig10_assembly_costs.png', dpi=300, bbox_inches='tight'); plt.close()
        return [f'{output_dir}/fig9_assembly_hierarchy.png', f'{output_dir}/fig10_assembly_costs.png']

def run_module5():
    print("\n"+"="*60); print("MODULE 5: Hierarchical Gibson Assembly Strategy"); print("="*60)
    d = GibsonAssemblyDesigner(530000, 42)
    h = d.design_fragment_hierarchy()
    print(f"  L0: {len(h['level0'])} fragments, L1: {len(h['level1'])} Gibson, L2: {len(h['level2'])} TAR")
    costs = d.estimate_assembly_costs()
    with open('results/module5_assembly.json','w') as f: json.dump(costs, f, indent=2)
    print(f"  Cost: ${costs['total']['total_cost_usd']:,.0f}, Time: {costs['total']['total_time_weeks']} weeks")
    qc = d.quality_control_checkpoints()
    with open('results/qc_checkpoints.json','w') as f: json.dump(qc, f, indent=2)
    hs = {'level0_count':len(h['level0']),'level1_count':len(h['level1']),'level2_count':len(h['level2']),
          'level0_avg_size':int(np.mean([f['size_bp'] for f in h['level0']])),
          'level1_avg_size':int(np.mean([f['assembled_size_bp'] for f in h['level1']])),
          'level2_avg_size':int(np.mean([f['assembled_size_bp'] for f in h['level2']])),
          'final_genome_size':h['level3']['assembled_size_bp']}
    with open('results/assembly_hierarchy.json','w') as f: json.dump(hs, f, indent=2)
    plots = d.generate_plots()
    for p in plots: print(f"  Saved: {p}")
    return d, costs

if __name__ == '__main__': run_module5()
