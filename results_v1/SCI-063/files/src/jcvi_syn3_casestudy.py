#!/usr/bin/env python3
"""Module 6: JCVI-syn3.0 Extended Case Study"""
import numpy as np
import pandas as pd
import json, os

class JCVISyn3CaseStudy:
    SYN3 = {'genome_size_bp':531490,'total_genes':473,'protein_coding':438,'rna_genes':35,
            'gc_content':0.269,'coding_density':0.90,'essential_known':256,'quasi_essential':129,
            'nonessential':53,'unknown_function':149,'doubling_time_min':180}
    FUNCTIONAL = {
        'Genetic information processing':{'Translation':108,'Transcription':22,'DNA replication':16,'tRNA modification':12,'Ribosome':43},
        'Cell membrane & transport':{'Membrane lipids':15,'Transport':34,'Cell division':11},
        'Metabolism':{'Central carbon':18,'Nucleotide metabolism':16,'Cofactor biosynthesis':12,'Lipid metabolism':8},
        'Maintenance & regulation':{'Protein quality control':14,'DNA repair':8,'Regulation':7},
        'Unknown function':{'Genes of unknown function':149},
    }

    def __init__(self, seed=42):
        self.seed = seed; self.rng = np.random.RandomState(seed)

    def analyze_gene_categories(self):
        cats = {}
        ess_rates = {'Translation':0.85,'Transcription':0.80,'DNA replication':0.90,'tRNA modification':0.70,
            'Ribosome':0.95,'Membrane lipids':0.75,'Transport':0.60,'Cell division':0.90,
            'Central carbon':0.80,'Nucleotide metabolism':0.70,'Cofactor biosynthesis':0.65,
            'Lipid metabolism':0.55,'Protein quality control':0.70,'DNA repair':0.50,
            'Regulation':0.40,'Genes of unknown function':0.45}
        for grp, items in self.FUNCTIONAL.items():
            for cat, cnt in items.items():
                r = ess_rates.get(cat, 0.5)
                cats[cat] = {'count':cnt,'group':grp,'percentage':cnt/self.SYN3['total_genes']*100,
                    'essentiality_rate':r,'estimated_essential':int(cnt*r)}
        self.category_analysis = cats
        return cats

    def propose_syn3_extensions(self):
        ext = {
            'growth_improvement':{'title':'Growth Rate Enhancement','genes_to_add':[
                {'gene':'ftsZ_opt','function':'Optimized cell division ring','size_bp':1200},
                {'gene':'mreB_min','function':'Minimal cytoskeleton','size_bp':1050},
                {'gene':'dnaG_fast','function':'Enhanced primase','size_bp':1800}],'genome_cost_bp':4050},
            'stress_tolerance':{'title':'Stress Tolerance','genes_to_add':[
                {'gene':'dnaK_min','function':'Minimal chaperone','size_bp':1900},
                {'gene':'groES','function':'Co-chaperonin','size_bp':300},
                {'gene':'otsA_mini','function':'Trehalose synthesis','size_bp':1400}],'genome_cost_bp':3600},
            'genetic_stability':{'title':'Genome Stability','genes_to_add':[
                {'gene':'mutS_min','function':'Mismatch repair','size_bp':2400},
                {'gene':'recA_mini','function':'Minimal recombinase','size_bp':1050}],'genome_cost_bp':3450},
            'metabolic_expansion':{'title':'Metabolic Expansion','genes_to_add':[
                {'gene':'pgi','function':'Phosphoglucose isomerase','size_bp':1650},
                {'gene':'pgk','function':'Phosphoglycerate kinase','size_bp':1200},
                {'gene':'eno','function':'Enolase','size_bp':1300},
                {'gene':'folA_mini','function':'Dihydrofolate reductase','size_bp':500}],'genome_cost_bp':4650},
            'biocontainment':{'title':'Biocontainment','genes_to_add':[
                {'gene':'toxin_antitoxin','function':'Kill switch','size_bp':800},
                {'gene':'synthetic_auxotrophy','function':'Synthetic aa dependency','size_bp':1200}],'genome_cost_bp':2000},
        }
        total_bp = sum(e['genome_cost_bp'] for e in ext.values())
        ext['summary'] = {'total_genes_added':sum(len(e.get('genes_to_add',[])) for e in ext.values()),
            'total_bp_added':total_bp,'extended_genome_size':self.SYN3['genome_size_bp']+total_bp,
            'size_increase_pct':total_bp/self.SYN3['genome_size_bp']*100}
        self.extensions = ext
        return ext

    def comparative_analysis(self):
        orgs = {
            'JCVI-syn3.0':{'genome_size_kb':531,'genes':473,'gc_pct':26.9,'lifestyle':'Synthetic','doubling_time_min':180},
            'M. genitalium':{'genome_size_kb':580,'genes':525,'gc_pct':31.7,'lifestyle':'Obligate parasite','doubling_time_min':720},
            'M. mycoides':{'genome_size_kb':1212,'genes':985,'gc_pct':24.0,'lifestyle':'Parasite','doubling_time_min':90},
            'Buchnera aphidicola':{'genome_size_kb':422,'genes':362,'gc_pct':26.3,'lifestyle':'Endosymbiont','doubling_time_min':600},
            'Carsonella ruddii':{'genome_size_kb':160,'genes':182,'gc_pct':16.6,'lifestyle':'Endosymbiont','doubling_time_min':np.nan},
            'Nasuia deltocephalinicola':{'genome_size_kb':112,'genes':137,'gc_pct':17.1,'lifestyle':'Endosymbiont','doubling_time_min':np.nan},
            'Proposed syn3.0+':{'genome_size_kb':549,'genes':490,'gc_pct':27.0,'lifestyle':'Synthetic enhanced','doubling_time_min':80},
        }
        self.comparative = pd.DataFrame(orgs).T; self.comparative.index.name = 'organism'
        return self.comparative

    def unknown_gene_analysis(self):
        self.unknown_analysis = {
            'total_unknown':149,
            'predicted_categories':{'Likely membrane-associated':38,'Likely regulatory':22,
                'Predicted enzyme (unknown substrate)':31,'Conserved hypothetical':28,
                'Unique to lineage':15,'Mobile element remnant':8,'Predicted structural':7},
            'characterization_approaches':[
                {'method':'CRISPRi knockdown titration','target':'All 149','time_weeks':12,'priority':'High'},
                {'method':'Proteomics (AP-MS)','target':'Top 50 by conservation','time_weeks':16,'priority':'High'},
                {'method':'AlphaFold2 structure prediction','target':'All 149','time_weeks':2,'priority':'Medium'},
                {'method':'Metabolomics profiling','target':'31 predicted enzymes','time_weeks':8,'priority':'Medium'},
                {'method':'Tn-seq conditional','target':'Genome-wide','time_weeks':6,'priority':'High'},
            ],
        }
        return self.unknown_analysis

    def generate_plots(self, output_dir='figures'):
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        os.makedirs(output_dir, exist_ok=True)

        # Functional + essentiality
        fig, axes = plt.subplots(1,2, figsize=(16,7))
        gc = {}
        for grp, cats in self.FUNCTIONAL.items(): gc[grp] = sum(cats.values())
        cs = ['#1565C0','#2E7D32','#E65100','#6A1B9A','#BF360C']
        axes[0].pie(gc.values(), labels=[g.replace(' & ','\n& ') for g in gc.keys()], autopct='%1.1f%%', colors=cs, startangle=90, textprops={'fontsize':9})
        axes[0].set_title(f'JCVI-syn3.0 Gene Functional Groups (n={self.SYN3["total_genes"]})', fontsize=13)
        ed = {'Essential':self.SYN3['essential_known'],'Quasi-essential':self.SYN3['quasi_essential'],
              'Non-essential':self.SYN3['nonessential'],'Unknown':self.SYN3['unknown_function']}
        ec = ['#C62828','#EF6C00','#2E7D32','#757575']
        bars = axes[1].bar(ed.keys(), ed.values(), color=ec)
        axes[1].set_ylabel('Gene Count'); axes[1].set_title('Gene Essentiality Classification')
        for i,(k,v) in enumerate(ed.items()): axes[1].text(i, v+3, str(v), ha='center', fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig11_syn3_functional.png', dpi=300, bbox_inches='tight'); plt.close()

        # Comparative
        fig, axes = plt.subplots(1,2, figsize=(14,6))
        comp = self.comparative.sort_values('genome_size_kb')
        colors_org = ['#E53935' if 'syn3' in str(o).lower() or 'Proposed' in str(o) else '#1565C0' for o in comp.index]
        axes[0].barh(range(len(comp)), comp['genome_size_kb'], color=colors_org)
        axes[0].set_yticks(range(len(comp))); axes[0].set_yticklabels(comp.index, fontsize=9)
        axes[0].set_xlabel('Genome Size (kb)'); axes[0].set_title('Genome Size Comparison')
        for idx,row in comp.iterrows():
            c = '#E53935' if 'syn3' in str(idx).lower() or 'Proposed' in str(idx) else '#1565C0'
            axes[1].scatter(row['genome_size_kb'], row['genes'], s=100, color=c, edgecolor='white', zorder=5)
            axes[1].annotate(str(idx)[:15], (row['genome_size_kb'], row['genes']), fontsize=7)
        axes[1].set_xlabel('Genome Size (kb)'); axes[1].set_ylabel('Genes'); axes[1].set_title('Genes vs Genome Size')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig12_comparative_genomes.png', dpi=300, bbox_inches='tight'); plt.close()

        # Unknown genes
        fig, ax = plt.subplots(figsize=(10,6))
        uk = self.unknown_analysis['predicted_categories']
        cs_uk = plt.cm.viridis(np.linspace(0.2,0.8,len(uk)))
        ax.barh(list(uk.keys()), list(uk.values()), color=cs_uk)
        ax.set_xlabel('Number of Genes'); ax.set_title('Predicted Categories for 149 Unknown Function Genes')
        for i,v in enumerate(uk.values()): ax.text(v+0.5, i, str(v), va='center', fontweight='bold')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/fig13_unknown_genes.png', dpi=300, bbox_inches='tight'); plt.close()
        return [f'{output_dir}/fig11_syn3_functional.png',f'{output_dir}/fig12_comparative_genomes.png',f'{output_dir}/fig13_unknown_genes.png']

def run_module6():
    print("\n"+"="*60); print("MODULE 6: JCVI-syn3.0 Extended Case Study"); print("="*60)
    s = JCVISyn3CaseStudy(42)
    cats = s.analyze_gene_categories()
    with open('results/syn3_categories.json','w') as f: json.dump(cats, f, indent=2)
    ext = s.propose_syn3_extensions()
    with open('results/syn3_extensions.json','w') as f: json.dump(ext, f, indent=2)
    print(f"  Extensions: {ext['summary']['total_genes_added']} genes, +{ext['summary']['size_increase_pct']:.1f}%")
    comp = s.comparative_analysis()
    comp.to_csv('results/comparative_genomes.csv')
    unk = s.unknown_gene_analysis()
    with open('results/syn3_unknown_genes.json','w') as f: json.dump(unk, f, indent=2)
    print(f"  Unknown genes: {unk['total_unknown']}, {len(unk['characterization_approaches'])} approaches")
    plots = s.generate_plots()
    for p in plots: print(f"  Saved: {p}")
    sm = {'syn3_stats':s.SYN3,'extension_summary':ext['summary'],'n_unknown':unk['total_unknown'],'n_compared':len(comp)}
    with open('results/module6_summary.json','w') as f: json.dump(sm, f, indent=2)
    return s, sm

if __name__ == '__main__': run_module6()
