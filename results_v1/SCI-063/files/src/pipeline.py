#!/usr/bin/env python3
"""
Minimal Genome Design Pipeline — Main Entry Point
Integrates all 6 modules for rational minimal genome design and synthesis.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))

from essential_gene_predictor import run_module1
from codon_optimizer import run_module2
from gene_arrangement import run_module3
from refactoring_strategy import run_module4
from assembly_strategy import run_module5
from jcvi_syn3_casestudy import run_module6

def main():
    start = time.time()
    print("╔" + "═"*58 + "╗")
    print("║  MINIMAL GENOME RATIONAL DESIGN & SYNTHESIS PIPELINE   ║")
    print("╚" + "═"*58 + "╝")
    print()

    os.makedirs('data', exist_ok=True)
    os.makedirs('results', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    os.makedirs('logs', exist_ok=True)

    all_results = {}

    # Module 1
    _, s1 = run_module1()
    all_results['module1'] = s1

    # Module 2
    _, s2 = run_module2()
    all_results['module2'] = s2

    # Module 3
    _, s3 = run_module3()
    all_results['module3'] = s3

    # Module 4
    _, s4 = run_module4()
    all_results['module4'] = s4

    # Module 5
    _, s5 = run_module5()
    all_results['module5'] = s5

    # Module 6
    _, s6 = run_module6()
    all_results['module6'] = s6

    elapsed = time.time() - start

    # Save complete results
    with open('results/pipeline_summary.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    # Process log
    log_entry = {
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
        'phase': 'complete',
        'event_type': 'run_completed',
        'actor': 'co-scientist',
        'skill_or_tool': 'minimal-genome-pipeline',
        'elapsed_seconds': round(elapsed, 1),
        'files_written': [],
        'status': 'ok',
    }
    with open('logs/process-log.jsonl', 'a') as f:
        f.write(json.dumps(log_entry) + '\n')

    print("\n" + "="*60)
    print("PIPELINE COMPLETE")
    print("="*60)
    print(f"  Time elapsed: {elapsed:.1f}s")
    print(f"  Results saved to: results/")
    print(f"  Figures saved to: figures/")
    print(f"  Data saved to: data/")

if __name__ == '__main__':
    main()
