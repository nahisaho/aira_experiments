"""
Task 6: GFP Fluorescence Optimization Case Study
- ESM-2 based scoring of GFP variants
- In silico directed evolution
- Fluorescence-activity landscape analysis
"""
import json
import os
import time
import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from transformers import AutoTokenizer, EsmForMaskedLM
from scipy.stats import spearmanr

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

# avGFP chromophore-region fragment (positions around 60-110 of real GFP)
GFP_WT = "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"

# Known important positions (0-indexed within our fragment)
CHROMOPHORE_POSITIONS = [5, 6, 7]  # "ELF" -> chromophore region
BETA_BARREL = list(range(10, 20)) + list(range(30, 45))

def simulate_gfp_fitness(wt_seq, mut_seq):
    """Simulate GFP fluorescence based on sequence features."""
    L = len(wt_seq)
    fitness = 1.0  # WT = 1.0

    for i, (wt, mt) in enumerate(zip(wt_seq, mut_seq)):
        if wt == mt:
            continue
        # Chromophore mutations are very deleterious
        if i in CHROMOPHORE_POSITIONS:
            fitness *= 0.1
        # Beta-barrel mutations are moderately deleterious
        elif i in BETA_BARREL:
            fitness *= 0.7
        # Surface mutations are neutral or slightly positive
        else:
            fitness *= np.random.choice([0.9, 1.0, 1.05], p=[0.3, 0.5, 0.2])

    return float(np.clip(fitness + np.random.normal(0, 0.05), 0, 2))

def generate_gfp_variants(wt_seq, n_variants=200, max_mutations=4):
    """Generate GFP variants with simulated fluorescence."""
    variants = []
    L = len(wt_seq)

    for i in range(n_variants):
        np.random.seed(SEED + i)
        n_mut = np.random.randint(1, max_mutations + 1)
        positions = np.random.choice(L, n_mut, replace=False)
        mut_seq = list(wt_seq)
        mutations = []
        for pos in positions:
            new_aa = np.random.choice([a for a in AA_LIST if a != wt_seq[pos]])
            mutations.append(f"{wt_seq[pos]}{pos+1}{new_aa}")
            mut_seq[pos] = new_aa
        mut_seq_str = ''.join(mut_seq)
        fitness = simulate_gfp_fitness(wt_seq, mut_seq_str)
        variants.append({
            'sequence': mut_seq_str,
            'mutations': mutations,
            'n_mutations': n_mut,
            'fluorescence': fitness
        })

    return variants

def score_variants_esm(tokenizer, model, wt_seq, variants):
    """Score variants using ESM-2 masked marginal approach."""
    model.eval()
    scores = []

    for v in variants:
        total_llr = 0
        for mut_str in v['mutations']:
            wt_aa = mut_str[0]
            mut_aa = mut_str[-1]
            pos = int(mut_str[1:-1]) - 1

            seq_list = list(wt_seq)
            seq_list[pos] = tokenizer.mask_token
            masked_str = ''.join(seq_list)
            inputs = tokenizer(masked_str, return_tensors='pt')

            with torch.no_grad():
                logits = model(**inputs).logits[0]
            log_probs = torch.log_softmax(logits[pos + 1], dim=-1)

            wt_lp = log_probs[tokenizer.convert_tokens_to_ids(wt_aa)].item()
            mt_lp = log_probs[tokenizer.convert_tokens_to_ids(mut_aa)].item()
            total_llr += (mt_lp - wt_lp)

        scores.append(total_llr)

    return scores

def in_silico_directed_evolution(tokenizer, model, wt_seq, n_rounds=5, pool_size=20, top_k=5, temperature=0.8):
    """Simulate directed evolution guided by ESM-2 scoring."""
    model.eval()
    current_seq = wt_seq
    evolution_history = [{'round': 0, 'sequence': current_seq,
                          'fitness': 1.0, 'pll': 0, 'mutations': []}]
    L = len(wt_seq)

    for round_num in range(1, n_rounds + 1):
        candidates = []
        for _ in range(pool_size):
            # Random single mutation
            pos = np.random.randint(0, L)
            seq_list = list(current_seq)
            seq_list[pos] = tokenizer.mask_token
            masked_str = ''.join(seq_list)
            inputs = tokenizer(masked_str, return_tensors='pt')

            with torch.no_grad():
                logits = model(**inputs).logits[0]
            probs = torch.softmax(logits[pos + 1] / temperature, dim=-1)

            # Sample from model distribution
            top_probs, top_idx = torch.topk(probs, 5)
            top_probs = top_probs / top_probs.sum()
            chosen = top_idx[torch.multinomial(top_probs, 1).item()]
            new_aa = tokenizer.decode(chosen.item()).strip()

            if new_aa not in AA_LIST or new_aa == current_seq[pos]:
                continue

            new_seq = list(current_seq)
            new_seq[pos] = new_aa
            new_seq_str = ''.join(new_seq)

            # Score the candidate
            wt_lp = torch.log_softmax(logits[pos + 1], dim=-1)[
                tokenizer.convert_tokens_to_ids(current_seq[pos])].item()
            mt_lp = torch.log_softmax(logits[pos + 1], dim=-1)[chosen].item()

            fitness = simulate_gfp_fitness(wt_seq, new_seq_str)
            candidates.append({
                'sequence': new_seq_str,
                'mutation': f"{current_seq[pos]}{pos+1}{new_aa}",
                'llr': float(mt_lp - wt_lp),
                'fitness': fitness
            })

        if not candidates:
            continue

        # Select best by ESM score
        candidates.sort(key=lambda x: x['llr'], reverse=True)
        best = candidates[0]
        current_seq = best['sequence']
        evolution_history.append({
            'round': round_num,
            'sequence': current_seq,
            'fitness': best['fitness'],
            'mutation': best['mutation'],
            'llr': best['llr']
        })

    return evolution_history

def main():
    start = time.time()
    print("=" * 60)
    print("Task 6: GFP Fluorescence Optimization Case Study")
    print("=" * 60)

    print(f"\nModel: {MODEL_NAME}")
    print(f"GFP WT: {GFP_WT[:30]}... (length={len(GFP_WT)})")

    print("\n[1/4] Loading model...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME)

    print("[2/4] Generating & scoring GFP variants...")
    variants = generate_gfp_variants(GFP_WT, n_variants=150, max_mutations=3)
    esm_scores = score_variants_esm(tokenizer, model, GFP_WT, variants)
    for v, s in zip(variants, esm_scores):
        v['esm_score'] = float(s)

    fluorescence_vals = [v['fluorescence'] for v in variants]
    sp_rho, sp_p = spearmanr(esm_scores, fluorescence_vals)
    print(f"  ESM score vs Fluorescence: Spearman ρ = {sp_rho:.3f} (p={sp_p:.2e})")

    print("[3/4] Running in silico directed evolution...")
    evolution = in_silico_directed_evolution(tokenizer, model, GFP_WT,
                                             n_rounds=5, pool_size=15, top_k=3)
    for e in evolution:
        print(f"  Round {e['round']}: fitness={e['fitness']:.3f}"
              + (f" mutation={e.get('mutation', 'N/A')}" if e['round'] > 0 else ""))

    print("[4/4] Creating visualizations...")

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. ESM score vs fluorescence
    n_muts = [v['n_mutations'] for v in variants]
    scatter = axes[0, 0].scatter(esm_scores, fluorescence_vals, c=n_muts,
                                  cmap='viridis', alpha=0.6, s=40, edgecolors='white', linewidths=0.5)
    plt.colorbar(scatter, ax=axes[0, 0], label='# Mutations')
    z = np.polyfit(esm_scores, fluorescence_vals, 1)
    p = np.poly1d(z)
    x_range = np.linspace(min(esm_scores), max(esm_scores), 100)
    axes[0, 0].plot(x_range, p(x_range), 'r--', linewidth=2)
    axes[0, 0].set_xlabel('ESM-2 Score (sum LLR)')
    axes[0, 0].set_ylabel('Fluorescence (simulated)')
    axes[0, 0].set_title(f'ESM-2 Score vs GFP Fluorescence\nSpearman ρ = {sp_rho:.3f}')

    # 2. Directed evolution trajectory
    rounds = [e['round'] for e in evolution]
    fitnesses = [e['fitness'] for e in evolution]
    axes[0, 1].plot(rounds, fitnesses, 'go-', linewidth=2, markersize=10)
    for e in evolution:
        if 'mutation' in e and e['round'] > 0:
            axes[0, 1].annotate(e['mutation'], (e['round'], e['fitness']),
                                textcoords="offset points", xytext=(5, 10), fontsize=8)
    axes[0, 1].set_xlabel('Evolution Round')
    axes[0, 1].set_ylabel('Fitness (Fluorescence)')
    axes[0, 1].set_title('In Silico Directed Evolution\n(ESM-2 guided)')
    axes[0, 1].set_xticks(rounds)

    # 3. Fluorescence by mutation count
    mut_groups = {}
    for v in variants:
        n = v['n_mutations']
        if n not in mut_groups:
            mut_groups[n] = []
        mut_groups[n].append(v['fluorescence'])
    positions = sorted(mut_groups.keys())
    bp_data = [mut_groups[p] for p in positions]
    bp = axes[1, 0].boxplot(bp_data, positions=positions, widths=0.6, patch_artist=True)
    colors_bp = plt.cm.viridis(np.linspace(0.2, 0.8, len(positions)))
    for patch, color in zip(bp['boxes'], colors_bp):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    axes[1, 0].set_xlabel('Number of Mutations')
    axes[1, 0].set_ylabel('Fluorescence')
    axes[1, 0].set_title('Fluorescence vs Mutation Count')
    axes[1, 0].axhline(y=1.0, color='red', linestyle='--', label='WT level', alpha=0.5)
    axes[1, 0].legend()

    # 4. Position importance (avg effect per position)
    pos_effects = {}
    for v in variants:
        for mut in v['mutations']:
            pos = int(mut[1:-1]) - 1
            if pos not in pos_effects:
                pos_effects[pos] = []
            pos_effects[pos].append(v['fluorescence'])

    L = len(GFP_WT)
    avg_effects = np.ones(L)
    for pos, effects in pos_effects.items():
        if pos < L:
            avg_effects[pos] = np.mean(effects)

    bar_colors = ['#d62728' if i in CHROMOPHORE_POSITIONS else
                  '#ff7f0e' if i in BETA_BARREL else '#2ca02c' for i in range(L)]
    axes[1, 1].bar(range(L), avg_effects, color=bar_colors, alpha=0.7)
    axes[1, 1].axhline(y=1.0, color='black', linestyle='--', alpha=0.5)
    axes[1, 1].set_xlabel('Residue Position')
    axes[1, 1].set_ylabel('Mean Fluorescence of Mutants')
    axes[1, 1].set_title('Position-wise Mutation Tolerance\n(Red=chromophore, Orange=β-barrel, Green=surface)')
    axes[1, 1].tick_params(axis='x', labelsize=6)

    plt.suptitle('GFP Fluorescence Optimization with ESM-2', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task6_gfp_optimization.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Summary table
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis('off')
    summary_data = [
        ['Metric', 'Value'],
        ['Model', MODEL_NAME],
        ['GFP sequence length', str(len(GFP_WT))],
        ['N variants scored', str(len(variants))],
        ['Spearman ρ (score vs fluorescence)', f'{sp_rho:.3f}'],
        ['Directed evolution rounds', str(len(evolution) - 1)],
        ['Final evolved fitness', f'{evolution[-1]["fitness"]:.3f}'],
        ['Mean fluorescence (1-mut)', f'{np.mean(mut_groups.get(1, [0])):.3f}'],
        ['Mean fluorescence (3-mut)', f'{np.mean(mut_groups.get(3, [0])):.3f}'],
    ]
    table = ax.table(cellText=summary_data[1:], colLabels=summary_data[0],
                     loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    plt.title('GFP Optimization Summary', fontsize=14, pad=20)
    plt.savefig(f'{FIGURES_DIR}/task6_summary_table.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save results
    results = {
        'model': MODEL_NAME,
        'gfp_sequence_length': len(GFP_WT),
        'n_variants': len(variants),
        'correlation': {
            'spearman_rho': float(sp_rho),
            'spearman_pvalue': float(sp_p)
        },
        'variant_statistics': {
            'mean_fluorescence': float(np.mean(fluorescence_vals)),
            'std_fluorescence': float(np.std(fluorescence_vals)),
            'fraction_functional': float(np.mean([1 for f in fluorescence_vals if f > 0.5]))
        },
        'directed_evolution': evolution,
        'fluorescence_by_n_mutations': {
            str(k): {'mean': float(np.mean(v)), 'std': float(np.std(v)), 'n': len(v)}
            for k, v in mut_groups.items()
        },
        'elapsed_seconds': round(time.time() - start, 2)
    }

    with open(f'{RESULTS_DIR}/task6_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved. Elapsed: {results['elapsed_seconds']:.1f}s")

if __name__ == '__main__':
    main()
