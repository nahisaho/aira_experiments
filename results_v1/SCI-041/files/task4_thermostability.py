"""
Task 4: Zero-shot Thermostability Prediction
- Pseudo-log-likelihood based stability scoring
- Mutation scanning for stability-enhancing variants
- No training data required (zero-shot)
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

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

# Example: a thermophilic-like protein sequence
WT_SEQUENCE = "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPATKI"
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def compute_sequence_pll(tokenizer, model, sequence):
    """Compute pseudo-log-likelihood by masking each position."""
    model.eval()
    position_logprobs = []
    for pos in range(len(sequence)):
        seq_list = list(sequence)
        seq_list[pos] = tokenizer.mask_token
        masked_str = ''.join(seq_list)
        inputs = tokenizer(masked_str, return_tensors='pt')
        with torch.no_grad():
            logits = model(**inputs).logits[0]
        log_probs = torch.log_softmax(logits[pos + 1], dim=-1)
        true_token = tokenizer.convert_tokens_to_ids(sequence[pos])
        position_logprobs.append(log_probs[true_token].item())
    return sum(position_logprobs), position_logprobs

def scan_all_mutations(tokenizer, model, sequence, positions=None):
    """Scan all possible single mutations and score with masked marginal."""
    model.eval()
    if positions is None:
        positions = range(len(sequence))

    mutation_scores = {}
    for pos in positions:
        wt_aa = sequence[pos]
        seq_list = list(sequence)
        seq_list[pos] = tokenizer.mask_token
        masked_str = ''.join(seq_list)
        inputs = tokenizer(masked_str, return_tensors='pt')

        with torch.no_grad():
            logits = model(**inputs).logits[0]
        log_probs = torch.log_softmax(logits[pos + 1], dim=-1)

        wt_logprob = log_probs[tokenizer.convert_tokens_to_ids(wt_aa)].item()
        for mut_aa in AA_LIST:
            if mut_aa == wt_aa:
                continue
            mut_logprob = log_probs[tokenizer.convert_tokens_to_ids(mut_aa)].item()
            llr = mut_logprob - wt_logprob
            mutation_scores[f"{wt_aa}{pos+1}{mut_aa}"] = {
                'position': pos,
                'wt_aa': wt_aa,
                'mut_aa': mut_aa,
                'llr_score': float(llr),
                'mut_logprob': float(mut_logprob),
                'wt_logprob': float(wt_logprob)
            }
    return mutation_scores

def main():
    start = time.time()
    print("=" * 60)
    print("Task 4: Zero-shot Thermostability Prediction")
    print("=" * 60)

    print(f"\nModel: {MODEL_NAME}")
    print(f"WT Sequence: {WT_SEQUENCE[:30]}... (length={len(WT_SEQUENCE)})")

    print("\n[1/3] Loading model...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME)

    print("[2/3] Computing WT pseudo-log-likelihood...")
    wt_pll, position_logprobs = compute_sequence_pll(tokenizer, model, WT_SEQUENCE)
    print(f"  WT PLL: {wt_pll:.2f}")

    print("[3/3] Scanning all single mutations...")
    mutation_scores = scan_all_mutations(tokenizer, model, WT_SEQUENCE)
    print(f"  Scanned {len(mutation_scores)} mutations")

    # Identify top stabilizing mutations (highest LLR = most favored by model)
    sorted_mutations = sorted(mutation_scores.items(), key=lambda x: x[1]['llr_score'], reverse=True)
    top_stabilizing = sorted_mutations[:20]
    top_destabilizing = sorted_mutations[-20:]

    print("\n  Top 10 predicted stabilizing mutations:")
    for name, info in top_stabilizing[:10]:
        print(f"    {name}: LLR = {info['llr_score']:+.3f}")

    print("\n  Top 10 predicted destabilizing mutations:")
    for name, info in top_destabilizing[:10]:
        print(f"    {name}: LLR = {info['llr_score']:+.3f}")

    # --- Visualizations ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. Position-wise vulnerability
    axes[0, 0].bar(range(len(WT_SEQUENCE)), position_logprobs, color='steelblue', alpha=0.8)
    axes[0, 0].set_xlabel('Residue Position')
    axes[0, 0].set_ylabel('Log P(wt | context)')
    axes[0, 0].set_title('Position-wise WT Log-Probability\n(lower = more conserved/important)')
    axes[0, 0].axhline(y=np.mean(position_logprobs), color='red', linestyle='--', label='Mean')
    axes[0, 0].legend()

    # 2. Mutation landscape heatmap
    L = len(WT_SEQUENCE)
    heatmap_data = np.full((20, L), np.nan)
    for name, info in mutation_scores.items():
        aa_idx = AA_LIST.index(info['mut_aa'])
        heatmap_data[aa_idx, info['position']] = info['llr_score']

    sns.heatmap(heatmap_data, cmap='RdBu_r', center=0, ax=axes[0, 1],
                yticklabels=AA_LIST, xticklabels=[f"{WT_SEQUENCE[i]}{i+1}" for i in range(L)],
                cbar_kws={'label': 'LLR (positive = stabilizing)'})
    axes[0, 1].set_title('Mutation Landscape (Zero-shot Prediction)')
    axes[0, 1].set_xlabel('Position (WT residue)')
    axes[0, 1].set_ylabel('Mutant Amino Acid')
    axes[0, 1].tick_params(axis='x', rotation=90, labelsize=6)

    # 3. LLR score distribution
    all_llrs = [info['llr_score'] for info in mutation_scores.values()]
    axes[1, 0].hist(all_llrs, bins=50, color='steelblue', alpha=0.7, edgecolor='white')
    axes[1, 0].axvline(x=0, color='red', linestyle='--', linewidth=2, label='Neutral')
    axes[1, 0].set_xlabel('Log-Likelihood Ratio')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Distribution of Mutation Scores')
    n_stabilizing = sum(1 for s in all_llrs if s > 0)
    n_destabilizing = sum(1 for s in all_llrs if s <= 0)
    axes[1, 0].legend([f'Neutral line', f'Stabilizing: {n_stabilizing}', f'Destabilizing: {n_destabilizing}'])

    # 4. Top mutations bar chart
    top_names = [name for name, _ in top_stabilizing[:15]]
    top_scores = [info['llr_score'] for _, info in top_stabilizing[:15]]
    bottom_names = [name for name, _ in top_destabilizing[:15]]
    bottom_scores = [info['llr_score'] for _, info in top_destabilizing[:15]]

    combined_names = top_names + bottom_names
    combined_scores = top_scores + bottom_scores
    colors = ['#2ca02c'] * len(top_names) + ['#d62728'] * len(bottom_names)
    sorted_idx = np.argsort(combined_scores)[::-1]
    axes[1, 1].barh(range(len(combined_names)),
                     [combined_scores[i] for i in sorted_idx],
                     color=[colors[i] for i in sorted_idx], alpha=0.8)
    axes[1, 1].set_yticks(range(len(combined_names)))
    axes[1, 1].set_yticklabels([combined_names[i] for i in sorted_idx], fontsize=7)
    axes[1, 1].set_xlabel('LLR Score')
    axes[1, 1].set_title('Top Stabilizing (green) & Destabilizing (red) Mutations')
    axes[1, 1].axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    axes[1, 1].invert_yaxis()

    plt.suptitle('Zero-shot Thermostability Prediction with ESM-2', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task4_thermostability.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save results
    results = {
        'model': MODEL_NAME,
        'wt_sequence_length': len(WT_SEQUENCE),
        'wt_pll': float(wt_pll),
        'n_mutations_scanned': len(mutation_scores),
        'n_stabilizing_predicted': n_stabilizing,
        'n_destabilizing_predicted': n_destabilizing,
        'top_stabilizing': [{name: info} for name, info in top_stabilizing[:10]],
        'top_destabilizing': [{name: info} for name, info in top_destabilizing[:10]],
        'score_statistics': {
            'mean': float(np.mean(all_llrs)),
            'std': float(np.std(all_llrs)),
            'min': float(np.min(all_llrs)),
            'max': float(np.max(all_llrs))
        },
        'elapsed_seconds': round(time.time() - start, 2)
    }

    with open(f'{RESULTS_DIR}/task4_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved. Elapsed: {results['elapsed_seconds']:.1f}s")

if __name__ == '__main__':
    main()
