"""
Task 3: Mutation Effect Prediction using Deep Mutational Scanning (DMS) data
- ESM-2 log-likelihood ratio scoring
- Marginal / masked marginal approach
- Correlation with fitness
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
from scipy.stats import spearmanr, pearsonr

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

# Simulated wildtype sequence (GB1-like, 56 residues)
WT_SEQUENCE = "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

def generate_dms_data(wt_seq, n_mutations=300):
    """Simulate a DMS dataset with single-point mutations and fitness scores."""
    mutations = []
    L = len(wt_seq)
    positions_seen = set()
    for _ in range(n_mutations):
        pos = np.random.randint(0, L)
        wt_aa = wt_seq[pos]
        mut_aas = [a for a in AA_LIST if a != wt_aa]
        mut_aa = np.random.choice(mut_aas)
        key = f"{wt_aa}{pos+1}{mut_aa}"
        if key in positions_seen:
            continue
        positions_seen.add(key)

        # Simulate fitness: conservative mutations = less effect
        blosum_like_penalty = abs(AA_LIST.index(wt_aa) - AA_LIST.index(mut_aa)) * 0.05
        position_importance = np.sin(pos / L * np.pi) * 0.5  # Central residues more important
        fitness = 1.0 - blosum_like_penalty - position_importance + np.random.normal(0, 0.2)
        mutations.append({
            'mutation': key,
            'position': pos,
            'wt_aa': wt_aa,
            'mut_aa': mut_aa,
            'fitness': float(np.clip(fitness, -2, 2))
        })
    return mutations

def compute_masked_marginal_scores(tokenizer, model, wt_seq, mutations):
    """Compute masked marginal log-likelihood ratio for each mutation."""
    model.eval()
    scores = []

    for mut in mutations:
        pos = mut['position']
        # Mask the position
        seq_list = list(wt_seq)
        masked_seq = seq_list.copy()
        masked_seq[pos] = '<mask>'
        masked_str = ''.join(masked_seq).replace('<mask>', tokenizer.mask_token)

        inputs = tokenizer(masked_str, return_tensors='pt')
        with torch.no_grad():
            outputs = model(**inputs)
        logits = outputs.logits[0]

        # Token position (account for CLS token)
        token_pos = pos + 1
        log_probs = torch.log_softmax(logits[token_pos], dim=-1)

        wt_token_id = tokenizer.convert_tokens_to_ids(mut['wt_aa'])
        mut_token_id = tokenizer.convert_tokens_to_ids(mut['mut_aa'])

        wt_logprob = log_probs[wt_token_id].item()
        mut_logprob = log_probs[mut_token_id].item()

        score = mut_logprob - wt_logprob  # Positive = mutation favored
        scores.append(score)

    return scores

def compute_wildtype_loglikelihood(tokenizer, model, wt_seq):
    """Compute pseudo-log-likelihood of the wild-type sequence."""
    model.eval()
    total_ll = 0
    for pos in range(len(wt_seq)):
        seq_list = list(wt_seq)
        seq_list[pos] = tokenizer.mask_token
        masked_str = ''.join(seq_list)
        inputs = tokenizer(masked_str, return_tensors='pt')
        with torch.no_grad():
            logits = model(**inputs).logits[0]
        log_probs = torch.log_softmax(logits[pos + 1], dim=-1)
        true_token_id = tokenizer.convert_tokens_to_ids(wt_seq[pos])
        total_ll += log_probs[true_token_id].item()
    return total_ll

def main():
    start = time.time()
    print("=" * 60)
    print("Task 3: Mutation Effect Prediction (DMS)")
    print("=" * 60)

    print(f"\nModel: {MODEL_NAME}")
    print(f"WT Sequence length: {len(WT_SEQUENCE)}")

    print("\n[1/4] Loading model...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME)
    model.eval()

    print("[2/4] Generating synthetic DMS data...")
    mutations = generate_dms_data(WT_SEQUENCE, n_mutations=150)
    print(f"  Generated {len(mutations)} unique mutations")

    print("[3/4] Computing masked marginal scores...")
    scores = compute_masked_marginal_scores(tokenizer, model, WT_SEQUENCE, mutations)

    # Compute WT pseudo-log-likelihood
    wt_pll = compute_wildtype_loglikelihood(tokenizer, model, WT_SEQUENCE)
    print(f"  WT pseudo-log-likelihood: {wt_pll:.2f}")

    print("[4/4] Analyzing correlations...")
    fitness_values = [m['fitness'] for m in mutations]

    spearman_r, spearman_p = spearmanr(scores, fitness_values)
    pearson_r, pearson_p = pearsonr(scores, fitness_values)

    print(f"  Spearman ρ: {spearman_r:.3f} (p={spearman_p:.2e})")
    print(f"  Pearson r:  {pearson_r:.3f} (p={pearson_p:.2e})")

    # --- Visualizations ---
    # 1. Score vs Fitness scatter
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].scatter(scores, fitness_values, alpha=0.5, s=25, c='steelblue', edgecolors='white', linewidths=0.5)
    z = np.polyfit(scores, fitness_values, 1)
    p = np.poly1d(z)
    x_range = np.linspace(min(scores), max(scores), 100)
    axes[0].plot(x_range, p(x_range), 'r--', linewidth=2, label=f'ρ={spearman_r:.3f}')
    axes[0].set_xlabel('ESM-2 Masked Marginal Score')
    axes[0].set_ylabel('Fitness (simulated)')
    axes[0].set_title('Mutation Effect Prediction')
    axes[0].legend()

    # 2. Position-wise average scores
    positions = [m['position'] for m in mutations]
    pos_scores = {}
    pos_fitness = {}
    for m, s in zip(mutations, scores):
        p = m['position']
        if p not in pos_scores:
            pos_scores[p] = []
            pos_fitness[p] = []
        pos_scores[p].append(s)
        pos_fitness[p].append(m['fitness'])

    pos_list = sorted(pos_scores.keys())
    avg_scores = [np.mean(pos_scores[p]) for p in pos_list]
    avg_fitness = [np.mean(pos_fitness[p]) for p in pos_list]

    axes[1].bar(pos_list, avg_scores, color='steelblue', alpha=0.7, label='ESM-2 Score')
    ax2 = axes[1].twinx()
    ax2.plot(pos_list, avg_fitness, 'ro-', markersize=3, alpha=0.7, label='Fitness')
    axes[1].set_xlabel('Residue Position')
    axes[1].set_ylabel('Avg ESM-2 Score', color='steelblue')
    ax2.set_ylabel('Avg Fitness', color='red')
    axes[1].set_title('Position-wise Mutation Tolerance')
    axes[1].legend(loc='upper left')
    ax2.legend(loc='upper right')

    # 3. Mutation type heatmap (AA substitution matrix)
    sub_matrix = np.full((20, 20), np.nan)
    count_matrix = np.zeros((20, 20))
    for m, s in zip(mutations, scores):
        wi = AA_LIST.index(m['wt_aa'])
        mi = AA_LIST.index(m['mut_aa'])
        if np.isnan(sub_matrix[wi, mi]):
            sub_matrix[wi, mi] = s
        else:
            sub_matrix[wi, mi] = (sub_matrix[wi, mi] * count_matrix[wi, mi] + s) / (count_matrix[wi, mi] + 1)
        count_matrix[wi, mi] += 1

    sns.heatmap(sub_matrix, xticklabels=AA_LIST, yticklabels=AA_LIST,
                cmap='RdBu_r', center=0, ax=axes[2], square=True, linewidths=0.5,
                cbar_kws={'label': 'LLR Score'})
    axes[2].set_xlabel('Mutant AA')
    axes[2].set_ylabel('Wild-type AA')
    axes[2].set_title('AA Substitution Score Matrix')

    plt.suptitle('Deep Mutational Scanning Analysis with ESM-2', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task3_dms_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save results
    results = {
        'model': MODEL_NAME,
        'wt_sequence_length': len(WT_SEQUENCE),
        'n_mutations': len(mutations),
        'wt_pseudo_loglikelihood': float(wt_pll),
        'correlation': {
            'spearman_rho': float(spearman_r),
            'spearman_pvalue': float(spearman_p),
            'pearson_r': float(pearson_r),
            'pearson_pvalue': float(pearson_p)
        },
        'score_statistics': {
            'mean': float(np.mean(scores)),
            'std': float(np.std(scores)),
            'min': float(np.min(scores)),
            'max': float(np.max(scores))
        },
        'elapsed_seconds': round(time.time() - start, 2)
    }

    with open(f'{RESULTS_DIR}/task3_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved. Elapsed: {results['elapsed_seconds']:.1f}s")

if __name__ == '__main__':
    main()
