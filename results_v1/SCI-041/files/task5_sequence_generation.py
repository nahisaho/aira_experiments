"""
Task 5: Sequence Generation
- Conditional sequence generation using masked language model
- Iterative refinement strategy
- Sequence diversity analysis
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
from collections import Counter

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
TEMPLATE_SEQUENCE = "MRVLKFGGTSVANAERFLRVADILESNARQGQVATVLSAPATKI"

def generate_from_scratch(tokenizer, model, length=30, n_iter=5, temperature=1.0):
    """Generate sequence from fully masked template using iterative unmasking."""
    model.eval()
    # Start with all masks
    seq = [tokenizer.mask_token] * length
    history = []

    for iteration in range(n_iter):
        n_to_unmask = max(1, length // n_iter)
        masked_positions = [i for i, c in enumerate(seq) if c == tokenizer.mask_token]
        if not masked_positions:
            break

        # Score all masked positions
        seq_str = ''.join(seq)
        inputs = tokenizer(seq_str, return_tensors='pt')
        with torch.no_grad():
            logits = model(**inputs).logits[0]

        # Find positions with highest confidence
        confidences = []
        for pos in masked_positions:
            probs = torch.softmax(logits[pos + 1] / temperature, dim=-1)
            max_prob = probs.max().item()
            confidences.append((pos, max_prob))

        # Unmask the most confident positions
        confidences.sort(key=lambda x: x[1], reverse=True)
        positions_to_unmask = [p for p, _ in confidences[:n_to_unmask]]

        for pos in positions_to_unmask:
            probs = torch.softmax(logits[pos + 1] / temperature, dim=-1)
            # Sample from top-k
            top_k = 5
            top_probs, top_indices = torch.topk(probs, top_k)
            top_probs = top_probs / top_probs.sum()
            chosen = top_indices[torch.multinomial(top_probs, 1).item()]
            seq[pos] = tokenizer.decode(chosen.item()).strip()

        current_seq = ''.join([c if c != tokenizer.mask_token else '_' for c in seq])
        n_masked = sum(1 for c in seq if c == tokenizer.mask_token)
        history.append({'iteration': iteration + 1, 'sequence': current_seq, 'n_masked': n_masked})

    final_seq = ''.join([c if c != tokenizer.mask_token else 'X' for c in seq])
    return final_seq, history

def conditional_generation(tokenizer, model, template, mask_fraction=0.3, n_samples=20, temperature=0.8):
    """Conditionally generate sequences from a template with partial masking."""
    model.eval()
    L = len(template)
    n_mask = max(1, int(L * mask_fraction))
    generated = []

    for sample_idx in range(n_samples):
        np.random.seed(SEED + sample_idx)
        mask_positions = np.random.choice(L, n_mask, replace=False)
        seq_list = list(template)
        for pos in mask_positions:
            seq_list[pos] = tokenizer.mask_token
        seq_str = ''.join(seq_list)

        inputs = tokenizer(seq_str, return_tensors='pt')
        with torch.no_grad():
            logits = model(**inputs).logits[0]

        for pos in mask_positions:
            probs = torch.softmax(logits[pos + 1] / temperature, dim=-1)
            top_k = 3
            top_probs, top_indices = torch.topk(probs, top_k)
            top_probs = top_probs / top_probs.sum()
            chosen = top_indices[torch.multinomial(top_probs, 1).item()]
            seq_list[pos] = tokenizer.decode(chosen.item()).strip()

        gen_seq = ''.join(seq_list)
        # Sequence identity to template
        identity = sum(1 for a, b in zip(gen_seq, template) if a == b) / L
        generated.append({
            'sequence': gen_seq,
            'identity': float(identity),
            'n_mutations': int(L - sum(1 for a, b in zip(gen_seq, template) if a == b))
        })

    return generated

def analyze_diversity(sequences, template):
    """Analyze amino acid composition and diversity of generated sequences."""
    L = len(template)
    # Position-wise AA frequency
    pos_freq = {pos: Counter() for pos in range(L)}
    for seq_info in sequences:
        for pos, aa in enumerate(seq_info['sequence']):
            pos_freq[pos][aa] += 1

    # Diversity: Shannon entropy per position
    entropies = []
    for pos in range(L):
        total = sum(pos_freq[pos].values())
        if total == 0:
            entropies.append(0)
            continue
        h = 0
        for count in pos_freq[pos].values():
            p = count / total
            if p > 0:
                h -= p * np.log2(p)
        entropies.append(h)

    # Overall AA composition
    all_aas = ''.join([s['sequence'] for s in sequences])
    aa_counts = Counter(all_aas)
    total_aa = len(all_aas)

    return entropies, aa_counts, total_aa

def compute_sequence_pll(tokenizer, model, sequence):
    """Pseudo-log-likelihood for generated sequences."""
    model.eval()
    total_ll = 0
    for pos in range(len(sequence)):
        seq_list = list(sequence)
        seq_list[pos] = tokenizer.mask_token
        masked_str = ''.join(seq_list)
        inputs = tokenizer(masked_str, return_tensors='pt')
        with torch.no_grad():
            logits = model(**inputs).logits[0]
        log_probs = torch.log_softmax(logits[pos + 1], dim=-1)
        token_id = tokenizer.convert_tokens_to_ids(sequence[pos])
        total_ll += log_probs[token_id].item()
    return total_ll

def main():
    start = time.time()
    print("=" * 60)
    print("Task 5: Sequence Generation")
    print("=" * 60)

    print(f"\nModel: {MODEL_NAME}")
    print(f"Template: {TEMPLATE_SEQUENCE[:30]}... (length={len(TEMPLATE_SEQUENCE)})")

    print("\n[1/4] Loading model...")
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME)

    print("[2/4] De novo generation (iterative unmasking)...")
    denovo_seq, denovo_history = generate_from_scratch(tokenizer, model, length=30, n_iter=6, temperature=0.9)
    print(f"  Generated: {denovo_seq}")

    print("[3/4] Conditional generation (partial masking)...")
    generated_seqs = conditional_generation(tokenizer, model, TEMPLATE_SEQUENCE,
                                            mask_fraction=0.3, n_samples=20, temperature=0.8)
    identities = [s['identity'] for s in generated_seqs]
    print(f"  Generated {len(generated_seqs)} variants")
    print(f"  Mean identity: {np.mean(identities):.3f} ± {np.std(identities):.3f}")

    print("[4/4] Analyzing diversity & quality...")
    entropies, aa_counts, total_aa = analyze_diversity(generated_seqs, TEMPLATE_SEQUENCE)

    # Compute PLL for top 5 generated + template
    print("  Computing PLL for template and top-5 variants...")
    template_pll = compute_sequence_pll(tokenizer, model, TEMPLATE_SEQUENCE)
    gen_plls = []
    for seq_info in generated_seqs[:5]:
        pll = compute_sequence_pll(tokenizer, model, seq_info['sequence'])
        gen_plls.append({'sequence': seq_info['sequence'][:20] + '...', 'pll': float(pll),
                         'identity': seq_info['identity']})
    print(f"  Template PLL: {template_pll:.2f}")

    # --- Visualizations ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1. De novo generation process
    iters = [h['iteration'] for h in denovo_history]
    n_masked = [h['n_masked'] for h in denovo_history]
    axes[0, 0].plot(iters, n_masked, 'bo-', linewidth=2, markersize=8)
    axes[0, 0].set_xlabel('Iteration')
    axes[0, 0].set_ylabel('Masked Positions Remaining')
    axes[0, 0].set_title('De Novo Generation: Iterative Unmasking')
    axes[0, 0].fill_between(iters, n_masked, alpha=0.2, color='blue')

    # 2. Position-wise entropy
    axes[0, 1].bar(range(len(TEMPLATE_SEQUENCE)), entropies, color='teal', alpha=0.7)
    axes[0, 1].set_xlabel('Residue Position')
    axes[0, 1].set_ylabel('Shannon Entropy (bits)')
    axes[0, 1].set_title('Position-wise Diversity\n(Conditional Generation, 30% masking)')
    axes[0, 1].axhline(y=np.mean(entropies), color='red', linestyle='--', label=f'Mean={np.mean(entropies):.2f}')
    axes[0, 1].legend()

    # 3. AA composition comparison
    template_counts = Counter(TEMPLATE_SEQUENCE)
    gen_fracs = {aa: aa_counts.get(aa, 0) / total_aa for aa in AA_LIST}
    template_fracs = {aa: template_counts.get(aa, 0) / len(TEMPLATE_SEQUENCE) for aa in AA_LIST}
    x = np.arange(len(AA_LIST))
    w = 0.35
    axes[1, 0].bar(x - w/2, [template_fracs[aa] for aa in AA_LIST], w, label='Template', color='steelblue', alpha=0.7)
    axes[1, 0].bar(x + w/2, [gen_fracs[aa] for aa in AA_LIST], w, label='Generated', color='coral', alpha=0.7)
    axes[1, 0].set_xticks(x)
    axes[1, 0].set_xticklabels(AA_LIST)
    axes[1, 0].set_xlabel('Amino Acid')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('AA Composition: Template vs Generated')
    axes[1, 0].legend()

    # 4. Identity distribution
    axes[1, 1].hist(identities, bins=10, color='teal', alpha=0.7, edgecolor='white')
    axes[1, 1].axvline(x=np.mean(identities), color='red', linestyle='--',
                       label=f'Mean={np.mean(identities):.3f}')
    axes[1, 1].set_xlabel('Sequence Identity to Template')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title('Generated Sequence Identity Distribution')
    axes[1, 1].legend()

    plt.suptitle('Protein Sequence Generation with ESM-2 MLM', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task5_generation.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save results
    results = {
        'model': MODEL_NAME,
        'denovo_generation': {
            'sequence': denovo_seq,
            'length': len(denovo_seq),
            'n_iterations': len(denovo_history),
            'history': denovo_history
        },
        'conditional_generation': {
            'n_samples': len(generated_seqs),
            'mask_fraction': 0.3,
            'mean_identity': float(np.mean(identities)),
            'std_identity': float(np.std(identities)),
            'top_5_plls': gen_plls,
            'template_pll': float(template_pll)
        },
        'diversity': {
            'mean_positional_entropy': float(np.mean(entropies)),
            'max_positional_entropy': float(np.max(entropies))
        },
        'elapsed_seconds': round(time.time() - start, 2)
    }

    with open(f'{RESULTS_DIR}/task5_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved. Elapsed: {results['elapsed_seconds']:.1f}s")

if __name__ == '__main__':
    main()
