"""
Task 1: Pre-trained ESM-2 Internal Representation Analysis
- Attention pattern extraction and visualization
- Contact prediction from attention maps
- Layer-wise embedding analysis
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
from transformers import AutoTokenizer, AutoModel, EsmForMaskedLM
from scipy.spatial.distance import squareform, pdist
from sklearn.metrics import precision_score, average_precision_score

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"  # Small model for demo

# Example protein: lysozyme fragment (first 50 residues for speed)
SEQUENCE = "MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRA"

def load_model():
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForMaskedLM.from_pretrained(MODEL_NAME, output_attentions=True, output_hidden_states=True)
    model.eval()
    return tokenizer, model

def extract_representations(tokenizer, model, sequence):
    inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
    with torch.no_grad():
        outputs = model(**inputs)
    return outputs, inputs

def analyze_attention_patterns(outputs, sequence):
    """Extract and visualize attention patterns across layers and heads."""
    attentions = outputs.attentions  # tuple of (batch, heads, seq_len, seq_len)
    n_layers = len(attentions)
    n_heads = attentions[0].shape[1]
    seq_len = len(sequence)

    # Remove special tokens (CLS, EOS)
    attention_matrices = []
    for layer_att in attentions:
        att = layer_att[0, :, 1:seq_len+1, 1:seq_len+1].numpy()
        attention_matrices.append(att)
    attention_matrices = np.array(attention_matrices)  # (layers, heads, L, L)

    # 1. Average attention per layer
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    for i, ax in enumerate(axes.flat):
        if i < n_layers:
            avg_att = attention_matrices[i].mean(axis=0)
            im = ax.imshow(avg_att, cmap='viridis', aspect='auto')
            ax.set_title(f'Layer {i+1} (avg over heads)', fontsize=12)
            ax.set_xlabel('Position')
            ax.set_ylabel('Position')
            plt.colorbar(im, ax=ax, fraction=0.046)
    plt.suptitle('ESM-2 Attention Patterns by Layer', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task1_attention_layers.png', dpi=150, bbox_inches='tight')
    plt.close()

    # 2. Attention head diversity (entropy per head)
    head_entropies = []
    for l in range(n_layers):
        for h in range(n_heads):
            att = attention_matrices[l, h]
            att_clipped = np.clip(att, 1e-10, 1.0)
            entropy = -np.sum(att_clipped * np.log(att_clipped), axis=-1).mean()
            head_entropies.append({'layer': l+1, 'head': h+1, 'entropy': float(entropy)})

    fig, ax = plt.subplots(figsize=(10, 6))
    entropy_matrix = np.zeros((n_layers, n_heads))
    for entry in head_entropies:
        entropy_matrix[entry['layer']-1, entry['head']-1] = entry['entropy']
    sns.heatmap(entropy_matrix, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
                xticklabels=[f'H{i+1}' for i in range(n_heads)],
                yticklabels=[f'L{i+1}' for i in range(n_layers)])
    ax.set_title('Attention Head Entropy (higher = more distributed)', fontsize=14)
    ax.set_xlabel('Head')
    ax.set_ylabel('Layer')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task1_head_entropy.png', dpi=150, bbox_inches='tight')
    plt.close()

    return attention_matrices, head_entropies

def contact_prediction_from_attention(attention_matrices, sequence):
    """Predict residue contacts from attention (symmetrized, APC-corrected)."""
    n_layers, n_heads, L, _ = attention_matrices.shape

    # Symmetrize attention
    sym_att = (attention_matrices + attention_matrices.transpose(0, 1, 3, 2)) / 2

    # Average across all heads and layers
    contact_score = sym_att.mean(axis=(0, 1))

    # APC correction
    mean_i = contact_score.mean(axis=1, keepdims=True)
    mean_j = contact_score.mean(axis=0, keepdims=True)
    mean_all = contact_score.mean()
    apc = (mean_i * mean_j) / (mean_all + 1e-10)
    corrected = contact_score - apc
    np.fill_diagonal(corrected, 0)

    # Simulate "true" contacts using sequence distance heuristic
    # (In real work, use PDB structures)
    true_contacts = np.zeros((L, L))
    for i in range(L):
        for j in range(L):
            if abs(i - j) >= 6:
                # Simulate contacts: some residues are "close" in 3D
                np.random.seed(i * L + j)
                if np.random.random() < 0.05:
                    true_contacts[i, j] = 1
                    true_contacts[j, i] = 1

    # Evaluate top-L predictions
    mask = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :]) >= 6
    scores_filtered = corrected * mask
    true_filtered = true_contacts * mask

    upper_idx = np.triu_indices(L, k=6)
    pred_scores = scores_filtered[upper_idx]
    true_labels = true_filtered[upper_idx]

    top_L = min(L, len(pred_scores))
    top_indices = np.argsort(pred_scores)[::-1][:top_L]
    precision_at_L = true_labels[top_indices].mean() if len(top_indices) > 0 else 0

    # Visualize contact map
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].imshow(corrected, cmap='hot', aspect='auto')
    axes[0].set_title('Predicted Contact Scores\n(APC-corrected attention)', fontsize=12)

    axes[1].imshow(true_contacts, cmap='binary', aspect='auto')
    axes[1].set_title('Simulated True Contacts', fontsize=12)

    # Overlay
    axes[2].imshow(true_contacts, cmap='Greys', alpha=0.3, aspect='auto')
    top_pred = np.zeros_like(corrected)
    flat_idx = np.argsort(corrected.ravel())[::-1][:L]
    for idx in flat_idx:
        i, j = divmod(idx, L)
        if abs(i - j) >= 6:
            top_pred[i, j] = 1
    axes[2].imshow(top_pred, cmap='Reds', alpha=0.5, aspect='auto')
    axes[2].set_title(f'Top-L Predictions Overlay\nPrecision@L = {precision_at_L:.3f}', fontsize=12)

    for ax in axes:
        ax.set_xlabel('Residue Position')
        ax.set_ylabel('Residue Position')

    plt.suptitle('Contact Prediction from ESM-2 Attention Maps', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task1_contact_prediction.png', dpi=150, bbox_inches='tight')
    plt.close()

    return {'precision_at_L': float(precision_at_L), 'n_true_contacts': int(true_contacts.sum()//2)}

def analyze_hidden_states(outputs, sequence):
    """Analyze layer-wise hidden state representations."""
    hidden_states = outputs.hidden_states  # tuple of (batch, seq_len, hidden_dim)
    n_layers = len(hidden_states)
    seq_len = len(sequence)

    # Compute CLS token representation similarity across layers
    cls_embeddings = []
    for hs in hidden_states:
        cls_embeddings.append(hs[0, 0].numpy())
    cls_embeddings = np.array(cls_embeddings)

    # Layer-wise cosine similarity
    from sklearn.metrics.pairwise import cosine_similarity
    layer_sim = cosine_similarity(cls_embeddings)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    im = axes[0].imshow(layer_sim, cmap='coolwarm', vmin=0, vmax=1, aspect='auto')
    axes[0].set_title('Layer-wise CLS Embedding Similarity', fontsize=12)
    axes[0].set_xlabel('Layer')
    axes[0].set_ylabel('Layer')
    plt.colorbar(im, ax=axes[0])

    # Per-residue embedding norm across layers
    norms = []
    for hs in hidden_states:
        residue_norms = torch.norm(hs[0, 1:seq_len+1], dim=-1).numpy()
        norms.append(residue_norms)
    norms = np.array(norms)

    im2 = axes[1].imshow(norms, cmap='viridis', aspect='auto')
    axes[1].set_title('Per-residue Embedding Norm by Layer', fontsize=12)
    axes[1].set_xlabel('Residue Position')
    axes[1].set_ylabel('Layer')
    plt.colorbar(im2, ax=axes[1])

    plt.suptitle('ESM-2 Hidden State Analysis', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task1_hidden_states.png', dpi=150, bbox_inches='tight')
    plt.close()

    return {
        'n_layers': n_layers,
        'hidden_dim': int(hidden_states[0].shape[-1]),
        'mean_cls_similarity': float(layer_sim.mean()),
        'mean_norm_final_layer': float(norms[-1].mean())
    }

def main():
    start = time.time()
    print("=" * 60)
    print("Task 1: ESM-2 Internal Representation Analysis")
    print("=" * 60)

    print(f"\nModel: {MODEL_NAME}")
    print(f"Sequence length: {len(SEQUENCE)}")

    print("\n[1/4] Loading model...")
    tokenizer, model = load_model()

    print("[2/4] Extracting representations...")
    outputs, inputs = extract_representations(tokenizer, model, SEQUENCE)

    print("[3/4] Analyzing attention patterns...")
    attention_matrices, head_entropies = analyze_attention_patterns(outputs, SEQUENCE)

    print("[4/4] Contact prediction & hidden state analysis...")
    contact_results = contact_prediction_from_attention(attention_matrices, SEQUENCE)
    hidden_results = analyze_hidden_states(outputs, SEQUENCE)

    results = {
        'model': MODEL_NAME,
        'sequence_length': len(SEQUENCE),
        'attention_analysis': {
            'n_layers': len(outputs.attentions),
            'n_heads': outputs.attentions[0].shape[1],
            'head_entropies': head_entropies
        },
        'contact_prediction': contact_results,
        'hidden_state_analysis': hidden_results,
        'elapsed_seconds': round(time.time() - start, 2)
    }

    with open(f'{RESULTS_DIR}/task1_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to {RESULTS_DIR}/task1_results.json")
    print(f"Figures saved to {FIGURES_DIR}/task1_*.png")
    print(f"Contact precision@L: {contact_results['precision_at_L']:.3f}")
    print(f"Hidden dim: {hidden_results['hidden_dim']}, Layers: {hidden_results['n_layers']}")
    print(f"Elapsed: {results['elapsed_seconds']:.1f}s")

if __name__ == '__main__':
    main()
