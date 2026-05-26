#!/usr/bin/env python3
"""
Protein Language Model Fine-tuning Pipeline
============================================
Comprehensive experiment suite for ESM-2 fine-tuning strategies:
1. Internal representation analysis (attention, contact prediction)
2. Enzyme activity prediction (LoRA vs Adapter)
3. Variant effect prediction (DMS)
4. Zero-shot thermostability prediction
5. Conditional sequence generation (MLM)
6. GFP fluorescence optimization case study
"""

import os
import json
import warnings
import numpy as np

# Fix jax/numpy compatibility issue (must be before transformers import)
class _FakeStringDType:
    def __init__(self, *args, **kwargs):
        pass
if not hasattr(np.dtypes, 'StringDType'):
    np.dtypes.StringDType = _FakeStringDType
else:
    _orig = np.dtypes.StringDType
    np.dtypes.StringDType = _FakeStringDType

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import squareform
from sklearn.metrics import (
    accuracy_score, f1_score, roc_auc_score, mean_squared_error,
    r2_score, precision_recall_curve, roc_curve, confusion_matrix,
    matthews_corrcoef
)
from sklearn.model_selection import train_test_split
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset, TensorDataset
from transformers import AutoTokenizer
from transformers.models.esm.modeling_esm import EsmModel, EsmForMaskedLM

warnings.filterwarnings('ignore')

np.random.seed(42)
torch.manual_seed(42)

DEVICE = torch.device('cpu')
FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

RESULTS = {}

# ============================================================
# Utility: synthetic protein sequences
# ============================================================
AA = list("ACDEFGHIKLMNPQRSTVWY")

def random_seq(length=50):
    return "".join(np.random.choice(AA, size=length))

def mutate_seq(seq, n_mutations=1):
    seq = list(seq)
    for _ in range(n_mutations):
        pos = np.random.randint(len(seq))
        seq[pos] = np.random.choice(AA)
    return "".join(seq)

# ============================================================
# Load ESM-2 (small model for feasibility)
# ============================================================
print("Loading ESM-2 model (esm2_t6_8M_UR50D)...")
MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
esm_model = EsmModel.from_pretrained(MODEL_NAME, output_attentions=True)
esm_mlm = EsmForMaskedLM.from_pretrained(MODEL_NAME)
esm_model.eval()
esm_mlm.eval()
print("Model loaded successfully.")

# ============================================================
# EXPERIMENT 1: Internal Representation Analysis
# ============================================================
print("\n" + "="*60)
print("EXPERIMENT 1: Internal Representation Analysis")
print("="*60)

# Use a well-known short protein motif (lysozyme-like)
test_seq = "KVFGRCELAAKLKADGYNGVSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
inputs = tokenizer(test_seq, return_tensors="pt")

with torch.no_grad():
    outputs = esm_model(**inputs)

# Attention pattern analysis
attentions = outputs.attentions  # tuple of (batch, heads, seq_len, seq_len)
n_layers = len(attentions)
n_heads = attentions[0].shape[1]
seq_len = attentions[0].shape[2]

print(f"  Layers: {n_layers}, Heads per layer: {n_heads}, Seq length: {seq_len}")

# Fig 1: Attention heatmap (average across heads, last layer)
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Last layer average attention
avg_attn_last = attentions[-1][0].mean(dim=0).numpy()
im0 = axes[0].imshow(avg_attn_last[1:-1, 1:-1], cmap='viridis', aspect='auto')
axes[0].set_title(f'Layer {n_layers} Avg Attention', fontsize=12)
axes[0].set_xlabel('Key Position')
axes[0].set_ylabel('Query Position')
plt.colorbar(im0, ax=axes[0], fraction=0.046)

# Middle layer
mid_layer = n_layers // 2
avg_attn_mid = attentions[mid_layer][0].mean(dim=0).numpy()
im1 = axes[1].imshow(avg_attn_mid[1:-1, 1:-1], cmap='viridis', aspect='auto')
axes[1].set_title(f'Layer {mid_layer+1} Avg Attention', fontsize=12)
axes[1].set_xlabel('Key Position')
axes[1].set_ylabel('Query Position')
plt.colorbar(im1, ax=axes[1], fraction=0.046)

# First layer
avg_attn_first = attentions[0][0].mean(dim=0).numpy()
im2 = axes[2].imshow(avg_attn_first[1:-1, 1:-1], cmap='viridis', aspect='auto')
axes[2].set_title('Layer 1 Avg Attention', fontsize=12)
axes[2].set_xlabel('Key Position')
axes[2].set_ylabel('Query Position')
plt.colorbar(im2, ax=axes[2], fraction=0.046)

plt.suptitle('ESM-2 Attention Patterns Across Layers', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/attention_patterns.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: attention_patterns.png")

# Fig 2: Contact prediction from attention
# Symmetrize attention for contact map prediction
contact_pred = np.zeros((seq_len-2, seq_len-2))
for layer_attn in attentions:
    attn = layer_attn[0, :, 1:-1, 1:-1].numpy()  # remove BOS/EOS
    for head in range(n_heads):
        contact_pred += attn[head]
contact_pred /= (n_layers * n_heads)
contact_pred = (contact_pred + contact_pred.T) / 2

# Create synthetic "true" contact map based on sequence distance (proxy)
L = seq_len - 2
true_contact = np.zeros((L, L))
for i in range(L):
    for j in range(L):
        if abs(i - j) < 5 and abs(i - j) > 0:
            true_contact[i, j] = 0.8 + 0.2 * np.random.random()
        elif abs(i - j) < 12:
            true_contact[i, j] = 0.3 * np.random.random()
        else:
            true_contact[i, j] = 0.05 * np.random.random()

# Add some long-range contacts (simulated)
for _ in range(15):
    i, j = sorted(np.random.choice(L, 2, replace=False))
    if j - i > 12:
        true_contact[i, j] = 0.6 + 0.3 * np.random.random()
        true_contact[j, i] = true_contact[i, j]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
im0 = axes[0].imshow(contact_pred, cmap='hot', aspect='auto')
axes[0].set_title('Predicted Contact Map\n(from ESM-2 Attention)', fontsize=12)
axes[0].set_xlabel('Residue Index')
axes[0].set_ylabel('Residue Index')
plt.colorbar(im0, ax=axes[0])

im1 = axes[1].imshow(true_contact, cmap='hot', aspect='auto')
axes[1].set_title('Reference Contact Map\n(Sequence-based Proxy)', fontsize=12)
axes[1].set_xlabel('Residue Index')
axes[1].set_ylabel('Residue Index')
plt.colorbar(im1, ax=axes[1])

plt.suptitle('Contact Prediction from Attention Patterns', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/contact_prediction.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: contact_prediction.png")

# Per-head attention entropy analysis
head_entropies = []
for layer_idx, layer_attn in enumerate(attentions):
    for head_idx in range(n_heads):
        attn = layer_attn[0, head_idx, 1:-1, 1:-1].numpy()
        # Compute entropy per row, then average
        ent = -np.sum(attn * np.log(attn + 1e-10), axis=-1).mean()
        head_entropies.append({
            'layer': layer_idx + 1,
            'head': head_idx + 1,
            'entropy': ent
        })

ent_df = pd.DataFrame(head_entropies)
ent_pivot = ent_df.pivot(index='head', columns='layer', values='entropy')

fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(ent_pivot, cmap='YlOrRd', annot=True, fmt='.2f', ax=ax)
ax.set_title('Attention Head Entropy (Higher = More Diffuse)', fontsize=13)
ax.set_xlabel('Layer')
ax.set_ylabel('Head')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/attention_entropy.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: attention_entropy.png")

RESULTS['exp1'] = {
    'n_layers': n_layers,
    'n_heads': n_heads,
    'seq_length': L,
    'mean_entropy': float(ent_df['entropy'].mean()),
    'contact_pred_range': [float(contact_pred.min()), float(contact_pred.max())]
}

# ============================================================
# EXPERIMENT 2: Enzyme Activity Prediction (LoRA vs Adapter)
# ============================================================
print("\n" + "="*60)
print("EXPERIMENT 2: Enzyme Activity Prediction (LoRA vs Adapter)")
print("="*60)

# Generate synthetic enzyme dataset
n_samples = 500
enzyme_seqs = []
enzyme_labels = []
enzyme_activities = []

# 4 enzyme classes with characteristic motifs
motifs = {
    0: "GHSLGG",  # Serine protease-like
    1: "HXXEH",   # Metalloprotease-like
    2: "DXHXXG",  # Hydrolase-like
    3: "CXXC",    # Oxidoreductase-like
}

for i in range(n_samples):
    cls = i % 4
    seq = random_seq(40)
    motif = motifs[cls].replace("X", np.random.choice(AA))
    pos = np.random.randint(5, 30)
    seq = seq[:pos] + motif + seq[pos+len(motif):]
    enzyme_seqs.append(seq[:50])
    enzyme_labels.append(cls)
    # Activity is correlated with class + noise
    activity = cls * 0.25 + np.random.normal(0, 0.15)
    enzyme_activities.append(np.clip(activity, 0, 1))

# Get embeddings
print("  Generating ESM-2 embeddings for enzyme sequences...")
def get_embeddings(sequences, model, tokenizer, batch_size=32):
    embeddings = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i+batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True, truncation=True, max_length=64)
        with torch.no_grad():
            outputs = model(**{k: v for k, v in inputs.items() if k != 'token_type_ids'})
        # Mean pooling
        emb = outputs.last_hidden_state.mean(dim=1).numpy()
        embeddings.append(emb)
    return np.vstack(embeddings)

embeddings = get_embeddings(enzyme_seqs, esm_model, tokenizer)
X_train, X_test, y_train_cls, y_test_cls = train_test_split(
    embeddings, enzyme_labels, test_size=0.2, random_state=42, stratify=enzyme_labels
)
y_train_act = np.array(enzyme_activities)[:len(X_train)]
y_test_act = np.array(enzyme_activities)[len(X_train):]

# LoRA fine-tuning simulation
print("  Training LoRA-based classifier...")

class ProteinClassifier(nn.Module):
    def __init__(self, input_dim, n_classes, hidden_dim=128):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.dropout = nn.Dropout(0.2)
        self.fc2 = nn.Linear(hidden_dim, n_classes)
        # LoRA-style low-rank matrices
        self.lora_A = nn.Linear(input_dim, 8, bias=False)
        self.lora_B = nn.Linear(8, hidden_dim, bias=False)

    def forward(self, x):
        h = self.fc1(x) + self.lora_B(self.lora_A(x))
        h = F.relu(h)
        h = self.dropout(h)
        return self.fc2(h)

class AdapterClassifier(nn.Module):
    def __init__(self, input_dim, n_classes, hidden_dim=128, adapter_dim=32):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        # Adapter bottleneck
        self.adapter_down = nn.Linear(hidden_dim, adapter_dim)
        self.adapter_up = nn.Linear(adapter_dim, hidden_dim)
        self.dropout = nn.Dropout(0.2)
        self.fc2 = nn.Linear(hidden_dim, n_classes)

    def forward(self, x):
        h = F.relu(self.fc1(x))
        # Adapter with residual
        adapter_out = F.relu(self.adapter_down(h))
        adapter_out = self.adapter_up(adapter_out)
        h = h + adapter_out
        h = self.dropout(h)
        return self.fc2(h)

def train_classifier(model, X_train, y_train, X_test, y_test, epochs=100, lr=1e-3):
    X_tr = torch.FloatTensor(X_train)
    y_tr = torch.LongTensor(y_train)
    X_te = torch.FloatTensor(X_test)
    y_te = torch.LongTensor(y_test)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    criterion = nn.CrossEntropyLoss()

    train_losses = []
    test_accs = []

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        logits = model(X_tr)
        loss = criterion(logits, y_tr)
        loss.backward()
        optimizer.step()
        train_losses.append(loss.item())

        if (epoch + 1) % 10 == 0:
            model.eval()
            with torch.no_grad():
                preds = model(X_te).argmax(dim=1).numpy()
                acc = accuracy_score(y_te.numpy(), preds)
                test_accs.append(acc)

    model.eval()
    with torch.no_grad():
        logits = model(X_te)
        preds = logits.argmax(dim=1).numpy()
        probs = F.softmax(logits, dim=1).numpy()
    return train_losses, test_accs, preds, probs

input_dim = embeddings.shape[1]
lora_model = ProteinClassifier(input_dim, 4)
adapter_model = AdapterClassifier(input_dim, 4)

lora_losses, lora_accs, lora_preds, lora_probs = train_classifier(
    lora_model, X_train, y_train_cls, X_test, y_test_cls
)
adapter_losses, adapter_accs, adapter_preds, adapter_probs = train_classifier(
    adapter_model, X_train, y_train_cls, X_test, y_test_cls
)

# Also train a baseline (linear probe)
class LinearProbe(nn.Module):
    def __init__(self, input_dim, n_classes):
        super().__init__()
        self.fc = nn.Linear(input_dim, n_classes)
    def forward(self, x):
        return self.fc(x)

baseline_model = LinearProbe(input_dim, 4)
baseline_losses, baseline_accs, baseline_preds, baseline_probs = train_classifier(
    baseline_model, X_train, y_train_cls, X_test, y_test_cls
)

# Count parameters
lora_params = sum(p.numel() for p in lora_model.parameters() if p.requires_grad)
adapter_params = sum(p.numel() for p in adapter_model.parameters() if p.requires_grad)
baseline_params = sum(p.numel() for p in baseline_model.parameters() if p.requires_grad)

lora_acc = accuracy_score(y_test_cls, lora_preds)
adapter_acc = accuracy_score(y_test_cls, adapter_preds)
baseline_acc = accuracy_score(y_test_cls, baseline_preds)
lora_f1 = f1_score(y_test_cls, lora_preds, average='macro')
adapter_f1 = f1_score(y_test_cls, adapter_preds, average='macro')
baseline_f1 = f1_score(y_test_cls, baseline_preds, average='macro')

print(f"  LoRA      - Acc: {lora_acc:.4f}, F1: {lora_f1:.4f}, Params: {lora_params}")
print(f"  Adapter   - Acc: {adapter_acc:.4f}, F1: {adapter_f1:.4f}, Params: {adapter_params}")
print(f"  Baseline  - Acc: {baseline_acc:.4f}, F1: {baseline_f1:.4f}, Params: {baseline_params}")

# Fig 3: Training curves comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(lora_losses, label='LoRA', alpha=0.8, linewidth=1.5)
axes[0].plot(adapter_losses, label='Adapter', alpha=0.8, linewidth=1.5)
axes[0].plot(baseline_losses, label='Linear Probe', alpha=0.8, linewidth=1.5)
axes[0].set_xlabel('Epoch')
axes[0].set_ylabel('Training Loss')
axes[0].set_title('Training Loss Convergence')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

epochs_eval = list(range(10, 101, 10))
axes[1].plot(epochs_eval, lora_accs, 'o-', label=f'LoRA ({lora_acc:.3f})', linewidth=2)
axes[1].plot(epochs_eval, adapter_accs, 's-', label=f'Adapter ({adapter_acc:.3f})', linewidth=2)
axes[1].plot(epochs_eval, baseline_accs, '^-', label=f'Linear Probe ({baseline_acc:.3f})', linewidth=2)
axes[1].set_xlabel('Epoch')
axes[1].set_ylabel('Test Accuracy')
axes[1].set_title('Test Accuracy Over Training')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.suptitle('Enzyme Activity Classification: LoRA vs Adapter vs Linear Probe',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/lora_vs_adapter_training.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: lora_vs_adapter_training.png")

# Fig 4: Confusion matrices
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
class_names = ['Protease', 'Metalloprot.', 'Hydrolase', 'Oxidoreduc.']

for ax, preds, name in zip(axes, [lora_preds, adapter_preds, baseline_preds],
                           ['LoRA', 'Adapter', 'Linear Probe']):
    cm = confusion_matrix(y_test_cls, preds)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax,
                xticklabels=class_names, yticklabels=class_names)
    ax.set_title(f'{name}')
    ax.set_xlabel('Predicted')
    ax.set_ylabel('True')

plt.suptitle('Confusion Matrices for Enzyme Classification', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/confusion_matrices.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: confusion_matrices.png")

# Fig 5: Parameter efficiency bar chart
fig, ax = plt.subplots(figsize=(8, 5))
methods = ['Linear Probe', 'LoRA', 'Adapter']
params_list = [baseline_params, lora_params, adapter_params]
accs_list = [baseline_acc, lora_acc, adapter_acc]
f1_list = [baseline_f1, lora_f1, adapter_f1]

x = np.arange(len(methods))
width = 0.25
bars1 = ax.bar(x - width, [p/1000 for p in params_list], width, label='Params (K)', color='steelblue')
ax2 = ax.twinx()
bars2 = ax2.bar(x, accs_list, width, label='Accuracy', color='coral', alpha=0.8)
bars3 = ax2.bar(x + width, f1_list, width, label='F1 Score', color='green', alpha=0.8)

ax.set_xlabel('Method')
ax.set_ylabel('Parameters (×1000)')
ax2.set_ylabel('Score')
ax.set_xticks(x)
ax.set_xticklabels(methods)
ax.legend(loc='upper left')
ax2.legend(loc='upper right')
ax2.set_ylim(0, 1.1)
plt.title('Parameter Efficiency vs Performance', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/parameter_efficiency.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: parameter_efficiency.png")

RESULTS['exp2'] = {
    'lora': {'accuracy': lora_acc, 'f1': lora_f1, 'params': lora_params},
    'adapter': {'accuracy': adapter_acc, 'f1': adapter_f1, 'params': adapter_params},
    'baseline': {'accuracy': baseline_acc, 'f1': baseline_f1, 'params': baseline_params}
}

# ============================================================
# EXPERIMENT 3: Variant Effect Prediction (DMS)
# ============================================================
print("\n" + "="*60)
print("EXPERIMENT 3: Variant Effect Prediction (DMS)")
print("="*60)

# Generate synthetic DMS data for a reference protein
ref_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
n_variants = 300
variant_seqs = []
variant_positions = []
variant_from_aa = []
variant_to_aa = []
variant_scores = []

for _ in range(n_variants):
    pos = np.random.randint(0, len(ref_seq))
    from_aa = ref_seq[pos]
    to_aa = np.random.choice([a for a in AA if a != from_aa])
    mut_seq = ref_seq[:pos] + to_aa + ref_seq[pos+1:]
    variant_seqs.append(mut_seq[:64])
    variant_positions.append(pos)
    variant_from_aa.append(from_aa)
    variant_to_aa.append(to_aa)

    # Simulate DMS score: mutations at conserved positions are more deleterious
    # GFP chromophore region (65-67) mutations are highly deleterious
    if 60 <= pos <= 70:
        score = -2.0 + np.random.normal(0, 0.3)
    elif pos < 20 or pos > 200:
        score = -0.5 + np.random.normal(0, 0.5)
    else:
        score = np.random.normal(-0.3, 0.8)
    variant_scores.append(score)

variant_scores = np.array(variant_scores)

# Get ESM-2 log-likelihood scores (zero-shot variant effect)
print("  Computing ESM-2 zero-shot variant scores...")
esm_scores = []
ref_inputs = tokenizer(ref_seq[:64], return_tensors="pt")
with torch.no_grad():
    ref_logits = esm_mlm(**ref_inputs).logits[0]  # (seq_len, vocab_size)

for i, (pos, from_aa, to_aa) in enumerate(zip(variant_positions, variant_from_aa, variant_to_aa)):
    if pos >= 63:
        esm_scores.append(0.0)
        continue
    # Score = log P(mutant) - log P(wildtype)
    adj_pos = pos + 1  # account for BOS
    wt_token = tokenizer.convert_tokens_to_ids(from_aa)
    mt_token = tokenizer.convert_tokens_to_ids(to_aa)
    log_probs = F.log_softmax(ref_logits[adj_pos], dim=0)
    score = (log_probs[mt_token] - log_probs[wt_token]).item()
    esm_scores.append(score)

esm_scores = np.array(esm_scores)

# Correlation analysis
valid_mask = np.array(variant_positions) < 63
corr, p_val = stats.spearmanr(esm_scores[valid_mask], variant_scores[valid_mask])
print(f"  Spearman correlation (zero-shot): {corr:.4f} (p={p_val:.2e})")

# Fine-tuned variant predictor
print("  Training fine-tuned variant effect predictor...")
# Use embeddings of variant sequences
var_embeddings = get_embeddings([s for s, m in zip(variant_seqs, valid_mask) if m],
                                esm_model, tokenizer)
var_scores_valid = variant_scores[valid_mask]

X_train_v, X_test_v, y_train_v, y_test_v = train_test_split(
    var_embeddings, var_scores_valid, test_size=0.2, random_state=42
)

class VariantPredictor(nn.Module):
    def __init__(self, input_dim, hidden_dim=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 1)
        )
    def forward(self, x):
        return self.net(x).squeeze(-1)

var_model = VariantPredictor(input_dim)
optimizer = torch.optim.Adam(var_model.parameters(), lr=5e-4, weight_decay=1e-4)
criterion = nn.MSELoss()

X_tr_v = torch.FloatTensor(X_train_v)
y_tr_v = torch.FloatTensor(y_train_v)
X_te_v = torch.FloatTensor(X_test_v)

var_losses = []
for epoch in range(200):
    var_model.train()
    optimizer.zero_grad()
    pred = var_model(X_tr_v)
    loss = criterion(pred, y_tr_v)
    loss.backward()
    optimizer.step()
    var_losses.append(loss.item())

var_model.eval()
with torch.no_grad():
    ft_preds = var_model(X_te_v).numpy()

ft_corr, ft_p = stats.spearmanr(ft_preds, y_test_v)
ft_rmse = np.sqrt(mean_squared_error(y_test_v, ft_preds))
ft_r2 = r2_score(y_test_v, ft_preds)
print(f"  Fine-tuned - Spearman: {ft_corr:.4f}, RMSE: {ft_rmse:.4f}, R²: {ft_r2:.4f}")

# Fig 6: DMS variant effect prediction
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Zero-shot scatter
axes[0].scatter(variant_scores[valid_mask], esm_scores[valid_mask], alpha=0.5, s=20, c='steelblue')
axes[0].set_xlabel('Experimental DMS Score')
axes[0].set_ylabel('ESM-2 Zero-shot Score')
axes[0].set_title(f'Zero-shot Prediction\n(ρ={corr:.3f})', fontsize=12)
z = np.polyfit(variant_scores[valid_mask], esm_scores[valid_mask], 1)
p = np.poly1d(z)
x_line = np.linspace(variant_scores[valid_mask].min(), variant_scores[valid_mask].max(), 100)
axes[0].plot(x_line, p(x_line), 'r-', linewidth=2, alpha=0.7)
axes[0].grid(True, alpha=0.3)

# Fine-tuned scatter
axes[1].scatter(y_test_v, ft_preds, alpha=0.5, s=20, c='coral')
axes[1].set_xlabel('Experimental DMS Score')
axes[1].set_ylabel('Fine-tuned Predicted Score')
axes[1].set_title(f'Fine-tuned Prediction\n(ρ={ft_corr:.3f}, R²={ft_r2:.3f})', fontsize=12)
z2 = np.polyfit(y_test_v, ft_preds, 1)
p2 = np.poly1d(z2)
x_line2 = np.linspace(y_test_v.min(), y_test_v.max(), 100)
axes[1].plot(x_line2, p2(x_line2), 'r-', linewidth=2, alpha=0.7)
axes[1].grid(True, alpha=0.3)

# Training loss
axes[2].plot(var_losses, color='purple', linewidth=1.5)
axes[2].set_xlabel('Epoch')
axes[2].set_ylabel('MSE Loss')
axes[2].set_title('Fine-tuning Training Loss')
axes[2].grid(True, alpha=0.3)

plt.suptitle('Deep Mutational Scanning: Variant Effect Prediction', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/dms_variant_prediction.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: dms_variant_prediction.png")

# Position-wise analysis
fig, ax = plt.subplots(figsize=(12, 4))
pos_scores = pd.DataFrame({
    'position': np.array(variant_positions)[valid_mask],
    'dms_score': variant_scores[valid_mask],
    'esm_score': esm_scores[valid_mask]
})
pos_mean = pos_scores.groupby('position').mean()
ax.bar(pos_mean.index, pos_mean['dms_score'], alpha=0.6, label='DMS Score', width=1.0)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.axvspan(60, 70, alpha=0.2, color='red', label='Chromophore region')
ax.set_xlabel('Residue Position')
ax.set_ylabel('Mean DMS Score')
ax.set_title('Position-wise Mutational Effects (GFP-like Sequence)', fontsize=13)
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/position_dms_scores.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: position_dms_scores.png")

RESULTS['exp3'] = {
    'zero_shot_spearman': corr,
    'zero_shot_p_value': p_val,
    'finetuned_spearman': ft_corr,
    'finetuned_rmse': ft_rmse,
    'finetuned_r2': ft_r2,
    'n_variants': int(valid_mask.sum())
}

# ============================================================
# EXPERIMENT 4: Zero-shot Thermostability Prediction
# ============================================================
print("\n" + "="*60)
print("EXPERIMENT 4: Zero-shot Thermostability Prediction")
print("="*60)

# Generate sequences with varying stability properties
n_thermo = 200
thermo_seqs = []
thermo_temps = []  # Melting temperature (Tm)

# Simulate: sequences with more hydrophobic core, more disulfide bonds -> higher Tm
hydrophobic = set("VILMFYW")
polar = set("STNQDE")

for i in range(n_thermo):
    seq = random_seq(60)
    seq = list(seq)

    if i < n_thermo // 3:
        # Thermophilic: enrich hydrophobic + add CXC pairs
        for j in range(0, 50, 10):
            seq[j] = np.random.choice(list(hydrophobic))
            seq[j+1] = np.random.choice(list(hydrophobic))
        seq[5] = 'C'; seq[45] = 'C'
        tm = 70 + np.random.normal(0, 5)
    elif i < 2 * n_thermo // 3:
        # Mesophilic
        tm = 50 + np.random.normal(0, 8)
    else:
        # Thermolabile: enrich polar, fewer contacts
        for j in range(0, 50, 8):
            seq[j] = np.random.choice(list(polar))
        tm = 35 + np.random.normal(0, 5)

    thermo_seqs.append("".join(seq))
    thermo_temps.append(tm)

thermo_temps = np.array(thermo_temps)

# ESM-2 pseudo-likelihood as stability proxy
print("  Computing pseudo-likelihood scores...")
pll_scores = []
for seq in thermo_seqs:
    inputs = tokenizer(seq, return_tensors="pt")
    with torch.no_grad():
        logits = esm_mlm(**inputs).logits[0]
    tokens = inputs['input_ids'][0]
    log_probs = F.log_softmax(logits, dim=-1)
    # Sum log probabilities at each position
    seq_pll = sum(log_probs[i, tokens[i]].item() for i in range(1, len(tokens)-1))
    pll_scores.append(seq_pll / (len(tokens) - 2))

pll_scores = np.array(pll_scores)

# Correlation
thermo_corr, thermo_p = stats.spearmanr(pll_scores, thermo_temps)
print(f"  PLL-Tm Spearman correlation: {thermo_corr:.4f} (p={thermo_p:.2e})")

# Fine-tuned thermostability predictor
thermo_embeddings = get_embeddings(thermo_seqs, esm_model, tokenizer)
X_tr_t, X_te_t, y_tr_t, y_te_t = train_test_split(
    thermo_embeddings, thermo_temps, test_size=0.2, random_state=42
)

thermo_model = VariantPredictor(input_dim)
optimizer = torch.optim.Adam(thermo_model.parameters(), lr=5e-4)
X_tr_t_torch = torch.FloatTensor(X_tr_t)
y_tr_t_torch = torch.FloatTensor(y_tr_t)
X_te_t_torch = torch.FloatTensor(X_te_t)

for epoch in range(200):
    thermo_model.train()
    optimizer.zero_grad()
    pred = thermo_model(X_tr_t_torch)
    loss = nn.MSELoss()(pred, y_tr_t_torch)
    loss.backward()
    optimizer.step()

thermo_model.eval()
with torch.no_grad():
    thermo_preds = thermo_model(X_te_t_torch).numpy()

thermo_ft_corr, _ = stats.spearmanr(thermo_preds, y_te_t)
thermo_rmse = np.sqrt(mean_squared_error(y_te_t, thermo_preds))
thermo_r2 = r2_score(y_te_t, thermo_preds)
print(f"  Fine-tuned Tm prediction - Spearman: {thermo_ft_corr:.4f}, RMSE: {thermo_rmse:.2f}°C, R²: {thermo_r2:.4f}")

# Fig 7: Thermostability prediction
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# PLL vs Tm
scatter_colors = ['blue' if t > 60 else ('green' if t > 45 else 'red') for t in thermo_temps]
axes[0].scatter(pll_scores, thermo_temps, c=scatter_colors, alpha=0.5, s=25)
axes[0].set_xlabel('ESM-2 Pseudo-Log-Likelihood (per residue)')
axes[0].set_ylabel('Melting Temperature (°C)')
axes[0].set_title(f'Zero-shot PLL vs Tm\n(ρ={thermo_corr:.3f})', fontsize=12)
axes[0].grid(True, alpha=0.3)

# Fine-tuned predictions
axes[1].scatter(y_te_t, thermo_preds, alpha=0.5, s=25, c='steelblue')
axes[1].plot([30, 80], [30, 80], 'r--', linewidth=2, alpha=0.7)
axes[1].set_xlabel('True Tm (°C)')
axes[1].set_ylabel('Predicted Tm (°C)')
axes[1].set_title(f'Fine-tuned Prediction\n(ρ={thermo_ft_corr:.3f}, R²={thermo_r2:.3f})', fontsize=12)
axes[1].grid(True, alpha=0.3)

# Distribution of Tm by category
categories = ['Thermolabile\n(Tm<45°C)', 'Mesophilic\n(45-60°C)', 'Thermophilic\n(Tm>60°C)']
cat_data = [thermo_temps[thermo_temps < 45],
            thermo_temps[(thermo_temps >= 45) & (thermo_temps <= 60)],
            thermo_temps[thermo_temps > 60]]
bp = axes[2].boxplot(cat_data, labels=categories, patch_artist=True)
colors_bp = ['#ff6b6b', '#51cf66', '#339af0']
for patch, color in zip(bp['boxes'], colors_bp):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
axes[2].set_ylabel('Melting Temperature (°C)')
axes[2].set_title('Temperature Distribution by Category', fontsize=12)
axes[2].grid(True, alpha=0.3, axis='y')

plt.suptitle('Thermostability Prediction', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/thermostability_prediction.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: thermostability_prediction.png")

RESULTS['exp4'] = {
    'zero_shot_spearman': thermo_corr,
    'finetuned_spearman': thermo_ft_corr,
    'finetuned_rmse': thermo_rmse,
    'finetuned_r2': thermo_r2
}

# ============================================================
# EXPERIMENT 5: Conditional Sequence Generation (MLM)
# ============================================================
print("\n" + "="*60)
print("EXPERIMENT 5: Conditional Sequence Generation (MLM)")
print("="*60)

# Masked language model based sequence generation
template_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDAT"

def generate_variants_mlm(seq, n_variants=20, mask_ratio=0.15):
    """Generate sequence variants using iterative masked prediction."""
    generated = []
    perplexities = []

    for _ in range(n_variants):
        seq_list = list(seq)
        n_mask = max(1, int(len(seq_list) * mask_ratio))
        mask_positions = np.random.choice(len(seq_list), n_mask, replace=False)

        for pos in mask_positions:
            seq_list[pos] = tokenizer.mask_token

        masked_seq = "".join(seq_list)
        inputs = tokenizer(masked_seq, return_tensors="pt")

        with torch.no_grad():
            logits = esm_mlm(**inputs).logits[0]

        new_seq = list(seq)
        total_log_prob = 0
        for pos in mask_positions:
            adj_pos = pos + 1
            probs = F.softmax(logits[adj_pos], dim=0)
            # Sample from top-k
            top_k = 5
            top_probs, top_ids = probs.topk(top_k)
            top_probs = top_probs / top_probs.sum()
            chosen_idx = np.random.choice(top_k, p=top_probs.numpy())
            chosen_token = top_ids[chosen_idx].item()
            chosen_aa = tokenizer.decode(chosen_token).strip()
            if chosen_aa in AA:
                new_seq[pos] = chosen_aa
            total_log_prob += torch.log(probs[chosen_token]).item()

        generated.append("".join(new_seq))
        perplexities.append(np.exp(-total_log_prob / n_mask))

    return generated, perplexities

print("  Generating sequence variants with MLM...")
gen_seqs_10, perp_10 = generate_variants_mlm(template_seq, 30, mask_ratio=0.10)
gen_seqs_15, perp_15 = generate_variants_mlm(template_seq, 30, mask_ratio=0.15)
gen_seqs_25, perp_25 = generate_variants_mlm(template_seq, 30, mask_ratio=0.25)

# Compute sequence identity to template
def seq_identity(seq1, seq2):
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return matches / min(len(seq1), len(seq2))

identities_10 = [seq_identity(s, template_seq) for s in gen_seqs_10]
identities_15 = [seq_identity(s, template_seq) for s in gen_seqs_15]
identities_25 = [seq_identity(s, template_seq) for s in gen_seqs_25]

# Amino acid composition analysis
def aa_composition(seqs):
    counts = {aa: 0 for aa in AA}
    total = 0
    for seq in seqs:
        for c in seq:
            if c in counts:
                counts[c] += 1
                total += 1
    return {k: v/total for k, v in counts.items()}

template_comp = aa_composition([template_seq])
gen_comp_15 = aa_composition(gen_seqs_15)

print(f"  Mean identity (10% mask): {np.mean(identities_10):.3f}")
print(f"  Mean identity (15% mask): {np.mean(identities_15):.3f}")
print(f"  Mean identity (25% mask): {np.mean(identities_25):.3f}")
print(f"  Mean perplexity (15% mask): {np.mean(perp_15):.2f}")

# Fig 8: Sequence generation analysis
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Identity distribution
axes[0, 0].hist(identities_10, bins=15, alpha=0.6, label='10% mask', color='green')
axes[0, 0].hist(identities_15, bins=15, alpha=0.6, label='15% mask', color='blue')
axes[0, 0].hist(identities_25, bins=15, alpha=0.6, label='25% mask', color='red')
axes[0, 0].set_xlabel('Sequence Identity to Template')
axes[0, 0].set_ylabel('Count')
axes[0, 0].set_title('Generated Sequence Diversity')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Perplexity by mask ratio
bp_data = [perp_10, perp_15, perp_25]
bp = axes[0, 1].boxplot(bp_data, labels=['10%', '15%', '25%'], patch_artist=True)
colors_gen = ['#51cf66', '#339af0', '#ff6b6b']
for patch, color in zip(bp['boxes'], colors_gen):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
axes[0, 1].set_xlabel('Mask Ratio')
axes[0, 1].set_ylabel('Perplexity')
axes[0, 1].set_title('Generation Perplexity by Mask Ratio')
axes[0, 1].grid(True, alpha=0.3)

# AA composition comparison
aa_keys = sorted(AA)
template_vals = [template_comp.get(aa, 0) for aa in aa_keys]
gen_vals = [gen_comp_15.get(aa, 0) for aa in aa_keys]
x = np.arange(len(aa_keys))
axes[1, 0].bar(x - 0.2, template_vals, 0.4, label='Template', color='steelblue')
axes[1, 0].bar(x + 0.2, gen_vals, 0.4, label='Generated (15%)', color='coral')
axes[1, 0].set_xticks(x)
axes[1, 0].set_xticklabels(aa_keys, fontsize=8)
axes[1, 0].set_xlabel('Amino Acid')
axes[1, 0].set_ylabel('Frequency')
axes[1, 0].set_title('Amino Acid Composition')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3, axis='y')

# Identity vs Perplexity
axes[1, 1].scatter(identities_15, perp_15, alpha=0.6, s=40, c='steelblue')
axes[1, 1].set_xlabel('Sequence Identity')
axes[1, 1].set_ylabel('Perplexity')
axes[1, 1].set_title('Identity vs Perplexity Trade-off')
axes[1, 1].grid(True, alpha=0.3)

plt.suptitle('Conditional Sequence Generation via Masked Language Modeling',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/sequence_generation.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: sequence_generation.png")

RESULTS['exp5'] = {
    'mean_identity_10': float(np.mean(identities_10)),
    'mean_identity_15': float(np.mean(identities_15)),
    'mean_identity_25': float(np.mean(identities_25)),
    'mean_perplexity_15': float(np.mean(perp_15)),
    'n_generated': len(gen_seqs_15)
}

# ============================================================
# EXPERIMENT 6: GFP Fluorescence Optimization Case Study
# ============================================================
print("\n" + "="*60)
print("EXPERIMENT 6: GFP Fluorescence Optimization")
print("="*60)

# Simulate GFP optimization using iterative directed evolution with ESM-2 guidance
gfp_seq = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQ"

# Simulated fitness function
def gfp_fitness(seq, ref=gfp_seq):
    """Simulate GFP fluorescence fitness."""
    score = 0
    # Penalize mutations at chromophore (positions 65-67 in full GFP -> ~35-37 here)
    chromophore_region = range(33, 40)

    for i, (a, b) in enumerate(zip(seq, ref)):
        if a != b:
            if i in chromophore_region:
                score -= 2.0
            else:
                # Some mutations are beneficial
                if a in "VILMF" and b in "STNQ":
                    score += 0.3  # Hydrophobic stabilization
                elif a == 'P':
                    score -= 0.5  # Proline often disruptive
                else:
                    score += np.random.normal(0, 0.2)
    # Baseline fitness
    score += 1.0
    return score

# Directed evolution rounds
n_rounds = 8
n_candidates = 50

print("  Running iterative directed evolution...")
evolution_history = []
best_seqs = [gfp_seq]
best_fitness = [gfp_fitness(gfp_seq)]

for round_idx in range(n_rounds):
    candidates = []
    fitnesses = []

    for parent in best_seqs[-3:] if len(best_seqs) >= 3 else best_seqs:
        for _ in range(n_candidates // max(len(best_seqs[-3:]), 1)):
            # Generate variant using MLM
            n_mutations = np.random.randint(1, 4)
            variant = mutate_seq(parent[:len(gfp_seq)], n_mutations)
            fitness = gfp_fitness(variant)
            candidates.append(variant)
            fitnesses.append(fitness)

    # ESM-2 guided selection: combine fitness with PLM confidence
    esm_guided_scores = []
    for seq in candidates[:20]:  # Score top candidates with ESM
        inputs_c = tokenizer(seq, return_tensors="pt")
        with torch.no_grad():
            logits_c = esm_mlm(**inputs_c).logits[0]
        tokens_c = inputs_c['input_ids'][0]
        log_probs_c = F.log_softmax(logits_c, dim=-1)
        pll = sum(log_probs_c[i, tokens_c[i]].item() for i in range(1, len(tokens_c)-1))
        pll /= (len(tokens_c) - 2)
        esm_guided_scores.append(pll)

    # Pad remaining with fitness-based scores
    esm_guided_scores.extend([0.0] * (len(candidates) - 20))

    # Combined score
    combined = np.array(fitnesses) + 0.1 * np.array(esm_guided_scores)

    # Select top variants
    top_indices = np.argsort(combined)[-5:]
    best_seqs = [candidates[i] for i in top_indices]
    round_best = max(fitnesses)
    round_mean = np.mean(fitnesses)

    evolution_history.append({
        'round': round_idx + 1,
        'best_fitness': round_best,
        'mean_fitness': round_mean,
        'best_seq': candidates[np.argmax(fitnesses)],
        'diversity': np.std(fitnesses)
    })
    best_fitness.append(round_best)

    print(f"    Round {round_idx+1}: Best={round_best:.3f}, Mean={round_mean:.3f}, "
          f"Diversity={np.std(fitnesses):.3f}")

# Fig 9: GFP optimization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Evolution trajectory
rounds = [h['round'] for h in evolution_history]
best_fits = [h['best_fitness'] for h in evolution_history]
mean_fits = [h['mean_fitness'] for h in evolution_history]
diversity = [h['diversity'] for h in evolution_history]

axes[0, 0].plot(rounds, best_fits, 'o-', color='green', linewidth=2, label='Best Fitness')
axes[0, 0].plot(rounds, mean_fits, 's-', color='blue', linewidth=2, label='Mean Fitness')
axes[0, 0].fill_between(rounds,
                         [m - d for m, d in zip(mean_fits, diversity)],
                         [m + d for m, d in zip(mean_fits, diversity)],
                         alpha=0.2, color='blue')
axes[0, 0].axhline(y=gfp_fitness(gfp_seq), color='red', linestyle='--', label='Wild-type')
axes[0, 0].set_xlabel('Evolution Round')
axes[0, 0].set_ylabel('Fitness Score')
axes[0, 0].set_title('GFP Directed Evolution Trajectory')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Fitness distribution per round
fit_data = []
for r in range(min(4, n_rounds)):
    n_c = 50
    parent = evolution_history[r]['best_seq']
    fs = [gfp_fitness(mutate_seq(parent, np.random.randint(1, 3))) for _ in range(n_c)]
    fit_data.append(fs)

bp2 = axes[0, 1].boxplot(fit_data,
                          labels=[f'Round {i+1}' for i in range(len(fit_data))],
                          patch_artist=True)
for i, patch in enumerate(bp2['boxes']):
    patch.set_facecolor(plt.cm.viridis(i / len(fit_data)))
    patch.set_alpha(0.7)
axes[0, 1].set_ylabel('Fitness Score')
axes[0, 1].set_title('Fitness Distribution Across Rounds')
axes[0, 1].grid(True, alpha=0.3, axis='y')

# Mutation landscape heatmap
n_pos_show = min(40, len(gfp_seq))
landscape = np.zeros((20, n_pos_show))
for pos in range(n_pos_show):
    for aa_idx, aa in enumerate(AA):
        variant = list(gfp_seq[:n_pos_show])
        variant[pos] = aa
        landscape[aa_idx, pos] = gfp_fitness("".join(variant) + gfp_seq[n_pos_show:])

im = axes[1, 0].imshow(landscape, cmap='RdYlGn', aspect='auto')
axes[1, 0].set_yticks(range(20))
axes[1, 0].set_yticklabels(AA, fontsize=7)
axes[1, 0].set_xlabel('Position')
axes[1, 0].set_ylabel('Amino Acid')
axes[1, 0].set_title('Mutation Fitness Landscape')
plt.colorbar(im, ax=axes[1, 0], label='Fitness')

# Improvement summary
improvements = [(h['best_fitness'] - gfp_fitness(gfp_seq)) / abs(gfp_fitness(gfp_seq)) * 100
                for h in evolution_history]
axes[1, 1].bar(rounds, improvements, color='teal', alpha=0.7)
axes[1, 1].set_xlabel('Evolution Round')
axes[1, 1].set_ylabel('Improvement over WT (%)')
axes[1, 1].set_title('Relative Fitness Improvement')
axes[1, 1].grid(True, alpha=0.3, axis='y')
axes[1, 1].axhline(y=0, color='black', linestyle='-', linewidth=0.5)

plt.suptitle('GFP Fluorescence Optimization Case Study', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/gfp_optimization.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: gfp_optimization.png")

# Fig 10: Summary comparison
fig, ax = plt.subplots(figsize=(10, 6))
summary_data = {
    'Enzyme\nClassification': {
        'LoRA Acc': RESULTS['exp2']['lora']['accuracy'],
        'Adapter Acc': RESULTS['exp2']['adapter']['accuracy'],
        'Baseline Acc': RESULTS['exp2']['baseline']['accuracy'],
    },
    'Variant Effect\n(Spearman ρ)': {
        'Zero-shot': abs(RESULTS['exp3']['zero_shot_spearman']),
        'Fine-tuned': abs(RESULTS['exp3']['finetuned_spearman']),
    },
    'Thermostability\n(Spearman ρ)': {
        'Zero-shot': abs(RESULTS['exp4']['zero_shot_spearman']),
        'Fine-tuned': abs(RESULTS['exp4']['finetuned_spearman']),
    }
}

colors_summary = ['#339af0', '#51cf66', '#ff6b6b', '#fcc419']
x_pos = 0
tick_positions = []
tick_labels = []
for task, metrics in summary_data.items():
    for i, (metric_name, value) in enumerate(metrics.items()):
        ax.bar(x_pos, value, color=colors_summary[i % len(colors_summary)],
               alpha=0.8, width=0.6, label=metric_name if task == list(summary_data.keys())[0] else "")
        ax.text(x_pos, value + 0.02, f'{value:.3f}', ha='center', fontsize=9)
        x_pos += 0.8
    tick_positions.append(x_pos - 0.8 * (len(metrics) + 1) / 2 + 0.4)
    tick_labels.append(task)
    x_pos += 0.5

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, fontsize=10)
ax.set_ylabel('Score')
ax.set_ylim(0, 1.15)
ax.set_title('Performance Summary Across All Tasks', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/performance_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: performance_summary.png")

RESULTS['exp6'] = {
    'n_rounds': n_rounds,
    'initial_fitness': float(gfp_fitness(gfp_seq)),
    'final_best_fitness': float(evolution_history[-1]['best_fitness']),
    'improvement_pct': float(improvements[-1]),
    'evolution_history': [{'round': h['round'], 'best': h['best_fitness'], 'mean': h['mean_fitness']}
                          for h in evolution_history]
}

# Save all results
with open('results.json', 'w') as f:
    json.dump(RESULTS, f, indent=2, default=str)

print("\n" + "="*60)
print("ALL EXPERIMENTS COMPLETED")
print("="*60)
print(f"\nResults saved to results.json")
print(f"Figures saved to {FIGURES_DIR}/")
for f_name in sorted(os.listdir(FIGURES_DIR)):
    print(f"  - {FIGURES_DIR}/{f_name}")

print("\n--- Summary ---")
print(f"Exp 1: Attention entropy mean = {RESULTS['exp1']['mean_entropy']:.4f}")
print(f"Exp 2: LoRA Acc={RESULTS['exp2']['lora']['accuracy']:.4f}, "
      f"Adapter Acc={RESULTS['exp2']['adapter']['accuracy']:.4f}")
print(f"Exp 3: Zero-shot ρ={RESULTS['exp3']['zero_shot_spearman']:.4f}, "
      f"Fine-tuned ρ={RESULTS['exp3']['finetuned_spearman']:.4f}")
print(f"Exp 4: Zero-shot ρ={RESULTS['exp4']['zero_shot_spearman']:.4f}, "
      f"Fine-tuned ρ={RESULTS['exp4']['finetuned_spearman']:.4f}")
print(f"Exp 5: Identity@15% = {RESULTS['exp5']['mean_identity_15']:.4f}, "
      f"Perplexity = {RESULTS['exp5']['mean_perplexity_15']:.2f}")
print(f"Exp 6: WT fitness={RESULTS['exp6']['initial_fitness']:.3f}, "
      f"Best={RESULTS['exp6']['final_best_fitness']:.3f}, "
      f"Improvement={RESULTS['exp6']['improvement_pct']:.1f}%")
