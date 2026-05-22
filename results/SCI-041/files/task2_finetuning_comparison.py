"""
Task 2: Enzyme Activity Prediction Fine-tuning
- LoRA vs Adapter vs Full Fine-tuning comparison
- HuggingFace Transformers + PEFT pipeline
"""
import json
import os
import time
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from transformers import AutoTokenizer, EsmModel, get_linear_schedule_with_warmup
from peft import LoraConfig, get_peft_model, TaskType
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import spearmanr

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

MODEL_NAME = "facebook/esm2_t6_8M_UR50D"
SEED = 42
N_SAMPLES = 200
SEQ_LEN_RANGE = (30, 80)
N_EPOCHS = 5
BATCH_SIZE = 8
LR = 1e-4

torch.manual_seed(SEED)
np.random.seed(SEED)

# --- Synthetic enzyme activity dataset ---
AA = list("ACDEFGHIKLMNPQRSTVWY")

def generate_synthetic_enzyme_data(n_samples=N_SAMPLES):
    """Generate synthetic enzyme sequences with activity labels."""
    data = []
    for i in range(n_samples):
        length = np.random.randint(*SEQ_LEN_RANGE)
        seq = ''.join(np.random.choice(AA, length))
        # Activity correlated with hydrophobic content + some noise
        hydrophobic = sum(1 for c in seq if c in 'AILMFWVP') / length
        charged = sum(1 for c in seq if c in 'DEKRH') / length
        activity = 2.0 * hydrophobic - 1.5 * charged + np.random.normal(0, 0.15)
        data.append({'sequence': seq, 'activity': float(activity)})
    return data

class EnzymeDataset(Dataset):
    def __init__(self, data, tokenizer, max_len=82):
        self.data = data
        self.tokenizer = tokenizer
        self.max_len = max_len

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        item = self.data[idx]
        enc = self.tokenizer(item['sequence'], return_tensors='pt',
                             padding='max_length', truncation=True,
                             max_length=self.max_len)
        return {
            'input_ids': enc['input_ids'].squeeze(0),
            'attention_mask': enc['attention_mask'].squeeze(0),
            'label': torch.tensor(item['activity'], dtype=torch.float32)
        }

class ESMRegressor(nn.Module):
    def __init__(self, base_model):
        super().__init__()
        self.esm = base_model
        hidden_size = base_model.config.hidden_size
        self.regressor = nn.Sequential(
            nn.Linear(hidden_size, 128),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(128, 1)
        )

    def forward(self, input_ids, attention_mask):
        outputs = self.esm(input_ids=input_ids, attention_mask=attention_mask)
        cls_output = outputs.last_hidden_state[:, 0, :]
        return self.regressor(cls_output).squeeze(-1)

def count_trainable_params(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def count_total_params(model):
    return sum(p.numel() for p in model.parameters())

def train_model(model, train_loader, val_loader, n_epochs=N_EPOCHS, lr=LR, method_name=""):
    optimizer = torch.optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()), lr=lr)
    total_steps = len(train_loader) * n_epochs
    scheduler = get_linear_schedule_with_warmup(optimizer, num_warmup_steps=total_steps//10,
                                                 num_training_steps=total_steps)
    criterion = nn.MSELoss()

    history = {'train_loss': [], 'val_loss': [], 'val_r2': [], 'val_spearman': []}

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0
        for batch in train_loader:
            optimizer.zero_grad()
            preds = model(batch['input_ids'], batch['attention_mask'])
            loss = criterion(preds, batch['label'])
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            scheduler.step()
            epoch_loss += loss.item()
        avg_train_loss = epoch_loss / len(train_loader)
        history['train_loss'].append(avg_train_loss)

        # Validation
        model.eval()
        all_preds, all_labels = [], []
        val_loss_sum = 0
        with torch.no_grad():
            for batch in val_loader:
                preds = model(batch['input_ids'], batch['attention_mask'])
                loss = criterion(preds, batch['label'])
                val_loss_sum += loss.item()
                all_preds.extend(preds.numpy().tolist())
                all_labels.extend(batch['label'].numpy().tolist())

        avg_val_loss = val_loss_sum / len(val_loader)
        r2 = r2_score(all_labels, all_preds)
        sp, _ = spearmanr(all_labels, all_preds)
        history['val_loss'].append(avg_val_loss)
        history['val_r2'].append(r2)
        history['val_spearman'].append(sp)

        print(f"  [{method_name}] Epoch {epoch+1}/{n_epochs} | "
              f"Train Loss: {avg_train_loss:.4f} | Val Loss: {avg_val_loss:.4f} | "
              f"R²: {r2:.3f} | Spearman: {sp:.3f}")

    return history, all_preds, all_labels

def run_finetuning_experiment():
    print("=" * 60)
    print("Task 2: Enzyme Activity Prediction Fine-tuning")
    print("=" * 60)
    start = time.time()

    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)

    # Generate data
    print("\n[1/5] Generating synthetic enzyme dataset...")
    data = generate_synthetic_enzyme_data()
    split = int(0.8 * len(data))
    train_data, val_data = data[:split], data[split:]
    train_ds = EnzymeDataset(train_data, tokenizer)
    val_ds = EnzymeDataset(val_data, tokenizer)
    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=BATCH_SIZE)

    results = {}

    # --- Method 1: Frozen backbone (linear probe) ---
    print("\n[2/5] Training: Frozen backbone (linear probe)...")
    base_model = EsmModel.from_pretrained(MODEL_NAME)
    model_frozen = ESMRegressor(base_model)
    for param in model_frozen.esm.parameters():
        param.requires_grad = False
    trainable = count_trainable_params(model_frozen)
    total = count_total_params(model_frozen)
    print(f"  Trainable params: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")
    hist_frozen, preds_frozen, labels_frozen = train_model(
        model_frozen, train_loader, val_loader, method_name="Frozen")
    results['frozen'] = {
        'history': hist_frozen, 'trainable_params': trainable, 'total_params': total,
        'final_r2': hist_frozen['val_r2'][-1], 'final_spearman': hist_frozen['val_spearman'][-1]
    }

    # --- Method 2: LoRA ---
    print("\n[3/5] Training: LoRA fine-tuning...")
    base_model_lora = EsmModel.from_pretrained(MODEL_NAME)
    lora_config = LoraConfig(
        task_type=TaskType.FEATURE_EXTRACTION,
        r=8, lora_alpha=16, lora_dropout=0.1,
        target_modules=["query", "value"]
    )
    lora_model = get_peft_model(base_model_lora, lora_config)
    model_lora = ESMRegressor(lora_model)
    trainable = count_trainable_params(model_lora)
    total = count_total_params(model_lora)
    print(f"  Trainable params: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")
    hist_lora, preds_lora, labels_lora = train_model(
        model_lora, train_loader, val_loader, method_name="LoRA")
    results['lora'] = {
        'history': hist_lora, 'trainable_params': trainable, 'total_params': total,
        'final_r2': hist_lora['val_r2'][-1], 'final_spearman': hist_lora['val_spearman'][-1]
    }

    # --- Method 3: Adapter (bottleneck via manual injection) ---
    print("\n[4/5] Training: Last-2-layers fine-tuning (adapter-like)...")
    base_model_adapter = EsmModel.from_pretrained(MODEL_NAME)
    model_adapter = ESMRegressor(base_model_adapter)
    # Freeze all but last 2 encoder layers + regressor
    for param in model_adapter.esm.parameters():
        param.requires_grad = False
    for param in model_adapter.esm.encoder.layer[-2:].parameters():
        param.requires_grad = True
    trainable = count_trainable_params(model_adapter)
    total = count_total_params(model_adapter)
    print(f"  Trainable params: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")
    hist_adapter, preds_adapter, labels_adapter = train_model(
        model_adapter, train_loader, val_loader, method_name="Last2Layers")
    results['last2layers'] = {
        'history': hist_adapter, 'trainable_params': trainable, 'total_params': total,
        'final_r2': hist_adapter['val_r2'][-1], 'final_spearman': hist_adapter['val_spearman'][-1]
    }

    # --- Method 4: Full fine-tuning ---
    print("\n[5/5] Training: Full fine-tuning...")
    base_model_full = EsmModel.from_pretrained(MODEL_NAME)
    model_full = ESMRegressor(base_model_full)
    trainable = count_trainable_params(model_full)
    total = count_total_params(model_full)
    print(f"  Trainable params: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")
    hist_full, preds_full, labels_full = train_model(
        model_full, train_loader, val_loader, lr=5e-5, method_name="FullFT")
    results['full_ft'] = {
        'history': hist_full, 'trainable_params': trainable, 'total_params': total,
        'final_r2': hist_full['val_r2'][-1], 'final_spearman': hist_full['val_spearman'][-1]
    }

    # --- Visualization ---
    methods = ['frozen', 'lora', 'last2layers', 'full_ft']
    labels_map = {'frozen': 'Frozen (Linear Probe)', 'lora': 'LoRA (r=8)',
                  'last2layers': 'Last-2-Layers FT', 'full_ft': 'Full Fine-tuning'}
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Training loss
    for m, c in zip(methods, colors):
        axes[0].plot(results[m]['history']['train_loss'], label=labels_map[m], color=c, linewidth=2)
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Training Loss (MSE)')
    axes[0].set_title('Training Loss')
    axes[0].legend()

    # Validation R²
    for m, c in zip(methods, colors):
        axes[1].plot(results[m]['history']['val_r2'], label=labels_map[m], color=c, linewidth=2)
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('R²')
    axes[1].set_title('Validation R²')
    axes[1].legend()

    # Parameter efficiency
    param_data = [(labels_map[m], results[m]['trainable_params'], results[m]['final_r2']) for m in methods]
    x_pos = range(len(param_data))
    bars = axes[2].bar(x_pos, [d[1] for d in param_data], color=colors, alpha=0.7)
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels([d[0] for d in param_data], rotation=25, ha='right', fontsize=9)
    axes[2].set_ylabel('Trainable Parameters')
    axes[2].set_title('Parameter Efficiency')
    for bar, d in zip(bars, param_data):
        axes[2].text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                     f'R²={d[2]:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.suptitle('Fine-tuning Strategy Comparison for Enzyme Activity Prediction', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task2_finetuning_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Scatter plot: predicted vs actual
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
    all_preds_list = [preds_frozen, preds_lora, preds_adapter, preds_full]
    for ax, m, preds, c in zip(axes, methods, all_preds_list, colors):
        ax.scatter(labels_frozen, preds, alpha=0.6, color=c, s=30)
        ax.plot([min(labels_frozen), max(labels_frozen)],
                [min(labels_frozen), max(labels_frozen)], 'k--', alpha=0.5)
        ax.set_xlabel('True Activity')
        ax.set_ylabel('Predicted Activity')
        ax.set_title(f'{labels_map[m]}\nR²={results[m]["final_r2"]:.3f}')
    plt.suptitle('Predicted vs True Enzyme Activity', fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/task2_scatter_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Save results
    save_results = {}
    for m in methods:
        save_results[m] = {
            'trainable_params': results[m]['trainable_params'],
            'total_params': results[m]['total_params'],
            'param_ratio': results[m]['trainable_params'] / results[m]['total_params'],
            'final_r2': results[m]['final_r2'],
            'final_spearman': results[m]['final_spearman'],
            'final_train_loss': results[m]['history']['train_loss'][-1],
            'final_val_loss': results[m]['history']['val_loss'][-1]
        }

    save_results['metadata'] = {
        'model': MODEL_NAME, 'n_samples': N_SAMPLES, 'n_epochs': N_EPOCHS,
        'batch_size': BATCH_SIZE, 'seed': SEED,
        'elapsed_seconds': round(time.time() - start, 2)
    }

    with open(f'{RESULTS_DIR}/task2_results.json', 'w') as f:
        json.dump(save_results, f, indent=2)

    print(f"\n{'='*60}")
    print("Summary:")
    for m in methods:
        print(f"  {labels_map[m]:25s} | Params: {save_results[m]['trainable_params']:>8,} "
              f"({100*save_results[m]['param_ratio']:.2f}%) | R²: {save_results[m]['final_r2']:.3f} "
              f"| ρ: {save_results[m]['final_spearman']:.3f}")
    print(f"Elapsed: {save_results['metadata']['elapsed_seconds']:.1f}s")

if __name__ == '__main__':
    run_finetuning_experiment()
