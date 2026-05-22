"""
Deep Learning Epigenetic Clock using PyTorch.
"""
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import mean_absolute_error, r2_score
from scipy import stats


class CpGAttention(nn.Module):
    def __init__(self, d_model, n_heads=4):
        super().__init__()
        self.attn = nn.MultiheadAttention(d_model, n_heads, batch_first=True)
        self.norm = nn.LayerNorm(d_model)

    def forward(self, x):
        attn_out, weights = self.attn(x, x, x)
        return self.norm(x + attn_out), weights


class ResBlock(nn.Module):
    def __init__(self, in_dim, out_dim, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, out_dim), nn.BatchNorm1d(out_dim),
            nn.GELU(), nn.Dropout(dropout),
        )
        self.skip = nn.Linear(in_dim, out_dim) if in_dim != out_dim else nn.Identity()

    def forward(self, x):
        return self.net(x) + self.skip(x)


class DeepEpigeneticClock(nn.Module):
    def __init__(self, n_cpg, n_tissues=5, embed_dim=16, hidden_dims=None,
                 dropout=0.3, use_attention=True):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128, 64]
        self.use_attention = use_attention
        self.tissue_embedding = nn.Embedding(n_tissues, embed_dim)
        self.cpg_proj = nn.Sequential(
            nn.Linear(n_cpg, hidden_dims[0]), nn.BatchNorm1d(hidden_dims[0]),
            nn.GELU(), nn.Dropout(dropout),
        )
        if use_attention:
            self.reshape_dim = 32
            self.n_tokens = hidden_dims[0] // self.reshape_dim
            self.attention = CpGAttention(self.reshape_dim, n_heads=4)
        feat_dim = hidden_dims[0] + embed_dim + 1
        layers = []
        in_dim = feat_dim
        for h in hidden_dims[1:]:
            layers.append(ResBlock(in_dim, h, dropout))
            in_dim = h
        self.mlp = nn.Sequential(*layers)
        self.head = nn.Linear(hidden_dims[-1], 1)

    def forward(self, cpg, tissue_id, sex):
        h = self.cpg_proj(cpg)
        if self.use_attention:
            h_r = h.view(h.size(0), self.n_tokens, self.reshape_dim)
            h_a, _ = self.attention(h_r)
            h = h_a.view(h.size(0), -1)
        t_emb = self.tissue_embedding(tissue_id)
        combined = torch.cat([h, t_emb, sex.unsqueeze(1)], dim=1)
        return self.head(self.mlp(combined)).squeeze(1)


class DeepClockTrainer:
    def __init__(self, n_cpg, n_tissues=5, lr=1e-3, weight_decay=1e-4,
                 epochs=100, batch_size=64, patience=15, device=None):
        self.n_cpg = n_cpg
        self.n_tissues = n_tissues
        self.lr = lr
        self.weight_decay = weight_decay
        self.epochs = epochs
        self.batch_size = batch_size
        self.patience = patience
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.scaler = StandardScaler()
        self.tissue_encoder = LabelEncoder()
        self.history = {"train_loss": [], "val_loss": []}

    def _prepare_data(self, df, fit=True):
        cpg_cols = [c for c in df.columns if c.startswith("cg")]
        X = df[cpg_cols].values
        if fit:
            X = self.scaler.fit_transform(X)
            self.tissue_encoder.fit(df["tissue"].values)
        else:
            X = self.scaler.transform(X)
        tissue = self.tissue_encoder.transform(df["tissue"].values)
        sex = df["sex"].values.astype(np.float32)
        age = df["true_bio_age"].values.astype(np.float32)
        return (torch.tensor(X, dtype=torch.float32),
                torch.tensor(tissue, dtype=torch.long),
                torch.tensor(sex, dtype=torch.float32),
                torch.tensor(age, dtype=torch.float32))

    def fit(self, train_df, val_df=None):
        X_t, tissue_t, sex_t, y_t = self._prepare_data(train_df, fit=True)
        train_loader = DataLoader(TensorDataset(X_t, tissue_t, sex_t, y_t),
                                  batch_size=self.batch_size, shuffle=True)
        val_loader = None
        if val_df is not None:
            X_v, tissue_v, sex_v, y_v = self._prepare_data(val_df, fit=False)
            val_loader = DataLoader(TensorDataset(X_v, tissue_v, sex_v, y_v),
                                    batch_size=self.batch_size)

        self.model = DeepEpigeneticClock(self.n_cpg, self.n_tissues).to(self.device)
        optimizer = optim.AdamW(self.model.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=self.epochs)
        criterion = nn.HuberLoss(delta=5.0)
        best_val_loss, patience_counter = float("inf"), 0

        for epoch in range(self.epochs):
            self.model.train()
            losses = []
            for X_b, t_b, s_b, y_b in train_loader:
                X_b, t_b, s_b, y_b = [v.to(self.device) for v in [X_b, t_b, s_b, y_b]]
                optimizer.zero_grad()
                loss = criterion(self.model(X_b, t_b, s_b), y_b)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                optimizer.step()
                losses.append(loss.item())
            scheduler.step()
            self.history["train_loss"].append(np.mean(losses))

            if val_loader:
                val_loss = self._eval_loss(val_loader, criterion)
                self.history["val_loss"].append(val_loss)
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    patience_counter = 0
                    self.best_state = {k: v.cpu().clone() for k, v in self.model.state_dict().items()}
                else:
                    patience_counter += 1
                    if patience_counter >= self.patience:
                        print(f"  Early stopping at epoch {epoch+1}")
                        break
                if (epoch + 1) % 20 == 0:
                    print(f"  Epoch {epoch+1}: train={np.mean(losses):.4f}, val={val_loss:.4f}")

        if hasattr(self, "best_state"):
            self.model.load_state_dict(self.best_state)
        return self

    def _eval_loss(self, loader, criterion):
        self.model.eval()
        losses = []
        with torch.no_grad():
            for X_b, t_b, s_b, y_b in loader:
                X_b, t_b, s_b, y_b = [v.to(self.device) for v in [X_b, t_b, s_b, y_b]]
                losses.append(criterion(self.model(X_b, t_b, s_b), y_b).item())
        return np.mean(losses)

    def predict(self, df):
        X, tissue, sex, _ = self._prepare_data(df, fit=False)
        loader = DataLoader(TensorDataset(X, tissue, sex), batch_size=self.batch_size)
        self.model.eval()
        preds = []
        with torch.no_grad():
            for X_b, t_b, s_b in loader:
                X_b, t_b, s_b = [v.to(self.device) for v in [X_b, t_b, s_b]]
                preds.append(self.model(X_b, t_b, s_b).cpu().numpy())
        return np.concatenate(preds)
