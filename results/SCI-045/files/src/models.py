"""
Epigenetic Clock Models: Baseline and Deep Learning approaches.
Implements ElasticNet (Horvath-like), tissue-aware models, and neural network clocks.
"""
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from torch.utils.data import DataLoader, TensorDataset


class HorvathBaseline:
    """ElasticNet-based clock (Horvath-like baseline)."""
    def __init__(self, alpha=0.1, l1_ratio=0.5):
        self.model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=5000)
        self.scaler = StandardScaler()
    
    def fit(self, X, y):
        X_scaled = self.scaler.fit_transform(X)
        self.model.fit(X_scaled, y)
        return self
    
    def predict(self, X):
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)
    
    def get_selected_cpgs(self):
        return np.where(np.abs(self.model.coef_) > 0)[0]


class TissueAwareClock:
    """Tissue-specific clock with separate models per tissue type."""
    def __init__(self):
        self.models = {}
        self.scalers = {}
    
    def fit(self, X, y, tissues):
        unique_tissues = np.unique(tissues)
        for tissue in unique_tissues:
            mask = tissues == tissue
            self.scalers[tissue] = StandardScaler()
            X_tissue = self.scalers[tissue].fit_transform(X[mask])
            self.models[tissue] = ElasticNet(alpha=0.05, l1_ratio=0.5, max_iter=10000)
            self.models[tissue].fit(X_tissue, y[mask])
        return self
    
    def predict(self, X, tissues):
        predictions = np.zeros(len(X))
        for tissue in np.unique(tissues):
            mask = tissues == tissue
            X_tissue = self.scalers[tissue].transform(X[mask])
            predictions[mask] = self.models[tissue].predict(X_tissue)
        return predictions


class DeepEpiClock(nn.Module):
    """Deep neural network-based epigenetic clock with attention."""
    def __init__(self, n_features, n_tissues=5):
        super().__init__()
        self.tissue_embedding = nn.Embedding(n_tissues, 32)
        
        self.feature_encoder = nn.Sequential(
            nn.Linear(n_features, 512),
            nn.BatchNorm1d(512),
            nn.GELU(),
            nn.Dropout(0.2),
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.GELU(),
            nn.Dropout(0.15),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.GELU(),
        )
        
        # Attention mechanism for CpG importance
        self.attention = nn.Sequential(
            nn.Linear(128, 64),
            nn.Tanh(),
            nn.Linear(64, 128),
            nn.Softmax(dim=1),
        )
        
        self.predictor = nn.Sequential(
            nn.Linear(128 + 32, 64),
            nn.GELU(),
            nn.Dropout(0.1),
            nn.Linear(64, 32),
            nn.GELU(),
            nn.Linear(32, 1),
        )
    
    def forward(self, x_cpg, x_tissue):
        encoded = self.feature_encoder(x_cpg)
        attn_weights = self.attention(encoded)
        attended = encoded * attn_weights
        tissue_emb = self.tissue_embedding(x_tissue)
        combined = torch.cat([attended, tissue_emb], dim=1)
        return self.predictor(combined).squeeze(1)


class DeepClockTrainer:
    """Training pipeline for DeepEpiClock."""
    def __init__(self, n_features, n_tissues=5, lr=0.001, epochs=100, batch_size=64):
        self.model = DeepEpiClock(n_features, n_tissues)
        self.lr = lr
        self.epochs = epochs
        self.batch_size = batch_size
        self.scaler = StandardScaler()
        self.le = LabelEncoder()
        self.train_losses = []
        self.val_losses = []
    
    def fit(self, X, y, tissues, X_val=None, y_val=None, tissues_val=None):
        X_scaled = self.scaler.fit_transform(X)
        tissue_encoded = self.le.fit_transform(tissues)
        
        X_t = torch.FloatTensor(X_scaled)
        y_t = torch.FloatTensor(y)
        tissue_t = torch.LongTensor(tissue_encoded)
        
        dataset = TensorDataset(X_t, tissue_t, y_t)
        loader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True)
        
        optimizer = torch.optim.AdamW(self.model.parameters(), lr=self.lr, weight_decay=1e-3)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=self.epochs)
        criterion = nn.MSELoss()
        
        # Validation data
        has_val = X_val is not None
        if has_val:
            X_val_s = self.scaler.transform(X_val)
            tissue_val_enc = self.le.transform(tissues_val)
            X_val_t = torch.FloatTensor(X_val_s)
            y_val_t = torch.FloatTensor(y_val)
            tissue_val_t = torch.LongTensor(tissue_val_enc)
        
        self.model.train()
        for epoch in range(self.epochs):
            epoch_loss = 0
            for batch_x, batch_tissue, batch_y in loader:
                optimizer.zero_grad()
                pred = self.model(batch_x, batch_tissue)
                loss = criterion(pred, batch_y)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()
            
            avg_loss = epoch_loss / len(loader)
            self.train_losses.append(avg_loss)
            
            if has_val:
                self.model.eval()
                with torch.no_grad():
                    val_pred = self.model(X_val_t, tissue_val_t)
                    val_loss = criterion(val_pred, y_val_t).item()
                self.val_losses.append(val_loss)
                self.model.train()
            
            scheduler.step()
        
        return self
    
    def predict(self, X, tissues):
        self.model.eval()
        X_scaled = self.scaler.transform(X)
        tissue_encoded = self.le.transform(tissues)
        
        X_t = torch.FloatTensor(X_scaled)
        tissue_t = torch.LongTensor(tissue_encoded)
        
        with torch.no_grad():
            predictions = self.model(X_t, tissue_t).numpy()
        return predictions


class GradientBoostClock:
    """Gradient Boosting based clock for comparison."""
    def __init__(self):
        self.model = GradientBoostingRegressor(
            n_estimators=200, max_depth=5, learning_rate=0.05, subsample=0.8
        )
        self.scaler = StandardScaler()
    
    def fit(self, X, y):
        X_scaled = self.scaler.fit_transform(X)
        self.model.fit(X_scaled, y)
        return self
    
    def predict(self, X):
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)


def evaluate_model(y_true, y_pred):
    """Calculate evaluation metrics."""
    return {
        'MAE': mean_absolute_error(y_true, y_pred),
        'RMSE': np.sqrt(mean_squared_error(y_true, y_pred)),
        'R2': r2_score(y_true, y_pred),
        'Pearson_r': np.corrcoef(y_true, y_pred)[0, 1],
        'Median_AE': np.median(np.abs(y_true - y_pred)),
    }
