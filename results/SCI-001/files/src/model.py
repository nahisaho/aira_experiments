"""
CRISPR-Cas9 Off-Target Prediction Model: CNN + Attention Architecture
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Dict, Optional, Tuple


class MultiScaleCNNBlock(nn.Module):
    """Multi-scale 1D CNN for capturing sequence motifs at different resolutions.
    
    Uses parallel convolutions with kernel sizes 3, 5, 7 to capture
    local and extended sequence patterns.
    """
    
    def __init__(self, in_channels: int, out_channels: int, dropout: float = 0.2):
        super().__init__()
        # Three parallel conv branches
        self.conv3 = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=3, padding=1),
            nn.BatchNorm1d(out_channels),
            nn.ReLU(),
        )
        self.conv5 = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=5, padding=2),
            nn.BatchNorm1d(out_channels),
            nn.ReLU(),
        )
        self.conv7 = nn.Sequential(
            nn.Conv1d(in_channels, out_channels, kernel_size=7, padding=3),
            nn.BatchNorm1d(out_channels),
            nn.ReLU(),
        )
        # Fusion
        self.fusion = nn.Sequential(
            nn.Conv1d(out_channels * 3, out_channels, kernel_size=1),
            nn.BatchNorm1d(out_channels),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (batch, channels, seq_len)
        Returns:
            (batch, out_channels, seq_len)
        """
        h3 = self.conv3(x)
        h5 = self.conv5(x)
        h7 = self.conv7(x)
        h = torch.cat([h3, h5, h7], dim=1)
        return self.fusion(h)


class PositionalEncoding(nn.Module):
    """Sinusoidal positional encoding for sequence positions."""
    
    def __init__(self, d_model: int, max_len: int = 30):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-np.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term[:d_model // 2])
        pe = pe.unsqueeze(0)  # (1, max_len, d_model)
        self.register_buffer('pe', pe)
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Add positional encoding.
        Args: x of shape (batch, seq_len, d_model)
        """
        return x + self.pe[:, :x.size(1), :]


class MultiHeadSelfAttention(nn.Module):
    """Multi-head self-attention for capturing long-range dependencies
    in guide-target alignment.
    """
    
    def __init__(self, d_model: int, n_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        assert d_model % n_heads == 0
        self.d_model = d_model
        self.n_heads = n_heads
        self.d_k = d_model // n_heads
        
        self.W_q = nn.Linear(d_model, d_model)
        self.W_k = nn.Linear(d_model, d_model)
        self.W_v = nn.Linear(d_model, d_model)
        self.W_o = nn.Linear(d_model, d_model)
        
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(d_model)
        
        # Store attention weights for interpretability
        self.attention_weights = None
    
    def forward(self, x: torch.Tensor, 
                mask: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Args:
            x: (batch, seq_len, d_model)
            mask: optional attention mask
        Returns:
            (batch, seq_len, d_model)
        """
        batch_size, seq_len, _ = x.shape
        residual = x
        
        Q = self.W_q(x).view(batch_size, seq_len, self.n_heads, self.d_k).transpose(1, 2)
        K = self.W_k(x).view(batch_size, seq_len, self.n_heads, self.d_k).transpose(1, 2)
        V = self.W_v(x).view(batch_size, seq_len, self.n_heads, self.d_k).transpose(1, 2)
        
        # Scaled dot-product attention
        scores = torch.matmul(Q, K.transpose(-2, -1)) / np.sqrt(self.d_k)
        if mask is not None:
            scores = scores.masked_fill(mask == 0, -1e9)
        
        attn = F.softmax(scores, dim=-1)
        self.attention_weights = attn.detach()
        attn = self.dropout(attn)
        
        context = torch.matmul(attn, V)
        context = context.transpose(1, 2).contiguous().view(batch_size, seq_len, self.d_model)
        
        output = self.W_o(context)
        output = self.dropout(output)
        output = self.layer_norm(output + residual)
        
        return output


class GuideTargetCrossAttention(nn.Module):
    """Cross-attention between guide RNA and target DNA representations.
    
    Captures position-specific interactions between guide and target,
    emphasizing mismatch positions.
    """
    
    def __init__(self, d_model: int, n_heads: int = 4, dropout: float = 0.1):
        super().__init__()
        self.d_model = d_model
        self.n_heads = n_heads
        self.d_k = d_model // n_heads
        
        self.W_q = nn.Linear(d_model, d_model)
        self.W_k = nn.Linear(d_model, d_model)
        self.W_v = nn.Linear(d_model, d_model)
        self.W_o = nn.Linear(d_model, d_model)
        
        self.dropout = nn.Dropout(dropout)
        self.layer_norm = nn.LayerNorm(d_model)
        
        self.cross_attention_weights = None
    
    def forward(self, query: torch.Tensor, 
                key_value: torch.Tensor) -> torch.Tensor:
        """
        Args:
            query: guide features (batch, guide_len, d_model)
            key_value: target features (batch, target_len, d_model)
        Returns:
            (batch, guide_len, d_model)
        """
        batch_size = query.size(0)
        guide_len = query.size(1)
        target_len = key_value.size(1)
        
        Q = self.W_q(query).view(batch_size, guide_len, self.n_heads, self.d_k).transpose(1, 2)
        K = self.W_k(key_value).view(batch_size, target_len, self.n_heads, self.d_k).transpose(1, 2)
        V = self.W_v(key_value).view(batch_size, target_len, self.n_heads, self.d_k).transpose(1, 2)
        
        scores = torch.matmul(Q, K.transpose(-2, -1)) / np.sqrt(self.d_k)
        attn = F.softmax(scores, dim=-1)
        self.cross_attention_weights = attn.detach()
        attn = self.dropout(attn)
        
        context = torch.matmul(attn, V)
        context = context.transpose(1, 2).contiguous().view(batch_size, guide_len, self.d_model)
        
        output = self.W_o(context)
        output = self.dropout(output)
        output = self.layer_norm(output + query)
        
        return output


class EpigeneticEncoder(nn.Module):
    """MLP encoder for epigenetic features."""
    
    def __init__(self, input_dim: int = 7, hidden_dim: int = 32, 
                 output_dim: int = 64, dropout: float = 0.2):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim),
            nn.ReLU(),
        )
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args: x of shape (batch, input_dim)
        Returns: (batch, output_dim)
        """
        return self.mlp(x)


class CRISPROffTargetNet(nn.Module):
    """Main model: CNN + Multi-Head Attention for CRISPR off-target prediction.
    
    Architecture:
    1. Separate CNN encoders for guide RNA and target DNA
    2. Mismatch pattern CNN encoder
    3. Positional encoding
    4. Self-attention on concatenated features
    5. Cross-attention between guide and target
    6. Epigenetic feature integration via gated fusion
    7. Classification head
    
    Input:
        guide_onehot: (B, 4, 20)
        target_onehot: (B, 4, 23)
        mismatch_features: (B, 14, 20)
        pam_encoding: (B, 4, 3)
        epigenetic_features: (B, 7)
    
    Output:
        logits: (B, 1) off-target cleavage probability
    """
    
    def __init__(self, 
                 cnn_channels: int = 64,
                 d_model: int = 128,
                 n_attention_heads: int = 4,
                 n_attention_layers: int = 2,
                 epigenetic_dim: int = 7,
                 dropout: float = 0.2,
                 use_epigenetics: bool = True):
        super().__init__()
        self.use_epigenetics = use_epigenetics
        self.d_model = d_model
        
        # Guide RNA encoder
        self.guide_cnn = nn.Sequential(
            MultiScaleCNNBlock(4, cnn_channels, dropout),
            MultiScaleCNNBlock(cnn_channels, cnn_channels, dropout),
        )
        
        # Target DNA encoder
        self.target_cnn = nn.Sequential(
            MultiScaleCNNBlock(4, cnn_channels, dropout),
            MultiScaleCNNBlock(cnn_channels, cnn_channels, dropout),
        )
        
        # Mismatch pattern encoder
        self.mismatch_cnn = nn.Sequential(
            MultiScaleCNNBlock(14, cnn_channels, dropout),
            MultiScaleCNNBlock(cnn_channels, cnn_channels, dropout),
        )
        
        # PAM encoder
        self.pam_encoder = nn.Sequential(
            nn.Conv1d(4, cnn_channels, kernel_size=3, padding=1),
            nn.BatchNorm1d(cnn_channels),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
            nn.Linear(cnn_channels, d_model),
            nn.ReLU(),
        )
        
        # Project CNN outputs to d_model
        self.guide_proj = nn.Linear(cnn_channels, d_model)
        self.target_proj = nn.Linear(cnn_channels, d_model)
        self.mismatch_proj = nn.Linear(cnn_channels, d_model)
        
        # Positional encoding
        self.pos_encoding = PositionalEncoding(d_model, max_len=30)
        
        # Self-attention layers
        self.self_attention_layers = nn.ModuleList([
            MultiHeadSelfAttention(d_model, n_attention_heads, dropout)
            for _ in range(n_attention_layers)
        ])
        
        # Cross-attention: guide queries target
        self.cross_attention = GuideTargetCrossAttention(
            d_model, n_attention_heads, dropout)
        
        # Epigenetic encoder
        if use_epigenetics:
            self.epigenetic_encoder = EpigeneticEncoder(
                epigenetic_dim, 32, d_model, dropout)
            # Gated fusion for epigenetic features
            self.epigenetic_gate = nn.Sequential(
                nn.Linear(d_model * 2, d_model),
                nn.Sigmoid(),
            )
        
        # Classification head
        classifier_input = d_model * 2  # pooled features + PAM
        if use_epigenetics:
            classifier_input += d_model
        
        self.classifier = nn.Sequential(
            nn.Linear(classifier_input, d_model),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, d_model // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(d_model // 2, 1),
        )
        
        # Initialize weights
        self._init_weights()
    
    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
    
    def forward(self, guide_onehot: torch.Tensor,
                target_onehot: torch.Tensor,
                mismatch_features: torch.Tensor,
                pam_encoding: torch.Tensor,
                epigenetic_features: Optional[torch.Tensor] = None
                ) -> Tuple[torch.Tensor, Dict[str, torch.Tensor]]:
        """
        Forward pass.
        
        Returns:
            logits: (B, 1) prediction scores
            attention_info: dict of attention weights for interpretability
        """
        batch_size = guide_onehot.size(0)
        
        # === CNN Feature Extraction ===
        guide_feat = self.guide_cnn(guide_onehot)       # (B, C, 20)
        target_feat = self.target_cnn(target_onehot)     # (B, C, 23)
        mismatch_feat = self.mismatch_cnn(mismatch_features)  # (B, C, 20)
        
        # Transpose for attention: (B, seq_len, C)
        guide_feat = guide_feat.transpose(1, 2)
        target_feat = target_feat.transpose(1, 2)
        mismatch_feat = mismatch_feat.transpose(1, 2)
        
        # Project to d_model
        guide_feat = self.guide_proj(guide_feat)      # (B, 20, d_model)
        target_feat = self.target_proj(target_feat)    # (B, 23, d_model)
        mismatch_feat = self.mismatch_proj(mismatch_feat)  # (B, 20, d_model)
        
        # Add positional encoding
        guide_feat = self.pos_encoding(guide_feat)
        target_feat = self.pos_encoding(target_feat)
        
        # Combine guide and mismatch features
        combined = guide_feat + mismatch_feat  # (B, 20, d_model)
        
        # === Self-Attention ===
        for attn_layer in self.self_attention_layers:
            combined = attn_layer(combined)
        
        # === Cross-Attention (guide queries target) ===
        cross_out = self.cross_attention(combined, target_feat)  # (B, 20, d_model)
        
        # === Global Pooling ===
        # Average pooling + max pooling
        avg_pool = cross_out.mean(dim=1)   # (B, d_model)
        max_pool = cross_out.max(dim=1)[0]  # (B, d_model)
        pooled = avg_pool + max_pool        # (B, d_model)
        
        # PAM features
        pam_feat = self.pam_encoder(pam_encoding)  # (B, d_model)
        
        # === Feature Fusion ===
        features = [pooled, pam_feat]
        
        if self.use_epigenetics and epigenetic_features is not None:
            epi_feat = self.epigenetic_encoder(epigenetic_features)  # (B, d_model)
            # Gated fusion
            gate_input = torch.cat([pooled, epi_feat], dim=1)
            gate = self.epigenetic_gate(gate_input)
            gated_epi = gate * epi_feat
            features.append(gated_epi)
        
        # Concatenate all features
        final_features = torch.cat(features, dim=1)
        
        # === Classification ===
        logits = self.classifier(final_features)
        
        # Collect attention weights for interpretability
        attention_info = {
            'self_attention': [
                layer.attention_weights 
                for layer in self.self_attention_layers
                if layer.attention_weights is not None
            ],
            'cross_attention': self.cross_attention.cross_attention_weights,
        }
        
        return logits, attention_info
    
    def get_num_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters())
    
    def get_trainable_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


def build_model(config: Optional[Dict] = None) -> CRISPROffTargetNet:
    """Build model with configuration."""
    default_config = {
        'cnn_channels': 64,
        'd_model': 128,
        'n_attention_heads': 4,
        'n_attention_layers': 2,
        'epigenetic_dim': 7,
        'dropout': 0.2,
        'use_epigenetics': True,
    }
    if config:
        default_config.update(config)
    
    model = CRISPROffTargetNet(**default_config)
    return model


if __name__ == '__main__':
    print("=== CRISPROffTargetNet Architecture ===")
    
    model = build_model()
    print(f"Total parameters: {model.get_num_parameters():,}")
    print(f"Trainable parameters: {model.get_trainable_parameters():,}")
    
    # Test forward pass
    batch_size = 4
    guide = torch.randn(batch_size, 4, 20)
    target = torch.randn(batch_size, 4, 23)
    mismatch = torch.randn(batch_size, 14, 20)
    pam = torch.randn(batch_size, 4, 3)
    epi = torch.randn(batch_size, 7)
    
    logits, attn_info = model(guide, target, mismatch, pam, epi)
    print(f"\nOutput shape: {logits.shape}")
    print(f"Self-attention layers: {len(attn_info['self_attention'])}")
    
    probs = torch.sigmoid(logits)
    print(f"Prediction probabilities: {probs.detach().numpy().flatten()}")
    print("\n✓ Model forward pass successful.")
