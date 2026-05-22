"""
EEG Conformer: Time-Series Transformer for EEG Decoding.
Combines convolutional local feature extraction with Transformer self-attention
for global temporal context modeling.

Reference: Song et al. (2023), EEG Conformer: Convolutional Transformer for EEG
           Signal Decoding and Generation.
"""

import numpy as np
from scipy import linalg
from typing import Optional, Tuple, List, Dict
import time
import math


# ---------------------------------------------------------------------------
# Pure NumPy Transformer Components
# ---------------------------------------------------------------------------

def gelu(x: np.ndarray) -> np.ndarray:
    """Gaussian Error Linear Unit activation."""
    return 0.5 * x * (1.0 + np.tanh(math.sqrt(2.0 / math.pi) * (x + 0.044715 * x ** 3)))


def layer_norm(x: np.ndarray, weight: np.ndarray, bias: np.ndarray,
               eps: float = 1e-5) -> np.ndarray:
    """Layer normalization over last dimension."""
    mean = x.mean(axis=-1, keepdims=True)
    var = ((x - mean) ** 2).mean(axis=-1, keepdims=True)
    return weight * (x - mean) / np.sqrt(var + eps) + bias


def scaled_dot_product_attention(Q: np.ndarray, K: np.ndarray, V: np.ndarray,
                                  mask: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Scaled dot-product attention.
    Q, K, V: (..., seq_len, d_k)
    Returns: (..., seq_len, d_v)
    """
    d_k = Q.shape[-1]
    scores = Q @ K.swapaxes(-2, -1) / math.sqrt(d_k)  # (..., seq_len, seq_len)
    if mask is not None:
        scores = np.where(mask, scores, -1e9)
    attn_weights = np.exp(scores - scores.max(axis=-1, keepdims=True))
    attn_weights /= attn_weights.sum(axis=-1, keepdims=True) + 1e-8
    return attn_weights @ V, attn_weights


class MultiHeadAttention:
    """Multi-head self-attention (inference only, no backprop)."""

    def __init__(self, d_model: int, n_heads: int):
        self.d_model = d_model
        self.n_heads = n_heads
        self.d_k = d_model // n_heads
        # Parameters (initialized randomly; in practice loaded from PyTorch checkpoint)
        rng = np.random.RandomState(42)
        scale = math.sqrt(2.0 / (d_model + self.d_k))
        self.W_q = rng.randn(d_model, d_model).astype(np.float32) * scale
        self.W_k = rng.randn(d_model, d_model).astype(np.float32) * scale
        self.W_v = rng.randn(d_model, d_model).astype(np.float32) * scale
        self.W_o = rng.randn(d_model, d_model).astype(np.float32) * scale
        self.b_q = np.zeros(d_model, dtype=np.float32)
        self.b_k = np.zeros(d_model, dtype=np.float32)
        self.b_v = np.zeros(d_model, dtype=np.float32)
        self.b_o = np.zeros(d_model, dtype=np.float32)
        self.last_attn_weights: Optional[np.ndarray] = None

    def forward(self, x: np.ndarray,
                mask: Optional[np.ndarray] = None) -> np.ndarray:
        """x: (batch, seq_len, d_model) -> (batch, seq_len, d_model)"""
        B, T, D = x.shape
        Q = x @ self.W_q + self.b_q  # (B, T, D)
        K = x @ self.W_k + self.b_k
        V = x @ self.W_v + self.b_v

        # Reshape to (B, n_heads, T, d_k)
        Q = Q.reshape(B, T, self.n_heads, self.d_k).transpose(0, 2, 1, 3)
        K = K.reshape(B, T, self.n_heads, self.d_k).transpose(0, 2, 1, 3)
        V = V.reshape(B, T, self.n_heads, self.d_k).transpose(0, 2, 1, 3)

        out, self.last_attn_weights = scaled_dot_product_attention(Q, K, V, mask)
        # (B, n_heads, T, d_k) -> (B, T, D)
        out = out.transpose(0, 2, 1, 3).reshape(B, T, D)
        return out @ self.W_o + self.b_o


class FeedForward:
    """Two-layer feedforward network."""

    def __init__(self, d_model: int, d_ff: int):
        rng = np.random.RandomState(43)
        scale_1 = math.sqrt(2.0 / (d_model + d_ff))
        scale_2 = math.sqrt(2.0 / (d_ff + d_model))
        self.W1 = rng.randn(d_model, d_ff).astype(np.float32) * scale_1
        self.b1 = np.zeros(d_ff, dtype=np.float32)
        self.W2 = rng.randn(d_ff, d_model).astype(np.float32) * scale_2
        self.b2 = np.zeros(d_model, dtype=np.float32)

    def forward(self, x: np.ndarray) -> np.ndarray:
        return gelu(x @ self.W1 + self.b1) @ self.W2 + self.b2


class TransformerEncoderBlock:
    """Single Transformer encoder block with pre-LN."""

    def __init__(self, d_model: int, n_heads: int, d_ff: int):
        self.attn = MultiHeadAttention(d_model, n_heads)
        self.ff = FeedForward(d_model, d_ff)
        self.ln1_w = np.ones(d_model, dtype=np.float32)
        self.ln1_b = np.zeros(d_model, dtype=np.float32)
        self.ln2_w = np.ones(d_model, dtype=np.float32)
        self.ln2_b = np.zeros(d_model, dtype=np.float32)

    def forward(self, x: np.ndarray,
                mask: Optional[np.ndarray] = None) -> np.ndarray:
        """Pre-LN Transformer block."""
        # Self-attention sub-layer
        x_norm = layer_norm(x, self.ln1_w, self.ln1_b)
        x = x + self.attn.forward(x_norm, mask)
        # FFN sub-layer
        x_norm = layer_norm(x, self.ln2_w, self.ln2_b)
        x = x + self.ff.forward(x_norm)
        return x


# ---------------------------------------------------------------------------
# EEG Conformer Architecture
# ---------------------------------------------------------------------------

class EEGConformerConfig:
    """EEG Conformer hyperparameters."""
    n_channels: int = 64
    n_times: int = 256           # temporal samples per trial
    n_classes: int = 4
    # Convolutional frontend
    n_temporal_filters: int = 40
    temporal_kernel: int = 25
    n_spatial_filters: int = 40  # = n_temporal_filters
    pool_size: int = 75
    pool_stride: int = 15
    # Transformer
    d_model: int = 40            # = n_spatial_filters
    n_heads: int = 8
    n_transformer_layers: int = 6
    d_ff: int = 160              # 4 * d_model
    dropout: float = 0.5
    # Positional encoding
    max_seq_len: int = 512


def sinusoidal_positional_encoding(seq_len: int, d_model: int) -> np.ndarray:
    """Standard sinusoidal positional encoding (Vaswani et al., 2017)."""
    PE = np.zeros((seq_len, d_model), dtype=np.float32)
    positions = np.arange(seq_len)[:, None]
    div_term = np.exp(np.arange(0, d_model, 2) * (-math.log(10000.0) / d_model))
    PE[:, 0::2] = np.sin(positions * div_term)
    if d_model % 2 == 1:
        PE[:, 1::2] = np.cos(positions * div_term[:-1])
    else:
        PE[:, 1::2] = np.cos(positions * div_term)
    return PE


class EEGConformerNumpy:
    """
    EEG Conformer: Convolutional Transformer for EEG decoding.
    NumPy inference implementation of the full architecture.

    Architecture:
    1. Temporal Conv (local feature extraction)
    2. Depthwise Spatial Conv
    3. Temporal average pooling
    4. Positional encoding
    5. Transformer encoder (multi-head attention + FFN)
    6. Classification head
    """

    def __init__(self, config: EEGConformerConfig):
        self.config = config
        cfg = config

        # --- Convolutional frontend ---
        rng = np.random.RandomState(0)
        # Layer 1: Temporal conv (1, 1, 1, kern_temporal)
        self.W_temp = rng.randn(cfg.n_temporal_filters, 1, 1, cfg.temporal_kernel
                                ).astype(np.float32) * 0.1
        self.bn1_w = np.ones(cfg.n_temporal_filters, dtype=np.float32)
        self.bn1_b = np.zeros(cfg.n_temporal_filters, dtype=np.float32)
        self.bn1_mean = np.zeros(cfg.n_temporal_filters, dtype=np.float32)
        self.bn1_var = np.ones(cfg.n_temporal_filters, dtype=np.float32)

        # Layer 2: Depthwise spatial conv (n_ch, 1)
        self.W_spat = rng.randn(cfg.n_spatial_filters, cfg.n_channels, 1
                                ).astype(np.float32) * 0.1
        self.bn2_w = np.ones(cfg.n_spatial_filters, dtype=np.float32)
        self.bn2_b = np.zeros(cfg.n_spatial_filters, dtype=np.float32)
        self.bn2_mean = np.zeros(cfg.n_spatial_filters, dtype=np.float32)
        self.bn2_var = np.ones(cfg.n_spatial_filters, dtype=np.float32)

        # --- Transformer ---
        # After pooling: seq_len = (n_times - temporal_kernel + 1 - pool_size) // pool_stride + 1
        seq_after_temp_conv = cfg.n_times - cfg.temporal_kernel + 1
        seq_after_pool = (seq_after_temp_conv - cfg.pool_size) // cfg.pool_stride + 1
        self.seq_len = max(seq_after_pool, 1)

        self.pos_enc = sinusoidal_positional_encoding(cfg.max_seq_len, cfg.d_model)
        self.transformer_blocks = [
            TransformerEncoderBlock(cfg.d_model, cfg.n_heads, cfg.d_ff)
            for _ in range(cfg.n_transformer_layers)
        ]
        self.ln_final_w = np.ones(cfg.d_model, dtype=np.float32)
        self.ln_final_b = np.zeros(cfg.d_model, dtype=np.float32)

        # --- Classification head ---
        flat_dim = self.seq_len * cfg.d_model
        self.W_cls = rng.randn(flat_dim, cfg.n_classes).astype(np.float32) * 0.01
        self.b_cls = np.zeros(cfg.n_classes, dtype=np.float32)

        self._count_parameters()

    def _count_parameters(self) -> int:
        cfg = self.config
        p = 0
        p += cfg.n_temporal_filters * cfg.temporal_kernel          # temporal conv
        p += cfg.n_spatial_filters * cfg.n_channels                # spatial conv
        per_block = (4 * cfg.d_model ** 2 + 2 * cfg.d_model * cfg.d_ff
                     + 4 * cfg.d_model)
        p += cfg.n_transformer_layers * per_block
        p += self.seq_len * cfg.d_model * cfg.n_classes + cfg.n_classes
        self.n_parameters = p
        return p

    def _conv1d_same(self, x: np.ndarray, W: np.ndarray) -> np.ndarray:
        """
        1D convolution over temporal dimension with 'same' padding.
        x: (batch, in_ch, T)
        W: (out_ch, in_ch, kern)
        Returns: (batch, out_ch, T)
        """
        kern = W.shape[2]
        pad = kern // 2
        x_pad = np.pad(x, ((0, 0), (0, 0), (pad, pad)))
        B, _, T_out = x.shape
        out = np.zeros((B, W.shape[0], T_out), dtype=np.float32)
        for oc in range(W.shape[0]):
            for ic in range(W.shape[1]):
                # Sliding correlation (conv without flip for simplicity)
                for b in range(B):
                    out[b, oc] += np.convolve(x_pad[b, ic], W[oc, ic, ::-1], mode='valid')[:T_out]
        return out

    def _batch_norm_inference(self, x: np.ndarray, w: np.ndarray, b: np.ndarray,
                               mean: np.ndarray, var: np.ndarray,
                               eps: float = 1e-5) -> np.ndarray:
        """BN in inference mode: x is (batch, channels, T)."""
        x_norm = (x - mean[None, :, None]) / np.sqrt(var[None, :, None] + eps)
        return w[None, :, None] * x_norm + b[None, :, None]

    def forward(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Forward pass through EEG Conformer.
        X: (batch, n_ch, n_times) -> logits (batch, n_classes)
        Returns: (logits, attention_weights_last_layer)
        """
        B = X.shape[0]
        cfg = self.config

        # --- Temporal convolution ---
        # Treat EEG as (batch, n_ch, n_times); apply temporal filter per channel
        # Simplified: collapse channels after spatial conv, then apply temporal
        x = X.astype(np.float32)   # (B, n_ch, T)

        # Temporal conv across all channels: (B, n_temp_filt, T-kern+1)
        x_flat = x[:, np.newaxis, :, :]  # (B, 1, n_ch, T)
        # Apply temporal kernel along time dimension for each spatial position
        T = x.shape[2]
        kern = cfg.temporal_kernel
        T_conv = T - kern + 1
        x_temp = np.zeros((B, cfg.n_temporal_filters, cfg.n_channels, T_conv), dtype=np.float32)
        for f in range(cfg.n_temporal_filters):
            k = self.W_temp[f, 0, 0]  # (kern,) — shared across channels
            for b in range(B):
                for ch in range(cfg.n_channels):
                    x_temp[b, f, ch] = np.convolve(x[b, ch], k[::-1], mode='valid')

        # Spatial conv: (B, n_spat_filt, 1, T_conv)
        x_spat = np.zeros((B, cfg.n_spatial_filters, T_conv), dtype=np.float32)
        for f in range(cfg.n_spatial_filters):
            w_sp = self.W_spat[f, :, 0]  # (n_ch,) — spatial weights
            for b in range(B):
                x_spat[b, f] = (x_temp[b, f] * w_sp[:, None]).sum(axis=0)

        # BN + ELU
        x_spat = self._batch_norm_inference(x_spat, self.bn2_w, self.bn2_b,
                                             self.bn2_mean, self.bn2_var)
        x_spat = np.where(x_spat > 0, x_spat, 0.01 * x_spat)  # leaky ReLU

        # Average pooling
        pool_size = cfg.pool_size
        pool_stride = cfg.pool_stride
        n_pool_steps = (T_conv - pool_size) // pool_stride + 1
        x_pooled = np.zeros((B, cfg.n_spatial_filters, n_pool_steps), dtype=np.float32)
        for i in range(n_pool_steps):
            start = i * pool_stride
            x_pooled[:, :, i] = x_spat[:, :, start:start + pool_size].mean(axis=2)

        # --- Transformer ---
        # (B, n_spat_filt, seq_len) -> (B, seq_len, d_model)
        seq_len = x_pooled.shape[2]
        x_trans = x_pooled.transpose(0, 2, 1)  # (B, seq_len, d_model)

        # Add positional encoding
        x_trans = x_trans + self.pos_enc[:seq_len][None]

        # Transformer blocks
        for block in self.transformer_blocks:
            x_trans = block.forward(x_trans)

        x_trans = layer_norm(x_trans, self.ln_final_w, self.ln_final_b)

        # Attention weights from last block
        attn_weights = self.transformer_blocks[-1].attn.last_attn_weights

        # Classification head
        x_flat = x_trans.reshape(B, -1)  # (B, seq_len * d_model)
        # Pad or truncate to expected flat_dim
        flat_dim_expected = self.W_cls.shape[0]
        if x_flat.shape[1] < flat_dim_expected:
            pad_size = flat_dim_expected - x_flat.shape[1]
            x_flat = np.pad(x_flat, ((0, 0), (0, pad_size)))
        elif x_flat.shape[1] > flat_dim_expected:
            x_flat = x_flat[:, :flat_dim_expected]

        logits = x_flat @ self.W_cls + self.b_cls  # (B, n_classes)
        return logits, attn_weights

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Returns predicted class indices."""
        logits, _ = self.forward(X)
        return np.argmax(logits, axis=-1)

    def get_architecture_summary(self) -> str:
        cfg = self.config
        return f"""
EEG Conformer Architecture Summary
====================================
Input        : ({cfg.n_channels}, {cfg.n_times})

[Conv Frontend]
Temporal Conv: (F={cfg.n_temporal_filters}, kern={cfg.temporal_kernel}) → ({cfg.n_temporal_filters}, {cfg.n_channels}, {cfg.n_times - cfg.temporal_kernel + 1})
Spatial Conv : (F={cfg.n_spatial_filters}, kern={cfg.n_channels}) → ({cfg.n_spatial_filters}, 1, {cfg.n_times - cfg.temporal_kernel + 1})
AvgPool      : (size={cfg.pool_size}, stride={cfg.pool_stride}) → ({cfg.n_spatial_filters}, {self.seq_len})

[Transformer Encoder]
Positional Encoding: sinusoidal ({self.seq_len}, {cfg.d_model})
{cfg.n_transformer_layers} × TransformerBlock:
  MultiHeadAttention(d_model={cfg.d_model}, n_heads={cfg.n_heads}, d_k={cfg.d_model // cfg.n_heads})
  FeedForward({cfg.d_model} → {cfg.d_ff} → {cfg.d_model})
  LayerNorm (pre-LN)

[Classification Head]
Flatten   : ({self.seq_len} × {cfg.d_model},) = ({self.seq_len * cfg.d_model},)
Linear    : {self.seq_len * cfg.d_model} → {cfg.n_classes}
Softmax

Total parameters: {self.n_parameters:,}
""".strip()


# ---------------------------------------------------------------------------
# Temporal Attention Visualization Helper
# ---------------------------------------------------------------------------

class AttentionAnalyzer:
    """Extract and analyze attention patterns from EEG Conformer."""

    @staticmethod
    def compute_attention_rollout(attn_weights_per_layer: List[np.ndarray]) -> np.ndarray:
        """
        Compute attention rollout (Abnar & Zuidema, 2020).
        attn_weights_per_layer: list of (batch, n_heads, T, T) per layer.
        Returns: (batch, T, T) total attention.
        """
        rollout = attn_weights_per_layer[0].mean(axis=1)  # avg over heads (B, T, T)
        for attn in attn_weights_per_layer[1:]:
            layer_attn = attn.mean(axis=1)
            # Residual connection: A = 0.5 * A + 0.5 * I
            eye = np.eye(rollout.shape[-1])[None]
            layer_attn = 0.5 * layer_attn + 0.5 * eye
            rollout = layer_attn @ rollout
        return rollout

    @staticmethod
    def temporal_importance(rollout: np.ndarray) -> np.ndarray:
        """Sum rollout over query dimension to get per-timestep importance."""
        return rollout.mean(axis=-2)  # (batch, T)


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def demo_eeg_conformer():
    print("=== EEG Conformer Architecture Demo ===\n")
    rng = np.random.RandomState(5)
    cfg = EEGConformerConfig()
    cfg.n_channels = 32
    cfg.n_times = 256
    cfg.n_classes = 4
    cfg.n_transformer_layers = 3
    cfg.n_heads = 4

    model = EEGConformerNumpy(cfg)
    print(model.get_architecture_summary())

    # Benchmark forward pass
    batch_size = 16
    X_test = rng.randn(batch_size, cfg.n_channels, cfg.n_times).astype(np.float32)

    print(f"\nRunning forward pass: batch={batch_size}, ch={cfg.n_channels}, T={cfg.n_times}")
    t0 = time.perf_counter()
    logits, attn = model.forward(X_test)
    elapsed = (time.perf_counter() - t0) * 1000

    print(f"  Output shape   : {logits.shape}")
    print(f"  Forward pass   : {elapsed:.1f} ms ({elapsed/batch_size:.1f} ms/sample)")
    print(f"  Predicted classes: {model.predict(X_test)}")

    # Temporal importance analysis
    if attn is not None:
        analyzer = AttentionAnalyzer()
        rollout = analyzer.compute_attention_rollout([attn])
        importance = analyzer.temporal_importance(rollout)
        print(f"\nAttention rollout shape : {rollout.shape}")
        print(f"Temporal importance shape: {importance.shape}")

    return {
        "n_parameters": model.n_parameters,
        "forward_ms_per_sample": elapsed / batch_size,
        "output_shape": list(logits.shape),
    }


if __name__ == "__main__":
    demo_eeg_conformer()
