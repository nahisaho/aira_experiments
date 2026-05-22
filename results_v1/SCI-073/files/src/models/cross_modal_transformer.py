"""
Module 3: CrossModalTransformer - Tactile-Visual Multimodal Fusion
Fuses tactile and visual modalities via cross-attention transformer
for robust object recognition and manipulation planning.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision import models
import math


class PatchEmbedding(nn.Module):
    """Convert feature maps to patch token sequences."""

    def __init__(self, in_channels: int, embed_dim: int, patch_size: int = 2):
        super().__init__()
        self.proj = nn.Conv2d(in_channels, embed_dim, patch_size, patch_size)
        self.norm = nn.LayerNorm(embed_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """(B, C, H, W) -> (B, N, D)"""
        x = self.proj(x)  # (B, D, H', W')
        B, D, H, W = x.shape
        x = x.flatten(2).transpose(1, 2)  # (B, N, D)
        return self.norm(x)


class CrossAttentionBlock(nn.Module):
    """Bidirectional cross-attention between two modalities."""

    def __init__(self, dim: int, num_heads: int = 8, dropout: float = 0.1):
        super().__init__()
        self.norm_q = nn.LayerNorm(dim)
        self.norm_kv = nn.LayerNorm(dim)
        self.cross_attn = nn.MultiheadAttention(
            dim, num_heads, dropout=dropout, batch_first=True
        )
        self.ffn = nn.Sequential(
            nn.Linear(dim, dim * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim * 4, dim),
            nn.Dropout(dropout),
        )
        self.norm_ffn = nn.LayerNorm(dim)

    def forward(self, query: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        # Cross-attention
        q = self.norm_q(query)
        kv = self.norm_kv(context)
        attn_out, _ = self.cross_attn(q, kv, kv)
        x = query + attn_out
        # FFN
        x = x + self.ffn(self.norm_ffn(x))
        return x


class CrossModalTransformer(nn.Module):
    """
    Tactile-Visual fusion via cross-modal transformer.

    Architecture:
    1. Modality-specific encoders (ResNet-based)
    2. Patch tokenization
    3. Modality-specific positional encoding
    4. Cross-attention layers (tactile ↔ visual)
    5. Fusion & task heads
    """

    def __init__(
        self,
        tactile_backbone: str = "resnet18",
        visual_backbone: str = "resnet34",
        fusion_dim: int = 256,
        num_layers: int = 4,
        num_heads: int = 8,
        dropout: float = 0.1,
        num_object_classes: int = 50,
        num_material_classes: int = 20,
    ):
        super().__init__()
        self.fusion_dim = fusion_dim

        # Tactile encoder
        tact_enc = getattr(models, tactile_backbone)(pretrained=False)
        self.tactile_encoder = nn.Sequential(*list(tact_enc.children())[:-2])
        tact_channels = 512 if "18" in tactile_backbone else 2048

        # Visual encoder
        vis_enc = getattr(models, visual_backbone)(pretrained=False)
        self.visual_encoder = nn.Sequential(*list(vis_enc.children())[:-2])
        vis_channels = 512 if "34" in visual_backbone else 2048

        # Patch embedding
        self.tact_patch = PatchEmbedding(tact_channels, fusion_dim)
        self.vis_patch = PatchEmbedding(vis_channels, fusion_dim)

        # Learnable modality tokens
        self.tact_token = nn.Parameter(torch.randn(1, 1, fusion_dim))
        self.vis_token = nn.Parameter(torch.randn(1, 1, fusion_dim))

        # Positional encodings
        max_tokens = 256
        self.pos_encoding = nn.Parameter(
            self._sinusoidal_encoding(max_tokens, fusion_dim), requires_grad=False
        )

        # Cross-attention layers
        self.cross_layers_t2v = nn.ModuleList([
            CrossAttentionBlock(fusion_dim, num_heads, dropout)
            for _ in range(num_layers)
        ])
        self.cross_layers_v2t = nn.ModuleList([
            CrossAttentionBlock(fusion_dim, num_heads, dropout)
            for _ in range(num_layers)
        ])

        # Self-attention on fused tokens
        self.self_attn = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(
                d_model=fusion_dim, nhead=num_heads,
                dim_feedforward=fusion_dim * 4, dropout=dropout,
                activation='gelu', batch_first=True
            ),
            num_layers=2
        )

        # Task heads
        self.object_head = nn.Sequential(
            nn.LayerNorm(fusion_dim * 2),
            nn.Linear(fusion_dim * 2, fusion_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim, num_object_classes),
        )

        self.material_head = nn.Sequential(
            nn.LayerNorm(fusion_dim * 2),
            nn.Linear(fusion_dim * 2, fusion_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(fusion_dim, num_material_classes),
        )

        self.grasp_quality_head = nn.Sequential(
            nn.LayerNorm(fusion_dim * 2),
            nn.Linear(fusion_dim * 2, 128),
            nn.GELU(),
            nn.Linear(128, 1),
            nn.Sigmoid(),
        )

    @staticmethod
    def _sinusoidal_encoding(max_len: int, dim: int) -> torch.Tensor:
        pe = torch.zeros(max_len, dim)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div_term = torch.exp(
            torch.arange(0, dim, 2).float() * (-math.log(10000.0) / dim)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        return pe.unsqueeze(0)

    def encode_tactile(self, tactile_img: torch.Tensor) -> torch.Tensor:
        features = self.tactile_encoder(tactile_img)
        tokens = self.tact_patch(features)
        B, N, _ = tokens.shape
        cls_token = self.tact_token.expand(B, -1, -1)
        tokens = torch.cat([cls_token, tokens], dim=1)
        tokens = tokens + self.pos_encoding[:, :N + 1, :]
        return tokens

    def encode_visual(self, visual_img: torch.Tensor) -> torch.Tensor:
        features = self.visual_encoder(visual_img)
        tokens = self.vis_patch(features)
        B, N, _ = tokens.shape
        cls_token = self.vis_token.expand(B, -1, -1)
        tokens = torch.cat([cls_token, tokens], dim=1)
        tokens = tokens + self.pos_encoding[:, :N + 1, :]
        return tokens

    def forward(self, tactile_img: torch.Tensor,
                visual_img: torch.Tensor) -> dict:
        """
        Args:
            tactile_img: (B, 3, 224, 224) tactile image
            visual_img: (B, 3, 224, 224) visual/RGB image
        Returns:
            dict with object_logits, material_logits, grasp_quality,
                 tactile_features, visual_features
        """
        tact_tokens = self.encode_tactile(tactile_img)
        vis_tokens = self.encode_visual(visual_img)

        # Bidirectional cross-attention
        for t2v_layer, v2t_layer in zip(self.cross_layers_t2v, self.cross_layers_v2t):
            tact_tokens_new = t2v_layer(tact_tokens, vis_tokens)
            vis_tokens_new = v2t_layer(vis_tokens, tact_tokens)
            tact_tokens = tact_tokens_new
            vis_tokens = vis_tokens_new

        # Extract CLS tokens
        tact_cls = tact_tokens[:, 0]
        vis_cls = vis_tokens[:, 0]

        # Fused representation
        fused = torch.cat([tact_cls, vis_cls], dim=1)

        # Additional self-attention on all tokens
        all_tokens = torch.cat([tact_tokens, vis_tokens], dim=1)
        refined = self.self_attn(all_tokens)
        global_feat = refined.mean(dim=1)

        return {
            "object_logits": self.object_head(fused),
            "material_logits": self.material_head(fused),
            "grasp_quality": self.grasp_quality_head(fused).squeeze(-1),
            "tactile_features": tact_cls,
            "visual_features": vis_cls,
            "fused_features": global_feat,
        }


class MultiModalLoss(nn.Module):
    """Combined loss with contrastive alignment between modalities."""

    def __init__(self, contrastive_weight: float = 0.3, temperature: float = 0.07):
        super().__init__()
        self.ce = nn.CrossEntropyLoss()
        self.bce = nn.BCELoss()
        self.contrastive_weight = contrastive_weight
        self.temperature = temperature

    def info_nce_loss(self, feat_a: torch.Tensor, feat_b: torch.Tensor) -> torch.Tensor:
        """InfoNCE contrastive loss for cross-modal alignment."""
        feat_a = F.normalize(feat_a, dim=1)
        feat_b = F.normalize(feat_b, dim=1)
        sim = torch.matmul(feat_a, feat_b.T) / self.temperature
        labels = torch.arange(sim.shape[0], device=sim.device)
        loss_a2b = F.cross_entropy(sim, labels)
        loss_b2a = F.cross_entropy(sim.T, labels)
        return (loss_a2b + loss_b2a) / 2

    def forward(self, pred: dict, targets: dict) -> dict:
        obj_loss = self.ce(pred["object_logits"], targets["object_label"])
        mat_loss = self.ce(pred["material_logits"], targets["material_label"])
        grasp_loss = self.bce(pred["grasp_quality"], targets["grasp_quality"].float())
        contrastive = self.info_nce_loss(
            pred["tactile_features"], pred["visual_features"]
        )

        total = obj_loss + mat_loss + grasp_loss + self.contrastive_weight * contrastive

        return {
            "total": total,
            "object_loss": obj_loss,
            "material_loss": mat_loss,
            "grasp_loss": grasp_loss,
            "contrastive_loss": contrastive,
        }
