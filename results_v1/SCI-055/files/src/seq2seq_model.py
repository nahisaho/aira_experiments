"""Template-free Seq2Seq retrosynthesis model with Transformer architecture."""

import math
import torch
import torch.nn as nn
from typing import Optional


class PositionalEncoding(nn.Module):
    def __init__(self, d_model: int, max_len: int = 512, dropout: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model))
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)
        self.register_buffer("pe", pe)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.pe[:, :x.size(1)]
        return self.dropout(x)


class RetroSynthTransformer(nn.Module):
    """Transformer-based seq2seq model for single-step retrosynthesis.

    Input: product SMILES (tokenized)
    Output: reactant SMILES (tokenized)
    """

    def __init__(
        self,
        vocab_size: int,
        d_model: int = 256,
        nhead: int = 8,
        num_encoder_layers: int = 6,
        num_decoder_layers: int = 6,
        dim_feedforward: int = 1024,
        dropout: float = 0.1,
        max_len: int = 256,
    ):
        super().__init__()
        self.d_model = d_model
        self.embedding = nn.Embedding(vocab_size, d_model, padding_idx=0)
        self.pos_encoder = PositionalEncoding(d_model, max_len, dropout)

        self.transformer = nn.Transformer(
            d_model=d_model,
            nhead=nhead,
            num_encoder_layers=num_encoder_layers,
            num_decoder_layers=num_decoder_layers,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True,
        )
        self.fc_out = nn.Linear(d_model, vocab_size)

    def _generate_square_subsequent_mask(self, sz: int, device: torch.device) -> torch.Tensor:
        mask = torch.triu(torch.ones(sz, sz, device=device), diagonal=1).bool()
        return mask

    def forward(
        self,
        src: torch.Tensor,
        tgt: torch.Tensor,
        src_key_padding_mask: Optional[torch.Tensor] = None,
        tgt_key_padding_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        tgt_mask = self._generate_square_subsequent_mask(tgt.size(1), tgt.device)

        src_emb = self.pos_encoder(self.embedding(src) * math.sqrt(self.d_model))
        tgt_emb = self.pos_encoder(self.embedding(tgt) * math.sqrt(self.d_model))

        output = self.transformer(
            src_emb, tgt_emb,
            tgt_mask=tgt_mask,
            src_key_padding_mask=src_key_padding_mask,
            tgt_key_padding_mask=tgt_key_padding_mask,
        )
        return self.fc_out(output)

    def greedy_decode(
        self, src: torch.Tensor, max_len: int, sos_idx: int, eos_idx: int
    ) -> torch.Tensor:
        self.eval()
        with torch.no_grad():
            src_emb = self.pos_encoder(self.embedding(src) * math.sqrt(self.d_model))
            memory = self.transformer.encoder(src_emb)

            batch_size = src.size(0)
            ys = torch.full((batch_size, 1), sos_idx, dtype=torch.long, device=src.device)

            for _ in range(max_len - 1):
                tgt_emb = self.pos_encoder(self.embedding(ys) * math.sqrt(self.d_model))
                tgt_mask = self._generate_square_subsequent_mask(ys.size(1), ys.device)
                out = self.transformer.decoder(tgt_emb, memory, tgt_mask=tgt_mask)
                logits = self.fc_out(out[:, -1:, :])
                next_token = logits.argmax(dim=-1)
                ys = torch.cat([ys, next_token], dim=1)
                if (next_token == eos_idx).all():
                    break
        return ys

    def beam_search_decode(
        self, src: torch.Tensor, max_len: int, sos_idx: int, eos_idx: int, beam_width: int = 5
    ):
        """Beam search decoding for diverse retrosynthesis predictions."""
        self.eval()
        with torch.no_grad():
            src_emb = self.pos_encoder(self.embedding(src) * math.sqrt(self.d_model))
            memory = self.transformer.encoder(src_emb)

            beams = [(torch.tensor([[sos_idx]], device=src.device), 0.0)]
            completed = []

            for step in range(max_len):
                candidates = []
                for seq, score in beams:
                    if seq[0, -1].item() == eos_idx:
                        completed.append((seq, score))
                        continue
                    tgt_emb = self.pos_encoder(self.embedding(seq) * math.sqrt(self.d_model))
                    tgt_mask = self._generate_square_subsequent_mask(seq.size(1), seq.device)
                    out = self.transformer.decoder(tgt_emb, memory, tgt_mask=tgt_mask)
                    logits = self.fc_out(out[:, -1, :])
                    log_probs = torch.log_softmax(logits, dim=-1)
                    topk_probs, topk_ids = log_probs.topk(beam_width, dim=-1)
                    for k in range(beam_width):
                        new_seq = torch.cat([seq, topk_ids[:, k:k+1]], dim=1)
                        new_score = score + topk_probs[0, k].item()
                        candidates.append((new_seq, new_score))
                candidates.sort(key=lambda x: x[1], reverse=True)
                beams = candidates[:beam_width]
                if len(completed) >= beam_width:
                    break

            completed.extend(beams)
            completed.sort(key=lambda x: x[1] / max(x[0].size(1), 1), reverse=True)
            return completed[:beam_width]


class Graph2SMILESEncoder(nn.Module):
    """Graph Neural Network encoder for molecular graphs.

    Converts molecular graph to a sequence of node embeddings
    for use with the Transformer decoder.
    """

    def __init__(self, node_feat_dim: int = 64, edge_feat_dim: int = 16,
                 d_model: int = 256, num_layers: int = 4, dropout: float = 0.1):
        super().__init__()
        self.node_embed = nn.Linear(node_feat_dim, d_model)
        self.edge_embed = nn.Linear(edge_feat_dim, d_model)

        self.gnn_layers = nn.ModuleList()
        for _ in range(num_layers):
            self.gnn_layers.append(
                nn.ModuleDict({
                    "msg_mlp": nn.Sequential(
                        nn.Linear(d_model * 3, d_model),
                        nn.ReLU(),
                        nn.Dropout(dropout),
                        nn.Linear(d_model, d_model),
                    ),
                    "update_mlp": nn.Sequential(
                        nn.Linear(d_model * 2, d_model),
                        nn.ReLU(),
                        nn.Dropout(dropout),
                        nn.Linear(d_model, d_model),
                    ),
                    "norm": nn.LayerNorm(d_model),
                })
            )

    def forward(self, node_feats: torch.Tensor, edge_feats: torch.Tensor,
                adj: torch.Tensor) -> torch.Tensor:
        """
        Args:
            node_feats: (batch, num_nodes, node_feat_dim)
            edge_feats: (batch, num_nodes, num_nodes, edge_feat_dim)
            adj: (batch, num_nodes, num_nodes) adjacency matrix
        Returns:
            node_embeddings: (batch, num_nodes, d_model)
        """
        h = self.node_embed(node_feats)
        e = self.edge_embed(edge_feats)

        for layer in self.gnn_layers:
            batch_size, n_nodes, d = h.shape
            h_i = h.unsqueeze(2).expand(-1, -1, n_nodes, -1)
            h_j = h.unsqueeze(1).expand(-1, n_nodes, -1, -1)
            msg_input = torch.cat([h_i, h_j, e], dim=-1)
            messages = layer["msg_mlp"](msg_input)
            messages = messages * adj.unsqueeze(-1)
            agg = messages.sum(dim=2)
            h_new = layer["update_mlp"](torch.cat([h, agg], dim=-1))
            h = layer["norm"](h + h_new)
        return h
