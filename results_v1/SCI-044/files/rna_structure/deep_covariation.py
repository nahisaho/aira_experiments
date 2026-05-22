"""Deep learning-inspired covariation analysis for RNA MSAs.

This module is fully self-contained and uses only NumPy/SciPy-compatible
numerical operations. It provides:

* MSA parsing and preprocessing
* Classical covariation scores (MI, APC-MI, mean-field DCA DI)
* A small NumPy residual network for contact scoring
* Lightweight self-attention features over MSAs
* Thermodynamic/covariation score integration hooks
* Synthetic MSA generation for testing
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from typing import Callable, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np


ALPHABET: Tuple[str, ...] = ("A", "C", "G", "U", "-")
ALPHABET_INDEX = {symbol: idx for idx, symbol in enumerate(ALPHABET)}
GAP_INDEX = ALPHABET_INDEX["-"]
CANONICAL_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
)
EPS = 1e-8


ArrayLike = Union[np.ndarray, Sequence[Sequence[str]]]


@dataclass
class ResidualBlockParams:
    """Inference-only parameters for a residual block."""

    conv1_w: np.ndarray
    conv1_b: np.ndarray
    bn1_gamma: np.ndarray
    bn1_beta: np.ndarray
    bn1_mean: np.ndarray
    bn1_var: np.ndarray
    conv2_w: np.ndarray
    conv2_b: np.ndarray
    bn2_gamma: np.ndarray
    bn2_beta: np.ndarray
    bn2_mean: np.ndarray
    bn2_var: np.ndarray


class MSAProcessor:
    """Parser and preprocessing helper for RNA multiple sequence alignments."""

    def __init__(self, gap_threshold: float = 0.5) -> None:
        self.gap_threshold = float(gap_threshold)
        self.sequence_ids: List[str] = []
        self.msa: Optional[np.ndarray] = None
        self.original_msa: Optional[np.ndarray] = None
        self.kept_columns: Optional[np.ndarray] = None
        self.weights: Optional[np.ndarray] = None

    @staticmethod
    def _coerce_text(data: Union[str, os.PathLike[str]]) -> str:
        text = os.fspath(data)
        if os.path.exists(text):
            with open(text, "r", encoding="utf-8") as handle:
                return handle.read()
        return text

    @staticmethod
    def _normalize_sequence(sequence: str) -> str:
        cleaned = []
        for char in sequence.upper().replace("T", "U"):
            cleaned.append(char if char in ALPHABET_INDEX else "-")
        return "".join(cleaned)

    def _parse_fasta(self, text: str) -> Tuple[List[str], List[str]]:
        identifiers: List[str] = []
        sequences: List[str] = []
        header: Optional[str] = None
        chunks: List[str] = []
        for raw_line in text.splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    identifiers.append(header)
                    sequences.append(self._normalize_sequence("".join(chunks)))
                header = line[1:].strip() or f"seq_{len(identifiers)}"
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            identifiers.append(header)
            sequences.append(self._normalize_sequence("".join(chunks)))
        if not sequences:
            raise ValueError("No FASTA records found in MSA input.")
        return identifiers, sequences

    def _parse_stockholm(self, text: str) -> Tuple[List[str], List[str]]:
        seq_map: dict[str, List[str]] = {}
        order: List[str] = []
        for raw_line in text.splitlines():
            line = raw_line.rstrip()
            if not line or line.startswith("# STOCKHOLM") or line == "//":
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, segment = parts[0], parts[1]
            if name not in seq_map:
                seq_map[name] = []
                order.append(name)
            seq_map[name].append(segment)
        sequences = [self._normalize_sequence("".join(seq_map[name])) for name in order]
        if not sequences:
            raise ValueError("No Stockholm alignment rows found in MSA input.")
        return order, sequences

    @staticmethod
    def _to_array(sequences: Sequence[str]) -> np.ndarray:
        lengths = {len(seq) for seq in sequences}
        if len(lengths) != 1:
            raise ValueError("All MSA sequences must have identical aligned lengths.")
        return np.array([list(seq) for seq in sequences], dtype="<U1")

    def _remove_gappy_columns(self, msa: np.ndarray) -> np.ndarray:
        gap_fraction = np.mean(msa == "-", axis=0)
        keep = gap_fraction <= self.gap_threshold
        if not np.any(keep):
            raise ValueError("All columns were removed by the gap filter.")
        self.kept_columns = keep
        return msa[:, keep]

    def load_msa(self, data: Union[str, os.PathLike[str]], format: str = "fasta") -> np.ndarray:
        """Load an MSA from FASTA or Stockholm text or file path.

        Parameters
        ----------
        data:
            Raw alignment text or a file path.
        format:
            Either ``'fasta'`` or ``'stockholm'``.

        Returns
        -------
        np.ndarray
            Character array with shape ``(n_sequences, n_columns)``.
        """

        text = self._coerce_text(data)
        fmt = format.lower()
        if fmt == "fasta":
            ids, sequences = self._parse_fasta(text)
        elif fmt in {"stockholm", "sto"}:
            ids, sequences = self._parse_stockholm(text)
        else:
            raise ValueError("format must be 'fasta' or 'stockholm'.")
        msa = self._to_array(sequences)
        self.sequence_ids = ids
        self.original_msa = msa.copy()
        self.msa = self._remove_gappy_columns(msa)
        self.weights = None
        return self.msa

    def compute_weights(self, threshold: float = 0.8) -> np.ndarray:
        """Compute sequence weights using Henikoff weighting with redundancy tempering.

        The base weight follows the Henikoff position-based scheme. The optional
        ``threshold`` then discounts highly similar sequences using a simple
        sequence identity redundancy factor. This keeps the API compatible with
        common Neff-based workflows while retaining Henikoff-style weighting.
        """

        if self.msa is None:
            raise ValueError("Load an MSA before computing weights.")

        msa = self.msa
        n_seq, n_col = msa.shape
        henikoff = np.zeros(n_seq, dtype=float)
        supported = np.zeros(n_seq, dtype=float)

        for col_idx in range(n_col):
            column = msa[:, col_idx]
            non_gap = column != "-"
            if not np.any(non_gap):
                continue
            residues, counts = np.unique(column[non_gap], return_counts=True)
            n_types = float(len(residues))
            inv_count = {residue: 1.0 / (n_types * count) for residue, count in zip(residues, counts)}
            idx = np.where(non_gap)[0]
            henikoff[idx] += np.array([inv_count[res] for res in column[idx]], dtype=float)
            supported[idx] += 1.0

        henikoff /= np.maximum(supported, 1.0)
        if henikoff.sum() <= 0.0:
            henikoff[:] = 1.0
        henikoff *= n_seq / henikoff.sum()

        valid = (msa[:, None, :] != "-") & (msa[None, :, :] != "-")
        matches = (msa[:, None, :] == msa[None, :, :]) & valid
        denom = valid.sum(axis=2)
        identity = np.divide(
            matches.sum(axis=2),
            denom,
            out=np.zeros((n_seq, n_seq), dtype=float),
            where=denom > 0,
        )
        redundancy = np.maximum((identity >= float(threshold)).sum(axis=1), 1)
        self.weights = henikoff / redundancy
        return self.weights

    def get_neff(self) -> float:
        """Return the effective number of sequences."""

        if self.weights is None:
            self.compute_weights()
        assert self.weights is not None
        return float(self.weights.sum())

    def one_hot_encode(self) -> np.ndarray:
        """One-hot encode the processed MSA with alphabet A/C/G/U/gap."""

        if self.msa is None:
            raise ValueError("Load an MSA before one-hot encoding.")
        indices = np.vectorize(ALPHABET_INDEX.get)(self.msa)
        eye = np.eye(len(ALPHABET), dtype=float)
        return eye[indices]


class CovariationAnalyzer:
    """Classical covariation metrics for aligned RNA sequence families."""

    def __init__(self, msa: Union[MSAProcessor, ArrayLike], weights: Optional[np.ndarray] = None) -> None:
        if isinstance(msa, MSAProcessor):
            if msa.msa is None:
                raise ValueError("The supplied MSAProcessor has no loaded alignment.")
            self.msa = msa.msa
            if weights is None:
                weights = msa.weights if msa.weights is not None else msa.compute_weights()
        else:
            arr = np.asarray(msa)
            if arr.ndim != 2:
                raise ValueError("msa must be a 2D array-like object.")
            self.msa = arr.astype("<U1")
        self.n_seq, self.length = self.msa.shape
        self.q = len(ALPHABET)
        self.weights = np.ones(self.n_seq, dtype=float) if weights is None else np.asarray(weights, dtype=float)
        if self.weights.shape != (self.n_seq,):
            raise ValueError("weights must have shape (n_sequences,).")
        self.weights = np.clip(self.weights, EPS, None)
        self.weights /= self.weights.sum()
        self._one_hot = self._encode(self.msa)
        self._fi_cache: dict[float, np.ndarray] = {}
        self._fij_cache: dict[float, np.ndarray] = {}

    @staticmethod
    def _encode(msa: np.ndarray) -> np.ndarray:
        indices = np.vectorize(ALPHABET_INDEX.get)(msa)
        return np.eye(len(ALPHABET), dtype=float)[indices]

    def _compute_frequencies(self, pseudocount: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
        pc = float(pseudocount)
        if pc in self._fi_cache and pc in self._fij_cache:
            return self._fi_cache[pc], self._fij_cache[pc]

        fi = np.einsum("n,nla->la", self.weights, self._one_hot, optimize=True)
        fij = np.einsum("n,nla,nmb->lmab", self.weights, self._one_hot, self._one_hot, optimize=True)

        if pc > 0.0:
            fi = (1.0 - pc) * fi + pc / self.q
            fij = (1.0 - pc) * fij + pc / (self.q * self.q)

        for idx in range(self.length):
            fij[idx, idx, :, :] = 0.0
            fij[idx, idx, np.arange(self.q), np.arange(self.q)] = fi[idx]

        self._fi_cache[pc] = fi
        self._fij_cache[pc] = fij
        return fi, fij

    def compute_mi(self) -> np.ndarray:
        """Compute mutual information for all alignment column pairs."""

        fi, fij = self._compute_frequencies(pseudocount=0.0)
        denom = fi[:, None, :, None] * fi[None, :, None, :]
        ratio = np.log((fij + EPS) / (denom + EPS))
        mi = np.sum(fij * ratio, axis=(2, 3))
        np.fill_diagonal(mi, 0.0)
        return mi

    def compute_apc_mi(self) -> np.ndarray:
        """Apply average product correction to mutual information."""

        mi = self.compute_mi()
        row_mean = mi.mean(axis=1)
        global_mean = float(mi.mean()) + EPS
        apc = np.outer(row_mean, row_mean) / global_mean
        corrected = mi - apc
        np.fill_diagonal(corrected, 0.0)
        return corrected

    def _covariance_matrix(self, pseudocount: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        fi, fij = self._compute_frequencies(pseudocount=pseudocount)
        q_eff = self.q - 1
        fi_red = fi[:, :q_eff]
        fij_red = fij[:, :, :q_eff, :q_eff]
        cov_blocks = fij_red - fi_red[:, None, :, None] * fi_red[None, :, None, :]
        covariance = cov_blocks.transpose(0, 2, 1, 3).reshape(self.length * q_eff, self.length * q_eff)
        covariance += np.eye(covariance.shape[0]) * 1e-4
        precision = np.linalg.pinv(covariance)
        return fi, covariance, precision

    @staticmethod
    def _direct_probability(couplings: np.ndarray, fi: np.ndarray, fj: np.ndarray) -> np.ndarray:
        weights = np.exp(couplings - np.max(couplings))
        weights = np.clip(weights, EPS, None)
        fi = fi / np.clip(fi.sum(), EPS, None)
        fj = fj / np.clip(fj.sum(), EPS, None)
        mu_i = np.ones_like(fi) / len(fi)
        mu_j = np.ones_like(fj) / len(fj)

        for _ in range(100):
            new_mu_i = fi / np.clip(weights @ mu_j, EPS, None)
            new_mu_i /= np.clip(new_mu_i.sum(), EPS, None)
            new_mu_j = fj / np.clip(weights.T @ new_mu_i, EPS, None)
            new_mu_j /= np.clip(new_mu_j.sum(), EPS, None)
            delta = max(np.max(np.abs(new_mu_i - mu_i)), np.max(np.abs(new_mu_j - mu_j)))
            mu_i, mu_j = new_mu_i, new_mu_j
            if delta < 1e-6:
                break

        pdir = weights * np.outer(mu_i, mu_j)
        pdir /= np.clip(pdir.sum(), EPS, None)
        for _ in range(20):
            pdir *= fi[:, None] / np.clip(pdir.sum(axis=1, keepdims=True), EPS, None)
            pdir *= fj[None, :] / np.clip(pdir.sum(axis=0, keepdims=True), EPS, None)
        pdir /= np.clip(pdir.sum(), EPS, None)
        return pdir

    def compute_di_mfdca(self, pseudocount: float = 0.5) -> np.ndarray:
        """Compute mean-field DCA direct information."""

        fi, _, precision = self._covariance_matrix(pseudocount=pseudocount)
        q_eff = self.q - 1
        di = np.zeros((self.length, self.length), dtype=float)

        for i in range(self.length):
            start_i = i * q_eff
            end_i = start_i + q_eff
            for j in range(i + 1, self.length):
                start_j = j * q_eff
                end_j = start_j + q_eff
                couplings = np.zeros((self.q, self.q), dtype=float)
                couplings[:q_eff, :q_eff] = -precision[start_i:end_i, start_j:end_j]
                pdir = self._direct_probability(couplings, fi[i], fi[j])
                dij = float(np.sum(pdir * np.log((pdir + EPS) / (np.outer(fi[i], fi[j]) + EPS))))
                di[i, j] = di[j, i] = dij
        return di

    def get_pair_frequencies(self, pseudocount: float = 0.0) -> np.ndarray:
        """Return pair frequencies with shape (L, L, 5, 5)."""

        _, fij = self._compute_frequencies(pseudocount=pseudocount)
        return fij

    def get_single_frequencies(self, pseudocount: float = 0.0) -> np.ndarray:
        """Return single-site frequencies with shape (L, 5)."""

        fi, _ = self._compute_frequencies(pseudocount=pseudocount)
        return fi

    @staticmethod
    def get_top_pairs(matrix: np.ndarray, n: Optional[int] = None, min_separation: int = 4) -> List[Tuple[int, int, float]]:
        """Return ranked upper-triangular residue pairs as ``(i, j, score)``."""

        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError("matrix must be a square 2D array.")
        length = matrix.shape[0]
        pairs = [
            (i, j, float(matrix[i, j]))
            for i in range(length)
            for j in range(i + min_separation, length)
            if matrix[i, j] > 0.0
        ]
        pairs.sort(key=lambda item: item[2], reverse=True)
        return pairs if n is None else pairs[:n]


class AttentionCovariation:
    """Lightweight NumPy self-attention features over MSA rows and columns."""

    def __init__(self, latent_dim: int = 16, seed: Optional[int] = 0) -> None:
        self.latent_dim = int(latent_dim)
        rng = np.random.default_rng(seed)
        scale = 1.0 / math.sqrt(len(ALPHABET))
        self.wq = rng.normal(scale=scale, size=(len(ALPHABET), self.latent_dim))
        self.wk = rng.normal(scale=scale, size=(len(ALPHABET), self.latent_dim))
        self.wv = rng.normal(scale=scale, size=(len(ALPHABET), self.latent_dim))

    @staticmethod
    def _softmax(x: np.ndarray, axis: int = -1) -> np.ndarray:
        shifted = x - np.max(x, axis=axis, keepdims=True)
        exp_x = np.exp(shifted)
        return exp_x / np.clip(exp_x.sum(axis=axis, keepdims=True), EPS, None)

    def _self_attention(self, tokens: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        q = tokens @ self.wq
        k = tokens @ self.wk
        v = tokens @ self.wv
        scores = q @ k.T / math.sqrt(self.latent_dim)
        attn = self._softmax(scores, axis=-1)
        context = attn @ v
        return attn, context

    def compute_attention_features(self, msa_onehot: np.ndarray) -> np.ndarray:
        """Compute pairwise attention-derived features from one-hot encoded MSA.

        Parameters
        ----------
        msa_onehot:
            Array of shape ``(n_sequences, n_columns, 5)``.

        Returns
        -------
        np.ndarray
            Feature tensor with shape ``(L, L, 4)``.
        """

        if msa_onehot.ndim != 3 or msa_onehot.shape[2] != len(ALPHABET):
            raise ValueError("msa_onehot must have shape (N, L, 5).")

        n_seq, length, _ = msa_onehot.shape
        row_tokens = msa_onehot.mean(axis=1)
        row_attn, _ = self._self_attention(row_tokens)
        row_weights = row_attn.mean(axis=0)
        row_weights /= np.clip(row_weights.sum(), EPS, None)

        weighted_single = np.einsum("n,nla->la", row_weights, msa_onehot, optimize=True)
        weighted_pair = np.einsum("n,nla,nmb->lmab", row_weights, msa_onehot, msa_onehot, optimize=True)
        pair_strength = np.sqrt(np.sum(weighted_pair**2, axis=(2, 3)))

        col_tokens = msa_onehot.mean(axis=0)
        col_attn, col_context = self._self_attention(col_tokens)
        col_scores = 0.5 * (col_attn + col_attn.T)

        normalized_context = col_context / np.clip(np.linalg.norm(col_context, axis=1, keepdims=True), EPS, None)
        context_similarity = normalized_context @ normalized_context.T

        entropy = -np.sum(weighted_single * np.log(weighted_single + EPS), axis=1)
        conservation = 1.0 - entropy / math.log(len(ALPHABET))
        conservation_pair = 0.5 * (conservation[:, None] + conservation[None, :])

        return np.stack((col_scores, pair_strength, context_similarity, conservation_pair), axis=-1)


class CovariationNet:
    """Inference-only 2D residual network implemented entirely in NumPy."""

    def __init__(
        self,
        input_channels: int = 30,
        hidden_channels: int = 64,
        n_blocks: int = 4,
        seed: Optional[int] = 0,
        pretrained: Optional[dict[str, np.ndarray]] = None,
    ) -> None:
        self.input_channels = int(input_channels)
        self.hidden_channels = int(hidden_channels)
        self.n_blocks = int(n_blocks)
        self.rng = np.random.default_rng(seed)
        if pretrained is None:
            self._init_random_weights()
        else:
            self._load_pretrained(pretrained)

    def _he_conv(self, kernel_size: int, in_channels: int, out_channels: int) -> np.ndarray:
        scale = math.sqrt(2.0 / (kernel_size * kernel_size * in_channels))
        return self.rng.normal(scale=scale, size=(kernel_size, kernel_size, in_channels, out_channels))

    def _bn_params(self, channels: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return (
            np.ones(channels, dtype=float),
            np.zeros(channels, dtype=float),
            np.zeros(channels, dtype=float),
            np.ones(channels, dtype=float),
        )

    def _init_random_weights(self) -> None:
        self.input_w = self._he_conv(1, self.input_channels, self.hidden_channels)
        self.input_b = np.zeros(self.hidden_channels, dtype=float)
        self.blocks: List[ResidualBlockParams] = []
        for _ in range(self.n_blocks):
            bn1 = self._bn_params(self.hidden_channels)
            bn2 = self._bn_params(self.hidden_channels)
            self.blocks.append(
                ResidualBlockParams(
                    conv1_w=self._he_conv(3, self.hidden_channels, self.hidden_channels),
                    conv1_b=np.zeros(self.hidden_channels, dtype=float),
                    bn1_gamma=bn1[0],
                    bn1_beta=bn1[1],
                    bn1_mean=bn1[2],
                    bn1_var=bn1[3],
                    conv2_w=self._he_conv(3, self.hidden_channels, self.hidden_channels),
                    conv2_b=np.zeros(self.hidden_channels, dtype=float),
                    bn2_gamma=bn2[0],
                    bn2_beta=bn2[1],
                    bn2_mean=bn2[2],
                    bn2_var=bn2[3],
                )
            )
        self.output_w = self._he_conv(1, self.hidden_channels, 1)
        self.output_b = np.zeros(1, dtype=float)

    def _load_pretrained(self, weights: dict[str, np.ndarray]) -> None:
        self.input_w = weights["input_w"]
        self.input_b = weights["input_b"]
        self.blocks = weights["blocks"]
        self.output_w = weights["output_w"]
        self.output_b = weights["output_b"]

    @staticmethod
    def conv2d(x: np.ndarray, weight: np.ndarray, bias: np.ndarray, padding: int = 0) -> np.ndarray:
        """2D convolution for HWC tensors in inference mode."""

        if x.ndim != 3:
            raise ValueError("x must have shape (H, W, C).")
        kh, kw, in_channels, out_channels = weight.shape
        if x.shape[2] != in_channels:
            raise ValueError("Input channel count does not match convolution weight.")
        padded = np.pad(x, ((padding, padding), (padding, padding), (0, 0)), mode="constant")
        height, width = x.shape[:2]
        out = np.zeros((height, width, out_channels), dtype=float)
        for di in range(kh):
            for dj in range(kw):
                patch = padded[di : di + height, dj : dj + width, :]
                out += patch @ weight[di, dj]
        out += bias.reshape(1, 1, -1)
        return out

    @staticmethod
    def batch_norm(
        x: np.ndarray,
        gamma: np.ndarray,
        beta: np.ndarray,
        mean: np.ndarray,
        var: np.ndarray,
        eps: float = 1e-5,
    ) -> np.ndarray:
        """Batch normalization in inference mode."""

        return gamma.reshape(1, 1, -1) * (x - mean.reshape(1, 1, -1)) / np.sqrt(var.reshape(1, 1, -1) + eps) + beta.reshape(1, 1, -1)

    @staticmethod
    def relu(x: np.ndarray) -> np.ndarray:
        """ReLU activation."""

        return np.maximum(x, 0.0)

    @staticmethod
    def sigmoid(x: np.ndarray) -> np.ndarray:
        """Sigmoid activation."""

        return 1.0 / (1.0 + np.exp(-np.clip(x, -50.0, 50.0)))

    def _residual_block(self, x: np.ndarray, block: ResidualBlockParams) -> np.ndarray:
        y = self.conv2d(x, block.conv1_w, block.conv1_b, padding=1)
        y = self.batch_norm(y, block.bn1_gamma, block.bn1_beta, block.bn1_mean, block.bn1_var)
        y = self.relu(y)
        y = self.conv2d(y, block.conv2_w, block.conv2_b, padding=1)
        y = self.batch_norm(y, block.bn2_gamma, block.bn2_beta, block.bn2_mean, block.bn2_var)
        return self.relu(x + y)

    def forward(self, features: np.ndarray) -> np.ndarray:
        """Run an inference pass and return an ``L x L`` contact probability map."""

        if features.ndim != 3:
            raise ValueError("features must have shape (L, L, C).")
        if features.shape[2] != self.input_channels:
            raise ValueError(f"Expected {self.input_channels} input channels, got {features.shape[2]}.")
        x = self.conv2d(features, self.input_w, self.input_b, padding=0)
        x = self.relu(x)
        for block in self.blocks:
            x = self._residual_block(x, block)
        x = self.conv2d(x, self.output_w, self.output_b, padding=0)[..., 0]
        probs = self.sigmoid(x)
        probs = 0.5 * (probs + probs.T)
        np.fill_diagonal(probs, 0.0)
        return probs

    @staticmethod
    def build_input_features(
        msa: Union[MSAProcessor, np.ndarray],
        attention_features: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Construct pairwise input features for the covariation network."""

        if isinstance(msa, MSAProcessor):
            if msa.msa is None:
                raise ValueError("MSAProcessor must contain a loaded alignment.")
            processor = msa
        else:
            processor = MSAProcessor()
            processor.msa = np.asarray(msa, dtype="<U1")
            processor.weights = processor.compute_weights() if processor.msa.size else np.array([], dtype=float)
        if processor.weights is None:
            processor.compute_weights()

        analyzer = CovariationAnalyzer(processor, processor.weights)
        mi = analyzer.compute_mi()
        apc = analyzer.compute_apc_mi()
        di = analyzer.compute_di_mfdca()
        fi = analyzer.get_single_frequencies()
        fij = analyzer.get_pair_frequencies()

        entropy = -np.sum(fi * np.log(fi + EPS), axis=1)
        conservation = 1.0 - entropy / math.log(len(ALPHABET))
        gap_freq = fi[:, GAP_INDEX]
        cons_pair = 0.5 * (conservation[:, None] + conservation[None, :])
        gap_pair = 0.5 * (gap_freq[:, None] + gap_freq[None, :])
        pair_freq_flat = fij.reshape(fij.shape[0], fij.shape[1], -1)

        features = np.concatenate(
            [
                mi[..., None],
                apc[..., None],
                di[..., None],
                cons_pair[..., None],
                gap_pair[..., None],
                pair_freq_flat,
            ],
            axis=-1,
        )
        if attention_features is not None:
            if attention_features.shape[:2] != features.shape[:2]:
                raise ValueError("attention_features must share the same L x L dimensions.")
            features = np.concatenate([features, attention_features], axis=-1)
        return features

    def predict_contacts(
        self,
        msa: Union[MSAProcessor, np.ndarray],
        threshold: float = 0.5,
    ) -> List[Tuple[int, int, float]]:
        """Predict base-pair contacts from an MSA."""

        features = self.build_input_features(msa)
        probabilities = self.forward(features)
        return CovariationAnalyzer.get_top_pairs(np.where(probabilities >= threshold, probabilities, 0.0), min_separation=4)


class CovariationIntegrator:
    """Combine thermodynamic and covariation-derived evidence."""

    def __init__(self, lambda_weight: float = 1.0, min_loop_length: int = 3) -> None:
        self.lambda_weight = float(lambda_weight)
        self.min_loop_length = int(min_loop_length)

    def fit_lambda(
        self,
        thermo_scores: Sequence[np.ndarray],
        covariation_scores: Sequence[np.ndarray],
        targets: Sequence[np.ndarray],
        lambda_grid: Optional[Iterable[float]] = None,
    ) -> float:
        """Learn a scalar λ by minimizing squared error on training contact maps."""

        grid = list(lambda_grid) if lambda_grid is not None else list(np.linspace(-2.0, 2.0, 41))
        best_lambda = self.lambda_weight
        best_loss = float("inf")
        for candidate in grid:
            loss = 0.0
            for thermo, cov, target in zip(thermo_scores, covariation_scores, targets):
                pred = self.integrate(thermo, cov, lambda_weight=float(candidate))
                loss += float(np.mean((pred - target) ** 2))
            if loss < best_loss:
                best_loss = loss
                best_lambda = float(candidate)
        self.lambda_weight = best_lambda
        return best_lambda

    def integrate(self, thermo_scores: np.ndarray, covariation_scores: np.ndarray, lambda_weight: float = 1.0) -> np.ndarray:
        """Combine thermodynamic and covariation scores by weighted summation."""

        thermo = np.asarray(thermo_scores, dtype=float)
        cov = np.asarray(covariation_scores, dtype=float)
        if thermo.shape != cov.shape:
            raise ValueError("thermo_scores and covariation_scores must have the same shape.")
        return thermo + float(lambda_weight) * cov

    @staticmethod
    def _is_pair(a: str, b: str) -> bool:
        return (a, b) in CANONICAL_PAIRS

    @staticmethod
    def _pair_energy_from_model(sequence: str, i: int, j: int, model: object) -> float:
        if hasattr(model, "pair_energy"):
            return float(model.pair_energy(i, j, sequence))
        if hasattr(model, "get_pair_energy"):
            return float(model.get_pair_energy(i, j, sequence))
        if callable(model):
            return float(model(i, j, sequence))
        pair = (sequence[i], sequence[j])
        heuristic = {
            ("G", "C"): -3.0,
            ("C", "G"): -3.0,
            ("A", "U"): -2.0,
            ("U", "A"): -2.0,
            ("G", "U"): -1.0,
            ("U", "G"): -1.0,
        }
        return heuristic.get(pair, 2.0)

    def constrained_fold(
        self,
        sequence: str,
        covariation_matrix: np.ndarray,
        zuker_model: object,
    ) -> Tuple[str, float]:
        """Fold an RNA using a simple soft-constrained dynamic program.

        The input ``covariation_matrix`` is interpreted as an additive energy-like
        soft constraint through ``integrate``. To reward supported base pairs,
        provide negative values for strongly supported pairs.
        """

        seq = sequence.upper().replace("T", "U")
        length = len(seq)
        if covariation_matrix.shape != (length, length):
            raise ValueError("covariation_matrix must have shape (len(sequence), len(sequence)).")

        dp = np.zeros((length, length), dtype=float)
        trace = np.full((length, length), -1, dtype=int)

        for span in range(1, length):
            for i in range(length - span):
                j = i + span
                best = dp[i + 1, j] if i + 1 <= j else 0.0
                trace[i, j] = -2
                left_skip = dp[i, j - 1] if i <= j - 1 else 0.0
                if left_skip < best:
                    best = left_skip
                    trace[i, j] = -3
                for k in range(i, j):
                    split = (dp[i, k] if i <= k else 0.0) + (dp[k + 1, j] if k + 1 <= j else 0.0)
                    if split < best:
                        best = split
                        trace[i, j] = k
                if j - i > self.min_loop_length and self._is_pair(seq[i], seq[j]):
                    thermo = self._pair_energy_from_model(seq, i, j, zuker_model)
                    pair_energy = self.integrate(
                        np.array([[thermo]], dtype=float),
                        np.array([[covariation_matrix[i, j]]], dtype=float),
                        lambda_weight=self.lambda_weight,
                    )[0, 0]
                    candidate = (dp[i + 1, j - 1] if i + 1 <= j - 1 else 0.0) + pair_energy
                    if candidate < best:
                        best = candidate
                        trace[i, j] = -4
                dp[i, j] = best

        structure = ["."] * length

        def traceback(i: int, j: int) -> None:
            if i >= j:
                return
            code = trace[i, j]
            if code == -2:
                traceback(i + 1, j)
            elif code == -3:
                traceback(i, j - 1)
            elif code == -4:
                structure[i] = "("
                structure[j] = ")"
                traceback(i + 1, j - 1)
            elif code >= i:
                traceback(i, code)
                traceback(code + 1, j)

        if length > 1:
            traceback(0, length - 1)
        return "".join(structure), float(dp[0, length - 1] if length else 0.0)


class MSAGenerator:
    """Synthetic RNA MSA generator with paired covarying mutations."""

    def __init__(self, diversity: float = 1.0, conservation: float = 0.5, seed: Optional[int] = 0) -> None:
        self.diversity = float(diversity)
        self.conservation = float(np.clip(conservation, 0.0, 1.0))
        self.rng = np.random.default_rng(seed)

    @staticmethod
    def _parse_pairs(structure: str) -> List[Tuple[int, int]]:
        stack: List[int] = []
        pairs: List[Tuple[int, int]] = []
        for idx, char in enumerate(structure):
            if char == "(":
                stack.append(idx)
            elif char == ")":
                if not stack:
                    raise ValueError("Unbalanced dot-bracket structure.")
                pairs.append((stack.pop(), idx))
        if stack:
            raise ValueError("Unbalanced dot-bracket structure.")
        return pairs

    def _mutate_unpaired(self, residue: str, mutation_rate: float) -> str:
        if self.rng.random() >= mutation_rate * self.diversity * (1.0 - self.conservation):
            return residue
        choices = [x for x in ALPHABET[:-1] if x != residue]
        return str(self.rng.choice(choices))

    def _mutate_pair(self, left: str, right: str, mutation_rate: float) -> Tuple[str, str]:
        if self.rng.random() >= mutation_rate * self.diversity:
            return left, right
        compatible = list(CANONICAL_PAIRS)
        if (left, right) in compatible and self.rng.random() < self.conservation:
            return left, right
        new_left, new_right = compatible[int(self.rng.integers(0, len(compatible)))]
        return new_left, new_right

    def generate(self, sequence: str, structure: str, n_seqs: int = 100, mutation_rate: float = 0.3) -> List[str]:
        """Generate a synthetic alignment preserving paired covariation patterns."""

        seq = list(sequence.upper().replace("T", "U"))
        if len(seq) != len(structure):
            raise ValueError("sequence and structure must have identical lengths.")
        pairs = self._parse_pairs(structure)
        paired_positions = {i for pair in pairs for i in pair}

        msa: List[str] = []
        for _ in range(int(n_seqs)):
            mutated = seq.copy()
            for i, j in pairs:
                mutated[i], mutated[j] = self._mutate_pair(mutated[i], mutated[j], mutation_rate)
            for idx in range(len(mutated)):
                if idx not in paired_positions:
                    mutated[idx] = self._mutate_unpaired(mutated[idx], mutation_rate)
            msa.append("".join(mutated))
        return msa


__all__ = [
    "MSAProcessor",
    "CovariationAnalyzer",
    "CovariationNet",
    "AttentionCovariation",
    "CovariationIntegrator",
    "MSAGenerator",
]
