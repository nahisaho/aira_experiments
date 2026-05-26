#!/usr/bin/env python3
"""
Signal-level basecalling improvement module using Recurrent Neural Networks.
Part of the LongSV-Integra pipeline for structural variant detection.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
from enum import Enum
import math
import logging

logger = logging.getLogger(__name__)


@dataclass
class SignalConfig:
    """Configuration for signal processing and basecalling."""
    sample_rate: int = 4000
    window_size: int = 300
    stride: int = 30
    hidden_size: int = 256
    num_layers: int = 5
    dropout: float = 0.1
    num_classes: int = 5  # A, C, G, T, blank
    learning_rate: float = 0.001
    batch_size: int = 64
    normalization: str = "mad"  # median absolute deviation


class SignalNormalizer:
    """Normalize raw nanopore signals using MAD or Z-score normalization."""

    @staticmethod
    def mad_normalize(signal: np.ndarray) -> np.ndarray:
        median = np.median(signal)
        mad = np.median(np.abs(signal - median))
        if mad == 0:
            mad = 1.0
        return (signal - median) / (mad * 1.4826)

    @staticmethod
    def zscore_normalize(signal: np.ndarray) -> np.ndarray:
        mean = np.mean(signal)
        std = np.std(signal)
        if std == 0:
            std = 1.0
        return (signal - mean) / std

    @staticmethod
    def segment_signal(signal: np.ndarray, window_size: int, stride: int) -> np.ndarray:
        n_windows = max(1, (len(signal) - window_size) // stride + 1)
        segments = np.zeros((n_windows, window_size))
        for i in range(n_windows):
            start = i * stride
            end = min(start + window_size, len(signal))
            segments[i, :end - start] = signal[start:end]
        return segments


class GRUCell:
    """Manual GRU cell implementation for basecalling."""

    def __init__(self, input_size: int, hidden_size: int):
        self.input_size = input_size
        self.hidden_size = hidden_size
        scale = 1.0 / math.sqrt(hidden_size)
        # Initialize weight matrices
        self.W_z = np.random.uniform(-scale, scale, (hidden_size, input_size + hidden_size))
        self.W_r = np.random.uniform(-scale, scale, (hidden_size, input_size + hidden_size))
        self.W_h = np.random.uniform(-scale, scale, (hidden_size, input_size + hidden_size))
        self.b_z = np.zeros(hidden_size)
        self.b_r = np.zeros(hidden_size)
        self.b_h = np.zeros(hidden_size)

    def forward(self, x: np.ndarray, h_prev: np.ndarray) -> np.ndarray:
        combined = np.concatenate([x, h_prev])
        z = self._sigmoid(self.W_z @ combined + self.b_z)
        r = self._sigmoid(self.W_r @ combined + self.b_r)
        combined_r = np.concatenate([x, r * h_prev])
        h_tilde = np.tanh(self.W_h @ combined_r + self.b_h)
        h_new = (1 - z) * h_prev + z * h_tilde
        return h_new

    @staticmethod
    def _sigmoid(x: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


class BidirectionalGRU:
    """Bidirectional GRU for sequence modeling in basecalling."""

    def __init__(self, input_size: int, hidden_size: int, num_layers: int):
        self.layers_fwd = []
        self.layers_bwd = []
        for i in range(num_layers):
            in_sz = input_size if i == 0 else hidden_size * 2
            self.layers_fwd.append(GRUCell(in_sz, hidden_size))
            self.layers_bwd.append(GRUCell(in_sz, hidden_size))
        self.hidden_size = hidden_size
        self.num_layers = num_layers

    def forward(self, x_seq: np.ndarray) -> np.ndarray:
        """Process sequence through bidirectional GRU stack.
        Args:
            x_seq: (seq_len, input_size)
        Returns:
            output: (seq_len, hidden_size * 2)
        """
        seq_len = x_seq.shape[0]
        current_input = x_seq

        for layer_idx in range(self.num_layers):
            fwd_cell = self.layers_fwd[layer_idx]
            bwd_cell = self.layers_bwd[layer_idx]

            # Forward pass
            h_fwd = np.zeros(self.hidden_size)
            fwd_outputs = []
            for t in range(seq_len):
                h_fwd = fwd_cell.forward(current_input[t], h_fwd)
                fwd_outputs.append(h_fwd.copy())

            # Backward pass
            h_bwd = np.zeros(self.hidden_size)
            bwd_outputs = [None] * seq_len
            for t in range(seq_len - 1, -1, -1):
                h_bwd = bwd_cell.forward(current_input[t], h_bwd)
                bwd_outputs[t] = h_bwd.copy()

            current_input = np.array([
                np.concatenate([fwd_outputs[t], bwd_outputs[t]])
                for t in range(seq_len)
            ])

        return current_input


class CTCDecoder:
    """CTC beam search decoder for basecalling output."""

    BASES = ['A', 'C', 'G', 'T']
    BLANK = '-'

    @staticmethod
    def greedy_decode(log_probs: np.ndarray) -> str:
        labels = np.argmax(log_probs, axis=1)
        decoded = []
        prev = -1
        for label in labels:
            if label != 4 and label != prev:  # 4 = blank
                decoded.append(CTCDecoder.BASES[label])
            prev = label
        return ''.join(decoded)

    @staticmethod
    def beam_search_decode(log_probs: np.ndarray, beam_width: int = 10) -> str:
        seq_len, num_classes = log_probs.shape
        beams = [('', 0.0)]

        for t in range(seq_len):
            new_beams = {}
            for prefix, score in beams:
                for c in range(num_classes):
                    new_prefix = prefix
                    if c < 4:
                        base = CTCDecoder.BASES[c]
                        if not prefix or prefix[-1] != base:
                            new_prefix = prefix + base
                    new_score = score + log_probs[t, c]
                    if new_prefix not in new_beams or new_beams[new_prefix] < new_score:
                        new_beams[new_prefix] = new_score

            beams = sorted(new_beams.items(), key=lambda x: -x[1])[:beam_width]

        return beams[0][0] if beams else ''


class SignalBasecaller:
    """RNN-based basecaller for nanopore signal data."""

    def __init__(self, config: Optional[SignalConfig] = None):
        self.config = config or SignalConfig()
        self.normalizer = SignalNormalizer()
        self.rnn = BidirectionalGRU(
            input_size=self.config.window_size,
            hidden_size=self.config.hidden_size,
            num_layers=self.config.num_layers
        )
        self.decoder = CTCDecoder()
        # Output projection
        scale = 1.0 / math.sqrt(self.config.hidden_size * 2)
        self.W_out = np.random.uniform(
            -scale, scale,
            (self.config.num_classes, self.config.hidden_size * 2)
        )
        self.b_out = np.zeros(self.config.num_classes)

    def basecall(self, raw_signal: np.ndarray, use_beam_search: bool = True) -> str:
        normalized = self.normalizer.mad_normalize(raw_signal)
        segments = self.normalizer.segment_signal(
            normalized, self.config.window_size, self.config.stride
        )
        rnn_output = self.rnn.forward(segments)
        logits = rnn_output @ self.W_out.T + self.b_out
        log_probs = self._log_softmax(logits)

        if use_beam_search:
            return self.decoder.beam_search_decode(log_probs)
        return self.decoder.greedy_decode(log_probs)

    @staticmethod
    def _log_softmax(x: np.ndarray) -> np.ndarray:
        max_x = np.max(x, axis=1, keepdims=True)
        exp_x = np.exp(x - max_x)
        return x - max_x - np.log(np.sum(exp_x, axis=1, keepdims=True))

    def quality_score(self, log_probs: np.ndarray) -> np.ndarray:
        probs = np.exp(log_probs)
        max_probs = np.max(probs, axis=1)
        return -10 * np.log10(1 - max_probs + 1e-10)
