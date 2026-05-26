"""Communication support system for locked-in syndrome patients."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

try:
    from p300_classifier import AdaptiveP300Classifier, P300SpellerSimulation
except ImportError:  # pragma: no cover
    from .p300_classifier import AdaptiveP300Classifier, P300SpellerSimulation


@dataclass
class PredictiveTextEngine:
    """Simple n-gram predictive text engine for phrase completion."""

    corpus: Sequence[str] = (
        "HELLO HOW ARE YOU",
        "PLEASE HELP ME",
        "I NEED WATER",
        "THANK YOU",
        "GOOD MORNING",
        "YES NO MAYBE",
    )

    def __post_init__(self) -> None:
        self.unigrams = Counter()
        self.bigrams = defaultdict(Counter)
        for sentence in self.corpus:
            words = sentence.upper().split()
            for word in words:
                self.unigrams[word] += 1
            for previous, current in zip(["<s>"] + words[:-1], words):
                self.bigrams[previous][current] += 1

    def predict(self, context: str, prefix: str = "", top_k: int = 5) -> List[str]:
        previous = context.upper().split()[-1] if context.strip() else "<s>"
        candidates = self.bigrams.get(previous, self.unigrams)
        ranked = [word for word, _ in candidates.most_common() if word.startswith(prefix.upper())]
        if len(ranked) < top_k:
            extra = [word for word, _ in self.unigrams.most_common() if word.startswith(prefix.upper()) and word not in ranked]
            ranked.extend(extra)
        return ranked[:top_k]


@dataclass
class AdaptiveUI:
    """UI adaptation based on user performance."""

    flash_duration: float = 0.125
    inter_stimulus_interval: float = 0.075
    font_scale: float = 1.0

    def update(self, recent_accuracy: float) -> Dict[str, float]:
        if recent_accuracy < 0.7:
            self.flash_duration = min(0.2, self.flash_duration + 0.01)
            self.inter_stimulus_interval = min(0.12, self.inter_stimulus_interval + 0.005)
            self.font_scale = min(1.4, self.font_scale + 0.05)
        else:
            self.flash_duration = max(0.08, self.flash_duration - 0.005)
            self.inter_stimulus_interval = max(0.05, self.inter_stimulus_interval - 0.002)
            self.font_scale = max(1.0, self.font_scale - 0.02)
        return {
            "flash_duration": self.flash_duration,
            "inter_stimulus_interval": self.inter_stimulus_interval,
            "font_scale": self.font_scale,
        }


@dataclass
class PatientProfileManager:
    """Persistence for user-specific calibration data."""

    base_dir: Path = Path(__file__).resolve().parent.parent / "profiles"

    def __post_init__(self) -> None:
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def save_profile(self, patient_id: str, profile: Dict[str, object]) -> Path:
        path = self.base_dir / f"{patient_id}.json"
        with path.open("w", encoding="utf-8") as handle:
            json.dump(profile, handle, indent=2)
        return path

    def load_profile(self, patient_id: str) -> Dict[str, object]:
        path = self.base_dir / f"{patient_id}.json"
        if not path.exists():
            return {}
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)


@dataclass
class BCISpeller:
    """High-level P300 communication interface with suggestions and correction."""

    classifier: AdaptiveP300Classifier
    predictive_text: PredictiveTextEngine = field(default_factory=PredictiveTextEngine)
    ui: AdaptiveUI = field(default_factory=AdaptiveUI)
    profiles: PatientProfileManager = field(default_factory=PatientProfileManager)
    simulator: P300SpellerSimulation = field(default_factory=P300SpellerSimulation)

    def __post_init__(self) -> None:
        self.output_text = ""
        self.history: List[str] = []

    def append_character(self, char: str) -> str:
        self.history.append(self.output_text)
        self.output_text += char
        return self.output_text

    def undo(self) -> str:
        if self.history:
            self.output_text = self.history.pop()
        return self.output_text

    def correct_last(self, replacement: str) -> str:
        if self.output_text:
            self.output_text = self.output_text[:-1] + replacement
        return self.output_text

    def predict_words(self, prefix: str = "") -> List[str]:
        return self.predictive_text.predict(self.output_text, prefix=prefix)

    def spell_step(self, target_char: str, repetitions: int = 6) -> Dict[str, object]:
        trials, labels = self.simulator.generate_trials(target_char, repetitions=repetitions)
        probabilities = self.classifier.predict_proba(trials)[:, 1]
        flash_events = self.simulator.flash_sequence(target_char, repetitions=repetitions)
        row_scores = np.zeros(6)
        col_scores = np.zeros(6)
        for probability, (row, col, _) in zip(probabilities, flash_events):
            if row >= 0:
                row_scores[row] += probability
            if col >= 0:
                col_scores[col] += probability
        decoded_char = self.simulator.matrix[int(np.argmax(row_scores)), int(np.argmax(col_scores))]
        self.append_character(decoded_char)
        accuracy = float(np.mean((probabilities > 0.5) == labels))
        ui_state = self.ui.update(accuracy)
        suggestions = self.predict_words(prefix="")
        return {
            "decoded_char": str(decoded_char),
            "text": self.output_text,
            "accuracy": accuracy,
            "ui": ui_state,
            "suggestions": suggestions,
        }

    def save_patient_profile(self, patient_id: str, extra: Optional[Dict[str, object]] = None) -> Path:
        payload = {
            "text": self.output_text,
            "ui": self.ui.update(1.0),
            "suggestions": self.predict_words(),
        }
        if extra:
            payload.update(extra)
        return self.profiles.save_profile(patient_id, payload)

    def load_patient_profile(self, patient_id: str) -> Dict[str, object]:
        profile = self.profiles.load_profile(patient_id)
        self.output_text = str(profile.get("text", ""))
        return profile
