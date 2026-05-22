"""Variable Cartesian impedance control for deformable manipulation."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Mapping, Optional

import torch


class TaskPhase(str, Enum):
    """High-level manipulation phases used for gain scheduling."""

    APPROACH = "approach"
    CONTACT = "contact"
    MANIPULATION = "manipulation"
    RELEASE = "release"


class ComplianceMode(str, Enum):
    """Compliance presets for deformable-object handling."""

    PRECISE = "precise"
    COMPLIANT = "compliant"
    DEFORMABLE_SAFE = "deformable_safe"


@dataclass(frozen=True)
class ImpedanceGains:
    """Cartesian stiffness and damping gains."""

    stiffness: tuple[float, float, float, float, float, float]
    damping: tuple[float, float, float, float, float, float]

    def stiffness_tensor(self, device: Optional[torch.device] = None) -> torch.Tensor:
        return torch.tensor(self.stiffness, dtype=torch.float32, device=device)

    def damping_tensor(self, device: Optional[torch.device] = None) -> torch.Tensor:
        return torch.tensor(self.damping, dtype=torch.float32, device=device)


@dataclass
class CartesianState:
    """End-effector Cartesian state in 6-DoF pose/twist coordinates."""

    pose: torch.Tensor
    twist: torch.Tensor
    wrench: torch.Tensor = field(default_factory=lambda: torch.zeros(6, dtype=torch.float32))

    def __post_init__(self) -> None:
        self.pose = _ensure_vector(self.pose)
        self.twist = _ensure_vector(self.twist)
        self.wrench = _ensure_vector(self.wrench)


@dataclass(frozen=True)
class ControlCommand:
    """Controller output including safety and contact annotations."""

    wrench: torch.Tensor
    stiffness: torch.Tensor
    damping: torch.Tensor
    contact_detected: bool
    saturated: bool


class ImpedanceController:
    """6-DoF Cartesian impedance controller with variable compliance."""

    def __init__(
        self,
        phase_gains: Optional[Mapping[TaskPhase, ImpedanceGains]] = None,
        *,
        contact_force_threshold: float = 4.0,
        torque_threshold: float = 0.75,
        force_limit: float = 20.0,
        torque_limit: float = 2.0,
        force_feedback_gain: float = 0.1,
        device: Optional[torch.device] = None,
    ) -> None:
        self.device = device or torch.device("cpu")
        self.phase_gains = dict(phase_gains or _default_phase_gains())
        self.contact_force_threshold = contact_force_threshold
        self.torque_threshold = torque_threshold
        self.force_limit = force_limit
        self.torque_limit = torque_limit
        self.force_feedback_gain = force_feedback_gain
        self.mode_scaling = {
            ComplianceMode.PRECISE: 1.0,
            ComplianceMode.COMPLIANT: 0.6,
            ComplianceMode.DEFORMABLE_SAFE: 0.35,
        }

    def detect_contact(self, wrench: torch.Tensor) -> bool:
        wrench = _ensure_vector(wrench).to(device=self.device)
        return bool(
            torch.linalg.norm(wrench[:3]) >= self.contact_force_threshold
            or torch.linalg.norm(wrench[3:]) >= self.torque_threshold
        )

    def compute_command(
        self,
        current: CartesianState,
        target_pose: torch.Tensor,
        *,
        target_twist: Optional[torch.Tensor] = None,
        desired_wrench: Optional[torch.Tensor] = None,
        phase: TaskPhase = TaskPhase.MANIPULATION,
        mode: ComplianceMode = ComplianceMode.DEFORMABLE_SAFE,
    ) -> ControlCommand:
        target_pose = _ensure_vector(target_pose).to(device=self.device)
        target_twist = torch.zeros(6, device=self.device) if target_twist is None else _ensure_vector(target_twist).to(device=self.device)
        desired_wrench = torch.zeros(6, device=self.device) if desired_wrench is None else _ensure_vector(desired_wrench).to(device=self.device)

        pose_error = target_pose - current.pose.to(device=self.device)
        twist_error = target_twist - current.twist.to(device=self.device)
        measured_wrench = current.wrench.to(device=self.device)
        contact_detected = self.detect_contact(measured_wrench)

        stiffness, damping = self._scheduled_gains(phase, mode, contact_detected, measured_wrench)
        wrench = stiffness * pose_error + damping * twist_error + desired_wrench - self.force_feedback_gain * measured_wrench
        limited_wrench, saturated = self._limit_wrench(wrench)
        return ControlCommand(
            wrench=limited_wrench,
            stiffness=stiffness,
            damping=damping,
            contact_detected=contact_detected,
            saturated=saturated,
        )

    def _scheduled_gains(
        self,
        phase: TaskPhase,
        mode: ComplianceMode,
        contact_detected: bool,
        measured_wrench: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        gains = self.phase_gains[phase]
        scale = self.mode_scaling[mode]
        stiffness = gains.stiffness_tensor(self.device) * scale
        damping = gains.damping_tensor(self.device) * scale
        if contact_detected:
            force_ratio = torch.linalg.norm(measured_wrench[:3]).clamp_min(1.0) / max(self.contact_force_threshold, 1e-6)
            translational_scale = 1.0 / force_ratio
            stiffness[:3] = stiffness[:3] * translational_scale.clamp(min=0.25, max=1.0)
            damping[:3] = damping[:3] * translational_scale.clamp(min=0.5, max=1.0)
        return stiffness, damping

    def _limit_wrench(self, wrench: torch.Tensor) -> tuple[torch.Tensor, bool]:
        limited = wrench.clone()
        saturated = False
        force_norm = torch.linalg.norm(limited[:3])
        if force_norm > self.force_limit:
            limited[:3] *= self.force_limit / force_norm
            saturated = True
        torque_norm = torch.linalg.norm(limited[3:])
        if torque_norm > self.torque_limit:
            limited[3:] *= self.torque_limit / torque_norm
            saturated = True
        return limited, saturated


def _ensure_vector(vector: torch.Tensor) -> torch.Tensor:
    vector = vector.to(dtype=torch.float32)
    if vector.shape != (6,):
        raise ValueError("Expected a 6-D Cartesian vector.")
    return vector


def _default_phase_gains() -> Mapping[TaskPhase, ImpedanceGains]:
    return {
        TaskPhase.APPROACH: ImpedanceGains(
            stiffness=(900.0, 900.0, 900.0, 80.0, 80.0, 80.0),
            damping=(120.0, 120.0, 120.0, 12.0, 12.0, 12.0),
        ),
        TaskPhase.CONTACT: ImpedanceGains(
            stiffness=(450.0, 450.0, 350.0, 50.0, 50.0, 50.0),
            damping=(90.0, 90.0, 80.0, 10.0, 10.0, 10.0),
        ),
        TaskPhase.MANIPULATION: ImpedanceGains(
            stiffness=(300.0, 300.0, 250.0, 35.0, 35.0, 35.0),
            damping=(70.0, 70.0, 65.0, 8.0, 8.0, 8.0),
        ),
        TaskPhase.RELEASE: ImpedanceGains(
            stiffness=(200.0, 200.0, 160.0, 20.0, 20.0, 20.0),
            damping=(55.0, 55.0, 50.0, 5.0, 5.0, 5.0),
        ),
    }
