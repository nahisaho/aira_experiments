"""
Module 6: Exploratory Grasping - Safe Exploration Strategy for Unknown Objects
Bayesian exploration with safety-aware grasping for novel objects,
combining Gaussian process priors with learned grasp quality models.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import math
from typing import Optional


class GraspParameterEncoder(nn.Module):
    """Encode grasp parameters into a compact representation."""

    def __init__(self, grasp_dim: int = 7, embed_dim: int = 64):
        super().__init__()
        # Grasp params: position(3) + orientation(4, quaternion)
        self.encoder = nn.Sequential(
            nn.Linear(grasp_dim, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, embed_dim),
            nn.ReLU(inplace=True),
        )

    def forward(self, grasp_params: torch.Tensor) -> torch.Tensor:
        return self.encoder(grasp_params)


class VisualGraspProposer(nn.Module):
    """Generate candidate grasps from visual observation of the object."""

    def __init__(self, visual_dim: int = 512, num_candidates: int = 50):
        super().__init__()
        self.num_candidates = num_candidates
        self.backbone = nn.Sequential(
            nn.Conv2d(3, 32, 5, 2, 2), nn.ReLU(inplace=True),
            nn.Conv2d(32, 64, 3, 2, 1), nn.ReLU(inplace=True),
            nn.Conv2d(64, 128, 3, 2, 1), nn.ReLU(inplace=True),
            nn.Conv2d(128, 256, 3, 2, 1), nn.ReLU(inplace=True),
            nn.AdaptiveAvgPool2d(1), nn.Flatten(),
        )
        self.grasp_generator = nn.Sequential(
            nn.Linear(256, 512),
            nn.ReLU(inplace=True),
            nn.Linear(512, num_candidates * 7),  # 7 per grasp (pos + quat)
        )
        self.prior_quality = nn.Sequential(
            nn.Linear(256 + 7, 128),
            nn.ReLU(inplace=True),
            nn.Linear(128, 1),
            nn.Sigmoid(),
        )

    def forward(self, visual_img: torch.Tensor) -> dict:
        features = self.backbone(visual_img)
        B = features.shape[0]
        grasp_params = self.grasp_generator(features)
        grasp_params = grasp_params.view(B, self.num_candidates, 7)

        # Normalize quaternion part
        grasp_params_normed = grasp_params.clone()
        grasp_params_normed[:, :, 3:] = F.normalize(
            grasp_params[:, :, 3:], dim=2
        )

        # Prior quality estimate for each candidate
        features_expanded = features.unsqueeze(1).expand(-1, self.num_candidates, -1)
        quality_input = torch.cat([features_expanded, grasp_params_normed], dim=2)
        prior_quality = self.prior_quality(
            quality_input.view(-1, 256 + 7)
        ).view(B, self.num_candidates)

        return {
            "grasp_candidates": grasp_params_normed,
            "prior_quality": prior_quality,
            "visual_features": features,
        }


class GaussianProcessSurrogate(nn.Module):
    """
    Neural network approximation of Gaussian Process for grasp quality.
    Provides mean prediction and uncertainty estimation via MC dropout.
    """

    def __init__(self, input_dim: int = 64 + 256, hidden_dim: int = 256,
                 num_mc_samples: int = 20):
        super().__init__()
        self.num_mc_samples = num_mc_samples

        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, 128),
            nn.ReLU(inplace=True),
            nn.Dropout(0.2),
        )

        self.mean_head = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 1),
            nn.Sigmoid(),
        )

        self.log_var_head = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 1),
        )

    def forward(self, x: torch.Tensor, use_mc_dropout: bool = True) -> dict:
        """
        Args:
            x: (B, D) grasp + visual features
            use_mc_dropout: whether to use MC dropout for uncertainty
        Returns:
            dict with mean, variance, samples
        """
        if use_mc_dropout:
            self.train()  # Enable dropout
            samples = []
            for _ in range(self.num_mc_samples):
                h = self.network(x)
                mean = self.mean_head(h).squeeze(-1)
                samples.append(mean)
            samples = torch.stack(samples, dim=0)  # (S, B)
            pred_mean = samples.mean(dim=0)
            pred_var = samples.var(dim=0)
            self.eval()
        else:
            h = self.network(x)
            pred_mean = self.mean_head(h).squeeze(-1)
            log_var = self.log_var_head(h).squeeze(-1)
            pred_var = torch.exp(log_var)
            samples = None

        return {
            "mean": pred_mean,
            "variance": pred_var,
            "samples": samples,
        }


class AcquisitionFunction(nn.Module):
    """Acquisition functions for Bayesian optimization of grasp selection."""

    def __init__(self, method: str = "expected_improvement", xi: float = 0.01):
        super().__init__()
        self.method = method
        self.xi = xi

    def expected_improvement(self, mean: torch.Tensor, std: torch.Tensor,
                              best_so_far: float) -> torch.Tensor:
        improvement = mean - best_so_far - self.xi
        Z = improvement / std.clamp(min=1e-8)
        # Approximate CDF and PDF of standard normal
        cdf = 0.5 * (1 + torch.erf(Z / math.sqrt(2)))
        pdf = torch.exp(-0.5 * Z ** 2) / math.sqrt(2 * math.pi)
        ei = improvement * cdf + std * pdf
        return ei

    def upper_confidence_bound(self, mean: torch.Tensor, std: torch.Tensor,
                                beta: float = 2.0) -> torch.Tensor:
        return mean + beta * std

    def thompson_sampling(self, mean: torch.Tensor,
                           std: torch.Tensor) -> torch.Tensor:
        return mean + std * torch.randn_like(std)

    def forward(self, mean: torch.Tensor, variance: torch.Tensor,
                best_so_far: float = 0.0) -> torch.Tensor:
        std = variance.sqrt().clamp(min=1e-8)

        if self.method == "expected_improvement":
            return self.expected_improvement(mean, std, best_so_far)
        elif self.method == "ucb":
            return self.upper_confidence_bound(mean, std)
        elif self.method == "thompson":
            return self.thompson_sampling(mean, std)
        else:
            raise ValueError(f"Unknown acquisition function: {self.method}")


class SafetyConstraint(nn.Module):
    """Safety checker for grasp execution feasibility and force limits."""

    def __init__(self, max_force: float = 10.0,
                 max_velocity: float = 0.05,
                 emergency_threshold: float = 15.0):
        super().__init__()
        self.max_force = max_force
        self.max_velocity = max_velocity
        self.emergency_threshold = emergency_threshold

        # Learned safety classifier
        self.safety_net = nn.Sequential(
            nn.Linear(7 + 256, 128),  # grasp_params + visual_features
            nn.ReLU(inplace=True),
            nn.Linear(128, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 2),  # safe / unsafe
        )

    def check_force_limit(self, force: torch.Tensor) -> torch.Tensor:
        """Check if force is within safety limits."""
        force_mag = force.norm(dim=-1)
        is_safe = force_mag < self.max_force
        is_emergency = force_mag > self.emergency_threshold
        return is_safe, is_emergency

    def forward(self, grasp_params: torch.Tensor,
                visual_features: torch.Tensor) -> dict:
        x = torch.cat([grasp_params, visual_features], dim=-1)
        safety_logits = self.safety_net(x)
        safety_prob = F.softmax(safety_logits, dim=-1)[:, 0]  # P(safe)
        return {
            "safety_logits": safety_logits,
            "safety_probability": safety_prob,
            "is_safe": safety_prob > 0.8,
        }


class ExploratoryGraspingPolicy(nn.Module):
    """
    Complete exploratory grasping pipeline:
    1. Generate grasp candidates from visual input
    2. Evaluate candidates with GP surrogate
    3. Apply acquisition function for exploration/exploitation balance
    4. Safety filtering
    5. Execute best safe grasp
    6. Update GP with tactile feedback
    """

    def __init__(
        self,
        num_candidates: int = 50,
        num_refinement_steps: int = 5,
        exploration_budget: int = 20,
        acquisition: str = "expected_improvement",
    ):
        super().__init__()
        self.num_candidates = num_candidates
        self.num_refinement_steps = num_refinement_steps
        self.exploration_budget = exploration_budget

        # Components
        self.grasp_proposer = VisualGraspProposer(num_candidates=num_candidates)
        self.grasp_encoder = GraspParameterEncoder(grasp_dim=7, embed_dim=64)
        self.gp_surrogate = GaussianProcessSurrogate(input_dim=64 + 256)
        self.acquisition_fn = AcquisitionFunction(method=acquisition)
        self.safety_checker = SafetyConstraint()

        # Grasp refinement network
        self.refiner = nn.Sequential(
            nn.Linear(7 + 256 + 1, 128),  # grasp + visual + quality
            nn.ReLU(inplace=True),
            nn.Linear(128, 64),
            nn.ReLU(inplace=True),
            nn.Linear(64, 7),  # refined grasp delta
            nn.Tanh(),
        )
        self.refine_scale = nn.Parameter(torch.tensor(0.05))

        # Observation history encoder (for online learning)
        self.history_encoder = nn.GRU(
            input_size=64 + 1,  # grasp_embed + outcome
            hidden_size=128,
            num_layers=1,
            batch_first=True,
        )

    def propose_and_evaluate(self, visual_img: torch.Tensor,
                              best_so_far: float = 0.0) -> dict:
        """Generate and evaluate grasp candidates."""
        B = visual_img.shape[0]

        # Step 1: Propose candidates
        proposals = self.grasp_proposer(visual_img)
        candidates = proposals["grasp_candidates"]  # (B, K, 7)
        visual_feat = proposals["visual_features"]  # (B, 256)

        K = candidates.shape[1]

        # Step 2: Encode grasps
        grasp_embeds = self.grasp_encoder(candidates.view(-1, 7))
        grasp_embeds = grasp_embeds.view(B, K, -1)

        # Step 3: GP evaluation
        visual_expanded = visual_feat.unsqueeze(1).expand(-1, K, -1)
        gp_input = torch.cat([grasp_embeds, visual_expanded], dim=2)
        gp_input_flat = gp_input.view(B * K, -1)
        gp_result = self.gp_surrogate(gp_input_flat, use_mc_dropout=True)
        gp_mean = gp_result["mean"].view(B, K)
        gp_var = gp_result["variance"].view(B, K)

        # Step 4: Acquisition function
        acq_values = self.acquisition_fn(gp_mean, gp_var, best_so_far)

        # Step 5: Safety filtering
        candidates_flat = candidates.view(B * K, 7)
        vis_flat = visual_expanded.reshape(B * K, -1)
        safety = self.safety_checker(candidates_flat, vis_flat)
        safety_mask = safety["safety_probability"].view(B, K)

        # Combined score = acquisition * safety
        combined_score = acq_values * safety_mask

        # Step 6: Select best
        best_idx = combined_score.argmax(dim=1)
        best_grasp = candidates[torch.arange(B), best_idx]

        return {
            "best_grasp": best_grasp,
            "best_idx": best_idx,
            "all_candidates": candidates,
            "gp_mean": gp_mean,
            "gp_variance": gp_var,
            "acquisition_values": acq_values,
            "safety_scores": safety_mask,
            "combined_scores": combined_score,
            "prior_quality": proposals["prior_quality"],
        }

    def refine_grasp(self, grasp: torch.Tensor,
                     visual_features: torch.Tensor,
                     quality: torch.Tensor) -> torch.Tensor:
        """Iteratively refine grasp parameters."""
        for _ in range(self.num_refinement_steps):
            refine_input = torch.cat([
                grasp, visual_features, quality.unsqueeze(1)
            ], dim=1)
            delta = self.refiner(refine_input) * self.refine_scale
            grasp = grasp + delta
            # Re-normalize quaternion
            grasp_pos = grasp[:, :3]
            grasp_quat = F.normalize(grasp[:, 3:], dim=1)
            grasp = torch.cat([grasp_pos, grasp_quat], dim=1)
        return grasp

    def forward(self, visual_img: torch.Tensor,
                best_so_far: float = 0.0) -> dict:
        """
        Full exploratory grasping pipeline.
        Args:
            visual_img: (B, 3, H, W) visual observation
            best_so_far: best observed quality so far
        Returns:
            dict with selected grasp and exploration diagnostics
        """
        result = self.propose_and_evaluate(visual_img, best_so_far)

        # Refine the selected grasp
        proposals = self.grasp_proposer(visual_img)
        refined_grasp = self.refine_grasp(
            result["best_grasp"],
            proposals["visual_features"],
            result["gp_mean"][
                torch.arange(visual_img.shape[0]),
                result["best_idx"]
            ],
        )
        result["refined_grasp"] = refined_grasp
        return result


class ExplorationLoss(nn.Module):
    """Loss for training the exploratory grasping policy."""

    def __init__(self):
        super().__init__()

    def forward(self, pred: dict, target: dict) -> dict:
        # GP prediction loss
        gp_loss = F.mse_loss(
            pred["gp_mean"],
            target["true_quality"].unsqueeze(1).expand_as(pred["gp_mean"])
        )

        # Safety classification loss
        safety_loss = F.binary_cross_entropy(
            pred["safety_scores"],
            target["safety_labels"].float()
        )

        # Grasp quality maximization
        selected_quality = pred["gp_mean"][
            torch.arange(pred["best_idx"].shape[0]),
            pred["best_idx"]
        ]
        quality_loss = -selected_quality.mean()  # Maximize quality

        total = gp_loss + safety_loss + 0.5 * quality_loss

        return {
            "total": total,
            "gp_loss": gp_loss,
            "safety_loss": safety_loss,
            "quality_loss": quality_loss,
        }
