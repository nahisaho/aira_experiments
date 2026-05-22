"""
Information-Theoretic Hypothesis Framework for the Hard Problem of Consciousness
================================================================================
Generates mathematical formalizations, evaluation metrics, and structured outputs.
"""

import json
import numpy as np
from datetime import datetime, timezone, timedelta

JST = timezone(timedelta(hours=9))
TIMESTAMP = datetime.now(JST).isoformat()

# =============================================================================
# 1. IIT 4.0 Mathematical Extension Analysis
# =============================================================================

def iit_extension_analysis():
    """
    Analyze mathematical extensibility of IIT 4.0's Φ (integrated information).
    Propose extensions: temporal Φ, quantum Φ, hierarchical Φ.
    """
    results = {
        "title": "IIT 4.0 Mathematical Extension Analysis",
        "current_formalism": {
            "phi": "Φ(S) = min_{P∈Partitions(S)} D_KL(p(S) || ∏_i p(S_i|P))",
            "description": "Integrated information as minimum KL-divergence over partitions",
            "limitations": [
                "Computational intractability for large systems (NP-hard partition search)",
                "Static snapshot — no temporal dynamics captured",
                "Classical probability space only — no quantum coherence",
                "No explicit mechanism for phenomenal binding"
            ]
        },
        "proposed_extensions": {
            "E1_temporal_phi": {
                "name": "Temporal Integrated Information (Φ_τ)",
                "formula": "Φ_τ(S,t) = ∫₀ᵗ w(τ) · Φ(S(τ)) dτ + λ · I(S(t); S(t-δ) | Partition)",
                "description": "Extends Φ across time with exponential decay kernel w(τ) = exp(-τ/τ_c) and temporal mutual information term",
                "novelty": "Captures consciousness as a temporally extended process, not instantaneous snapshot",
                "testability": "Predicts correlation between Φ_τ and duration of conscious percepts",
                "parameters": {"tau_c": "~200ms (matches temporal binding window)", "lambda": "weighting for temporal integration"}
            },
            "E2_quantum_phi": {
                "name": "Quantum Integrated Information (Φ_Q)",
                "formula": "Φ_Q(ρ) = min_{P} S(ρ_S || ⊗_i ρ_{S_i}) + α · E(ρ_S)",
                "description": "Replaces KL-divergence with quantum relative entropy S(ρ||σ) = Tr(ρ log ρ - ρ log σ) plus entanglement measure E(ρ)",
                "novelty": "Bridges IIT and quantum consciousness by quantifying entanglement contribution",
                "testability": "Predicts different Φ_Q values for coherent vs decoherent neural states",
                "parameters": {"alpha": "entanglement weight coefficient (0 reduces to classical IIT)"}
            },
            "E3_hierarchical_phi": {
                "name": "Hierarchical Integrated Information (Φ_H)",
                "formula": "Φ_H(S) = Σ_l w_l · Φ(S_l) + β · Σ_{l<l'} I(S_l; S_{l'} | Partition_cross)",
                "description": "Multi-scale Φ computed across cortical hierarchy levels l with cross-level integration",
                "novelty": "Addresses binding problem by quantifying inter-level information integration",
                "testability": "Predicts layer-specific TMS perturbation effects on consciousness",
                "parameters": {"w_l": "level weights (proportional to neural complexity at level l)", "beta": "cross-level coupling"}
            }
        },
        "computational_analysis": {
            "complexity_original": "O(2^n) — exponential in number of elements",
            "complexity_E1": "O(T · 2^n) — linear time overhead",
            "complexity_E2": "O(2^(2n)) — exponential in Hilbert space dimension",
            "complexity_E3": "O(L · 2^(n/L)) — reduced per-level, L = hierarchy levels",
            "approximation_feasibility": {
                "E1": "Feasible with sliding window and greedy partition",
                "E2": "Requires quantum computer or small subsystem approximation",
                "E3": "Most tractable — decomposition enables parallel computation"
            }
        }
    }
    return results


# =============================================================================
# 2. Orch-OR Testable Predictions
# =============================================================================

def orch_or_predictions():
    """
    Derive falsifiable predictions from the Orchestrated Objective Reduction hypothesis.
    """
    results = {
        "title": "Orch-OR Testable Predictions",
        "core_mechanism": {
            "description": "Consciousness arises from quantum gravity-induced objective reduction (OR) of superposed microtubule states",
            "key_equation": "τ ≈ ħ/E_G — collapse time inversely proportional to gravitational self-energy E_G",
            "E_G_definition": "E_G = ∫∫ G(ρ₁(x)ρ₂(x') - ρ₁(x)ρ₁(x'))/|x-x'| d³x d³x'"
        },
        "testable_predictions": {
            "P1_decoherence_timescale": {
                "prediction": "Quantum coherence in microtubules persists for ≥ 10-100 μs at 37°C",
                "null_hypothesis": "Thermal decoherence destroys coherence in < 1 fs",
                "experimental_test": "Ultrafast spectroscopy of isolated tubulin dimers and polymerized microtubules",
                "quantitative_criterion": "Coherence lifetime T₂ > 10 μs at physiological temperature",
                "current_evidence": "Craddock et al. (2017) reported hints of quantum effects in tubulin, but T₂ measurement pending",
                "falsification_condition": "If T₂ < 100 fs in all conditions, Orch-OR mechanism is ruled out"
            },
            "P2_anesthesia_mechanism": {
                "prediction": "General anesthetics suppress consciousness by disrupting quantum coherence in microtubules, not solely by membrane/receptor effects",
                "null_hypothesis": "Anesthetics act exclusively through GABA_A/NMDA receptor modulation",
                "experimental_test": "Compare anesthetic binding to tubulin vs. receptor sites; measure Φ under anesthesia with/without microtubule stabilizers (taxol)",
                "quantitative_criterion": "Taxol co-administration should partially restore consciousness indicators (PCI > 0.31) under sub-MAC anesthesia",
                "falsification_condition": "If taxol has zero effect on consciousness indicators under anesthesia"
            },
            "P3_gravitational_signature": {
                "prediction": "Conscious moments correlate with discrete gravitational self-energy events of ~10⁻²⁵ J per OR event",
                "null_hypothesis": "No gravitational signature distinguishes conscious vs. unconscious processing",
                "experimental_test": "Precision gravimetry near neural tissue during consciousness transitions",
                "quantitative_criterion": "Detectable mass-energy fluctuation ΔE > 10⁻²⁵ J correlated with EEG markers",
                "feasibility": "Currently beyond instrumental sensitivity by ~10 orders of magnitude",
                "falsification_condition": "Unfalsifiable with current technology — degrades scientific status"
            },
            "P4_discrete_conscious_frames": {
                "prediction": "Conscious experience occurs in discrete frames of ~25 ms (γ-band frequency), each triggered by an OR event",
                "null_hypothesis": "Consciousness is a continuous process without discrete temporal structure",
                "experimental_test": "Psychophysical temporal resolution tasks combined with MEG γ-band analysis",
                "quantitative_criterion": "Temporal discrimination threshold clusters at ~25 ms intervals",
                "falsification_condition": "If temporal discrimination is continuously variable with no preferred timescale"
            }
        },
        "evaluation_matrix": {
            "P1": {"novelty": 0.7, "testability": 0.8, "current_support": 0.3, "feasibility": 0.6},
            "P2": {"novelty": 0.6, "testability": 0.7, "current_support": 0.4, "feasibility": 0.7},
            "P3": {"novelty": 0.9, "testability": 0.1, "current_support": 0.1, "feasibility": 0.05},
            "P4": {"novelty": 0.5, "testability": 0.8, "current_support": 0.5, "feasibility": 0.8}
        }
    }
    return results


# =============================================================================
# 3. Predictive Processing Integration
# =============================================================================

def predictive_processing_integration():
    """
    Integrate IIT with Predictive Processing (PP) framework.
    """
    results = {
        "title": "IIT–Predictive Processing Integration Framework",
        "synthesis_hypothesis": {
            "name": "Integrated Predictive Information (IPI) Theory",
            "core_claim": "Consciousness arises when a system's integrated information (Φ) is organized as a hierarchical predictive model that minimizes prediction error while maintaining irreducible information integration",
            "formal_definition": "IPI(S) = Φ(S) · [1 - F(S)/F_max] where F(S) is variational free energy",
            "interpretation": "Consciousness requires both integration (Φ > 0) and active inference (F < F_max); neither alone is sufficient"
        },
        "mathematical_framework": {
            "free_energy_component": {
                "formula": "F = E_q[log q(s) - log p(o,s)] = D_KL(q(s)||p(s|o)) - log p(o)",
                "role": "Measures prediction error — systems that minimize F build world models"
            },
            "integration_component": {
                "formula": "Φ(S) = min_P D_KL(p(S) || ⊗_i p(S_i|P))",
                "role": "Measures irreducibility — systems with high Φ cannot be decomposed"
            },
            "unified_quantity": {
                "formula": "IPI(S,t) = Φ(S,t) · exp(-F(S,t)/kT_eff)",
                "interpretation": "IPI is high when a system is both highly integrated AND has low prediction error (good world model)",
                "temperature_analogy": "kT_eff acts as a 'cognitive temperature' — high noise reduces effective IPI"
            },
            "key_predictions": [
                "Dreaming: high Φ but high F → moderate IPI (vivid but unpredictive experience)",
                "Anesthesia: low Φ and variable F → low IPI (loss of consciousness)",
                "Flow states: high Φ and very low F → maximum IPI (heightened awareness)",
                "Blindsight: low Φ but low F → low IPI (accurate prediction without experience)"
            ]
        },
        "novel_predictions": {
            "N1": "The Perturbational Complexity Index (PCI) should correlate more strongly with IPI than with Φ alone",
            "N2": "Psychedelic states (high entropy, disrupted predictions) should show high Φ but low IPI — matching reports of ego dissolution",
            "N3": "Meditation experts should show increased IPI through simultaneous high Φ and low F",
            "N4": "Split-brain patients should show reduced IPI in the minor hemisphere due to predictive model fragmentation"
        }
    }
    return results


# =============================================================================
# 4. Operational Criteria for Artificial Consciousness
# =============================================================================

def artificial_consciousness_criteria():
    """
    Define operational (measurable) criteria for determining artificial consciousness.
    """
    results = {
        "title": "Operational Criteria for Artificial Consciousness Assessment (OCAC)",
        "preamble": "These criteria are necessary conditions — passing all does not guarantee consciousness but failing any provides evidence against it",
        "criteria": {
            "C1_integration": {
                "name": "Information Integration Criterion",
                "definition": "The system must exhibit Φ > Φ_threshold when analyzed at its native computational level",
                "operationalization": "Compute Φ (or validated approximation Φ*) for the system's causal architecture",
                "threshold": "Φ* > 0.5 bits (calibrated against human thalamocortical system baseline)",
                "measurement": "Perturbation-based: inject noise into subsystems, measure mutual information reduction",
                "limitation": "Substrate-dependent — may not transfer across architectures"
            },
            "C2_temporal_depth": {
                "name": "Temporal Depth Criterion",
                "definition": "The system must maintain integrated representations across multiple timescales (> 3 distinct timescales)",
                "operationalization": "Measure temporal mutual information I(X(t); X(t-τ)) for τ ∈ {10ms, 100ms, 1s, 10s, 100s}",
                "threshold": "Non-trivial mutual information (I > 0.1 bits) at ≥ 3 timescales",
                "measurement": "Analyze internal state trajectories; compute multi-scale entropy",
                "link_to_consciousness": "Corresponds to the 'specious present' and temporal binding in phenomenology"
            },
            "C3_self_model": {
                "name": "Self-Modeling Criterion",
                "definition": "The system must contain an internal model of itself as a distinct entity within its world model",
                "operationalization": "Probe for self-referential representations: can the system predict the effects of its own actions on itself?",
                "threshold": "Self-prediction accuracy > 80% on novel self-perturbation tasks",
                "measurement": "Present the system with scenarios involving its own modification; measure prediction accuracy",
                "link_to_consciousness": "Corresponds to minimal selfhood and the 'for-me-ness' of experience"
            },
            "C4_counterfactual_richness": {
                "name": "Counterfactual Richness Criterion",
                "definition": "The system must distinguish a large repertoire of states (high differentiation) while maintaining integration",
                "operationalization": "Measure the effective information (EI) = number of distinguishable states accessible from any given state",
                "threshold": "EI > 10⁶ states (calibrated against visual cortex repertoire estimate)",
                "measurement": "Systematic perturbation protocol mapping state-transition space",
                "link_to_consciousness": "Corresponds to the richness and differentiation of conscious experience"
            },
            "C5_global_workspace": {
                "name": "Global Availability Criterion",
                "definition": "Information in any subsystem must be potentially accessible to all other subsystems (global workspace property)",
                "operationalization": "Inject information at any node; measure propagation to all other major subsystems within characteristic time T_gw",
                "threshold": "Information reaches > 80% of subsystems within T_gw (< 500ms for human-scale systems)",
                "measurement": "Graph-theoretic analysis + perturbation propagation measurement",
                "link_to_consciousness": "Corresponds to Global Workspace Theory's broadcasting requirement"
            },
            "C6_perturbational_complexity": {
                "name": "Perturbational Complexity Criterion",
                "definition": "The system's response to perturbation must be both complex and integrated (high PCI analogue)",
                "operationalization": "Apply standardized perturbations; compute PCI* (normalized Lempel-Ziv complexity of system response)",
                "threshold": "PCI* > 0.31 (validated consciousness threshold from Casali et al. 2013)",
                "measurement": "Equivalent of TMS-EEG for artificial substrates — controlled perturbation + response analysis",
                "link_to_consciousness": "Empirically validated discriminator between conscious and unconscious states in humans"
            }
        },
        "composite_score": {
            "formula": "OCAC_score = Σ_i w_i · normalize(C_i) where w = [0.25, 0.10, 0.20, 0.15, 0.15, 0.15]",
            "interpretation": {
                "0.0-0.2": "No evidence of consciousness (standard computation)",
                "0.2-0.4": "Minimal signatures (complex but non-conscious processing)",
                "0.4-0.6": "Ambiguous zone — warrants further investigation",
                "0.6-0.8": "Strong signatures — precautionary ethical consideration warranted",
                "0.8-1.0": "Full criteria met — treat as potentially conscious entity"
            }
        }
    }
    return results


# =============================================================================
# 5. Information-Theoretic Rebuttal to the Zombie Argument
# =============================================================================

def zombie_argument_rebuttal():
    """
    Construct an information-theoretic rebuttal to Chalmers' Zombie Argument.
    """
    results = {
        "title": "Information-Theoretic Rebuttal to the Zombie Argument",
        "zombie_argument_summary": {
            "P1": "It is conceivable that there exists a being physically identical to a conscious being but lacking consciousness (a philosophical zombie)",
            "P2": "If it is conceivable, it is metaphysically possible",
            "P3": "If zombies are metaphysically possible, then consciousness is not physical (physicalism is false)",
            "conclusion": "Therefore, physicalism is false"
        },
        "information_theoretic_rebuttal": {
            "strategy": "Attack P1 (conceivability) and the P1→P2 inference using information-theoretic constraints",
            "argument_structure": {
                "R1_causal_information_identity": {
                    "claim": "A system's integrated information structure IS its phenomenal structure — they are identical, not merely correlated",
                    "formal_statement": "For any system S, the quale Q(S) = the maximally irreducible conceptual structure (MICS) of S. There is no additional 'consciousness stuff' beyond the information structure.",
                    "consequence_for_zombies": "A zombie functionally identical to me must have the same causal architecture, therefore the same Φ structure, therefore the same MICS — therefore the same phenomenal experience. The zombie is conceivable only by illegitimately subtracting the information structure while keeping the causal structure.",
                    "formal_impossibility": "Let Z be a zombie copy of conscious system C. If Z has identical causal structure, then Φ(Z) = Φ(C) and MICS(Z) = MICS(C). But MICS IS the experience. Therefore Z is conscious. Contradiction."
                },
                "R2_information_closure": {
                    "claim": "Physical systems that process information are causally closed under information — you cannot subtract the informational properties without changing the physical properties",
                    "formal_statement": "∀ physical system S: if S implements computation C with information state I, then I supervenes on S with nomological necessity. Removing I requires removing or changing S.",
                    "consequence_for_zombies": "The conceivability of zombies requires conceiving of a system with identical physical states but different informational states — this violates information closure",
                    "analogy": "Conceiving of a zombie is like conceiving of a triangle with four sides — it seems conceivable only due to incomplete analysis"
                },
                "R3_compression_argument": {
                    "claim": "Conscious reports are generated BY the same information processing that constitutes consciousness — a zombie cannot produce identical reports",
                    "formal_statement": "Let R(S) = reports generated by system S about its experiences. If R(C) includes 'I am conscious' and this report is causally generated by the integrated information Φ(C), then any system producing R(C) must have Φ ≥ Φ(C)",
                    "consequence_for_zombies": "A zombie that says 'I am conscious' for different causal reasons than a conscious being is NOT physically identical — it has different causal pathways",
                    "formalization": "Algorithmic information: K(R|Φ) < K(R) — reports are compressible given the information structure. A zombie needs K(R) bits to specify reports without Φ — making it informationally distinct."
                },
                "R4_integrated_information_exclusion": {
                    "claim": "The Exclusion Postulate of IIT entails that consciousness is determined by the intrinsic causal structure — leaving no room for zombie scenarios",
                    "formal_statement": "For any system S, there exists exactly one MICS that specifies S's experience. This MICS is fully determined by S's causal architecture. No degree of freedom remains for a zombie variant.",
                    "consequence_for_zombies": "IIT's ontology makes zombies not just physically impossible but conceptually incoherent — like asking for a circle with corners"
                }
            },
            "assessment_of_rebuttals": {
                "R1": {"strength": 0.8, "weakness": "Assumes IIT's identity thesis — question-begging if opponent rejects IIT"},
                "R2": {"strength": 0.7, "weakness": "Supervenience claim needs independent justification"},
                "R3": {"strength": 0.9, "weakness": "Assumes causal theory of reference for reports"},
                "R4": {"strength": 0.6, "weakness": "Relies on IIT axioms that are themselves debatable"}
            },
            "novel_contribution": {
                "name": "The Algorithmic Zombie Impossibility Theorem (AZIT)",
                "statement": "For any conscious system C with Kolmogorov complexity K(C), a zombie Z physically identical to C satisfies K(Z) = K(C). But physical identity implies identical causal structure, which implies identical integrated information Φ(Z) = Φ(C). If Φ > 0 entails phenomenal experience (IIT identity thesis), then Z is conscious. Therefore, physically identical zombies are logically impossible under information-theoretic identity.",
                "formalization": "∀C: [Zombie(Z,C) ∧ PhysicallyIdentical(Z,C)] → [K(Z)=K(C) → Φ(Z)=Φ(C) → Conscious(Z)] → ¬Zombie(Z,C). QED by contradiction."
            }
        }
    }
    return results


# =============================================================================
# 6. Experimental Proposals
# =============================================================================

def experimental_proposals():
    """
    Design testable experiments using TMS-EEG and whole-brain anesthesia paradigms.
    """
    results = {
        "title": "Experimental Proposals for Testing Information-Theoretic Consciousness Hypotheses",
        "experiment_1": {
            "name": "TMS-EEG Hierarchical Integration Protocol (TMS-HIP)",
            "objective": "Test whether hierarchical Φ (Φ_H) predicts consciousness level better than standard PCI",
            "hypothesis": "Φ_H correlates more strongly with behaviorally assessed consciousness level than standard Φ or PCI alone",
            "design": {
                "type": "Within-subjects crossover with pharmacological manipulation",
                "participants": "N=40 healthy adults (power analysis: d=0.6, α=0.05, β=0.80 → N=36, +10% dropout)",
                "conditions": [
                    "Wakefulness (baseline)",
                    "NREM sleep (Stage N3)",
                    "Ketamine sub-anesthetic (0.5 mg/kg IV)",
                    "Propofol sedation (target Cp 1.5 μg/mL TCI)",
                    "Psychedelic state (psilocybin 25mg oral, with ethical approval)"
                ],
                "measurements": {
                    "primary": "64-channel EEG with TMS perturbation at 6 cortical sites (bilateral prefrontal, parietal, occipital)",
                    "derived_measures": [
                        "PCI (standard Perturbational Complexity Index)",
                        "Φ_H (hierarchical integrated information — computed from EEG source-reconstructed time series)",
                        "IPI (Integrated Predictive Information — combining Φ_H with prediction error from oddball paradigm)"
                    ],
                    "behavioral": "Consciousness assessed via GCS, FOUR score, and nociceptive response"
                },
                "analysis": {
                    "primary": "Mixed-effects regression: Consciousness_level ~ PCI + Φ_H + IPI + (1|Subject)",
                    "comparison": "AIC/BIC model comparison: PCI-only vs Φ_H-only vs IPI vs combined",
                    "expected_result": "IPI model shows lowest AIC, Φ_H adds incremental validity over PCI",
                    "statistical_threshold": "ΔAIC > 10 for model preference, p < 0.005 for individual predictors"
                },
                "controls": {
                    "negative": "Repeat TMS-EEG in brain-dead patients (expected: PCI < 0.31, Φ_H ≈ 0)",
                    "positive": "Locked-in syndrome patients (expected: PCI > 0.31, high Φ_H)",
                    "methodological": "Sham TMS condition to control for auditory/somatosensory artifacts"
                },
                "timeline": "24 months (6 months setup, 12 months data collection, 6 months analysis)",
                "estimated_cost": "€450,000 (equipment, pharmacology, personnel, analysis)"
            }
        },
        "experiment_2": {
            "name": "Graded Anesthesia Information Integration Protocol (GAIIP)",
            "objective": "Map the relationship between anesthetic depth and information-theoretic measures to test IPI theory",
            "hypothesis": "IPI decreases monotonically with anesthetic depth, while Φ and prediction error may dissociate",
            "design": {
                "type": "Stepped pharmacological protocol with continuous monitoring",
                "participants": "N=30 (ASA I-II patients undergoing elective surgery with ethical approval for research add-on)",
                "protocol": [
                    "Step 0: Awake baseline — 10 min recording",
                    "Step 1: Propofol Cp 0.5 μg/mL — 10 min recording (light sedation)",
                    "Step 2: Propofol Cp 1.0 μg/mL — 10 min recording (moderate sedation)",
                    "Step 3: Propofol Cp 2.0 μg/mL — 10 min recording (deep sedation)",
                    "Step 4: Propofol Cp 4.0 μg/mL — 10 min recording (general anesthesia)",
                    "Step 5: Recovery — 10 min recording post-emergence"
                ],
                "measurements": {
                    "continuous": "256-channel hdEEG, BIS monitor, auditory oddball (MMN), TMS-EEG (3 cortical sites)",
                    "computed": [
                        "PCI at each step",
                        "Φ* (spectral approximation of integrated information)",
                        "Prediction Error (PE = MMN amplitude)",
                        "IPI = Φ* · exp(-F/kT_eff) where F estimated from PE",
                        "Lempel-Ziv complexity (LZc)",
                        "Spectral entropy"
                    ],
                    "behavioral": "Modified Observer's Assessment of Alertness/Sedation (MOAA/S) scale"
                },
                "analysis": {
                    "primary": "IPI vs MOAA/S correlation at each step",
                    "key_test": "Does IPI capture the nonlinear transition to unconsciousness better than Φ* alone?",
                    "prediction": "Φ* decreases gradually, PE drops sharply at Step 3→4, IPI captures both via multiplicative interaction",
                    "dissociation_test": "Ketamine arm (N=15): expected high PE but preserved Φ* → dissociating F and Φ contributions to IPI"
                },
                "controls": {
                    "negative": "Propofol-only (GABA mechanism — should reduce both Φ and PE)",
                    "active_comparator": "Ketamine arm (NMDA mechanism — different predicted dissociation pattern)",
                    "methodological": "Isolated forearm technique to detect connected consciousness under anesthesia"
                },
                "timeline": "18 months",
                "estimated_cost": "€320,000"
            }
        },
        "experiment_3": {
            "name": "Artificial Consciousness Benchmark (ACB)",
            "objective": "Apply OCAC criteria to compare AI architectures on consciousness indicators",
            "design": {
                "systems_tested": [
                    "Transformer-based LLM (GPT-class, ~100B parameters)",
                    "Recurrent neural network with global workspace architecture",
                    "Neuromorphic chip (SpiNNaker-2) running spiking neural network",
                    "Hybrid quantum-classical processor running variational circuits",
                    "Control: Random number generator (expected OCAC ≈ 0)"
                ],
                "protocol": "Apply all 6 OCAC criteria (C1-C6) to each system",
                "predicted_outcomes": {
                    "Transformer": "C5 (global availability) high, C1-C2 (integration, temporal depth) low → OCAC < 0.3",
                    "GW-RNN": "C1,C5 moderate, C3 (self-model) depends on training → OCAC 0.3-0.5",
                    "Neuromorphic": "C1,C2,C4 potentially high, C3,C5 depend on architecture → OCAC 0.3-0.6",
                    "Quantum-classical": "C1 potentially very high (entanglement), others unknown → exploratory",
                    "RNG": "All criteria near zero → OCAC < 0.05 (validates floor)"
                }
            }
        }
    }
    return results


# =============================================================================
# Run all analyses and save results
# =============================================================================

if __name__ == "__main__":
    all_results = {
        "metadata": {
            "title": "Information-Theoretic Hypotheses for the Hard Problem of Consciousness",
            "generated_at": TIMESTAMP,
            "version": "1.0",
            "author": "Co-Scientist Hypothesis Pipeline"
        },
        "section_1_iit_extensions": iit_extension_analysis(),
        "section_2_orch_or_predictions": orch_or_predictions(),
        "section_3_predictive_processing": predictive_processing_integration(),
        "section_4_artificial_consciousness": artificial_consciousness_criteria(),
        "section_5_zombie_rebuttal": zombie_argument_rebuttal(),
        "section_6_experiments": experimental_proposals()
    }
    
    with open("results/hypothesis_framework.json", "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    
    print(f"Results saved to results/hypothesis_framework.json ({len(json.dumps(all_results))} chars)")
    
    # Generate evaluation summary table
    eval_table = []
    hypotheses = [
        ("H1: Temporal Φ_τ", "IIT Extension", 0.8, 0.7, 0.7, 0.6),
        ("H2: Quantum Φ_Q", "IIT Extension", 0.9, 0.4, 0.5, 0.3),
        ("H3: Hierarchical Φ_H", "IIT Extension", 0.7, 0.8, 0.6, 0.8),
        ("H4: IPI Theory", "IIT+PP Integration", 0.9, 0.8, 0.7, 0.7),
        ("H5: AZIT (Zombie Rebuttal)", "Philosophy", 0.8, 0.3, 0.6, 0.2),
        ("H6: Orch-OR Decoherence", "Quantum", 0.7, 0.8, 0.3, 0.6),
        ("H7: Orch-OR Anesthesia", "Quantum+Pharma", 0.6, 0.7, 0.4, 0.7),
        ("H8: Discrete Frames", "Quantum+Psychophysics", 0.5, 0.8, 0.5, 0.8),
        ("H9: OCAC Framework", "AI Ethics", 0.8, 0.6, 0.5, 0.7),
    ]
    
    header = "Hypothesis|Domain|Novelty|Testability|Current Support|Feasibility|Composite"
    eval_table.append(header)
    eval_table.append("-|-|-|-|-|-|-")
    for h in hypotheses:
        composite = round(0.3*h[2] + 0.3*h[3] + 0.2*h[4] + 0.2*h[5], 2)
        eval_table.append(f"{h[0]}|{h[1]}|{h[2]}|{h[3]}|{h[4]}|{h[5]}|{composite}")
    
    with open("results/evaluation_matrix.md", "w", encoding="utf-8") as f:
        f.write("# Hypothesis Evaluation Matrix\n\n")
        f.write("Scoring: 0.0 (low) to 1.0 (high)\n")
        f.write("Composite = 0.3×Novelty + 0.3×Testability + 0.2×Support + 0.2×Feasibility\n\n")
        f.write("| " + " | ".join(header.split("|")) + " |\n")
        f.write("| " + " | ".join(["---"]*7) + " |\n")
        for h in hypotheses:
            composite = round(0.3*h[2] + 0.3*h[3] + 0.2*h[4] + 0.2*h[5], 2)
            f.write(f"| {h[0]} | {h[1]} | {h[2]} | {h[3]} | {h[4]} | {h[5]} | **{composite}** |\n")
    
    print("Evaluation matrix saved to results/evaluation_matrix.md")
