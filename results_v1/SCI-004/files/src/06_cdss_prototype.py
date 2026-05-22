"""
Module 6: Clinical Decision Support System (CDSS) Prototype
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# CPIC/DPWG Guideline Rules Engine
# ─────────────────────────────────────────────────────────────
CPIC_RULES = {
    'Codeine': {
        'gene': 'CYP2D6',
        'rules': {
            'Poor Metabolizer (PM)':       {'recommendation': 'AVOID_USE', 'reason': 'Poor conversion to morphine; inadequate analgesia. Use alternative analgesic.', 'level': 'STRONG'},
            'Ultrarapid Metabolizer (UM)': {'recommendation': 'AVOID_USE', 'reason': 'Excess morphine production; risk of opioid toxicity/death.', 'level': 'STRONG'},
            'Intermediate Metabolizer (IM)': {'recommendation': 'USE_WITH_CAUTION', 'reason': 'Reduced but some CYP2D6 activity. Use lower dose or alternative.', 'level': 'MODERATE'},
            'Normal Metabolizer (NM)':     {'recommendation': 'STANDARD_DOSE', 'reason': 'Normal CYP2D6 activity. Standard dosing.', 'level': 'STANDARD'},
        }
    },
    'Clopidogrel': {
        'gene': 'CYP2C19',
        'rules': {
            'Poor Metabolizer (PM)':       {'recommendation': 'ALTERNATIVE_DRUG', 'reason': 'Markedly reduced active metabolite. Use prasugrel or ticagrelor.', 'level': 'STRONG'},
            'Intermediate Metabolizer (IM)': {'recommendation': 'CONSIDER_ALTERNATIVE', 'reason': 'Reduced active metabolite. Consider prasugrel/ticagrelor for ACS.', 'level': 'MODERATE'},
            'Rapid/Ultrarapid Metabolizer (RM/UM)': {'recommendation': 'STANDARD_DOSE', 'reason': 'Normal or increased active metabolite. Standard dosing.', 'level': 'STANDARD'},
            'Normal Metabolizer (NM)':     {'recommendation': 'STANDARD_DOSE', 'reason': 'Normal CYP2C19 activity.', 'level': 'STANDARD'},
        }
    },
    'Carbamazepine': {
        'gene': 'HLA-B',
        'rules': {
            'HLA-B*15:02_positive': {'recommendation': 'AVOID_USE', 'reason': 'High risk of SJS/TEN (OR=80). Use alternative AED.', 'level': 'STRONG'},
            'HLA-B*15:02_negative': {'recommendation': 'STANDARD_DOSE', 'reason': 'Low HLA-mediated SJS/TEN risk. Standard monitoring.', 'level': 'STANDARD'},
        }
    },
    'Warfarin': {
        'gene': 'CYP2C9_VKORC1',
        'rules': {
            'Sensitive':       {'recommendation': 'REDUCE_DOSE', 'reason': 'CYP2C9*2/*3 or VKORC1 -1639A: start with reduced dose, frequent INR monitoring.', 'level': 'STRONG'},
            'Normal':          {'recommendation': 'STANDARD_DOSE', 'reason': 'Standard warfarin dosing with INR monitoring.', 'level': 'STANDARD'},
            'Resistant':       {'recommendation': 'INCREASE_DOSE', 'reason': 'VKORC1 -1639G/G: may require higher doses.', 'level': 'MODERATE'},
        }
    },
    'Abacavir': {
        'gene': 'HLA-B',
        'rules': {
            'HLA-B*57:01_positive': {'recommendation': 'AVOID_USE', 'reason': 'High risk of hypersensitivity reaction (OR=117). Absolute contraindication.', 'level': 'STRONG'},
            'HLA-B*57:01_negative': {'recommendation': 'STANDARD_DOSE', 'reason': 'Low hypersensitivity risk. Standard monitoring.', 'level': 'STANDARD'},
        }
    },
    'Allopurinol': {
        'gene': 'HLA-B',
        'rules': {
            'HLA-B*58:01_positive': {'recommendation': 'AVOID_USE', 'reason': 'Extremely high SJS/TEN risk (OR=580). Use febuxostat.', 'level': 'STRONG'},
            'HLA-B*58:01_negative': {'recommendation': 'STANDARD_DOSE', 'reason': 'Low HLA-mediated risk.', 'level': 'STANDARD'},
        }
    },
}

def cdss_evaluate(patient_data: dict) -> dict:
    """Main CDSS evaluation function."""
    alerts = []
    recommendations = []
    risk_score = 0

    # Extract patient genomic profile
    cyp2d6_pheno  = patient_data.get('CYP2D6_phenotype', 'Normal Metabolizer (NM)')
    cyp2c19_pheno = patient_data.get('CYP2C19_phenotype', 'Normal Metabolizer (NM)')
    hla_b1502     = patient_data.get('HLA_B1502', 0)
    hla_b5701     = patient_data.get('HLA_B5701', 0)
    hla_b5801     = patient_data.get('HLA_B5801', 0)

    medications = patient_data.get('medications', [])

    for drug in medications:
        if drug in CPIC_RULES:
            rule_set = CPIC_RULES[drug]
            gene     = rule_set['gene']

            # Determine applicable phenotype key
            if gene == 'CYP2D6':
                pheno_key = cyp2d6_pheno
            elif gene == 'CYP2C19':
                pheno_key = cyp2c19_pheno
            elif gene == 'HLA-B':
                if drug == 'Carbamazepine':
                    pheno_key = 'HLA-B*15:02_positive' if hla_b1502 else 'HLA-B*15:02_negative'
                elif drug == 'Abacavir':
                    pheno_key = 'HLA-B*57:01_positive' if hla_b5701 else 'HLA-B*57:01_negative'
                elif drug == 'Allopurinol':
                    pheno_key = 'HLA-B*58:01_positive' if hla_b5801 else 'HLA-B*58:01_negative'
                else:
                    continue
            else:
                # Use patient-reported sensitivity
                pheno_key = patient_data.get('warfarin_sensitivity', 'Normal')

            if pheno_key in rule_set['rules']:
                rule  = rule_set['rules'][pheno_key]
                level = rule['level']
                if rule['recommendation'] == 'AVOID_USE':
                    risk_score += 10
                    alerts.append({
                        'type': 'CONTRAINDICATION',
                        'drug': drug,
                        'gene': gene,
                        'phenotype': pheno_key,
                        'message': rule['reason'],
                        'level': level,
                        'cpic_guideline': True,
                    })
                elif rule['recommendation'] in ['ALTERNATIVE_DRUG', 'CONSIDER_ALTERNATIVE', 'REDUCE_DOSE']:
                    risk_score += 5
                    alerts.append({
                        'type': 'DOSE_ADJUSTMENT',
                        'drug': drug,
                        'gene': gene,
                        'phenotype': pheno_key,
                        'message': rule['reason'],
                        'level': level,
                        'cpic_guideline': True,
                    })
                recommendations.append({
                    'drug': drug,
                    'recommendation': rule['recommendation'],
                    'reason': rule['reason'],
                    'level': level,
                })

    # Drug-drug interaction check (simplified)
    if 'Codeine' in medications and 'Fluoxetine' in medications:
        alerts.append({
            'type': 'DRUG_DRUG_INTERACTION',
            'drugs': ['Codeine', 'Fluoxetine'],
            'message': 'Fluoxetine inhibits CYP2D6, may increase codeine-related toxicity.',
            'level': 'MODERATE',
        })
        risk_score += 4

    overall_risk = 'HIGH' if risk_score >= 10 else 'MODERATE' if risk_score >= 5 else 'LOW'

    return {
        'patient_id': patient_data.get('patient_id', 'UNKNOWN'),
        'evaluation_timestamp': datetime.utcnow().isoformat(),
        'overall_risk': overall_risk,
        'risk_score': risk_score,
        'alerts': alerts,
        'recommendations': recommendations,
        'summary': f"{len(alerts)} alert(s) for {len(medications)} medication(s). Overall risk: {overall_risk}.",
    }

# ─────────────────────────────────────────────────────────────
# Test CDSS with 5 representative patient cases
# ─────────────────────────────────────────────────────────────
test_patients = [
    {
        'patient_id': 'PT001',
        'age': 45, 'sex': 'M', 'ethnicity': 'Asian',
        'CYP2D6_phenotype': 'Ultrarapid Metabolizer (UM)',
        'CYP2C19_phenotype': 'Poor Metabolizer (PM)',
        'HLA_B1502': 1, 'HLA_B5701': 0, 'HLA_B5801': 0,
        'medications': ['Codeine', 'Clopidogrel', 'Carbamazepine'],
        'diagnosis': 'Post-PCI + Epilepsy + Acute pain',
    },
    {
        'patient_id': 'PT002',
        'age': 62, 'sex': 'F', 'ethnicity': 'European',
        'CYP2D6_phenotype': 'Poor Metabolizer (PM)',
        'CYP2C19_phenotype': 'Normal Metabolizer (NM)',
        'HLA_B1502': 0, 'HLA_B5701': 1, 'HLA_B5801': 0,
        'medications': ['Abacavir', 'Codeine', 'Fluoxetine'],
        'diagnosis': 'HIV + Chronic pain + Depression',
    },
    {
        'patient_id': 'PT003',
        'age': 35, 'sex': 'M', 'ethnicity': 'European',
        'CYP2D6_phenotype': 'Normal Metabolizer (NM)',
        'CYP2C19_phenotype': 'Normal Metabolizer (NM)',
        'HLA_B1502': 0, 'HLA_B5701': 0, 'HLA_B5801': 0,
        'medications': ['Clopidogrel'],
        'diagnosis': 'ACS post-stent',
    },
    {
        'patient_id': 'PT004',
        'age': 70, 'sex': 'F', 'ethnicity': 'Asian',
        'CYP2D6_phenotype': 'Intermediate Metabolizer (IM)',
        'CYP2C19_phenotype': 'Intermediate Metabolizer (IM)',
        'HLA_B1502': 0, 'HLA_B5701': 0, 'HLA_B5801': 1,
        'medications': ['Allopurinol', 'Warfarin'],
        'warfarin_sensitivity': 'Sensitive',
        'diagnosis': 'Gout + Atrial Fibrillation',
    },
    {
        'patient_id': 'PT005',
        'age': 28, 'sex': 'M', 'ethnicity': 'African',
        'CYP2D6_phenotype': 'Normal Metabolizer (NM)',
        'CYP2C19_phenotype': 'Rapid/Ultrarapid Metabolizer (RM/UM)',
        'HLA_B1502': 0, 'HLA_B5701': 0, 'HLA_B5801': 0,
        'medications': ['Clopidogrel', 'Codeine'],
        'diagnosis': 'ACS + Acute pain',
    },
]

cdss_results = []
for patient in test_patients:
    result = cdss_evaluate(patient)
    cdss_results.append(result)
    print(f"  [{result['patient_id']}] Risk={result['overall_risk']} "
          f"({result['risk_score']}) | {len(result['alerts'])} alerts")

# ─────────────────────────────────────────────────────────────
# Figure 12: CDSS Risk Dashboard
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Risk scores by patient
patient_ids   = [r['patient_id'] for r in cdss_results]
risk_scores   = [r['risk_score'] for r in cdss_results]
n_alerts      = [len(r['alerts']) for r in cdss_results]
risk_colors   = ['#d73027' if s >= 10 else '#fc8d59' if s >= 5 else '#1a9850' for s in risk_scores]

ax = axes[0]
bars = ax.bar(patient_ids, risk_scores, color=risk_colors, alpha=0.85)
ax.set_title('CDSS Risk Score by Patient', fontsize=12)
ax.set_xlabel('Patient ID'); ax.set_ylabel('PGx Risk Score')
ax.set_ylim(0, 35)
for bar, score, n_al in zip(bars, risk_scores, n_alerts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f'{score}\n({n_al} alerts)', ha='center', va='bottom', fontsize=9)
patches_leg = [
    mpatches.Patch(color='#d73027', label='HIGH risk (≥10)'),
    mpatches.Patch(color='#fc8d59', label='MODERATE risk (5–9)'),
    mpatches.Patch(color='#1a9850', label='LOW risk (<5)'),
]
ax.legend(handles=patches_leg, loc='upper right', fontsize=9)

# Right: Alert type distribution
alert_types = {}
for r in cdss_results:
    for alert in r['alerts']:
        t = alert['type']
        alert_types[t] = alert_types.get(t, 0) + 1

ax2 = axes[1]
if alert_types:
    colors_alert = ['#d73027', '#fc8d59', '#4575b4', '#1a9850']
    wedges, texts, autotexts = ax2.pie(
        alert_types.values(), labels=alert_types.keys(),
        colors=colors_alert[:len(alert_types)], autopct='%1.1f%%', startangle=90
    )
    ax2.set_title('Alert Type Distribution Across Test Patients', fontsize=12)
else:
    ax2.text(0.5, 0.5, 'No alerts generated', ha='center', va='center')

plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig12_cdss_dashboard.png',
            dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# Figure 13: CDSS workflow architecture
# ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 7))
ax.axis('off')

def draw_box(ax, x, y, w, h, label, color, fontsize=9):
    rect = plt.Rectangle((x, y), w, h, linewidth=1.5, edgecolor='black',
                           facecolor=color, alpha=0.85)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, label, ha='center', va='center',
            fontsize=fontsize, fontweight='bold', wrap=True,
            multialignment='center')

def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

# Input layer
draw_box(ax, 0.05, 0.75, 0.18, 0.15, 'Patient EHR\n(Demographics\n+ Medications)', '#AED6F1')
draw_box(ax, 0.27, 0.75, 0.18, 0.15, 'Genomic\nProfile\n(SNP array/NGS)', '#A9DFBF')
draw_box(ax, 0.49, 0.75, 0.18, 0.15, 'Lab Values\n(Renal/Hepatic\nFunction)', '#F9E79F')

# Processing layer
draw_box(ax, 0.16, 0.50, 0.20, 0.15, 'PGx Phenotype\nClassifier\n(CYP/HLA)', '#85C1E9')
draw_box(ax, 0.40, 0.50, 0.20, 0.15, 'ML Prediction\nModels\n(Drug Response)', '#82E0AA')
draw_box(ax, 0.64, 0.50, 0.20, 0.15, 'CPIC/DPWG\nGuideline\nRules Engine', '#F8C471')

# Output layer
draw_box(ax, 0.10, 0.25, 0.22, 0.15, 'Risk Stratification\n(HIGH/MOD/LOW)', '#E8DAEF')
draw_box(ax, 0.36, 0.25, 0.28, 0.15, 'Drug Recommendations\n(AVOID/ADJUST/STANDARD)', '#FADBD8')
draw_box(ax, 0.68, 0.25, 0.22, 0.15, 'Alert &\nNotification\nSystem', '#D5F5E3')

# Bottom
draw_box(ax, 0.25, 0.05, 0.50, 0.12, 'Clinical Report (PDF/HL7 FHIR)\nAudit Trail & EHR Integration', '#D6EAF8')

# Arrows
for x_src in [0.14, 0.36, 0.58]:
    draw_arrow(ax, x_src, 0.75, 0.26+0.0, 0.65)
for x_src in [0.26, 0.50, 0.74]:
    draw_arrow(ax, x_src, 0.50, 0.50, 0.40)
for x_src in [0.21, 0.50, 0.79]:
    draw_arrow(ax, x_src, 0.25, 0.50, 0.17)

ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title('CDSS Architecture: Pharmacogenomics Clinical Decision Support System',
              fontsize=12, fontweight='bold', pad=10)
plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig13_cdss_architecture.png',
            dpi=150, bbox_inches='tight')
plt.close()

# Save CDSS results
with open('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/results/cdss_evaluation_results.json', 'w') as f:
    json.dump({'test_cases': cdss_results, 'cpic_rules': CPIC_RULES}, f, indent=2, ensure_ascii=False)

# Save CDSS patient report as readable table
rows = []
for p, r in zip(test_patients, cdss_results):
    for alert in r['alerts']:
        rows.append({
            'patient_id': p['patient_id'],
            'diagnosis': p['diagnosis'],
            'overall_risk': r['overall_risk'],
            'alert_type': alert['type'],
            'drug': alert.get('drug', str(alert.get('drugs',''))),
            'gene': alert.get('gene', '-'),
            'phenotype': alert.get('phenotype', '-'),
            'recommendation_message': alert['message'][:80],
            'cpic_level': alert.get('level','-'),
        })
df_report = pd.DataFrame(rows)
df_report.to_csv('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/results/cdss_patient_alerts.csv', index=False)

print("[CDSS Module] Done")
print(f"  Total alerts across {len(test_patients)} patients: {sum(len(r['alerts']) for r in cdss_results)}")
