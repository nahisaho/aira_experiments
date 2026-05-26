#!/usr/bin/env python3
"""
Food Supply Chain Safety Risk Prediction AI System
Integrated experiment: spatiotemporal prediction, NLP recall detection,
microbial growth modeling, HACCP scoring, blockchain traceability, Salmonella case study.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import (accuracy_score, precision_score, recall_score, f1_score,
                             mean_squared_error, mean_absolute_error, r2_score,
                             roc_auc_score, roc_curve, confusion_matrix, classification_report)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import warnings
import os
import json

warnings.filterwarnings('ignore')
np.random.seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

results = {}

# ============================================================
# Module 1: Spatiotemporal Foodborne Illness Prediction
# ============================================================
def generate_spatiotemporal_data(n_samples=2000):
    """Generate synthetic spatiotemporal foodborne illness data."""
    dates = pd.date_range('2018-01-01', periods=n_samples, freq='D')
    np.random.seed(42)
    
    temperature = 15 + 12 * np.sin(2 * np.pi * np.arange(n_samples) / 365) + np.random.normal(0, 3, n_samples)
    humidity = 60 + 15 * np.sin(2 * np.pi * np.arange(n_samples) / 365 + np.pi/4) + np.random.normal(0, 8, n_samples)
    month = np.array([d.month for d in dates])
    day_of_year = np.array([d.dayofyear for d in dates])
    
    # Seasonal risk: higher in summer
    seasonal_risk = 0.3 * np.sin(2 * np.pi * day_of_year / 365 - np.pi/6)
    
    lat = np.random.uniform(25, 45, n_samples)
    lon = np.random.uniform(-125, -70, n_samples)
    
    # Incident probability depends on temperature, humidity, season
    logit = (-3 + 0.08 * temperature + 0.02 * humidity + seasonal_risk +
             0.01 * (temperature - 25) * (humidity - 70) / 100 + np.random.normal(0, 0.5, n_samples))
    prob = 1 / (1 + np.exp(-logit))
    incidents = (np.random.random(n_samples) < prob).astype(int)
    
    df = pd.DataFrame({
        'date': dates, 'temperature': temperature, 'humidity': humidity,
        'month': month, 'day_of_year': day_of_year, 'latitude': lat,
        'longitude': lon, 'incident': incidents
    })
    return df


def run_spatiotemporal_model():
    print("=" * 60)
    print("Module 1: Spatiotemporal Foodborne Illness Prediction")
    print("=" * 60)
    
    df = generate_spatiotemporal_data()
    features = ['temperature', 'humidity', 'month', 'day_of_year', 'latitude', 'longitude']
    X = df[features].values
    y = df['incident'].values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)
    
    models = {
        'Logistic Regression': LogisticRegression(max_iter=1000),
        'Random Forest': RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42),
        'Gradient Boosting': GradientBoostingRegressor(n_estimators=200, max_depth=5, random_state=42),
        'MLP Neural Network': MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=500, random_state=42)
    }
    
    model_results = {}
    for name, model in models.items():
        if name == 'Gradient Boosting':
            model.fit(X_train, y_train)
            y_pred_prob = model.predict(X_test)
            y_pred = (y_pred_prob > 0.5).astype(int)
        else:
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            y_pred_prob = model.predict_proba(X_test)[:, 1] if hasattr(model, 'predict_proba') else y_pred.astype(float)
        
        acc = accuracy_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred, zero_division=0)
        rec = recall_score(y_test, y_pred, zero_division=0)
        f1 = f1_score(y_test, y_pred, zero_division=0)
        try:
            auc = roc_auc_score(y_test, y_pred_prob)
        except:
            auc = 0.0
        
        model_results[name] = {'accuracy': acc, 'precision': prec, 'recall': rec, 'f1': f1, 'auc': auc}
        print(f"  {name}: Acc={acc:.4f}, F1={f1:.4f}, AUC={auc:.4f}")
    
    # Feature importance from Random Forest
    rf = models['Random Forest']
    importances = rf.feature_importances_
    
    # Plot 1: Model comparison
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    metrics_df = pd.DataFrame(model_results).T
    metrics_df[['accuracy', 'precision', 'recall', 'f1']].plot(kind='bar', ax=axes[0])
    axes[0].set_title('Model Performance Comparison')
    axes[0].set_ylabel('Score')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
    axes[0].legend(loc='lower right')
    axes[0].set_ylim(0, 1)
    
    # Feature importance
    feat_imp = pd.Series(importances, index=features).sort_values(ascending=True)
    feat_imp.plot(kind='barh', ax=axes[1], color='steelblue')
    axes[1].set_title('Feature Importance (Random Forest)')
    axes[1].set_xlabel('Importance')
    
    # ROC curves
    for name, model in models.items():
        if name == 'Gradient Boosting':
            y_score = model.predict(X_test)
        elif hasattr(model, 'predict_proba'):
            y_score = model.predict_proba(X_test)[:, 1]
        else:
            continue
        fpr, tpr, _ = roc_curve(y_test, y_score)
        axes[2].plot(fpr, tpr, label=f"{name} (AUC={model_results[name]['auc']:.3f})")
    axes[2].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[2].set_title('ROC Curves')
    axes[2].set_xlabel('False Positive Rate')
    axes[2].set_ylabel('True Positive Rate')
    axes[2].legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'spatiotemporal_model_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Monthly incident pattern
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    monthly = df.groupby('month')['incident'].mean()
    axes[0].bar(monthly.index, monthly.values, color='coral', edgecolor='black')
    axes[0].set_title('Monthly Foodborne Illness Incidence Rate')
    axes[0].set_xlabel('Month')
    axes[0].set_ylabel('Incidence Rate')
    axes[0].set_xticks(range(1, 13))
    
    # Temperature vs incident scatter
    axes[1].scatter(df[df['incident']==0]['temperature'], df[df['incident']==0]['humidity'],
                    alpha=0.3, s=10, label='No Incident', c='blue')
    axes[1].scatter(df[df['incident']==1]['temperature'], df[df['incident']==1]['humidity'],
                    alpha=0.3, s=10, label='Incident', c='red')
    axes[1].set_title('Temperature vs Humidity by Incident')
    axes[1].set_xlabel('Temperature (°C)')
    axes[1].set_ylabel('Humidity (%)')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'spatiotemporal_patterns.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    results['spatiotemporal'] = model_results
    return model_results


# ============================================================
# Module 2: NLP-based Recall/Alert Early Detection
# ============================================================
def generate_recall_data(n_samples=1000):
    """Generate synthetic FDA/RASFF-style recall alert texts."""
    np.random.seed(42)
    
    categories = ['recall', 'alert', 'warning', 'information', 'normal']
    risk_levels = ['high', 'medium', 'low']
    
    contaminants = ['Salmonella', 'Listeria monocytogenes', 'E. coli O157:H7',
                    'Staphylococcus aureus', 'Clostridium botulinum', 'aflatoxin',
                    'pesticide residue', 'undeclared allergen', 'foreign body']
    
    products = ['chicken', 'beef', 'pork', 'lettuce', 'spinach', 'cheese',
                'milk', 'eggs', 'seafood', 'frozen vegetables', 'canned soup']
    
    templates_recall = [
        "RECALL: {product} contaminated with {contaminant}. Potential serious health risk. Remove from shelves immediately.",
        "FDA Class I Recall: {contaminant} detected in {product}. Risk of severe illness. Distribution across multiple states.",
        "URGENT recall of {product} due to {contaminant} contamination. Consumer advisory issued.",
        "Voluntary recall of {product} products after {contaminant} found during routine testing. Affected lot numbers listed.",
    ]
    templates_alert = [
        "RASFF Alert: Border rejection of {product} from Country X due to {contaminant} exceeding limits.",
        "Alert notification: {contaminant} in {product} detected by member state. Follow-up measures required.",
        "Serious risk: {contaminant} in {product}. Rapid alert triggered under RASFF notification system.",
    ]
    templates_warning = [
        "Warning: Elevated levels of {contaminant} found in imported {product}. Monitoring increased.",
        "Advisory: {product} from certain regions may contain {contaminant}. Testing recommended.",
    ]
    templates_info = [
        "Information: Routine sampling of {product} shows compliance with safety standards.",
        "Update: Previous alert regarding {product} resolved. No {contaminant} detected in follow-up.",
    ]
    templates_normal = [
        "Weekly market report: {product} prices stable. Quality inspection passed.",
        "Supply chain update: {product} shipment arrived on schedule. All documentation in order.",
        "Seasonal forecast: {product} production expected to increase next quarter.",
    ]
    
    texts = []
    labels = []
    risk_labels = []
    
    for _ in range(n_samples):
        cat = np.random.choice(categories, p=[0.25, 0.2, 0.15, 0.15, 0.25])
        product = np.random.choice(products)
        contaminant = np.random.choice(contaminants)
        
        if cat == 'recall':
            text = np.random.choice(templates_recall).format(product=product, contaminant=contaminant)
            risk = np.random.choice(['high', 'medium'], p=[0.7, 0.3])
        elif cat == 'alert':
            text = np.random.choice(templates_alert).format(product=product, contaminant=contaminant)
            risk = np.random.choice(['high', 'medium'], p=[0.5, 0.5])
        elif cat == 'warning':
            text = np.random.choice(templates_warning).format(product=product, contaminant=contaminant)
            risk = 'medium'
        elif cat == 'information':
            text = np.random.choice(templates_info).format(product=product, contaminant=contaminant)
            risk = 'low'
        else:
            text = np.random.choice(templates_normal).format(product=product, contaminant=contaminant)
            risk = 'low'
        
        texts.append(text)
        labels.append(cat)
        risk_labels.append(risk)
    
    return pd.DataFrame({'text': texts, 'category': labels, 'risk_level': risk_labels})


def run_nlp_recall_detection():
    print("\n" + "=" * 60)
    print("Module 2: NLP-based Recall/Alert Early Detection")
    print("=" * 60)
    
    df = generate_recall_data()
    
    # Binary classification: urgent (recall/alert) vs non-urgent
    df['urgent'] = df['category'].isin(['recall', 'alert']).astype(int)
    
    # TF-IDF features
    tfidf = TfidfVectorizer(max_features=500, ngram_range=(1, 2), stop_words='english')
    X_tfidf = tfidf.fit_transform(df['text'])
    y_binary = df['urgent'].values
    
    X_train, X_test, y_train, y_test = train_test_split(X_tfidf, y_binary, test_size=0.2, random_state=42)
    
    nlp_models = {
        'Logistic Regression': LogisticRegression(max_iter=1000),
        'SVM': SVC(kernel='rbf', probability=True),
        'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
        'MLP': MLPClassifier(hidden_layer_sizes=(128, 64), max_iter=500, random_state=42)
    }
    
    nlp_results = {}
    best_model = None
    best_f1 = 0
    
    for name, model in nlp_models.items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_prob = model.predict_proba(X_test)[:, 1]
        
        acc = accuracy_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred)
        rec = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        auc = roc_auc_score(y_test, y_prob)
        
        nlp_results[name] = {'accuracy': acc, 'precision': prec, 'recall': rec, 'f1': f1, 'auc': auc}
        print(f"  {name}: Acc={acc:.4f}, F1={f1:.4f}, AUC={auc:.4f}")
        
        if f1 > best_f1:
            best_f1 = f1
            best_model = model
            best_name = name
    
    # Multi-class classification
    le = LabelEncoder()
    y_multi = le.fit_transform(df['category'])
    X_train_m, X_test_m, y_train_m, y_test_m = train_test_split(X_tfidf, y_multi, test_size=0.2, random_state=42)
    
    rf_multi = RandomForestClassifier(n_estimators=200, random_state=42)
    rf_multi.fit(X_train_m, y_train_m)
    y_pred_m = rf_multi.predict(X_test_m)
    multi_acc = accuracy_score(y_test_m, y_pred_m)
    multi_f1 = f1_score(y_test_m, y_pred_m, average='weighted')
    print(f"\n  Multi-class (5 categories): Acc={multi_acc:.4f}, Weighted F1={multi_f1:.4f}")
    
    # Plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Binary classification comparison
    nlp_df = pd.DataFrame(nlp_results).T
    nlp_df[['accuracy', 'precision', 'recall', 'f1']].plot(kind='bar', ax=axes[0])
    axes[0].set_title('NLP Binary Classification Performance')
    axes[0].set_ylabel('Score')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
    axes[0].set_ylim(0, 1)
    axes[0].legend(loc='lower right')
    
    # Confusion matrix (multi-class)
    cm = confusion_matrix(y_test_m, y_pred_m)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1],
                xticklabels=le.classes_, yticklabels=le.classes_)
    axes[1].set_title('Multi-class Confusion Matrix')
    axes[1].set_xlabel('Predicted')
    axes[1].set_ylabel('True')
    
    # Top TF-IDF features
    feature_names = tfidf.get_feature_names_out()
    if hasattr(best_model, 'coef_'):
        top_idx = np.argsort(np.abs(best_model.coef_[0]))[-15:]
        axes[2].barh(range(len(top_idx)), np.abs(best_model.coef_[0][top_idx]), color='teal')
        axes[2].set_yticks(range(len(top_idx)))
        axes[2].set_yticklabels([feature_names[i] for i in top_idx])
    else:
        top_idx = np.argsort(best_model.feature_importances_)[-15:]
        axes[2].barh(range(len(top_idx)), best_model.feature_importances_[top_idx], color='teal')
        axes[2].set_yticks(range(len(top_idx)))
        axes[2].set_yticklabels([feature_names[i] for i in top_idx])
    axes[2].set_title(f'Top Features ({best_name})')
    axes[2].set_xlabel('Importance/Coefficient')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'nlp_recall_detection.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    nlp_results['multi_class'] = {'accuracy': multi_acc, 'weighted_f1': multi_f1}
    results['nlp'] = nlp_results
    return nlp_results


# ============================================================
# Module 3: Microbial Growth Prediction (Baranyi Model)
# ============================================================
def baranyi_model(t, y0, ymax, mu_max, lag):
    """Baranyi and Roberts bacterial growth model."""
    A = t + (1/mu_max) * np.log(np.exp(-mu_max * t) + np.exp(-mu_max * lag) - np.exp(-mu_max * (t + lag)))
    y = y0 + mu_max * A - np.log(1 + (np.exp(mu_max * A) - 1) / np.exp(ymax - y0))
    return y


def run_microbial_growth():
    print("\n" + "=" * 60)
    print("Module 3: Microbial Growth Prediction (Baranyi Model)")
    print("=" * 60)
    
    # Generate growth curves at different temperatures
    temperatures = [5, 10, 15, 20, 25, 30, 37]
    t = np.linspace(0, 48, 200)
    
    # Baranyi model parameters vary with temperature (based on Combase-like data)
    def get_params(temp):
        mu_max = 0.01 * np.exp(0.08 * temp)  # Ratkowsky-like relationship
        lag = max(1, 20 - 0.5 * temp + np.random.normal(0, 1))
        y0 = 2.0 + np.random.normal(0, 0.2)
        ymax = 9.0 + np.random.normal(0, 0.3)
        return y0, ymax, mu_max, lag
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot growth curves at different temperatures
    growth_data = {}
    for temp in temperatures:
        y0, ymax, mu_max, lag = get_params(temp)
        y = baranyi_model(t, y0, ymax, mu_max, lag)
        axes[0, 0].plot(t, y, label=f'{temp}°C (μ={mu_max:.3f})')
        growth_data[temp] = {'y0': y0, 'ymax': ymax, 'mu_max': mu_max, 'lag': lag}
    
    axes[0, 0].set_title('Baranyi Growth Model at Various Temperatures')
    axes[0, 0].set_xlabel('Time (hours)')
    axes[0, 0].set_ylabel('Log CFU/g')
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True, alpha=0.3)
    
    # Fit Baranyi model to noisy data (simulating real lab data)
    temp_fit = 25
    y0_true, ymax_true, mu_true, lag_true = 2.0, 9.0, 0.5, 3.0
    t_data = np.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30, 36, 42, 48])
    y_true = baranyi_model(t_data, y0_true, ymax_true, mu_true, lag_true)
    y_noisy = y_true + np.random.normal(0, 0.2, len(t_data))
    
    try:
        popt, pcov = curve_fit(baranyi_model, t_data, y_noisy, p0=[2.0, 9.0, 0.3, 2.0],
                               bounds=([0, 5, 0.01, 0.1], [5, 12, 2.0, 20]))
        y_fitted = baranyi_model(t, *popt)
        
        axes[0, 1].scatter(t_data, y_noisy, c='red', s=40, zorder=5, label='Observed')
        axes[0, 1].plot(t, y_fitted, 'b-', linewidth=2, label='Baranyi Fit')
        axes[0, 1].plot(t, baranyi_model(t, y0_true, ymax_true, mu_true, lag_true),
                        'g--', alpha=0.5, label='True Model')
        axes[0, 1].set_title(f'Baranyi Model Fitting (T={temp_fit}°C)')
        axes[0, 1].set_xlabel('Time (hours)')
        axes[0, 1].set_ylabel('Log CFU/g')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        fit_rmse = np.sqrt(mean_squared_error(y_noisy, baranyi_model(t_data, *popt)))
        fit_r2 = r2_score(y_noisy, baranyi_model(t_data, *popt))
        print(f"  Baranyi fit: RMSE={fit_rmse:.4f}, R²={fit_r2:.4f}")
        print(f"  Estimated params: y0={popt[0]:.2f}, ymax={popt[1]:.2f}, μ_max={popt[2]:.3f}, lag={popt[3]:.2f}")
    except Exception as e:
        print(f"  Fitting error: {e}")
        fit_rmse, fit_r2 = 0, 0
        popt = [2.0, 9.0, 0.5, 3.0]
    
    # ML-enhanced prediction: predict growth rate from environmental conditions
    n_ml = 500
    temps_ml = np.random.uniform(4, 40, n_ml)
    ph_ml = np.random.uniform(4.0, 8.0, n_ml)
    aw_ml = np.random.uniform(0.90, 1.0, n_ml)
    
    # Growth rate follows secondary model (simplified)
    mu_ml = 0.01 * np.exp(0.08 * temps_ml) * (1 - np.exp(-0.5 * (ph_ml - 3.5))) * (aw_ml - 0.88) / 0.12
    mu_ml = np.clip(mu_ml, 0, 3) + np.random.normal(0, 0.01, n_ml)
    mu_ml = np.clip(mu_ml, 0, None)
    
    X_ml = np.column_stack([temps_ml, ph_ml, aw_ml])
    X_train_ml, X_test_ml, y_train_ml, y_test_ml = train_test_split(X_ml, mu_ml, test_size=0.2, random_state=42)
    
    rf_growth = RandomForestRegressor(n_estimators=200, random_state=42)
    rf_growth.fit(X_train_ml, y_train_ml)
    y_pred_ml = rf_growth.predict(X_test_ml)
    
    ml_rmse = np.sqrt(mean_squared_error(y_test_ml, y_pred_ml))
    ml_r2 = r2_score(y_test_ml, y_pred_ml)
    ml_mae = mean_absolute_error(y_test_ml, y_pred_ml)
    print(f"  ML growth rate prediction: RMSE={ml_rmse:.4f}, R²={ml_r2:.4f}, MAE={ml_mae:.4f}")
    
    axes[1, 0].scatter(y_test_ml, y_pred_ml, alpha=0.5, s=15, c='purple')
    axes[1, 0].plot([0, max(y_test_ml)], [0, max(y_test_ml)], 'r--')
    axes[1, 0].set_title(f'ML Growth Rate Prediction (R²={ml_r2:.3f})')
    axes[1, 0].set_xlabel('Actual μ_max (h⁻¹)')
    axes[1, 0].set_ylabel('Predicted μ_max (h⁻¹)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Feature importance for growth rate prediction
    feat_names = ['Temperature', 'pH', 'Water Activity']
    imp = rf_growth.feature_importances_
    axes[1, 1].bar(feat_names, imp, color=['#e74c3c', '#3498db', '#2ecc71'], edgecolor='black')
    axes[1, 1].set_title('Feature Importance for Growth Rate')
    axes[1, 1].set_ylabel('Importance')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'microbial_growth.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    results['microbial'] = {
        'baranyi_fit': {'rmse': fit_rmse, 'r2': fit_r2, 'params': list(popt)},
        'ml_prediction': {'rmse': ml_rmse, 'r2': ml_r2, 'mae': ml_mae}
    }
    return results['microbial']


# ============================================================
# Module 4: HACCP Risk Scoring Automation
# ============================================================
def run_haccp_scoring():
    print("\n" + "=" * 60)
    print("Module 4: HACCP Risk Scoring Automation")
    print("=" * 60)
    
    np.random.seed(42)
    n_records = 1500
    
    # CCP monitoring data
    ccp_types = ['receiving', 'cooking', 'cooling', 'hot_holding', 'cold_storage', 'packaging', 'shipping']
    
    data = []
    for _ in range(n_records):
        ccp = np.random.choice(ccp_types)
        temp_deviation = np.random.exponential(2)
        time_deviation = np.random.exponential(5)
        humidity_deviation = np.random.exponential(3)
        equipment_age = np.random.uniform(0, 15)
        staff_training_score = np.random.uniform(50, 100)
        previous_violations = np.random.poisson(1)
        inspection_frequency = np.random.choice([1, 2, 4, 12])
        
        # Risk score calculation (simulated ground truth)
        risk = (0.25 * temp_deviation + 0.15 * time_deviation + 0.1 * humidity_deviation +
                0.15 * equipment_age / 15 * 10 + 0.1 * (100 - staff_training_score) / 10 +
                0.15 * previous_violations * 3 + 0.1 * (12 - inspection_frequency) / 3)
        risk = np.clip(risk + np.random.normal(0, 1), 0, 20)
        risk_category = 'low' if risk < 5 else ('medium' if risk < 10 else 'high')
        
        data.append({
            'ccp_type': ccp, 'temp_deviation': temp_deviation,
            'time_deviation': time_deviation, 'humidity_deviation': humidity_deviation,
            'equipment_age': equipment_age, 'staff_training_score': staff_training_score,
            'previous_violations': previous_violations, 'inspection_frequency': inspection_frequency,
            'risk_score': risk, 'risk_category': risk_category
        })
    
    df = pd.DataFrame(data)
    
    # Encode CCP type
    le_ccp = LabelEncoder()
    df['ccp_encoded'] = le_ccp.fit_transform(df['ccp_type'])
    
    features = ['ccp_encoded', 'temp_deviation', 'time_deviation', 'humidity_deviation',
                'equipment_age', 'staff_training_score', 'previous_violations', 'inspection_frequency']
    
    # Regression: predict risk score
    X = df[features].values
    y_score = df['risk_score'].values
    X_train, X_test, y_train, y_test = train_test_split(X, y_score, test_size=0.2, random_state=42)
    
    gbr = GradientBoostingRegressor(n_estimators=200, max_depth=5, random_state=42)
    gbr.fit(X_train, y_train)
    y_pred_score = gbr.predict(X_test)
    
    score_rmse = np.sqrt(mean_squared_error(y_test, y_pred_score))
    score_r2 = r2_score(y_test, y_pred_score)
    score_mae = mean_absolute_error(y_test, y_pred_score)
    print(f"  Risk Score Regression: RMSE={score_rmse:.4f}, R²={score_r2:.4f}, MAE={score_mae:.4f}")
    
    # Classification: predict risk category
    le_risk = LabelEncoder()
    y_cat = le_risk.fit_transform(df['risk_category'])
    X_train_c, X_test_c, y_train_c, y_test_c = train_test_split(X, y_cat, test_size=0.2, random_state=42)
    
    rf_cat = RandomForestClassifier(n_estimators=200, random_state=42)
    rf_cat.fit(X_train_c, y_train_c)
    y_pred_cat = rf_cat.predict(X_test_c)
    cat_acc = accuracy_score(y_test_c, y_pred_cat)
    cat_f1 = f1_score(y_test_c, y_pred_cat, average='weighted')
    print(f"  Risk Category Classification: Acc={cat_acc:.4f}, Weighted F1={cat_f1:.4f}")
    
    # Plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Actual vs Predicted risk score
    axes[0, 0].scatter(y_test, y_pred_score, alpha=0.4, s=15, c='darkorange')
    axes[0, 0].plot([0, 20], [0, 20], 'r--')
    axes[0, 0].set_title(f'HACCP Risk Score: Actual vs Predicted (R²={score_r2:.3f})')
    axes[0, 0].set_xlabel('Actual Risk Score')
    axes[0, 0].set_ylabel('Predicted Risk Score')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Risk distribution by CCP
    df.boxplot(column='risk_score', by='ccp_type', ax=axes[0, 1], rot=45)
    axes[0, 1].set_title('Risk Score Distribution by CCP Type')
    axes[0, 1].set_xlabel('CCP Type')
    axes[0, 1].set_ylabel('Risk Score')
    plt.suptitle('')
    
    # Feature importance
    feat_imp = pd.Series(gbr.feature_importances_, index=features).sort_values(ascending=True)
    feat_imp.plot(kind='barh', ax=axes[1, 0], color='forestgreen')
    axes[1, 0].set_title('Feature Importance for Risk Scoring')
    axes[1, 0].set_xlabel('Importance')
    
    # Risk category confusion matrix
    cm = confusion_matrix(y_test_c, y_pred_cat)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Oranges', ax=axes[1, 1],
                xticklabels=le_risk.classes_, yticklabels=le_risk.classes_)
    axes[1, 1].set_title('Risk Category Classification')
    axes[1, 1].set_xlabel('Predicted')
    axes[1, 1].set_ylabel('True')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'haccp_scoring.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    results['haccp'] = {
        'regression': {'rmse': score_rmse, 'r2': score_r2, 'mae': score_mae},
        'classification': {'accuracy': cat_acc, 'weighted_f1': cat_f1}
    }
    return results['haccp']


# ============================================================
# Module 5: Blockchain-integrated Traceability Simulation
# ============================================================
def run_blockchain_traceability():
    print("\n" + "=" * 60)
    print("Module 5: Blockchain-integrated Traceability Simulation")
    print("=" * 60)
    
    import hashlib
    import time as time_mod
    
    class Block:
        def __init__(self, index, data, previous_hash, timestamp=None):
            self.index = index
            self.timestamp = timestamp or time_mod.time()
            self.data = data
            self.previous_hash = previous_hash
            self.hash = self.compute_hash()
        
        def compute_hash(self):
            block_str = json.dumps({
                'index': self.index, 'timestamp': self.timestamp,
                'data': self.data, 'previous_hash': self.previous_hash
            }, sort_keys=True)
            return hashlib.sha256(block_str.encode()).hexdigest()
    
    class FoodChain:
        def __init__(self):
            self.chain = [self._genesis_block()]
        
        def _genesis_block(self):
            return Block(0, {'event': 'Genesis Block'}, '0')
        
        def add_block(self, data):
            prev = self.chain[-1]
            new_block = Block(len(self.chain), data, prev.hash)
            self.chain.append(new_block)
            return new_block
        
        def verify_chain(self):
            for i in range(1, len(self.chain)):
                current = self.chain[i]
                previous = self.chain[i - 1]
                if current.hash != current.compute_hash():
                    return False, i
                if current.previous_hash != previous.hash:
                    return False, i
            return True, -1
    
    # Simulate a chicken supply chain
    chain = FoodChain()
    
    supply_chain_events = [
        {'stage': 'Farm', 'event': 'Harvest', 'product': 'Chicken Batch #2024-001',
         'temperature': 4.2, 'location': 'Farm A, Iowa', 'certifications': ['USDA Organic']},
        {'stage': 'Processing', 'event': 'Slaughter & Processing',
         'temperature': 2.1, 'location': 'Plant B, Illinois', 'haccp_check': 'PASS'},
        {'stage': 'Cold Storage', 'event': 'Storage Entry',
         'temperature': -18.0, 'location': 'Warehouse C, Indiana', 'duration_hours': 48},
        {'stage': 'Transport', 'event': 'Refrigerated Transport',
         'temperature': -15.5, 'location': 'Route C→D', 'vehicle_id': 'TRK-4521'},
        {'stage': 'Distribution', 'event': 'Distribution Center',
         'temperature': -17.8, 'location': 'DC D, Ohio', 'quality_check': 'PASS'},
        {'stage': 'Retail', 'event': 'Store Delivery',
         'temperature': 3.5, 'location': 'Store E, Detroit', 'shelf_life_remaining': '5 days'},
    ]
    
    for event in supply_chain_events:
        block = chain.add_block(event)
    
    is_valid, _ = chain.verify_chain()
    print(f"  Chain length: {len(chain.chain)} blocks")
    print(f"  Chain integrity: {'VALID' if is_valid else 'INVALID'}")
    
    # Simulate anomaly detection in supply chain
    np.random.seed(42)
    n_shipments = 500
    
    normal_temps = np.random.normal(-18, 1.5, n_shipments)
    transit_times = np.random.exponential(12, n_shipments) + 4
    
    # Introduce anomalies (10%)
    anomaly_idx = np.random.choice(n_shipments, 50, replace=False)
    normal_temps[anomaly_idx] += np.random.uniform(10, 25, 50)
    transit_times[anomaly_idx] *= np.random.uniform(2, 4, 50)
    
    labels = np.zeros(n_shipments)
    labels[anomaly_idx] = 1
    
    X_bc = np.column_stack([normal_temps, transit_times])
    X_train_bc, X_test_bc, y_train_bc, y_test_bc = train_test_split(X_bc, labels, test_size=0.2, random_state=42)
    
    rf_anom = RandomForestClassifier(n_estimators=100, random_state=42)
    rf_anom.fit(X_train_bc, y_train_bc)
    y_pred_bc = rf_anom.predict(X_test_bc)
    y_prob_bc = rf_anom.predict_proba(X_test_bc)[:, 1]
    
    bc_acc = accuracy_score(y_test_bc, y_pred_bc)
    bc_f1 = f1_score(y_test_bc, y_pred_bc)
    bc_auc = roc_auc_score(y_test_bc, y_prob_bc)
    print(f"  Anomaly Detection: Acc={bc_acc:.4f}, F1={bc_f1:.4f}, AUC={bc_auc:.4f}")
    
    # Plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Supply chain temperature tracking
    stages = [e['stage'] for e in supply_chain_events]
    temps = [e['temperature'] for e in supply_chain_events]
    colors = ['green' if t <= 4 else ('blue' if t <= -10 else 'red') for t in temps]
    axes[0].bar(stages, temps, color=colors, edgecolor='black')
    axes[0].axhline(y=4, color='red', linestyle='--', alpha=0.5, label='Max Safe Temp (4°C)')
    axes[0].set_title('Supply Chain Temperature Tracking')
    axes[0].set_ylabel('Temperature (°C)')
    axes[0].set_xticklabels(stages, rotation=30, ha='right')
    axes[0].legend()
    
    # Anomaly scatter
    normal_mask = y_test_bc == 0
    axes[1].scatter(X_test_bc[normal_mask, 0], X_test_bc[normal_mask, 1],
                    alpha=0.5, s=20, c='blue', label='Normal')
    axes[1].scatter(X_test_bc[~normal_mask, 0], X_test_bc[~normal_mask, 1],
                    alpha=0.7, s=30, c='red', label='Anomaly')
    axes[1].set_title('Supply Chain Anomaly Detection')
    axes[1].set_xlabel('Temperature (°C)')
    axes[1].set_ylabel('Transit Time (hours)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Blockchain structure visualization
    block_indices = list(range(len(chain.chain)))
    block_hashes = [b.hash[:8] for b in chain.chain]
    axes[2].barh(block_indices, [1] * len(block_indices), color='steelblue', edgecolor='black')
    for i, (idx, h) in enumerate(zip(block_indices, block_hashes)):
        stage = 'Genesis' if i == 0 else supply_chain_events[i-1]['stage']
        axes[2].text(0.5, idx, f'{stage}\n{h}...', ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    axes[2].set_title('Blockchain Structure')
    axes[2].set_ylabel('Block Index')
    axes[2].set_xlabel('')
    axes[2].set_xticks([])
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'blockchain_traceability.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    results['blockchain'] = {
        'chain_valid': is_valid, 'chain_length': len(chain.chain),
        'anomaly_detection': {'accuracy': bc_acc, 'f1': bc_f1, 'auc': bc_auc}
    }
    return results['blockchain']


# ============================================================
# Module 6: Salmonella in Chicken - Case Study
# ============================================================
def run_salmonella_case_study():
    print("\n" + "=" * 60)
    print("Module 6: Salmonella in Chicken - Case Study")
    print("=" * 60)
    
    np.random.seed(42)
    n_samples = 1200
    
    # Generate realistic chicken processing data
    processing_temp = np.random.normal(4, 2, n_samples)
    cooking_temp = np.random.normal(74, 5, n_samples)
    storage_duration = np.random.exponential(24, n_samples)
    storage_temp = np.random.normal(3, 2, n_samples)
    humidity = np.random.normal(65, 10, n_samples)
    ambient_temp = np.random.normal(22, 8, n_samples)
    season = np.random.choice([0, 1, 2, 3], n_samples)  # 0=winter, 1=spring, 2=summer, 3=fall
    supplier_rating = np.random.uniform(60, 100, n_samples)
    chlorine_wash = np.random.choice([0, 1], n_samples, p=[0.3, 0.7])
    
    # Salmonella contamination probability
    logit = (0.5 + 0.25 * processing_temp - 0.02 * cooking_temp + 0.02 * storage_duration
             + 0.2 * storage_temp + 0.005 * humidity + 0.03 * ambient_temp
             + 0.6 * (season == 2).astype(float)  # summer risk
             - 0.02 * supplier_rating - 0.6 * chlorine_wash
             + np.random.normal(0, 0.8, n_samples))
    prob = 1 / (1 + np.exp(-logit))
    contaminated = (np.random.random(n_samples) < prob).astype(int)
    
    print(f"  Contamination rate: {contaminated.mean():.2%}")
    
    df = pd.DataFrame({
        'processing_temp': processing_temp, 'cooking_temp': cooking_temp,
        'storage_duration': storage_duration, 'storage_temp': storage_temp,
        'humidity': humidity, 'ambient_temp': ambient_temp, 'season': season,
        'supplier_rating': supplier_rating, 'chlorine_wash': chlorine_wash,
        'contaminated': contaminated
    })
    
    features = ['processing_temp', 'cooking_temp', 'storage_duration', 'storage_temp',
                'humidity', 'ambient_temp', 'season', 'supplier_rating', 'chlorine_wash']
    
    X = df[features].values
    y = df['contaminated'].values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)
    
    models = {
        'Logistic Regression': LogisticRegression(max_iter=1000),
        'Random Forest': RandomForestClassifier(n_estimators=300, max_depth=12, random_state=42),
        'Gradient Boosting': GradientBoostingRegressor(n_estimators=300, max_depth=5, random_state=42),
        'MLP': MLPClassifier(hidden_layer_sizes=(128, 64, 32), max_iter=500, random_state=42)
    }
    
    salmonella_results = {}
    for name, model in models.items():
        model.fit(X_train, y_train)
        if name == 'Gradient Boosting':
            y_pred_prob = np.clip(model.predict(X_test), 0, 1)
            y_pred = (y_pred_prob > 0.5).astype(int)
        else:
            y_pred = model.predict(X_test)
            y_pred_prob = model.predict_proba(X_test)[:, 1]
        
        acc = accuracy_score(y_test, y_pred)
        prec = precision_score(y_test, y_pred, zero_division=0)
        rec = recall_score(y_test, y_pred, zero_division=0)
        f1 = f1_score(y_test, y_pred, zero_division=0)
        try:
            auc = roc_auc_score(y_test, y_pred_prob)
        except:
            auc = 0
        
        salmonella_results[name] = {'accuracy': acc, 'precision': prec, 'recall': rec, 'f1': f1, 'auc': auc}
        print(f"  {name}: Acc={acc:.4f}, F1={f1:.4f}, AUC={auc:.4f}")
    
    # Cross-validation for best model
    rf = models['Random Forest']
    cv_scores = cross_val_score(rf, X_scaled, y, cv=5, scoring='f1')
    print(f"\n  Random Forest 5-fold CV F1: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    
    # Feature importance
    rf_importances = rf.feature_importances_
    
    # Plots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    
    # Model comparison
    sal_df = pd.DataFrame(salmonella_results).T
    sal_df[['accuracy', 'precision', 'recall', 'f1']].plot(kind='bar', ax=axes[0, 0])
    axes[0, 0].set_title('Salmonella Prediction: Model Comparison')
    axes[0, 0].set_ylabel('Score')
    axes[0, 0].set_xticklabels(axes[0, 0].get_xticklabels(), rotation=30, ha='right')
    axes[0, 0].set_ylim(0, 1)
    axes[0, 0].legend(loc='lower right', fontsize=8)
    
    # ROC curves
    for name, model in models.items():
        if name == 'Gradient Boosting':
            y_score = np.clip(model.predict(X_test), 0, 1)
        else:
            y_score = model.predict_proba(X_test)[:, 1]
        fpr, tpr, _ = roc_curve(y_test, y_score)
        axes[0, 1].plot(fpr, tpr, label=f"{name} (AUC={salmonella_results[name]['auc']:.3f})")
    axes[0, 1].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[0, 1].set_title('ROC Curves - Salmonella Detection')
    axes[0, 1].set_xlabel('FPR')
    axes[0, 1].set_ylabel('TPR')
    axes[0, 1].legend(fontsize=8)
    
    # Feature importance
    feat_imp = pd.Series(rf_importances, index=features).sort_values(ascending=True)
    feat_imp.plot(kind='barh', ax=axes[0, 2], color='crimson')
    axes[0, 2].set_title('Feature Importance (Random Forest)')
    axes[0, 2].set_xlabel('Importance')
    
    # Seasonal contamination rate
    season_names = ['Winter', 'Spring', 'Summer', 'Fall']
    season_rates = [df[df['season'] == i]['contaminated'].mean() for i in range(4)]
    axes[1, 0].bar(season_names, season_rates, color=['#3498db', '#2ecc71', '#e74c3c', '#f39c12'], edgecolor='black')
    axes[1, 0].set_title('Seasonal Contamination Rate')
    axes[1, 0].set_ylabel('Contamination Rate')
    
    # Temperature effect
    temp_bins = pd.cut(df['processing_temp'], bins=5)
    temp_contamination = df.groupby(temp_bins)['contaminated'].mean()
    axes[1, 1].bar(range(len(temp_contamination)), temp_contamination.values, color='salmon', edgecolor='black')
    axes[1, 1].set_xticks(range(len(temp_contamination)))
    axes[1, 1].set_xticklabels([f'{x.left:.0f}-{x.right:.0f}' for x in temp_contamination.index], rotation=30)
    axes[1, 1].set_title('Processing Temperature vs Contamination Rate')
    axes[1, 1].set_xlabel('Processing Temperature (°C)')
    axes[1, 1].set_ylabel('Contamination Rate')
    
    # Confusion matrix for best model
    y_pred_best = rf.predict(X_test)
    cm = confusion_matrix(y_test, y_pred_best)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Reds', ax=axes[1, 2],
                xticklabels=['Safe', 'Contaminated'], yticklabels=['Safe', 'Contaminated'])
    axes[1, 2].set_title('Confusion Matrix (Random Forest)')
    axes[1, 2].set_xlabel('Predicted')
    axes[1, 2].set_ylabel('True')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'salmonella_case_study.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    salmonella_results['cv_f1_mean'] = cv_scores.mean()
    salmonella_results['cv_f1_std'] = cv_scores.std()
    salmonella_results['contamination_rate'] = contaminated.mean()
    results['salmonella'] = salmonella_results
    return salmonella_results


# ============================================================
# Module 7: Integrated Risk Monitoring Dashboard
# ============================================================
def run_integrated_dashboard():
    print("\n" + "=" * 60)
    print("Module 7: Integrated Risk Monitoring System")
    print("=" * 60)
    
    # Create integrated summary visualization
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. System architecture overview (text-based)
    axes[0, 0].text(0.5, 0.95, 'Integrated Food Safety AI System', ha='center', va='top',
                    fontsize=14, fontweight='bold', transform=axes[0, 0].transAxes)
    
    components = [
        ('Spatiotemporal\nPrediction', 0.2, 0.7, '#3498db'),
        ('NLP Recall\nDetection', 0.5, 0.7, '#e74c3c'),
        ('Microbial\nGrowth Model', 0.8, 0.7, '#2ecc71'),
        ('HACCP\nScoring', 0.2, 0.3, '#f39c12'),
        ('Blockchain\nTraceability', 0.5, 0.3, '#9b59b6'),
        ('Salmonella\nCase Study', 0.8, 0.3, '#1abc9c'),
    ]
    
    for label, x, y, color in components:
        circle = plt.Circle((x, y), 0.12, color=color, alpha=0.7, transform=axes[0, 0].transAxes)
        axes[0, 0].add_patch(circle)
        axes[0, 0].text(x, y, label, ha='center', va='center', fontsize=8,
                        fontweight='bold', color='white', transform=axes[0, 0].transAxes)
    
    # Central hub
    center = plt.Circle((0.5, 0.5), 0.08, color='#34495e', alpha=0.9, transform=axes[0, 0].transAxes)
    axes[0, 0].add_patch(center)
    axes[0, 0].text(0.5, 0.5, 'Risk\nEngine', ha='center', va='center', fontsize=8,
                    fontweight='bold', color='white', transform=axes[0, 0].transAxes)
    
    axes[0, 0].set_xlim(0, 1)
    axes[0, 0].set_ylim(0, 1)
    axes[0, 0].set_xticks([])
    axes[0, 0].set_yticks([])
    axes[0, 0].set_title('System Architecture')
    
    # 2. Overall performance summary
    module_names = ['Spatiotemporal', 'NLP Detection', 'Microbial ML', 'HACCP Scoring', 'Blockchain', 'Salmonella']
    
    # Get best metrics from each module
    best_metrics = []
    if 'spatiotemporal' in results:
        best_f1 = max(v.get('f1', 0) for v in results['spatiotemporal'].values() if isinstance(v, dict))
        best_metrics.append(best_f1)
    else:
        best_metrics.append(0)
    
    if 'nlp' in results:
        nlp_vals = [v.get('f1', 0) for k, v in results['nlp'].items() if isinstance(v, dict) and 'f1' in v]
        best_metrics.append(max(nlp_vals) if nlp_vals else 0)
    else:
        best_metrics.append(0)
    
    if 'microbial' in results:
        best_metrics.append(results['microbial']['ml_prediction']['r2'])
    else:
        best_metrics.append(0)
    
    if 'haccp' in results:
        best_metrics.append(results['haccp']['regression']['r2'])
    else:
        best_metrics.append(0)
    
    if 'blockchain' in results:
        best_metrics.append(results['blockchain']['anomaly_detection']['f1'])
    else:
        best_metrics.append(0)
    
    if 'salmonella' in results:
        sal_f1s = [v.get('f1', 0) for v in results['salmonella'].values() if isinstance(v, dict)]
        best_metrics.append(max(sal_f1s) if sal_f1s else 0)
    else:
        best_metrics.append(0)
    
    colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
    axes[0, 1].bar(module_names, best_metrics, color=colors, edgecolor='black')
    axes[0, 1].set_title('Best Performance by Module (F1/R²)')
    axes[0, 1].set_ylabel('Score')
    axes[0, 1].set_xticklabels(module_names, rotation=30, ha='right')
    axes[0, 1].set_ylim(0, 1)
    axes[0, 1].axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='Target (0.8)')
    axes[0, 1].legend()
    
    # 3. Time series risk simulation
    np.random.seed(42)
    days = 365
    t_days = np.arange(days)
    
    # Simulated daily risk scores
    base_risk = 30 + 15 * np.sin(2 * np.pi * t_days / 365)
    noise = np.random.normal(0, 5, days)
    # Add some spikes (outbreaks)
    spikes = np.zeros(days)
    spike_days = [45, 120, 200, 280, 330]
    for sd in spike_days:
        spikes[max(0, sd-2):min(days, sd+5)] += np.random.uniform(20, 40)
    
    risk_timeline = np.clip(base_risk + noise + spikes, 0, 100)
    
    axes[1, 0].plot(t_days, risk_timeline, color='#e74c3c', alpha=0.8, linewidth=1)
    axes[1, 0].fill_between(t_days, risk_timeline, alpha=0.2, color='red')
    axes[1, 0].axhline(y=70, color='red', linestyle='--', label='High Risk Threshold')
    axes[1, 0].axhline(y=40, color='orange', linestyle='--', label='Medium Risk Threshold')
    for sd in spike_days:
        axes[1, 0].axvline(x=sd, color='gray', linestyle=':', alpha=0.3)
    axes[1, 0].set_title('Integrated Risk Score Timeline (1 Year)')
    axes[1, 0].set_xlabel('Day of Year')
    axes[1, 0].set_ylabel('Risk Score')
    axes[1, 0].legend(fontsize=8)
    
    # 4. Risk heatmap by category and month
    categories = ['Salmonella', 'Listeria', 'E. coli', 'Chemical', 'Allergen']
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
    risk_matrix = np.random.uniform(10, 90, (len(categories), 12))
    # Summer peaks for biological hazards
    risk_matrix[0, 5:9] += 30  # Salmonella
    risk_matrix[2, 5:9] += 25  # E. coli
    risk_matrix = np.clip(risk_matrix, 0, 100)
    
    sns.heatmap(risk_matrix, ax=axes[1, 1], cmap='YlOrRd', annot=True, fmt='.0f',
                xticklabels=months, yticklabels=categories, vmin=0, vmax=100)
    axes[1, 1].set_title('Monthly Risk Heatmap by Hazard Category')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'integrated_dashboard.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    print("  Integrated dashboard generated successfully.")
    return True


# ============================================================
# Main execution
# ============================================================
if __name__ == '__main__':
    print("Food Supply Chain Safety Risk Prediction AI System")
    print("=" * 60)
    
    r1 = run_spatiotemporal_model()
    r2 = run_nlp_recall_detection()
    r3 = run_microbial_growth()
    r4 = run_haccp_scoring()
    r5 = run_blockchain_traceability()
    r6 = run_salmonella_case_study()
    r7 = run_integrated_dashboard()
    
    # Save all results
    def convert_to_serializable(obj):
        if isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert_to_serializable(v) for k, v in obj.items()}
        return obj
    
    serializable_results = convert_to_serializable(results)
    with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results.json'), 'w') as f:
        json.dump(serializable_results, f, indent=2)
    
    print("\n" + "=" * 60)
    print("All experiments completed successfully!")
    print(f"Results saved to results.json")
    print(f"Figures saved to {FIGURES_DIR}/")
    print("=" * 60)
