# HEA ML Framework - Complete (fixed)
import sys, random, warnings, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from scipy.stats import norm
import xgboost as xgb
warnings.filterwarnings('ignore')

np.random.seed(42); random.seed(42)
FIGURES_DIR = '/app/projects/50b9b3f2-279d-4427-8ebf-e99b1e2beb9c/workspace/figures'
DATA_DIR    = '/app/projects/50b9b3f2-279d-4427-8ebf-e99b1e2beb9c/workspace/data/raw'

# ─── Element database ────────────────────────────────────────
ELEMENT_DATA = {
    'Cr': {'r':1.28,'chi':1.66,'Tm':2180,'G':115,'VEC':6},
    'Mn': {'r':1.27,'chi':1.55,'Tm':1519,'G': 80,'VEC':7},
    'Fe': {'r':1.26,'chi':1.83,'Tm':1811,'G': 82,'VEC':8},
    'Co': {'r':1.25,'chi':1.88,'Tm':1768,'G': 75,'VEC':9},
    'Ni': {'r':1.24,'chi':1.91,'Tm':1728,'G': 76,'VEC':10},
}
ELEMENTS = ['Cr','Mn','Fe','Co','Ni']
Omega_m = {('Cr','Mn'):-4.,('Cr','Fe'):-1.,('Cr','Co'):-4.,('Cr','Ni'):-7.,
            ('Mn','Fe'):0.,('Mn','Co'):-5.,('Mn','Ni'):-8.,
            ('Fe','Co'):-1.,('Fe','Ni'):-2.,('Co','Ni'):0.}

def calc_desc(x):
    """x: array of 5 fractions (Cr,Mn,Fe,Co,Ni)"""
    x = np.array(x, dtype=float); x /= x.sum()
    props = ELEMENT_DATA; el = ELEMENTS
    r_b  = sum(x[i]*props[e]['r']   for i,e in enumerate(el))
    chi_b= sum(x[i]*props[e]['chi'] for i,e in enumerate(el))
    Tm_b = sum(x[i]*props[e]['Tm']  for i,e in enumerate(el))
    G_b  = sum(x[i]*props[e]['G']   for i,e in enumerate(el))
    VEC_b= sum(x[i]*props[e]['VEC'] for i,e in enumerate(el))
    delt_r= np.sqrt(sum(x[i]*(1-props[e]['r']/r_b)**2 for i,e in enumerate(el)))*100
    delt_c= np.sqrt(sum(x[i]*(props[e]['chi']-chi_b)**2 for i,e in enumerate(el)))
    R=8.314; S=float(-R*sum(x[i]*np.log(x[i]+1e-12) for i in range(5)))
    H=sum(4*Omega_m.get((el[i],el[j]),Omega_m.get((el[j],el[i]),0.))*x[i]*x[j]
          for i in range(5) for j in range(i+1,5))
    Om= S*Tm_b/(abs(H)+1e-3)
    Ga= delt_r**2/(delt_c+1e-6)
    return np.array([x[0],x[1],x[2],x[3],x[4],
                     VEC_b,delt_r,delt_c,S,H,Om,Ga,Tm_b,G_b], dtype=float)

# ─── Cell 1: Generate dataset ────────────────────────────────
def gen_dataset(n=300, seed=42):
    rng = np.random.RandomState(seed)
    comps = rng.dirichlet(np.ones(5), n)
    rows=[]
    for c in comps:
        d=calc_desc(c)
        VEC=d[5]; dr=d[6]; Gb=d[13]; Hm=d[9]
        sy = 150 + 80*dr + 0.8*Gb + 10*abs(Hm) - 15*abs(VEC-8.0) + rng.normal(0,15)
        sy = max(100,sy)
        el = 40 - 4*dr - 0.05*Gb + 5*max(0,3-abs(VEC-8.5)) + 2*c[4] - 3*c[1] + rng.normal(0,3)
        el = float(np.clip(el,2,70))
        cr = 4.0 + 20*c[0] + 8*c[4] - 5*c[1] - 0.2*d[7] + rng.normal(0,0.5)
        cr = float(np.clip(cr,0,10))
        rows.append(list(d)+[sy,el,cr])
    cols = ['x_Cr','x_Mn','x_Fe','x_Co','x_Ni','VEC','delta_r','delta_chi',
            'S_mix','H_mix','Omega','Gamma','Tm','G_bar',
            'yield_strength','elongation','corrosion_resistance']
    return pd.DataFrame(rows, columns=cols)

df = gen_dataset(300)
print(f"[Cell 1] Dataset: {df.shape}")
df.to_csv(f'{DATA_DIR}/hea_dataset.csv', index=False)

feat_cols = ['x_Cr','x_Mn','x_Fe','x_Co','x_Ni','VEC','delta_r','delta_chi',
             'S_mix','H_mix','Omega','Gamma','Tm','G_bar']
tgt_cols  = ['yield_strength','elongation','corrosion_resistance']
X = df[feat_cols].values
Y = df[tgt_cols].values

# ─── Cell 2: 5-fold CV ───────────────────────────────────────
kf = KFold(n_splits=5, shuffle=True, random_state=42)
cv_results={}
for tgt in tgt_cols:
    y=df[tgt].values; cv_results[tgt]={}
    for name,mdl in [('RF', RandomForestRegressor(n_estimators=200,max_depth=8,random_state=42,n_jobs=-1)),
                     ('XGB',xgb.XGBRegressor(n_estimators=200,max_depth=5,learning_rate=0.05,
                                             subsample=0.8,random_state=42,verbosity=0))]:
        pipe=Pipeline([('sc',StandardScaler()),('m',mdl)])
        r2s = cross_val_score(pipe, X, y, cv=kf, scoring='r2')
        rmse_s = np.sqrt(-cross_val_score(pipe, X, y, cv=kf, scoring='neg_mean_squared_error'))
        cv_results[tgt][name]={'R2_mean':r2s.mean(),'R2_std':r2s.std(),
                               'RMSE_mean':rmse_s.mean(),'RMSE_std':rmse_s.std()}

print("\n[Cell 2] 5-fold CV Results:")
print(f"{'Target':22s} {'Model':5s} {'R2_val':>12s} {'RMSE_val':>14s}")
for tgt in tgt_cols:
    for nm,res in cv_results[tgt].items():
        print(f"{tgt:22s} {nm:5s} {res['R2_mean']:>7.3f}±{res['R2_std']:.3f}  {res['RMSE_mean']:>8.3f}±{res['RMSE_std']:.3f}")

# ─── Cell 3: Train final models ──────────────────────────────
final_models={}
for tgt in tgt_cols:
    pipe=Pipeline([('sc',StandardScaler()),
                   ('m',RandomForestRegressor(n_estimators=300,max_depth=8,random_state=42,n_jobs=-1))])
    pipe.fit(X, df[tgt].values)
    final_models[tgt]=pipe

# ─── Cell 4: Bayesian Optimization (GP surrogate) ────────────
print("\n[Cell 4] Bayesian Optimization with GP...")

sc_bo = StandardScaler().fit(X)
X_sc  = sc_bo.transform(X)

gp_surrogates={}
for tgt in ['yield_strength','elongation']:
    y=df[tgt].values; y_n=(y-y.mean())/y.std()
    gp=GaussianProcessRegressor(kernel=Matern(nu=2.5)*1.0, alpha=1e-3,
                                normalize_y=True, n_restarts_optimizer=3, random_state=42)
    gp.fit(X_sc, y_n)
    gp_surrogates[tgt]=(gp, y.mean(), y.std())

def EI(gp, X_cand, y_best, xi=0.01):
    mu,sigma=gp.predict(X_cand, return_std=True)
    sigma=np.maximum(sigma, 1e-9)
    z=(mu-y_best-xi)/sigma
    return float(np.sum((mu-y_best-xi)*norm.cdf(z)+sigma*norm.pdf(z)))

n_init=30
obs_idx=list(range(n_init)); pool=list(range(n_init,len(df)))
bo_log=[]
best_ys_init=df['yield_strength'].values[:n_init].max()

for it in range(20):
    gp_ys,ys_mu,ys_sd = gp_surrogates['yield_strength']
    y_obs_norm=(df['yield_strength'].values[obs_idx]-ys_mu)/ys_sd
    y_best_n=y_obs_norm.max()
    batch=pool[:50]
    X_cand=X_sc[batch]
    mu_all,sig_all=gp_ys.predict(X_cand, return_std=True)
    sig_all=np.maximum(sig_all,1e-9)
    z=(mu_all-y_best_n-0.01)/sig_all
    ei_vals=(mu_all-y_best_n-0.01)*norm.cdf(z)+sig_all*norm.pdf(z)
    best_local=int(np.argmax(ei_vals))
    chosen=batch[best_local]
    obs_idx.append(chosen); pool.remove(chosen)
    cur_best=df['yield_strength'].values[obs_idx].max()
    bo_log.append({'iter':it,'n_obs':len(obs_idx),'best_ys':cur_best})

bo_df=pd.DataFrame(bo_log)
bo_improvement=bo_df['best_ys'].max()-best_ys_init
print(f"[Cell 4] Initial best: {best_ys_init:.1f} MPa")
print(f"[Cell 4] BO best:      {bo_df['best_ys'].max():.1f} MPa")
print(f"[Cell 4] Improvement:  +{bo_improvement:.1f} MPa over {20} iterations")

# ─── Cell 5: Pareto front analysis ───────────────────────────
def pareto(costs):
    n=len(costs); eff=np.ones(n,bool)
    for i in range(n):
        if eff[i]:
            dominated = np.all(costs <= costs[i], axis=1) & np.any(costs < costs[i], axis=1)
            dominated[i]=False
            eff[dominated]=False
    return eff

costs2d=np.column_stack([-df['yield_strength'].values,-df['elongation'].values])
pareto_mask=pareto(costs2d)
print(f"\n[Cell 5] Pareto front: {pareto_mask.sum()} alloys")

# ─── Cell Figures 6: 
print("\n[Cell 6] Generating figures...")

# Fig1: Correlation + PCA
fig,axes=plt.subplots(1,2,figsize=(16,6))
corr_df=df[feat_cols+tgt_cols].corr()
tgt_corr=corr_df.loc[feat_cols,tgt_cols]
sns.heatmap(tgt_corr,annot=True,fmt='.2f',cmap='RdBu_r',center=0,
            vmin=-1,vmax=1,ax=axes[0],linewidths=0.5)
axes[0].set_title('Descriptor–Property Correlation\n(CrMnFeCoNi system)',
                   fontsize=13,fontweight='bold')

from sklearn.decomposition import PCA
X_sc2=StandardScaler().fit_transform(X)
pca=PCA(n_components=2,random_state=42); Xp=pca.fit_transform(X_sc2)
sc=axes[1].scatter(Xp[:,0],Xp[:,1],c=df['yield_strength'],cmap='plasma',s=30,alpha=0.7)
plt.colorbar(sc,ax=axes[1],label='Yield Strength (MPa)')
axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var.)")
axes[1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var.)")
axes[1].set_title('PCA projection (yield strength)',fontsize=13,fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig1_hea_descriptors_pca.png',dpi=150,bbox_inches='tight')
plt.close()
print("  fig1 saved")

# Fig2: Feature importance
fig,axes=plt.subplots(1,3,figsize=(18,5))
for ax,tgt in zip(axes,tgt_cols):
    m=final_models[tgt].named_steps['m']
    imp=m.feature_importances_; idx=np.argsort(imp)[::-1][:8]; idx_r=idx[::-1]
    ax.barh([feat_cols[i] for i in idx_r],[imp[i] for i in idx_r],
            color=plt.cm.viridis(np.linspace(0.2,0.9,len(idx_r))))
    ax.set_title(f'Feature Importance\n{tgt}',fontsize=11,fontweight='bold')
    ax.set_xlabel('Importance',fontsize=10)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_hea_feature_importance.png',dpi=150,bbox_inches='tight')
plt.close()
print("  fig2 saved")

# Fig3: Pareto + BO learning curve
fig,axes=plt.subplots(1,2,figsize=(14,5))
axes[0].scatter(df['yield_strength'][~pareto_mask],df['elongation'][~pareto_mask],
                c='lightgray',alpha=0.5,s=20,label='Non-Pareto')
axes[0].scatter(df['yield_strength'][pareto_mask],df['elongation'][pareto_mask],
                c='red',s=60,zorder=5,label=f'Pareto ({pareto_mask.sum()} pts)')
axes[0].set_xlabel('Yield Strength (MPa)'); axes[0].set_ylabel('Elongation (%)')
axes[0].set_title('Pareto Front: Strength vs Ductility',fontsize=12,fontweight='bold'); axes[0].legend()
axes[1].plot(bo_df['n_obs'],bo_df['best_ys'],'r-o',linewidth=2,label='BO (GP+EI)')
axes[1].axhline(best_ys_init,color='gray',linestyle='--',label='Random baseline')
axes[1].set_xlabel('# Observations'); axes[1].set_ylabel('Best Yield Strength (MPa)')
axes[1].set_title('BO Learning Curve',fontsize=12,fontweight='bold'); axes[1].legend()
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig3_hea_pareto_bo.png',dpi=150,bbox_inches='tight')
plt.close()
print("  fig3 saved")

# Fig4: Predicted vs actual + composition heatmap
X_sc3=StandardScaler().fit_transform(X)
rf_cv=RandomForestRegressor(n_estimators=200,random_state=42,n_jobs=-1)
y_pred_cv=cross_val_predict(rf_cv,X_sc3,df['yield_strength'].values,cv=5)
r2_cv=r2_score(df['yield_strength'].values,y_pred_cv)
rmse_cv=np.sqrt(mean_squared_error(df['yield_strength'].values,y_pred_cv))

fig,axes=plt.subplots(1,2,figsize=(14,5))
axes[0].scatter(df['yield_strength'],y_pred_cv,alpha=0.5,s=20,c='steelblue')
mn,mx=df['yield_strength'].min(),df['yield_strength'].max()
axes[0].plot([mn,mx],[mn,mx],'r--',linewidth=2)
axes[0].set_xlabel('True Yield Strength (MPa)'); axes[0].set_ylabel('CV-Predicted (MPa)')
axes[0].set_title(f'RF: CV Prediction\nR²={r2_cv:.3f}, RMSE={rmse_cv:.1f} MPa',
                   fontsize=12,fontweight='bold')
sc=axes[1].scatter(df['x_Cr'],df['x_Ni'],c=df['yield_strength'],cmap='hot_r',s=30,alpha=0.8)
plt.colorbar(sc,ax=axes[1],label='Yield Strength (MPa)')
axes[1].set_xlabel('x_Cr'); axes[1].set_ylabel('x_Ni')
axes[1].set_title('Composition Space (Cr-Ni)\ncolored by yield strength',
                   fontsize=12,fontweight='bold')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig4_hea_pred_composition.png',dpi=150,bbox_inches='tight')
plt.close()
print("  fig4 saved")
print(f"\n[Cell 6] CV R²={r2_cv:.3f}, RMSE={rmse_cv:.1f} MPa")

# ─── Cell 7: Case study - composition optimal ─────────────
print("\n[Cell 7] Optimal composition search...")
rng2=np.random.RandomState(123)
grid_comps=rng2.dirichlet(np.ones(5),1000)
grid_rows=[]
for c in grid_comps:
    feat=calc_desc(c).reshape(1,-1)
    ys=final_models['yield_strength'].predict(feat)[0]
    el_=final_models['elongation'].predict(feat)[0]
    cr_=final_models['corrosion_resistance'].predict(feat)[0]
    grid_rows.append({'x_Cr':c[0],'x_Mn':c[1],'x_Fe':c[2],'x_Co':c[3],'x_Ni':c[4],
                      'pred_ys':ys,'pred_el':el_,'pred_cr':cr_})
grd=pd.DataFrame(grid_rows)
grd['score']=(grd['pred_ys']/grd['pred_ys'].max()*0.4+
              grd['pred_el']/grd['pred_el'].max()*0.3+
              grd['pred_cr']/grd['pred_cr'].max()*0.3)
best5=grd.nlargest(5,'score')
print("[Cell 7] Top-5 optimized compositions:")
print(best5[['x_Cr','x_Mn','x_Fe','x_Co','x_Ni','pred_ys','pred_el','pred_cr']].round(3).to_string(index=False))

cantor_f=calc_desc([0.2]*5).reshape(1,-1)
c_ys=final_models['yield_strength'].predict(cantor_f)[0]
c_el=final_models['elongation'].predict(cantor_f)[0]
c_cr=final_models['corrosion_resistance'].predict(cantor_f)[0]
best1=best5.iloc[0]
print(f"\nCantor alloy: YS={c_ys:.1f} MPa, El={c_el:.1f}%, Corr={c_cr:.2f}")
print(f"Best optimized: YS={best1['pred_ys']:.1f} MPa, El={best1['pred_el']:.1f}%, Corr={best1['pred_cr']:.2f}")
print(f"  Cr{best1['x_Cr']:.3f} Mn{best1['x_Mn']:.3f} Fe{best1['x_Fe']:.3f} Co{best1['x_Co']:.3f} Ni{best1['x_Ni']:.3f}")

# ─── Cell 8: pip freeze ──────────────────────────────────────
res=subprocess.run([sys.executable,'-m','pip','freeze'],capture_output=True,text=True)
with open(f'{DATA_DIR}/pip_freeze.txt','w') as f: f.write(res.stdout)
print("\n[Cell 8] pip freeze saved")
for line in res.stdout.split('\n'):
    if any(p in line.lower() for p in ['numpy','pandas','scikit','xgboost','matplotlib','scipy','seaborn']):
        print(f"  {line}")

print("\n=== ALL CELLS COMPLETE ===")

