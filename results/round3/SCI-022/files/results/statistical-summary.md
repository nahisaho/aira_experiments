# Statistical Summary

Normality was checked with the Shapiro-Wilk test before comparing device efficiencies. Stable-group normality p = 0.023, non-stable-group normality p = 0.001, and Levene homoscedasticity p = 0.859. The mean PCE for geometrically stable candidates was 14.15 ± 4.93% versus 15.77 ± 5.22% for the remaining candidates. The mean difference was -1.62% (95% CI -4.37 to 1.21%), evaluated with Mann-Whitney U p = 0.1605; rank-biserial correlation = 0.24; Cohen's d = -0.32.

Model performance metrics are summarized below with mean ± SD across five folds:

- **bandgap_regressor**: r2_mean=0.945, r2_std=0.034, mae_mean=0.067, mae_std=0.014, rmse_mean=0.093, rmse_std=0.035, overfit_flag=0.000, prediction_r2=0.946, prediction_mae=0.068, prediction_rmse=0.099
- **pce_regressor**: r2_mean=0.958, r2_std=0.021, mae_mean=0.720, mae_std=0.104, rmse_mean=0.936, rmse_std=0.239, overfit_flag=0.000, prediction_r2=0.963, prediction_mae=0.723, prediction_rmse=0.968
- **stability_classifier**: accuracy_mean=0.980, accuracy_std=0.045, f1_mean=0.980, f1_std=0.045
