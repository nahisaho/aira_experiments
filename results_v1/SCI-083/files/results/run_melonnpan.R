
# ============================================================
# MelonnPan: Metabolite Prediction from Microbiome Data
# ============================================================
library(melonnpan)

# --- Data loading ---
taxa <- read.csv("../results/taxa_abundance.csv", row.names=1)
metabolites <- read.csv("../results/metabolite_abundance.csv", row.names=1)

# --- MelonnPan training ---
# Train the model on paired microbiome-metabolome data
melonnpan.train(
  metab = metabolites,
  taxa = taxa,
  output = "../results/melonnpan_output"
)

# --- Prediction on new samples ---
melonnpan.predict(
  taxa = taxa,
  weight.matrix = "../results/melonnpan_output/MelonnPan_Trained_Weights.txt",
  output = "../results/melonnpan_predictions"
)

# --- Evaluation ---
# Compare predicted vs measured metabolites
predicted <- read.table(
  "../results/melonnpan_predictions/MelonnPan_Predicted_Metabolites.txt",
  header=TRUE, row.names=1
)

# Correlation analysis
cors <- sapply(1:ncol(metabolites), function(i) {
  if(colnames(metabolites)[i] %in% colnames(predicted)) {
    cor(metabolites[,i], predicted[,colnames(metabolites)[i]], method="spearman")
  } else { NA }
})

cor.df <- data.frame(
  metabolite = colnames(metabolites),
  spearman_r = cors,
  well_predicted = cors > 0.3
)
write.csv(cor.df, "../results/melonnpan_prediction_quality.csv")

cat("MelonnPan pipeline completed.\n")
