
# ============================================================
# mixOmics DIABLO Integration Pipeline for IBD Case Study
# ============================================================
library(mixOmics)

# --- Data loading ---
taxa <- read.csv("../results/taxa_abundance.csv", row.names=1)
metabolites <- read.csv("../results/metabolite_abundance.csv", row.names=1)
metadata <- read.csv("../results/metadata.csv")

Y <- factor(metadata$group, levels = c("Control", "UC", "CD"))

# --- Pre-processing ---
# CLR transform for taxa
taxa.clr <- logratio.transfo(as.matrix(taxa), logratio = "CLR", offset = 1)

# Log transform for metabolites
metabolites.log <- log2(as.matrix(metabolites) + 1)

# --- Design matrix ---
# Specify expected correlation between blocks
design <- matrix(c(0,   1,   0.1,
                   1,   0,   0.1,
                   0.1, 0.1, 0  ), ncol=3, nrow=3,
                 dimnames = list(c("taxa", "metabolites", "clinical"),
                                 c("taxa", "metabolites", "clinical")))

# --- DIABLO model (block.splsda) ---
X <- list(taxa = taxa.clr, metabolites = metabolites.log)

# Tune keepX (number of features per component)
# tune.diablo <- tune.block.splsda(
#   X, Y, ncomp = 2,
#   test.keepX = list(taxa = c(5, 10, 15), metabolites = c(10, 15, 20)),
#   design = design[1:2, 1:2],
#   validation = "Mfold", folds = 5, nrepeat = 10
# )

# Fit final model
diablo.model <- block.splsda(
  X, Y, ncomp = 2,
  keepX = list(taxa = c(10, 10), metabolites = c(15, 15)),
  design = design[1:2, 1:2]
)

# --- Performance evaluation ---
perf.diablo <- perf(diablo.model, validation = "Mfold", folds = 5,
                     nrepeat = 10, progressBar = TRUE)

# --- Visualization ---
pdf("../results/../figures/diablo_plotIndiv.pdf")
plotIndiv(diablo.model, ind.names = FALSE, legend = TRUE,
          title = "DIABLO Sample Plot (IBD)")
dev.off()

pdf("../results/../figures/diablo_circosPlot.pdf")
circosPlot(diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks = c("steelblue", "darkorange"),
           color.cor = c("red", "blue"))
dev.off()

pdf("../results/../figures/diablo_loadings.pdf")
plotLoadings(diablo.model, comp = 1, contrib = "max",
             method = "median", legend.color = c("blue", "red", "green"))
dev.off()

pdf("../results/../figures/diablo_network.pdf")
network(diablo.model, blocks = c(1, 2), cutoff = 0.4,
        color.node = c("steelblue", "darkorange"))
dev.off()

# --- Export selected features ---
selected.taxa <- selectVar(diablo.model, block = "taxa", comp = 1)$taxa$name
selected.met <- selectVar(diablo.model, block = "metabolites", comp = 1)$metabolites$name

write.csv(data.frame(feature = selected.taxa), "../results/diablo_selected_taxa.csv")
write.csv(data.frame(feature = selected.met), "../results/diablo_selected_metabolites.csv")

# --- Save performance ---
sink("../results/diablo_performance.txt")
print(perf.diablo)
sink()

cat("DIABLO pipeline completed successfully.\n")
