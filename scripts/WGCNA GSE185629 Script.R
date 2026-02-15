# 2025 Nov 03  pol_LTMH WGCNA — cleaned & consolidated script
# ============================================================

# STEP 0: Packages & options
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ============================================================
# STEP 1: Load expression matrix (preserve names)
# ============================================================
expr_raw <- read.delim("expr_matrix_GSE185629.txt", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
cat("Expression matrix (genes x samples) dims (raw):", dim(expr_raw), "\n")
# Expression matrix (genes x samples) dims (raw): 1468 8 

# convert to numeric (columns are samples in GEO-style files)
expr_numeric <- as.data.frame(lapply(expr_raw, function(x) as.numeric(as.character(x))))
rownames(expr_numeric) <- rownames(expr_raw)   # restore gene names
expr <- expr_numeric

# Transpose for WGCNA: samples as rows, genes as columns
datExpr <- as.data.frame(t(expr))
cat("datExpr dimensions (samples x genes):", dim(datExpr), "\n")  # should be samples x genes
# datExpr dimensions (samples x genes): 8 1468 

# Standardize datExpr sample names: replace spaces/dots with underscore and trim
rownames(datExpr) <- gsub("\\.", "_", trimws(rownames(datExpr)))
rownames(datExpr) <- gsub(" ", "_", rownames(datExpr))

# ============================================================
# STEP 2: Load and prepare metadata (traits)
# ============================================================
traits <- read.table("metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# If FileName column exists, use it as rownames (and replace spaces with underscores)
if ("FileName" %in% colnames(traits)) {
  rownames(traits) <- gsub(" ", "_", trimws(traits$FileName))
  traits$FileName <- NULL
} else {
  # if first column is actually filenames, you may need to set it manually:
  # rownames(traits) <- gsub(" ", "_", trimws(traits[[1]])); traits <- traits[,-1]
  warning("No 'FileName' column found in traits; verify rownames manually.")
}

# Convert trait labels to numeric (Control = 0, LTMH = 1)
if ("sample" %in% colnames(traits)) {
  traits$sample <- as.numeric(factor(traits$sample, levels = c("Control", "LTMH"))) - 1
} else {
  warning("'sample' column not found in traits: check your metadata column names.")
}

# Standardize traits rownames as well
rownames(traits) <- gsub("\\.", "_", trimws(rownames(traits)))
rownames(traits) <- gsub(" ", "_", rownames(traits))

cat("Traits (after cleanup):\n")
print(traits)

# ============================================================
# STEP 3: QC - goodSamplesGenes and sample clustering
# ============================================================
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("Removing samples/genes with too many missing values...\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

# ============================================================
# STEP 4: Make sure traits & datExpr share same sample names (and order)
# ============================================================

# Standardize names again (defensive)
rownames(datExpr) <- gsub("\\.", "_", trimws(rownames(datExpr)))
rownames(traits)  <- gsub("\\.", "_", trimws(rownames(traits)))

# Find common samples
commonSamples <- intersect(rownames(datExpr), rownames(traits))
cat("Number of samples in datExpr:", nrow(datExpr), "\n")
# 8
cat("Number of samples in traits:", nrow(traits), "\n")
# 8
cat("Number of common samples:", length(commonSamples), "\n")
# 8

if (length(commonSamples) == 0) stop("No common samples between datExpr and traits — check sample names.")

# Subset and reorder both to the commonSamples
datExpr <- datExpr[commonSamples, , drop = FALSE]
traits  <- traits[commonSamples, , drop = FALSE]

# Final check
stopifnot(all(rownames(datExpr) == rownames(traits)))
cat("Sample names aligned between datExpr and traits. Samples:\n"); print(rownames(datExpr))

# ============================================================
# STEP 5: Pick soft threshold (optional visualization)
# ============================================================
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# (plotting code omitted here — leave as you prefer)

# Choose power (use your chosen softPower; e.g. 6)
softPower <- 6

# ============================================================
# STEP 6: Build network and detect modules (blockwiseModules)
# ============================================================
net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, saveTOMFileBase = "GSE185629_TOM",
                        verbose = 3)

# Get merged colors and plot dendrogram
mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# ============================================================
# STEP 7: Module eigengenes and alignment
# ============================================================
# Use net$MEs (module eigengenes computed by blockwiseModules)
MEs <- net$MEs

# Standardize MEs rownames
rownames(MEs) <- gsub("\\.", "_", trimws(rownames(MEs)))

# Ensure traits and datExpr names are standardized (again, defensive)
rownames(traits) <- gsub("\\.", "_", trimws(rownames(traits)))
rownames(datExpr) <- gsub("\\.", "_", trimws(rownames(datExpr)))

# Intersect to be safe (should already match)
commonSamples2 <- Reduce(intersect, list(rownames(MEs), rownames(traits), rownames(datExpr)))
cat("Common samples among MEs, traits, datExpr:", length(commonSamples2), "\n")
# Common samples among MEs, traits, datExpr: 8 

if (length(commonSamples2) == 0) stop("No common samples among MEs/traits/datExpr — check sample naming.")

MEs <- MEs[commonSamples2, , drop = FALSE]
traits <- traits[commonSamples2, , drop = FALSE]
datExpr <- datExpr[commonSamples2, , drop = FALSE]

# Final verification
stopifnot(all(rownames(MEs) == rownames(traits)))
stopifnot(all(rownames(MEs) == rownames(datExpr)))

# ============================================================
# STEP 8: Module–Trait Relationships
# ============================================================
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Heatmap visualization — saving PNG
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
png("Module_trait_heatmap.png", width = 1000, height = 800)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               main = "Module-trait relationships")
dev.off()

# Print numeric results to console
print(signif(moduleTraitCor, 3))
#     sample
ME1  0.985
ME2 -0.926
print(signif(moduleTraitPvalue, 3))
#       sample
ME1 7.77e-06
ME2 9.53e-04

# ============================================================
# STEP 9: Export Network for Cytoscape (choose module)
# ============================================================
# Note: probes are gene names (columns of datExpr)
probes <- colnames(datExpr)

# mergedColors corresponds to net$colors mapped to probes; ensure alignment
mergedColors <- labels2colors(net$colors)  # length should equal number of probes in that block
# If there are multiple blocks, you may need to handle net$colors appropriately (this script assumes one block)

module <- "turquoise"  # change as needed
inModule <- (mergedColors == module)
modGenes <- probes[inModule]

# Recompute TOM from expression (faster to reuse saved TOM if available)
TOM_full <- TOMsimilarityFromExpr(datExpr, power = softPower)

# Export for Cytoscape (use only rows/cols for inModule)
if (sum(inModule) > 1) {
  exportNetworkToCytoscape(TOM_full[inModule, inModule],
                           edgeFile = paste0("CytoscapeInput-edges-", module, ".txt"),
                           nodeFile = paste0("CytoscapeInput-nodes-", module, ".txt"),
                           weighted = TRUE, threshold = 0.05,
                           nodeNames = modGenes,
                           nodeAttr = mergedColors[inModule])
  cat("Cytoscape files created for module:", module, "\n")
} else {
  warning("Not enough genes in selected module to export network.")
}

# Done
cat("WGCNA pipeline completed. Samples used:", nrow(datExpr), "Modules detected:", length(unique(mergedColors)), "\n")

