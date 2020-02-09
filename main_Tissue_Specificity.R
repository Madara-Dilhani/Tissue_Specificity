########################################################################################################
## title: "Tissue Specific Gene Identification"
## author: "Madara Hetti-Arachchilage; mhettiar 'at' illinois.edu"
## date: "February 3, 2020"
########################################################################################################

##### Importing libraries and custom scripts 
# Loading R source code which has defined custom functions needed for later steps
# Loading R/Biconductor packages (checks if the requested package is already installed in 
# the environment. If not, the package is installed from the appropriate source and then 
# loaded into the workspace

source("get_Tau.R")
load_pkg("ggplot2", bioconductor = FALSE)
load_pkg("tidyverse", bioconductor = FALSE)
load_pkg("RColorBrewer", bioconductor = FALSE)
load_pkg("pheatmap", bioconductor = FALSE)
########################################################################################################

##### Importing the data
# Sample information/meta data
meta <- read.csv("Data/meta.csv", header = T, as.is = T)
# Normalized Counts; CPM
expr_norm <- read.csv(file = "Data/TPM.csv", header = T, as.is = T, row.names = 1)
########################################################################################################

##### Calculating Tissue-Specificity Index
# This step calls "mainf_Benchmark_Tau" function which calculates tissue specificity index (Tau) 
# and saves Tau indexes for all genes to a .csv file. This function takes 4 arguments, 
# normalized CPM values, sample information (there should be following colums; 'Tissue', 'Stage', 'Tissue_Stage'), 
# CPM cutoff for not expressed genes.

out <- mainf_Benchmark_Tau(expr_norm, meta, 5)

# Save output file for Tissue Specificity Scores calculated for all stages 
TS_by_AllStage <- "Output/Tissue_Specificity/Tau_by_AllStage.csv"
write.csv(out, file= TS_by_AllStage, quote=F)
########################################################################################################

##### Gene expression matrix for Tissue-specific genes
# Filter highly tissue-specific genes
out_TS <- subset(out, out$Tau >= 0.85 & out$Tissue == "tis1")
geneList <- unique(out_TS$Gene_ID)
# This step takes chosen tissue-specific genes and save it's CPM values for different samples
expr_norm_TS <- expr_norm[geneList,]
fltoutput <- "Output/Tissue_Specificity/Tau_0.85_CPM.csv"
write.csv(expr_norm_TS, file= fltoutput, quote=F, row.names = F)
########################################################################################################

##### Heatmap for tissue-specific genes
# This step generate heatmap for chosen tissue-specific genes.
log_expr <- log2(expr_norm_TS+1)
breaksList = seq(0, 10, by = 1)
pheatmap(log_expr, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 20, show_rownames = F, cluster_rows = T, legend = T, legend.cex = 0.9, 
         fontsize_row = 15, fontsize_col = 10,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList)
########################################################################################################

