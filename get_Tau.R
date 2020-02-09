########################################################################################################
## Title: R modules for Tissue-Specificity Index Calculations
## Author: Madara Hetti-Arachchilage
## Credit:Tissue-specificity score source code is adapted from  https://doi.org/10.1093/bib/bbw008
## date: "February 3, 2020"
########################################################################################################

# Defining 'load_pkg' function here
load_pkg <- function (pkg_name, bioconductor = FALSE) {
  # This function checks if the requested package is already installed in the
  # environment. If not, the package is installed from the appropriate source
  # and then loaded into the workspace
  
  if (!pkg_name %in% installed.packages()) {
    if (bioconductor == TRUE) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(pkg_name)
    } else {install.packages(pkg_name)}
  }
  if (pkg_name %in% installed.packages()) {
    library(pkg_name, character.only = TRUE)
  } else {stop(paste("Cannot install/load package", pkg_name, sep = ": "))}
}
########################################################################################################

##### FUNCTIONS NEEDED FOR TISSUE SPECIFICITY CALCULATIONS
########################################################################################################
# Function to get tissue names
fSampNames <- function(meta, column){	
  
  sampNames <- unique(meta[,column])
  print(sampNames)
  return(sampNames)
}	


########################################################################################################
# Function to get normalized expression matrix
# Takes three arguments; exprssion matrix, meta file with sample information
# column in meta file specifying the sample group
fNormalizeData <- function(expr, meta, column){
  
  y <- DGEList(expr, group= meta[,column])
  dge <- calcNormFactors(y, method="TMM")
  orgExpression <- cpm(dge, log = FALSE)
  return(orgExpression)
}

########################################################################################################

# Function requires data frame with expression values and set of strings that matches with sample names
# Mean values between replicates are calculated
fReplicateMean <- function(x, names){
  
  y <- data.frame(matrix(nrow = nrow(x), ncol = length(unique(names))))
  for (i in 1:length(names)){
    
    if (length(which(grepl(names[i], colnames(x)) >0)) > 1){
      y[,i] <- rowMeans(x[, grepl(names[i], colnames(x)) >0], dims = 1)
      names(y)[i] <- paste("Averaged.CPM.", names[i], sep="")
    } else if (length(which(grepl(names[i], colnames(x)) >0)) == 1){
      print(i)
      y[,i] <- x[,grepl(names[i], colnames(x)) >0]
      names(y)[i] <- paste("Averaged.CPM.", names[i], sep="")
    }
  } 
  row.names(y) <- row.names(x)
  return(y)
}
########################################################################################################

# This function quantile normalises the entire data frame of log normalised 
# counts, allowing comparisons across tissues. First, the row and column 
# names are collected and all zeros are changed to NA. After converting to a 
# matrix, the data frame is quantile normalised. Finally, all NAs are 
# reverted back to zero, the matrix is converted back to a data frame and 
# row and column names are re-attached
quantNorm <- function(x){
  x.cols <- names(x)
  x.rows <- row.names(x)
  x[x == 0] <- NA
  x_m <- as.matrix(x)
  x <- round(preprocessCore::normalize.quantiles(x_m), digits = 3)
  x[is.na(x)] <- 0
  x <-data.frame(x)
  names(x)[c(seq_along(x.cols))] <- x.cols
  row.names(x)[c(seq_along(x.rows))] <- x.rows
  return(x)
}

########################################################################################################

# Function require a vector with expression of one gene in different tissues.
# Mean value per gene is calculated
fmean <- function(x){
  if(!all(is.na(x))) {
    res <- mean(x, na.rm=TRUE)
  } else {res <- NA}
  return(res)
}
########################################################################################################

# Function require a vector with expression of one gene in different tissues.
# Max value per gene is calculated.
fmax <- function(x){
  if(!all(is.na(x))) {
    res <- max(x, na.rm=TRUE)
  } else {res <- NA}
  return(res)
}
########################################################################################################

# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
# Minimum 2 tissues

fTau <- function(x){
  if(all(!is.na(x))) {
    if(min(x, na.rm=TRUE) >= 0) {
      if(max(x)!=0) {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {res <- 0}
    } else {res <- NA} 
  } else {res <- NA} 
  return(res)
}
########################################################################################################

#Function require a data frame with expression data and 
#Function give back a vector with the tissue with highest expression
fTissue <- function(x)
{
  if(!all(is.na(x))) {
    x[,-1] <- t(apply(x[,-1], c(1), function(x){x <- ifelse(x==max(x),1,0)}))
    x$Organ <- apply(x[,-1], c(1), function(x){x <- which(x>0)})
    x$Organ.Number <- sapply(x$Organ, function(x){x <- as.numeric(x[1])})
    names <- gsub("Averaged.CPM.", "", colnames(x[,-1]))
    x$Organ.Name <- names[x$Organ.Number]
    res <- x$Organ.Name		
  } else {res <- NA}
  return(res)
}
########################################################################################################

#calling Tau
fTS <- function(expr, tNames){	
  
  nTissues <- length(tNames)
  temp <- c(paste("Averaged.CPM.", tNames[1:nTissues], sep=""))
  expr$Tau <- apply(expr[,temp], 1, fTau)
  return(expr)
}
########################################################################################################
# Function prepare the data for further analysis. 
# expr = data set, meta = sample info, CPM = cutt-off
# calculating mean for all tissue samples for each developmental stage, 
# Filter out genes that are not expressed, quantile normalization, calculates Tau for each gene

mainf_Benchmark_Tau <- function(expr_norm, meta, cutoff){
  
  # Pre-processing RNA-seq data
  tNames <- fSampNames (meta, "Tissue")
  expr_avg <- fReplicateMean(x = expr_norm, names = tNames)
  expr_avg[expr_avg == 0] <- 0.0000000001
  expr_avg <- log2(expr_avg)
  expr_avg[expr_avg < 0] <- 0
  keep_genes <- rowSums(expr_avg > log2(cutoff)) >= 1
  expr_avg <- expr_avg[keep_genes,]
  expr_quant <- quantNorm(expr_avg)
  
  # Calculating Tissue-Specificity
  Tissue_enriched <- fTS(expr_quant, tNames)
  Tissue_enriched_df <- as.data.frame(cbind(row.names(Tissue_enriched), Tissue_enriched[,1:(ncol(Tissue_enriched)-1)]))
  tissue <- as.data.frame(fTissue(Tissue_enriched_df))
  r_tissue <- cbind(tissue, row.names(Tissue_enriched), Tissue_enriched[,ncol(Tissue_enriched)])
  colnames(r_tissue) <- c("Tissue", "Gene_ID", "Tau")
  return(r_tissue)
  
}

###########################################################################################################
