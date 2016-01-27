# This script is written by Dr. Hemang Parikh as on January 26, 2016
# The Health Informatics Institute (HII) at the University of South Florida

# Median local background is performed
# BASH and HULK methods are used for beads artifact detection
# Summarization is performed with un-log transformation because neqc requires gene expression values to be un-logged
# Outliers are removed using the Illumina 3 M.A.D cut-off
# The mean and standard deviation for each bead type are reported

# In beadarray
# No background normalization
# Cubic spline normalization
# No transformation

# In lumi
# No function for cubic spline normalization

# To import beadarray, lumi and illuminaHumanv4.db libraries
suppressMessages(library(beadarray))
suppressMessages(library(illuminaHumanv4.db))
suppressMessages(library(lumi))
suppressMessages(library(limma))

# Versions of R packages: sessionInfo() or packageVersion("illuminaHumanv4.db")
# beadarray_2.20.0
# illuminaHumanv4.db_1.26.0
# org.Hs.eg.db_3.2.3
# RSQLite_1.0.0
# DBI_0.3.1
# AnnotationDbi_1.32.3
# IRanges_2.4.6
# S4Vectors_0.8.5
# ggplot2_2.0.0
# Biobase_2.30.0
# BiocGenerics_0.16.1
# Rcpp_0.12.2
# XVector_0.10.0
# magrittr_1.5
# GenomicRanges_1.22.3
# zlibbioc_1.16.0
# munsell_0.4.2
# colorspace_1.2-6
# stringr_1.0.0
# plyr_1.8.3
# GenomeInfoDb_1.6.1
# tools_3.2.3
# base64_1.2
# grid_3.2.3
# gtable_0.1.2
# reshape2_1.4.1
# limma_3.26.3
# stringi_1.0-1
# BeadDataPackR_1.22.0
# scales_0.3.0
# illuminaio_0.12.0
# lumi_2.22.0
# limma_3.26.3

#------------------------------------------------------------------------
# Get cli args and assign appropriate variables
#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# To use the combined filter data directory for storing qc data and results
filter_combined_result_dir <- toString(args[1])

# To change the directory
setwd(filter_combined_result_dir)

# dir.create(paste(getwd(), "/filter_results/", sep = ""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

# To read the gene expression data from the *.txt files
# Annotation from https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html
# To read data for beadarray
data.gene <- readBeadSummaryData(dataFile = "wo_control_expression.txt", ProbeID = "Probe_ID", skip = 0, columns = list(exprs = "AVG_Signal", se.exprs = "BEAD_STDEV", nObservations = "Avg_NBEADS", Detection = "Detection Pval"), annoCols = c("Entrez_ID"))
data.pheno <- read.table(file = "BeadArray_phenotype_details.txt", sep = "\t", header = TRUE)

# To apply the normaliseIllumina function to normalize and transform the intensities
# Perform no background normalization, no transformation and cubic spline normalization
data.gene.norm <- normaliseIllumina(data.gene, transform = "none", method = "qspline")

# Multidimensional scaling (MDS), assesses sample similarity based on pair-wise distances between samples. Ideally, sample should separate based on biological variables (RNA source, sex, treatment etc.) but often technical effects
# (such as samples processed together on the same BeadChip) may dominate
# To generate the MDS plot
pdf(file = "all_mds_plot_after_normalization.pdf", width = 11, height = 8.5)
plotMDS(exprs(data.gene.norm))
garbage <- dev.off()

# To check the effects of boxes in MDS plot
pdf(file = "all_mds_box_effects_after_normalization.pdf", width = 11, height = 8.5)
plotMDS(exprs(data.gene.norm), labels = data.pheno$Box)
garbage <- dev.off()

# To check the effects of BeadChips in MDS plot
pdf(file = "all_mds_beadchip_effects_after_normalization.pdf", width = 11, height = 8.5)
plotMDS(exprs(data.gene.norm), labels = data.pheno$Chip_Barcode)
garbage <- dev.off()

# To check the effects of case-controls in MDS plot
pdf(file = "all_mds_case_control_effects_after_normalization.pdf", width = 11, height = 8.5)
plotMDS(exprs(data.gene.norm), labels = data.pheno$casestatus)
garbage <- dev.off()

# To check the effects of sex in MDS plot
pdf(file = "all_mds_sex_effects_after_normalization.pdf", width = 11, height = 8.5)
plotMDS(exprs(data.gene.norm), labels = data.pheno$Sex)
garbage <- dev.off()

# To check the effects of donors in MDS plot
pdf(file = "all_mds_donors_effects_after_normalization.pdf", width = 11, height = 8.5)
plotMDS(exprs(data.gene.norm), labels = data.pheno$Donor_Number)
garbage <- dev.off()

# To retrieve quality information and verify that probes annotated as "Bad" or "No Match" generally have lower signal
ids.arrays <- as.character(rownames(exprs(data.gene.norm)))
qual <- unlist(mget(ids.arrays, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
arrays.qual = table(qual)
ave.signal.arrays = rowMeans(exprs(data.gene.norm))

# To make a box plot with array probe quality information
pdf(file = "all_probe_quality_after_normalization.pdf", width = 11, height = 8.5)
boxplot(ave.signal.arrays ~ qual)
garbage <- dev.off()

# To store normalized gene expression data for only high quality probes
rem <- qual == "No match" | qual == "Bad"
data.gene.norm.filt <- data.gene.norm[!rem, ]

# The "Bad" probes are hybridized to multiple places in the genome and often contain repetitive sequences
# queryIDs <- names(which(qual == "Bad" & ave.signal.arrays > 12))
# unlist(mget(queryIDs, illuminaHumanv4REPEATMASK))
# unlist(mget(queryIDs, illuminaHumanv4SECONDMATCHES))
# mget("ILMN_1692145", illuminaHumanv4PROBESEQUENCE)



# To make a set of highly-variable probes and cluster the samples
iqr <- apply(exprs(data.gene.norm.filt), 1, IQR, na.rm = TRUE)
top.var <- order(iqr, decreasing = TRUE)[1:500]

# To generate the cluster plot
pdf(file = "all_cluster_plot_after_normalization.pdf", width = 11, height = 8.5)
plot(hclust(dist(t(exprs(data.gene.norm.filt)[top.var,]))))
garbage <- dev.off()

# To check the effects of boxes in cluster plot
pdf(file = "all_cluster_box_effects_after_normalization.pdf", width = 11, height = 8.5)
plot(hclust(dist(t(exprs(data.gene.norm.filt)[top.var,]))), labels = data.pheno$Box)
garbage <- dev.off()

# To check the effects of BeadChips in cluster plot
pdf(file = "all_cluster_beadchip_effects_after_normalization.pdf", width = 11, height = 8.5)
plot(hclust(dist(t(exprs(data.gene.norm.filt)[top.var,]))), labels = data.pheno$Chip_Barcode)
garbage <- dev.off()

# To check the effects of case-controls in cluster plot
pdf(file = "all_cluster_case_control_effects_after_normalization.pdf", width = 11, height = 8.5)
plot(hclust(dist(t(exprs(data.gene.norm.filt)[top.var,]))), labels = data.pheno$casestatus)
garbage <- dev.off()

# To check the effects of sex in cluster plot
pdf(file = "all_cluster_sex_effects_after_normalization.pdf", width = 11, height = 8.5)
plot(hclust(dist(t(exprs(data.gene.norm.filt)[top.var,]))), labels = data.pheno$Sex)
garbage <- dev.off()

# To check the effects of donors in cluster plot
pdf(file = "all_cluster_donors_effects_after_normalization.pdf", width = 11, height = 8.5)
plot(hclust(dist(t(exprs(data.gene.norm.filt)[top.var,]))), labels = data.pheno$Donor_Number)
garbage <- dev.off()




# To generate the heatmap
pdf(file = "all_heatmap_plot_after_normalization.pdf", width = 11, height = 8.5)
heatmap(exprs(data.gene.norm.filt)[top.var,])
garbage <- dev.off()

# To check the effects of boxes in heatmap
pdf(file = "all_heatmap_box_effects_after_normalization.pdf", width = 11, height = 8.5)
heatmap(exprs(data.gene.norm.filt)[top.var,], labCol = data.pheno$Box)
garbage <- dev.off()

# To check the effects of BeadChips in heatmap
pdf(file = "all_heatmap_beadchip_effects_after_normalization.pdf", width = 11, height = 8.5)
heatmap(exprs(data.gene.norm.filt)[top.var,], labCol = data.pheno$Chip_Barcode)
garbage <- dev.off()

# To check the effects of case-controls in heatmap
pdf(file = "all_heatmap_case_control_effects_after_normalization.pdf", width = 11, height = 8.5)
heatmap(exprs(data.gene.norm.filt)[top.var,], labCol = data.pheno$casestatus)
garbage <- dev.off()

# To check the effects of sex in heatmap
pdf(file = "all_heatmap_sex_effects_after_normalization.pdf", width = 11, height = 8.5)
heatmap(exprs(data.gene.norm.filt)[top.var,], labCol = data.pheno$Sex)
garbage <- dev.off()

# To check the effects of donors in heatmap
pdf(file = "all_heatmap_donors_effects_after_normalization.pdf", width = 11, height = 8.5)
heatmap(exprs(data.gene.norm.filt)[top.var,], labCol = data.pheno$Donor_Number)
garbage <- dev.off()


# To annotate Illumina IDs
idsTosymbols = as.matrix(toTable(illuminaHumanv4ENTREZID))

# To store gene expression values for all the array
array.normalized <- matrix(, nrow = dim(as.matrix(exprs(data.gene.norm)))[1], ncol = (dim(as.matrix(exprs(data.gene.norm)))[2] + 2))

# To store Illumina IDs information
array.normalized[, 1] = row.names(exprs(data.gene.norm))

# To get annotation for each Illumina ID
for (j.id in 1:dim(as.matrix(array.normalized))[1]) {

  # To get the gene symbol for each Illumina ID
  j.symbol = which(as.character(idsTosymbols[, 1]) == as.character(as.matrix(array.normalized)[j.id, 1]))

  # To store Illumina IDs annotation
  if (length(j.symbol) == 0) {
    array.normalized[j.id, 2] = ""
  }

  else {
    array.normalized[j.id, 2] = as.character(idsTosymbols[j.symbol, 2])
  }
}

# To define column names
colnames(array.normalized) <- c("Probe_ID", "Entrez_ID", colnames(exprs(data.gene.norm)))


# To store gene expression values for all the array
array.normalized.filt <- matrix(, nrow = dim(as.matrix(exprs(data.gene.norm.filt)))[1], ncol = (dim(as.matrix(exprs(data.gene.norm.filt)))[2] + 2))

# To store Illumina IDs information
array.normalized.filt[, 1] = row.names(exprs(data.gene.norm.filt))

# To get annotation for each Illumina ID
for (k.id in 1:dim(as.matrix(array.normalized.filt))[1]) {

  # To get the gene symbol for each Illumina ID
  k.symbol = which(as.character(idsTosymbols[, 1]) == as.character(as.matrix(array.normalized.filt)[k.id, 1]))

  # To store Illumina IDs annotation
  if (length(k.symbol) == 0) {
    array.normalized.filt[k.id, 2] = ""
  }

  else {
    array.normalized.filt[k.id, 2] = as.character(idsTosymbols[k.symbol, 2])
  }
}

# To define column names
colnames(array.normalized.filt) <- c("Probe_ID", "Entrez_ID", colnames(exprs(data.gene.norm.filt)))

# To gene expression values
array.normalized[, 3:dim(as.matrix(array.normalized))[2]] = exprs(data.gene.norm)[, 1:dim(as.matrix(exprs(data.gene.norm)))[2]]
array.normalized.filt[, 3:dim(as.matrix(array.normalized.filt))[2]] = exprs(data.gene.norm.filt)[, 1:dim(as.matrix(exprs(data.gene.norm.filt)))[2]]

# To create different files with gene expression data
write.table(array.normalized, "expression_normalized.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt, "expression_normalized_filter.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# To filter out samples based on following
donor.1 <- data.frame(NULL)
n_donor.1 = 0

# To get the donors information
# To run a for loop for all the samples
for (donor.1.info in 1:dim(as.matrix(exprs(data.gene.norm)))[2]) {

  # To filter out samples based on donor information
  if ((!is.na(as.character(data.pheno$Donor_Number[donor.1.info]))) & ((as.character(data.pheno$Donor_Number[donor.1.info])) == as.character("Donor 1"))) {
      n_donor.1 = n_donor.1 + 1
      donor.1[n_donor.1, 1] = as.character(as.matrix(data.pheno)[donor.1.info, 1])
  }
}


donor.2 <- data.frame(NULL)
n_donor.2 = 0

# To get the donors information
# To run a for loop for all the samples
for (donor.2.info in 1:dim(as.matrix(exprs(data.gene.norm)))[2]) {

  # To filter out samples based on donor information
  if ((!is.na(as.character(data.pheno$Donor_Number[donor.2.info]))) & ((as.character(data.pheno$Donor_Number[donor.2.info])) == as.character("Donor 2"))) {
    n_donor.2 = n_donor.2 + 1
    donor.2[n_donor.2, 1] = as.character(as.matrix(data.pheno)[donor.2.info, 1])
  }
}

# To store gene expression values for donors and other samples
array.normalized.donor.1 <- matrix(, nrow = dim(as.matrix(array.normalized))[1], ncol = (dim(donor.1)[1]) + 2)
array.normalized.filt.donor.1 <- matrix(, nrow = dim(as.matrix(array.normalized.filt))[1], ncol = (dim(donor.1)[1]) + 2)
array.normalized.donor.2 <- matrix(, nrow = dim(as.matrix(array.normalized))[1], ncol = (dim(donor.2)[1]) + 2)
array.normalized.filt.donor.2 <- matrix(, nrow = dim(as.matrix(array.normalized.filt))[1], ncol = (dim(donor.2)[1]) + 2)
array.normalized.wo.donors <- matrix(, nrow = dim(as.matrix(array.normalized))[1], ncol = (dim(as.matrix(array.normalized))[2] - ((dim(donor.1)[1]) + (dim(donor.2)[1]))))
array.normalized.filt.wo.donors <- matrix(, nrow = dim(as.matrix(array.normalized.filt))[1], ncol = (dim(as.matrix(array.normalized.filt))[2] - ((dim(donor.1)[1]) + (dim(donor.2)[1]))))

array.normalized.donor.1[, 1:2] <- array.normalized[, 1:2]
array.normalized.filt.donor.1[, 1:2] <- array.normalized.filt[, 1:2]
array.normalized.donor.2[, 1:2] <- array.normalized[, 1:2]
array.normalized.filt.donor.2[, 1:2] <- array.normalized.filt[, 1:2]
array.normalized.wo.donors[, 1:2] <- array.normalized[, 1:2]
array.normalized.filt.wo.donors[, 1:2] <- array.normalized.filt[, 1:2]

# To give a column numbers
# To give a column numbers
n_donor.1.col = 2
n_donor.2.col = 2
n_wo.donors.col = 2
n_donor.1.col.filt = 2
n_donor.2.col.filt = 2
n_wo.donors.col.filt = 2

col.names.header = c("Probe_ID", "Entrez_ID")
n_donor.1.col.names = col.names.header
n_donor.2.col.names = col.names.header
n_wo.donors.col.names = col.names.header
n_donor.1.col.names.filt = col.names.header
n_donor.2.col.names.filt = col.names.header
n_wo.donors.col.names.filt = col.names.header

for (donor.data in 1:dim(as.matrix(exprs(data.gene.norm)))[2]) {

  # To filter out samples based on donor information
  if (as.character(colnames(exprs(data.gene.norm))[donor.data]) %in% as.character(as.matrix(donor.1))) {
    n_donor.1.col = n_donor.1.col + 1
    array.normalized.donor.1[, n_donor.1.col] = array.normalized[, (2 + donor.data)]
    n_donor.1.col.names = c(n_donor.1.col.names, as.character(colnames(exprs(data.gene.norm))[donor.data]))
  }

  else if (as.character(colnames(exprs(data.gene.norm))[donor.data]) %in% as.character(as.matrix(donor.2))) {
    n_donor.2.col = n_donor.2.col + 1
    array.normalized.donor.2[, n_donor.2.col] = array.normalized[, (2 + donor.data)]
    n_donor.2.col.names = c(n_donor.2.col.names, as.character(colnames(exprs(data.gene.norm))[donor.data]))
  }

  else {
    n_wo.donors.col = n_wo.donors.col + 1
    array.normalized.wo.donors[, n_wo.donors.col] = array.normalized[, (2 + donor.data)]
    n_wo.donors.col.names = c(n_wo.donors.col.names, as.character(colnames(exprs(data.gene.norm))[donor.data]))
  }

}

for (donor.data.filt in 1:dim(as.matrix(exprs(data.gene.norm.filt)))[2]) {

  # To filter out samples based on donor information
  if (as.character(colnames(exprs(data.gene.norm.filt))[donor.data.filt]) %in% as.character(as.matrix(donor.1))) {
    n_donor.1.col.filt = n_donor.1.col.filt + 1
    array.normalized.filt.donor.1[, n_donor.1.col.filt] = array.normalized.filt[, (2 + donor.data.filt)]
    n_donor.1.col.names.filt = c(n_donor.1.col.names.filt, as.character(colnames(exprs(data.gene.norm.filt))[donor.data.filt]))
  }

  else if (as.character(colnames(exprs(data.gene.norm.filt))[donor.data.filt]) %in% as.character(as.matrix(donor.2))) {
    n_donor.2.col.filt = n_donor.2.col.filt + 1
    array.normalized.filt.donor.2[, n_donor.2.col.filt] = array.normalized.filt[, (2 + donor.data.filt)]
    n_donor.2.col.names.filt = c(n_donor.2.col.names.filt, as.character(colnames(exprs(data.gene.norm.filt))[donor.data.filt]))
  }

  else {
    n_wo.donors.col.filt = n_wo.donors.col.filt + 1
    array.normalized.filt.wo.donors[, n_wo.donors.col.filt] = array.normalized.filt[, (2 + donor.data.filt)]
    n_wo.donors.col.names.filt = c(n_wo.donors.col.names.filt, as.character(colnames(exprs(data.gene.norm.filt))[donor.data.filt]))
  }

}

# To define column names
colnames(array.normalized.donor.1) <- n_donor.1.col.names
colnames(array.normalized.donor.2) <- n_donor.2.col.names
colnames(array.normalized.wo.donors) <- n_wo.donors.col.names
colnames(array.normalized.filt.donor.1) <- n_donor.1.col.names.filt
colnames(array.normalized.filt.donor.2) <- n_donor.2.col.names.filt
colnames(array.normalized.filt.wo.donors) <- n_wo.donors.col.names.filt

# To create different files with gene expression data
write.table(array.normalized.donor.1, "donor_1_expression_normalized.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.donor.1, "donor_1_expression_normalized_filter.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.donor.2, "donor_2_expression_normalized.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.donor.2, "donor_2_expression_normalized_filter.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.wo.donors, "wo_donors_expression_normalized.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.wo.donors, "wo_donors_expression_normalized_filter.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

###########################################################################################################################################
# To perform analysis using lumi package

# To input files for lumi
# To read data for lumi
data.gene.lumi <- lumiR("wo_control_expression_lumi.txt", sep = NULL, detectionTh = 0.01, na.rm = TRUE, convertNuID = FALSE, lib.mapping = NULL, dec = '.', parseColumnName = FALSE, checkDupId = TRUE, QC = TRUE, columnNameGrepPattern = list(exprs = 'AVG_Signal', se.exprs = 'BEAD_STDEV', detection = 'Detection Pval', beadNum = 'Avg_NBEADS'), inputAnnotation = TRUE, annotationColumn = c('Entrez_ID'))
controlFile.lumi <- "control_expression_lumi.txt"
data.control.lumi <- addControlData2lumi(controlFile.lumi, data.gene.lumi)

# To run a for loop to generate different plots
# To create a density plot
pdf(file = "all_density_before_normalization_lumi.pdf", width = 11, height = 8.5)

# To define i.lumi for looping
i.lumi = 1

while(i.lumi < dim(as.matrix(data.gene.lumi))[2]) {

  # To generate box plots for 10 arrays at a time
  # To get an integer value of number of arrays
  i.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

  # To compare if the arrays the last 10 arrays or not
  if (i.lumi < (i.arrays.lumi*10)) {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plot(data.gene.lumi[, i.lumi:(10*i.lumi)], what = "density")

  }

  else {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plot(data.gene.lumi[, i.lumi:(dim(as.matrix(data.gene.lumi))[2])], what = "density")

  }

  # To add 10 to a while loop
  i.lumi = i.lumi + 10

}

garbage <- dev.off()


# To run a for loop to generate different plots
# To create a Cumulative Distribution Function (CDF) plot
pdf(file = "all_cdf_before_normalization_lumi.pdf", width = 11, height = 8.5)

# To define j.lumi for looping
j.lumi = 1

while(j.lumi < dim(as.matrix(data.gene.lumi))[2]) {

  # To generate box plots for 10 arrays at a time
  # To get an integer value of number of arrays
  j.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

  # To compare if the arrays the last 10 arrays or not
  if (j.lumi < (j.arrays.lumi*10)) {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plotCDF(data.gene.lumi[, j.lumi:(10*j.lumi)], reverse = TRUE)

  }

  else {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plotCDF(data.gene.lumi[, j.lumi:(dim(as.matrix(data.gene.lumi))[2])], reverse = TRUE)

  }

  # To add 10 to a while loop
  j.lumi = j.lumi + 10

}

garbage <- dev.off()


# To run a for loop to generate different plots
# To create a housekeeping plot
pdf(file = "all_housekeeping_before_normalization_lumi.pdf", width = 11, height = 8.5)

# To define k.lumi for looping
k.lumi = 1

while(k.lumi < dim(as.matrix(data.gene.lumi))[2]) {

  # To generate box plots for 10 arrays at a time
  # To get an integer value of number of arrays
  k.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

  # To compare if the arrays the last 10 arrays or not
  if (k.lumi < (k.arrays.lumi*10)) {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plotHousekeepingGene(data.control.lumi[, k.lumi:(10*k.lumi)])

  }

  else {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plotHousekeepingGene(data.control.lumi[, k.lumi:(dim(as.matrix(data.gene.lumi))[2])])

  }

  # To add 10 to a while loop
  k.lumi = k.lumi + 10

}

garbage <- dev.off()



# To run a for loop to generate different plots
# To create a stringency gene plot
pdf(file = "all_stringency_gene_before_normalization_lumi.pdf", width = 11, height = 8.5)

# To define l.lumi for looping
l.lumi = 1

while(l.lumi < dim(as.matrix(data.gene.lumi))[2]) {

  # To generate box plots for 10 arrays at a time
  # To get an integer value of number of arrays
  l.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

  # To compare if the arrays the last 10 arrays or not
  if (l.lumi < (l.arrays.lumi*10)) {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plotStringencyGene(data.control.lumi[, l.lumi:(10*l.lumi)])

  }

  else {

    # To create a box plots for the regular probes and negative control probes before normalization and the regular probes after normalization
    plotStringencyGene(data.control.lumi[, l.lumi:(dim(as.matrix(data.gene.lumi))[2])])

  }

  # To add 10 to a while loop
  l.lumi = l.lumi + 10

}

garbage <- dev.off()

# To create a sample relations plot
pdf(file = "all_sample_relations_before_normalization_lumi.pdf", width = 11, height = 8.5)
plot(data.gene.lumi, what = 'sampleRelation')
garbage <- dev.off()

pdf(file = "all_sample_relations_mds_before_normalization_lumi.pdf", width = 11, height = 8.5)
plot(data.gene.lumi, what = 'sampleRelation', method = "mds")
garbage <- dev.off()
