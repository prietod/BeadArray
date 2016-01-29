# This script is written by Dr. Hemang Parikh as on January 26, 2016
# The Health Informatics Institute (HII) at the University of South Florida

# Median local background is performed
# BASH and HULK methods are used for beads artifact detection
# Summarization is performed with un-log transformation because neqc requires gene expression values to be un-logged
# Outliers are removed using the Illumina 3 M.A.D cut-off
# The mean and standard deviation for each bead type are reported

# In beadarray
# No function for variance stabilizing transformation

# In lumi
# No background normalization
# Variance stabilizing transformation
# Quantile normalization

# To import beadarray, lumi and illuminaHumanv4.db libraries
suppressMessages(library(beadarray))
suppressMessages(library(illuminaHumanv4.db))
suppressMessages(library(lumi))
suppressMessages(library(limma))

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

###########################################################################################################################################
# To perform analysis using lumi package

# To input files for lumi
# To read data for lumi
data.gene.lumi <- lumiR("wo_control_expression_lumi.txt", sep = "\t", detectionTh = 0.01, na.rm = TRUE, convertNuID = FALSE, lib.mapping = NULL, dec = '.', parseColumnName = FALSE, checkDupId = TRUE, QC = TRUE, columnNameGrepPattern = list(exprs = 'AVG_Signal', se.exprs = 'BEAD_STDEV', detection = 'Detection Pval', beadNum = 'Avg_NBEADS'), inputAnnotation = TRUE, annotationColumn = c('Entrez_ID'))
controlFile.lumi <- "control_expression_lumi.txt"
data.control.lumi <- addControlData2lumi(controlFile.lumi, data.gene.lumi)

# To perform variance stabilizing transformation
data.gene.lumi.t <- suppressWarnings(lumiT(data.gene.lumi, method = "vst", verbose = FALSE))

# Perform no background normalization and quantile normalization
data.gene.lumi.n <- lumiN(data.gene.lumi.t, method = "quantile", verbose = FALSE)

# To perform quality control estimation after normalization
data.gene.lumi.n.q <- lumiQ(data.gene.lumi.n)

# To annotate Illumina IDs
idsTosymbols.lumi = as.matrix(toTable(illuminaHumanv4ENTREZID))

# To store gene expression values for all the array
array.normalized.lumi <- matrix(, nrow = dim(as.matrix(data.gene.lumi.n.q))[1], ncol = (dim(as.matrix(data.gene.lumi.n.q))[2] + 2))

# To store Illumina IDs information
array.normalized.lumi[, 1] = row.names(data.gene.lumi.n.q)

# To get annotation for each Illumina ID
for (p.id.lumi in 1:dim(as.matrix(array.normalized.lumi))[1]) {

  # To get the gene symbol for each Illumina ID
  p.symbol.lumi = which(as.character(idsTosymbols.lumi[, 1]) == as.character(as.matrix(array.normalized.lumi)[p.id.lumi, 1]))

  # To store Illumina IDs annotation
  if (length(p.symbol.lumi) == 0) {
    array.normalized.lumi[p.id.lumi, 2] = ""
  }

  else {
    array.normalized.lumi[p.id.lumi, 2] = as.character(idsTosymbols.lumi[p.symbol.lumi, 2])
  }
}

# To define column names
colnames(array.normalized.lumi) <- c("Probe_ID", "Entrez_ID", colnames(data.gene.lumi.n.q))

# To retrieve quality information and verify that probes annotated as "Bad" or "No Match" generally have lower signal
ids.arrays.lumi <- as.character(rownames(data.gene.lumi.n.q))
qual.lumi <- unlist(mget(ids.arrays.lumi, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))

# To store normalized gene expression data for only high quality probes
rem.lumi <- qual.lumi == "No match" | qual.lumi == "Bad"
data.gene.lumi.n.q.filt <- data.gene.lumi.n.q[!rem.lumi, ]

# To store gene expression values for all the array
array.normalized.filt.lumi <- matrix(, nrow = dim(as.matrix(data.gene.lumi.n.q.filt))[1], ncol = (dim(as.matrix(data.gene.lumi.n.q.filt))[2] + 2))

# To store Illumina IDs information
array.normalized.filt.lumi[, 1] = row.names(data.gene.lumi.n.q.filt)

# To get annotation for each Illumina ID
for (q.id.lumi in 1:dim(as.matrix(array.normalized.filt.lumi))[1]) {

  # To get the gene symbol for each Illumina ID
  q.symbol.lumi = which(as.character(idsTosymbols.lumi[, 1]) == as.character(as.matrix(array.normalized.filt.lumi)[q.id.lumi, 1]))

  # To store Illumina IDs annotation
  if (length(q.symbol.lumi) == 0) {
    array.normalized.filt.lumi[q.id.lumi, 2] = ""
  }

  else {
    array.normalized.filt.lumi[q.id.lumi, 2] = as.character(idsTosymbols.lumi[q.symbol.lumi, 2])
  }
}

# To define column names
colnames(array.normalized.filt.lumi) <- c("Probe_ID", "Entrez_ID", colnames(data.gene.lumi.n.q.filt))

# To gene expression values
array.normalized.lumi[, 3:dim(as.matrix(array.normalized.lumi))[2]] = exprs(data.gene.lumi.n.q)[, 1:dim(as.matrix(data.gene.lumi.n.q))[2]]
array.normalized.filt.lumi[, 3:dim(as.matrix(array.normalized.filt.lumi))[2]] = exprs(data.gene.lumi.n.q.filt)[, 1:dim(as.matrix(data.gene.lumi.n.q.filt))[2]]

# To create different files with gene expression data
write.table(array.normalized.lumi, "expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.lumi, "expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# To filter out samples based on following
donor.1.lumi <- data.frame(NULL)
n_donor.1.lumi = 0

# To get the donors information
# To run a for loop for all the samples
for (donor.1.info.lumi in 1:dim(as.matrix(data.gene.lumi.n.q))[2]) {

  # To filter out samples based on donor information
  if ((!is.na(as.character(data.pheno$Donor_Number[donor.1.info.lumi]))) & ((as.character(data.pheno$Donor_Number[donor.1.info.lumi])) == as.character("Donor 1"))) {
      n_donor.1.lumi = n_donor.1.lumi + 1
      donor.1.lumi[n_donor.1.lumi, 1] = as.character(as.matrix(data.pheno)[donor.1.info.lumi, 1])
  }
}


donor.2.lumi <- data.frame(NULL)
n_donor.2.lumi = 0

# To get the donors information
# To run a for loop for all the samples
for (donor.2.info.lumi in 1:dim(as.matrix(data.gene.lumi.n.q))[2]) {

  # To filter out samples based on donor information
  if ((!is.na(as.character(data.pheno$Donor_Number[donor.2.info.lumi]))) & ((as.character(data.pheno$Donor_Number[donor.2.info.lumi])) == as.character("Donor 2"))) {
    n_donor.2.lumi = n_donor.2.lumi + 1
    donor.2.lumi[n_donor.2.lumi, 1] = as.character(as.matrix(data.pheno)[donor.2.info.lumi, 1])
  }
}

# To store gene expression values for donors and other samples
array.normalized.lumi.donor.1 <- matrix(, nrow = dim(as.matrix(array.normalized.lumi))[1], ncol = (dim(donor.1.lumi)[1]) + 2)
array.normalized.filt.lumi.donor.1 <- matrix(, nrow = dim(as.matrix(array.normalized.filt.lumi))[1], ncol = (dim(donor.1.lumi)[1]) + 2)
array.normalized.lumi.donor.2 <- matrix(, nrow = dim(as.matrix(array.normalized.lumi))[1], ncol = (dim(donor.2.lumi)[1]) + 2)
array.normalized.filt.lumi.donor.2 <- matrix(, nrow = dim(as.matrix(array.normalized.filt.lumi))[1], ncol = (dim(donor.2.lumi)[1]) + 2)
array.normalized.lumi.wo.donors <- matrix(, nrow = dim(as.matrix(array.normalized.lumi))[1], ncol = (dim(as.matrix(array.normalized.lumi))[2] - ((dim(donor.1.lumi)[1]) + (dim(donor.2.lumi)[1]))))
array.normalized.filt.lumi.wo.donors <- matrix(, nrow = dim(as.matrix(array.normalized.filt.lumi))[1], ncol = (dim(as.matrix(array.normalized.filt.lumi))[2] - ((dim(donor.1.lumi)[1]) + (dim(donor.2.lumi)[1]))))

array.normalized.lumi.donor.1[, 1:2] <- array.normalized.lumi[, 1:2]
array.normalized.filt.lumi.donor.1[, 1:2] <- array.normalized.filt.lumi[, 1:2]
array.normalized.lumi.donor.2[, 1:2] <- array.normalized.lumi[, 1:2]
array.normalized.filt.lumi.donor.2[, 1:2] <- array.normalized.filt.lumi[, 1:2]
array.normalized.lumi.wo.donors[, 1:2] <- array.normalized.lumi[, 1:2]
array.normalized.filt.lumi.wo.donors[, 1:2] <- array.normalized.filt.lumi[, 1:2]

# To give a column numbers
n_donor.1.col.lumi = 2
n_donor.2.col.lumi = 2
n_wo.donors.col.lumi = 2
n_donor.1.col.filt.lumi = 2
n_donor.2.col.filt.lumi = 2
n_wo.donors.col.filt.lumi = 2

col.names.header.lumi = c("Probe_ID", "Entrez_ID")
n_donor.1.col.names.lumi = col.names.header.lumi
n_donor.2.col.names.lumi = col.names.header.lumi
n_wo.donors.col.names.lumi = col.names.header.lumi
n_donor.1.col.names.filt.lumi = col.names.header.lumi
n_donor.2.col.names.filt.lumi = col.names.header.lumi
n_wo.donors.col.names.filt.lumi = col.names.header.lumi


for (donor.data.lumi in 1:dim(as.matrix(data.gene.lumi.n.q))[2]) {

  # To filter out samples based on donor information
  if (as.character(colnames(data.gene.lumi.n.q)[donor.data.lumi]) %in% as.character(as.matrix(donor.1.lumi))) {
    n_donor.1.col.lumi = n_donor.1.col.lumi + 1
    array.normalized.lumi.donor.1[, n_donor.1.col.lumi] = array.normalized.lumi[, (2 + donor.data.lumi)]
    n_donor.1.col.names.lumi = c(n_donor.1.col.names.lumi, as.character(colnames(data.gene.lumi.n.q)[donor.data.lumi]))
  }

  else if (as.character(colnames(data.gene.lumi.n.q)[donor.data.lumi]) %in% as.character(as.matrix(donor.2.lumi))) {
    n_donor.2.col.lumi = n_donor.2.col.lumi + 1
    array.normalized.lumi.donor.2[, n_donor.2.col.lumi] = array.normalized.lumi[, (2 + donor.data.lumi)]
    n_donor.2.col.names.lumi = c(n_donor.2.col.names.lumi, as.character(colnames(data.gene.lumi.n.q)[donor.data.lumi]))
  }

  else {
    n_wo.donors.col.lumi = n_wo.donors.col.lumi + 1
    array.normalized.lumi.wo.donors[, n_wo.donors.col.lumi] = array.normalized.lumi[, (2 + donor.data.lumi)]
    n_wo.donors.col.names.lumi = c(n_wo.donors.col.names.lumi, as.character(colnames(data.gene.lumi.n.q)[donor.data.lumi]))
  }

}

for (donor.data.filt.lumi in 1:dim(as.matrix(data.gene.lumi.n.q.filt))[2]) {

  # To filter out samples based on donor information
  if (as.character(colnames(data.gene.lumi.n.q.filt)[donor.data.filt.lumi]) %in% as.character(as.matrix(donor.1.lumi))) {
    n_donor.1.col.filt.lumi = n_donor.1.col.filt.lumi + 1
    array.normalized.filt.lumi.donor.1[, n_donor.1.col.filt.lumi] = array.normalized.filt.lumi[, (2 + donor.data.filt.lumi)]
    n_donor.1.col.names.filt.lumi = c(n_donor.1.col.names.filt.lumi, as.character(colnames(data.gene.lumi.n.q.filt)[donor.data.filt.lumi]))
  }

  else if (as.character(colnames(data.gene.lumi.n.q.filt)[donor.data.filt.lumi]) %in% as.character(as.matrix(donor.2.lumi))) {
    n_donor.2.col.filt.lumi = n_donor.2.col.filt.lumi + 1
    array.normalized.filt.lumi.donor.2[, n_donor.2.col.filt.lumi] = array.normalized.filt.lumi[, (2 + donor.data.filt.lumi)]
    n_donor.2.col.names.filt.lumi = c(n_donor.2.col.names.filt.lumi, as.character(colnames(data.gene.lumi.n.q.filt)[donor.data.filt.lumi]))
  }

  else {
    n_wo.donors.col.filt.lumi = n_wo.donors.col.filt.lumi + 1
    array.normalized.filt.lumi.wo.donors[, n_wo.donors.col.filt.lumi] = array.normalized.filt.lumi[, (2 + donor.data.filt.lumi)]
    n_wo.donors.col.names.filt.lumi = c(n_wo.donors.col.names.filt.lumi, as.character(colnames(data.gene.lumi.n.q.filt)[donor.data.filt.lumi]))
  }

}

# To define column names
colnames(array.normalized.lumi.donor.1) <- n_donor.1.col.names.lumi
colnames(array.normalized.lumi.donor.2) <- n_donor.2.col.names.lumi
colnames(array.normalized.lumi.wo.donors) <- n_wo.donors.col.names.lumi
colnames(array.normalized.filt.lumi.donor.1) <- n_donor.1.col.names.filt.lumi
colnames(array.normalized.filt.lumi.donor.2) <- n_donor.2.col.names.filt.lumi
colnames(array.normalized.filt.lumi.wo.donors) <- n_wo.donors.col.names.filt.lumi


# To create different files with gene expression data
write.table(array.normalized.lumi.donor.1, "donor_1_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.lumi.donor.1, "donor_1_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.lumi.donor.2, "donor_2_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.lumi.donor.2, "donor_2_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.lumi.wo.donors, "wo_donors_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.normalized.filt.lumi.wo.donors, "wo_donors_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
