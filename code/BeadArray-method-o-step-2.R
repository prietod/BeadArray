# This script is written by Dr. Hemang Parikh as on February 04, 2016
# The Health Informatics Institute (HII) at the University of South Florida, Tampa, FL

# Illumina average local background is performed
# BASH and HULK methods are used for beads artifact detection
# Summarization is performed with un-log transformation
# Outliers are removed using the Illumina 3 M.A.D cut-off
# The mean and standard deviation for each bead type are reported for summarization

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

# To enter input parameter for generating figures
generate.figures <- toString(args[2])

# To change the directory
setwd(filter_combined_result_dir)

# dir.create(paste(getwd(), "/filter_results/", sep = ""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

# To read the gene expression data from the *.txt files
# Annotation from https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html
# To read data for beadarray
data.gene <- readBeadSummaryData(dataFile = "wo_control_expression.txt", ProbeID = "Probe_ID", skip = 0, columns = list(exprs = "AVG_Signal", se.exprs = "BEAD_STDEV", nObservations = "Avg_NBEADS", Detection = "Detection Pval"), annoCols = c("Entrez_ID"))
data.pheno <- read.table(file = "BeadArray_phenotype_details.txt", sep = "\t", header = TRUE)

# To get status of genes distributions
# table(data.gene$genes$Status)

# To estimate the proportion of probes which are expressed above the level of the negative controls
# proportion.expr <- propexpr(data.gene)
# write.table(as.matrix(proportion.expr), "probe_proportion_expressed.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = c("Proportion of probes"))

###########################################################################################################################################
# To perform analysis using lumi package

# To input files for lumi
# To read data for lumi
data.gene.lumi <- lumiR("wo_control_expression_lumi.txt", sep = "\t", detectionTh = 0.01, na.rm = TRUE, convertNuID = FALSE, lib.mapping = NULL, dec = '.', parseColumnName = FALSE, checkDupId = TRUE, QC = TRUE, columnNameGrepPattern = list(exprs = 'AVG_Signal', se.exprs = 'BEAD_STDEV', detection = 'Detection Pval', beadNum = 'Avg_NBEADS'), inputAnnotation = TRUE, annotationColumn = c('Entrez_ID'))
controlFile.lumi <- "control_expression_lumi.txt"
data.gene.lumi <- addControlData2lumi((getControlData(controlFile.lumi, sep = "\t")), data.gene.lumi)

# To store control data for plotting
data.control.lumi <- addControlData2lumi((getControlData(controlFile.lumi, sep = "\t")), data.gene.lumi)

# To check the input parameter for generating figures
if (generate.figures == "yes")
{

  # To run a while loop to generate different plots
  # To create a density plot
  pdf(file = "all_density_before_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define i.lumi for looping
  i.lumi = 1

  while(i.lumi <= dim(as.matrix(data.gene.lumi))[2]) {

    # To generate density plots for 10 arrays at a time
    # To get an integer value of number of arrays
    i.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (i.lumi < (i.arrays.lumi*10)) {

      # To create a density plots before normalization
      plot(data.gene.lumi[, i.lumi:(i.lumi + 9)], what = "density")

    }

    else {

      # To create a density plots before normalization
      plot(data.gene.lumi[, i.lumi:(dim(as.matrix(data.gene.lumi))[2])], what = "density")

    }

    # To add 10 to a while loop
    i.lumi = i.lumi + 10

  }

  garbage <- dev.off()


  # To run a while loop to generate different plots
  # To create a Cumulative Distribution Function (CDF) plot
  pdf(file = "all_cdf_before_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define j.lumi for looping
  j.lumi = 1

  while(j.lumi <= dim(as.matrix(data.gene.lumi))[2]) {

    # To generate CDF plots for 10 arrays at a time
    # To get an integer value of number of arrays
    j.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (j.lumi < (j.arrays.lumi*10)) {

      # To create a CDF plots before normalization
      plotCDF(data.gene.lumi[, j.lumi:(j.lumi + 9)], reverse = TRUE)

    }

    else {

      # To create a CDF plots before normalization
      plotCDF(data.gene.lumi[, j.lumi:(dim(as.matrix(data.gene.lumi))[2])], reverse = TRUE)

    }

    # To add 10 to a while loop
    j.lumi = j.lumi + 10

  }

  garbage <- dev.off()


  # To run a while loop to generate different plots
  # To create a housekeeping plot
  pdf(file = "all_housekeeping_before_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define k.lumi for looping
  k.lumi = 1

  while(k.lumi < dim(as.matrix(data.gene.lumi))[2]) {

    # To generate housekeeping plots for 10 arrays at a time
    # To get an integer value of number of arrays
    k.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (k.lumi < (k.arrays.lumi*10)) {

      # To create a housekeeping plots before normalization
      plotHousekeepingGene(data.control.lumi[, k.lumi:(k.lumi + 9)])

    }

    else {

      # To create a housekeeping plots before normalization
      plotHousekeepingGene(data.control.lumi[, k.lumi:(dim(as.matrix(data.gene.lumi))[2])])

    }

    # To add 10 to a while loop
    k.lumi = k.lumi + 10

  }

  garbage <- dev.off()



  # To run a while loop to generate different plots
  # To create a stringency gene plot
  pdf(file = "all_stringency_gene_before_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define l.lumi for looping
  l.lumi = 1

  while(l.lumi < dim(as.matrix(data.gene.lumi))[2]) {

    # To generate stringency gene plots for 10 arrays at a time
    # To get an integer value of number of arrays
    l.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (l.lumi < (l.arrays.lumi*10)) {

      # To create a stringency gene plots before normalization
      plotStringencyGene(data.control.lumi[, l.lumi:(l.lumi + 9)])

    }

    else {

      # To create a stringency gene plots before normalization
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
}

# To perform variance stabilizing transform (VST)
data.gene.lumi.t <- suppressWarnings(lumiT(data.gene.lumi, method = "vst", verbose = FALSE))

# To check the input parameter for generating figures
if (generate.figures == "yes")
{

  # To run a while loop to generate different plots
  # To create a vst transformation plot
  pdf(file = "all_vst_transformation_before_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define m.lumi for looping
  m.lumi = 1

  while(m.lumi <= dim(as.matrix(data.gene.lumi.t))[2]) {

    # To generate a vst transformation plot for 10 arrays at a time
    # To get an integer value of number of arrays
    m.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi.t))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (m.lumi < (m.arrays.lumi*10)) {

      # To create a vst transformation plot before normalization
      plot.trans <- plotVST(data.gene.lumi.t[, m.lumi:(m.lumi + 9)])
      matplot(log2(plot.trans$untransformed), plot.trans$transformed, main = "compare VST and log2 transform")

    }

    else {

      # To create a vst transformation plot before normalization
      plot.trans <- plotVST(data.gene.lumi.t[, m.lumi:(dim(as.matrix(data.gene.lumi.t))[2])])
      matplot(log2(plot.trans$untransformed), plot.trans$transformed, main = "compare VST and log2 transform")

    }

    # To add 10 to a while loop
    m.lumi = m.lumi + 10

  }

  garbage <- dev.off()
}

# Perform no background normalization and quantile normalization
data.gene.lumi.n <- lumiN(data.gene.lumi.t, method = "quantile", verbose = FALSE)

# To perform quality control estimation after normalization
data.gene.lumi.n.q <- lumiQ(data.gene.lumi.n)

# To check the input parameter for generating figures
if (generate.figures == "yes")
{

  # To run a while loop to generate different plots
  # To create a density plot
  pdf(file = "all_density_after_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define n.lumi for looping
  n.lumi = 1

  while(n.lumi <= dim(as.matrix(data.gene.lumi.n.q))[2]) {

    # To generate density plots for 10 arrays at a time
    # To get an integer value of number of arrays
    n.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi.n.q))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (n.lumi < (n.arrays.lumi*10)) {

      # To create a density plots after normalization
      plot(data.gene.lumi.n.q[, n.lumi:(n.lumi + 9)], what = "density")

    }

    else {

      # To create a density plots after normalization
      plot(data.gene.lumi.n.q[, n.lumi:(dim(as.matrix(data.gene.lumi.n.q))[2])], what = "density")

    }

    # To add 10 to a while loop
    n.lumi = n.lumi + 10

  }

  garbage <- dev.off()


  # To run a while loop to generate different plots
  # To create a boxplot plot
  pdf(file = "all_boxplot_after_normalization_lumi.pdf", width = 11, height = 8.5)

  # To define o.lumi for looping
  o.lumi = 1

  while(o.lumi < dim(as.matrix(data.gene.lumi.n.q))[2]) {

    # To generate box plots for 10 arrays at a time
    # To get an integer value of number of arrays
    o.arrays.lumi = as.integer((dim(as.matrix(data.gene.lumi.n.q))[2])/10)

    # To compare if the arrays are the last 10 arrays or not
    if (o.lumi < (o.arrays.lumi*10)) {

      # To create a box plots after normalization
      plot(data.gene.lumi.n.q[, o.lumi:(o.lumi + 9)], what = "boxplot")

    }

    else {

      # To create a box plots after normalization
      plot(data.gene.lumi.n.q[, o.lumi:(dim(as.matrix(data.gene.lumi.n.q))[2])], what = "boxplot")

    }

    # To add 10 to a while loop
    o.lumi = o.lumi + 10

  }

  garbage <- dev.off()

  # To create a sample relations plot
  pdf(file = "all_sample_relations_after_normalization_lumi.pdf", width = 11, height = 8.5)
  plot(data.gene.lumi.n.q, what = 'sampleRelation')
  garbage <- dev.off()

  pdf(file = "all_sample_relations_mds_after_normalization_lumi.pdf", width = 11, height = 8.5)
  plot(data.gene.lumi.n.q, what = 'sampleRelation', method = "mds")
  garbage <- dev.off()
}

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

# To filter out samples based on following
donor.1.lumi <- data.frame(NULL)
n_donor.1.lumi = 0

# To get the donors information
# To run a for loop for all the samples
for (donor.1.info.lumi in 1:dim(as.matrix(data.gene.lumi.n.q))[2]) {

  # To filter out samples based on donor information
  if ((!is.na(as.character(data.pheno$Donor_Number[donor.1.info.lumi]))) & ((as.character(data.pheno$Donor_Number[donor.1.info.lumi])) == as.character("Donor-1"))) {
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
  if ((!is.na(as.character(data.pheno$Donor_Number[donor.2.info.lumi]))) & ((as.character(data.pheno$Donor_Number[donor.2.info.lumi])) == as.character("Donor-2"))) {
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

# To give column numbers
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

# To create files for non-donors' gene expression data
if ((n_donor.1.col.lumi == 2) & (n_donor.2.col.lumi == 2) & (n_donor.1.col.filt.lumi == 2) & (n_donor.2.col.filt.lumi == 2) & (n_wo.donors.col.lumi > 2) & (n_wo.donors.col.filt.lumi > 2)) {
  write.table(array.normalized.lumi.wo.donors, "wo_donors_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi.wo.donors, "wo_donors_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# To create files for donors' gene expression data
if ((n_donor.1.col.lumi > 2) & (n_donor.2.col.lumi > 2) & (n_donor.1.col.filt.lumi > 2) & (n_donor.2.col.filt.lumi > 2) & (n_wo.donors.col.lumi == 2) & (n_wo.donors.col.filt.lumi == 2)) {
  write.table(array.normalized.lumi.donor.1, "donor_1_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi.donor.1, "donor_1_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.lumi.donor.2, "donor_2_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi.donor.2, "donor_2_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# To create all files for gene expression data
if ((n_donor.1.col.lumi > 2) & (n_donor.2.col.lumi > 2) & (n_donor.1.col.filt.lumi > 2) & (n_donor.2.col.filt.lumi > 2) & (n_wo.donors.col.lumi > 2) & (n_wo.donors.col.filt.lumi > 2)) {
  write.table(array.normalized.lumi, "expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi, "expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.lumi.donor.1, "donor_1_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi.donor.1, "donor_1_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.lumi.donor.2, "donor_2_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi.donor.2, "donor_2_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.lumi.wo.donors, "wo_donors_expression_normalized_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(array.normalized.filt.lumi.wo.donors, "wo_donors_expression_normalized_filter_lumi.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
