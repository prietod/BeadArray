# This script is written by Dr. Hemang Parikh as on January 12, 2016.
# No local background is performed
# Summarization is performed with un-log transformation
# Outliers are removed using the Illumina 3 M.A.D cut-off
# The mean and standard deviation for each bead type are reported


# To import beadarray and illuminaHumanv4.db libraries
suppressMessages(library(beadarray))
suppressMessages(library(illuminaHumanv4.db))

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

# To create two directories within for storing results
filter_data_dir <- toString(args[1])
filter_result_dir <- toString(args[2])


# To read the data using readIllumina(), which will extract the intensities from the .txt, .locs and also TIFF files
# Annotation from https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html
data <- readIllumina(dir = filter_data_dir, useImages = TRUE, illuminaAnnotation = "Humanv4")

# To store array names and numbers of beads
array.names = as.matrix(sectionNames(data))
array.numBeads = as.matrix(numBeads(data))

# To list all the TIFF files
tiff.Files = list.files(path = filter_data_dir, pattern = "*.tif")

# To read each of the TIFF file separately
for (i in 1:length(array.names)) {

  # To read each of the TIFF file
  TIFF <- readTIFF(file.path(filter_data_dir, tiff.Files[i]))

  # To find numbers of pixels of value zero in the image that must be an imaging artifact rather than a true measure of intensity for that location
  # cbind (col(TIFF)[which(TIFF == 0)], row(TIFF)[which(TIFF == 0)])
  xcoords <- getBeadData(data, array = i, what = "GrnX")
  ycoords <- getBeadData(data, array = i, what = "GrnY")

  # To calculate foreground values
  IllF <- illuminaForeground(TIFF, cbind(xcoords, ycoords))
  data <- insertBeadData(data, array = i, what = "GrnF", IllF)

  # To change the directory to result
  setwd(filter_result_dir)

  # To move up one directory and out of result directory
  setwd("..")

}

# To change the directory to result
setwd(filter_result_dir)

# All observations are extracted, transformed and then grouped together according to their ArrayAddressID. Outliers are removed and the mean and standard deviation of the remaining beads are calculated
# The default options of summarize apply a un-log transformation, remove outliers using the Illumina 3 M.A.D cut-off and report the mean and standard deviation for each bead type
grnchannel.unlogged <- new("illuminaChannel", transFun = greenChannelTransform, outlierFun = illuminaOutlierMethod, exprFun = function(x) mean(x, na.rm = TRUE),  varFun = function(x) sd(x, na.rm = TRUE), channelName = "G")
datasumm.unlogged <- summarize(BLData = data, useSampleFac=FALSE, channelList = list(grnchannel.unlogged))


# The detection score, or detection p-value is a standard measure for Illumina expression experiments, and can be viewed as an empirical estimate of the p-value for the null hypothesis that a particular probe in not expressed
# These can be calculated for summarized data provided that the identity of the negative controls on the array is known using the function calculateDetection
det <- calculateDetection(datasumm.unlogged)

# To store detection counts information with probes having P-value < 0.05
Detection(datasumm.unlogged) <- det

# To store gene expression values for all the array
array.datasumm <- matrix(, nrow = dim(as.matrix(exprs(datasumm.unlogged)))[1], ncol = ((4*dim(as.matrix(exprs(datasumm.unlogged)))[2]) + 1))

# To define column names
col.datasumm = "Probe_ID"

# To store Probe ID information
array.datasumm[, 1] = row.names(exprs(datasumm.unlogged))

# To write expression values similar to Genome Studio for each BeadChip
# To read each of the array file separately
for (j in 1:dim(as.matrix(exprs(datasumm.unlogged)))[2]) {

  col.datasumm = c(col.datasumm, paste(sapply(strsplit(tiff.Files[j], "_Grn.tif"), "[", 1), "AVG_Signal", sep = "."), paste(sapply(strsplit(tiff.Files[j], "_Grn.tif"), "[", 1), "BEAD_STDERR", sep = "."), paste(sapply(strsplit(tiff.Files[j], "_Grn.tif"), "[", 1), "Avg_NBEADS", sep = "."), paste(sapply(strsplit(tiff.Files[j], "_Grn.tif"), "[", 1), "Detection Pval", sep = "."))
  array.datasumm[, ((4*j) - 2)] <- exprs(datasumm.unlogged)[, j]
  array.datasumm[, ((4*j) - 1)] <- se.exprs(datasumm.unlogged)[, j]
  array.datasumm[, (4*j)] <- nObservations(datasumm.unlogged)[, j]
  array.datasumm[, ((4*j) + 1)] <- Detection(datasumm.unlogged)[, j]

}

# To define column names
colnames(array.datasumm) <- col.datasumm

# To convert Probe IDs to Illumina IDs
adrToIllumina = as.matrix(toTable(illuminaHumanv4ARRAYADDRESS))

# To get a list of control Probe IDs
control.probe.ids = as.matrix(makeControlProfile("Humanv4", excludeERCC = TRUE))

# To store control Probe IDs gene expression data
array.datasumm.controls <- matrix(, nrow = dim(control.probe.ids)[1], ncol = ((4*dim(as.matrix(exprs(datasumm.unlogged)))[2]) + 1))

# To count the row number for control gene expression array
n.controls = 0

# To loop over each of the control Probe IDs
for (k in 1:dim(control.probe.ids)[1]) {

  # To loop over each of arrayaddress ID to Illumina ID
  for (l in 1:dim(adrToIllumina)[1]) {

    # To match each of the control Probe ID to Illumina ID
    if (as.character(control.probe.ids[k, 1]) == as.character(adrToIllumina[l, 2])) {

      # To match the control Probe ID with gene expression data
      for (m in 1:dim(as.matrix(exprs(datasumm.unlogged)))[1]) {

        if (as.character(adrToIllumina[l, 1]) == as.character(array.datasumm[m, 1])) {

          # To store gene expression data from control Probe IDs
          n.controls = n.controls + 1
          array.datasumm.controls[n.controls, ] = array.datasumm[m, ]

        }
      }
    }
  }
}

# To define column names
colnames(array.datasumm.controls) <- col.datasumm

# To remove duplicated control Probe_ID rows
array.datasumm.controls.unique = as.matrix(unique(array.datasumm.controls))

# To store non-control Probe IDs gene expression data
array.datasumm.wocontrols <- matrix(, nrow = (dim(as.matrix(exprs(datasumm.unlogged)))[1] - dim(array.datasumm.controls.unique)[1]), ncol = ((4*dim(as.matrix(exprs(datasumm.unlogged)))[2]) + 1))

# To count the row number for without control gene expression array
n.wocontrols = 0

# To loop over gene expression array to find non-control probes
for (n in 1:dim(as.matrix(exprs(datasumm.unlogged)))[1]) {

  #To check if the Probe ID is a control Probe ID or not
  if(!is.element(as.character(array.datasumm[n, 1]), as.character(array.datasumm.controls.unique[, 1]))) {

    # To store non-control Probe ID
    n.wocontrols = n.wocontrols + 1
    array.datasumm.wocontrols[n.wocontrols, ] = array.datasumm[n, ]

  }
}

# To define column names
colnames(array.datasumm.controls.unique) <- col.datasumm
colnames(array.datasumm.wocontrols) <- col.datasumm

# To create different files with gene expression data
write.table(array.datasumm, paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "expression.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.datasumm.controls.unique, paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "control_expression.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(array.datasumm.wocontrols, paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "wo_control_expression.txt", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(exprs(datasumm.unlogged), paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "Avg_Signal.txt", sep = "_"), sep="\t", quote = FALSE, col.names = paste(array.names, "AVG_Signal", sep = "."))
write.table(se.exprs(datasumm.unlogged), paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "BEAD_STDERR.txt", sep = "_"), sep="\t", quote = FALSE, col.names = paste(array.names, "BEAD_STDERR", sep = "."))
write.table(nObservations(datasumm.unlogged), paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "Avg_NBEADS.txt", sep = "_"), sep="\t", quote = FALSE, col.names = paste(array.names, "Avg_NBEADS", sep = "."))
write.table(Detection(datasumm.unlogged), paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "Detection_Pval.txt", sep = "_"), sep="\t", quote = FALSE, col.names = paste(array.names, "Detection Pval", sep = "."))

# To move up one directory and out of result directory
setwd("..")
