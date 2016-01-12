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

#------------------------------------------------------------------------
# Get cli args and assign appropriate variables
#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# To create two directories within for storing qc data
raw_data_dir <- toString(args[1])
raw_qc_dir <- toString(args[2])

# dir.create(paste(getwd(), "/raw_qc/", sep = ""), showWarnings = TRUE, recursive = FALSE, mode = "0777")

# To read the data using readIllumina(), which will extract the intensities from the .txt, .locs and also TIFF files
# Annotation from https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html
data <- readIllumina(dir = raw_data_dir, useImages = TRUE, illuminaAnnotation = "Humanv4")

# To store array names and numbers of beads
array.names = as.matrix(sectionNames(data))
array.numBeads = as.matrix(numBeads(data))

# To list all the TIFF files
tiff.Files = list.files(path=raw_data_dir, pattern = "*.tif")

# To store green intensities for all the array which will be used for the box plot
array.GreenInt <- matrix(, nrow = max(array.numBeads), ncol = length(array.names))

# To read each of the TIFF file separately
for (i in 1:length(array.names)) {

  # To read each of the TIFF file
  TIFF <- readTIFF(file.path(raw_data_dir, tiff.Files[i]))

  # To find numbers of pixels of value zero in the image that must be an imaging artifact rather than a true measure of intensity for that location
  # cbind (col(TIFF)[which(TIFF == 0)], row(TIFF)[which(TIFF == 0)])
  xcoords <- getBeadData(data, array = i, what = "GrnX")
  ycoords <- getBeadData(data, array = i, what = "GrnY")

  # To calculate a robust measure of background for each array using median of the five lowest pixel values
  Brob <- medianBackground(TIFF, cbind(xcoords, ycoords))
  data <- insertBeadData(data, array = i, what = "GrnRB", Brob)

  # To calculate foreground values in the normal way
  TIFF2 <- illuminaSharpen(TIFF)

  # To calculate foreground values
  IllF <- illuminaForeground(TIFF2, cbind(xcoords, ycoords))
  data <- insertBeadData(data, array = i, what = "GrnF", IllF)

  # To subtract the median background values to get locally background corrected intensities
  data <- backgroundCorrectSingleSection(data, array = i, fg = "GrnF", bg = "GrnRB", newName = "GrnR")

  # To compare the Illumina intensity with the robust intensity, plot the locations of beads whose expressions change substantially and overlay the locations of the implausibly low-intensity pixels in red
  oldG <- getBeadData(data, array = i, "Grn")
  newG <- getBeadData(data, array = i, "GrnR")

  # To change the directory to QC
  setwd(raw_qc_dir)

  # To save image intensities and outliers in PDF format
  pdf(file = paste(sapply(strsplit(tiff.Files[i], ".tif"), "[", 1), "image_low_intensity.pdf", sep = "_"), width = 11, height = 8.5)
  par(mfrow = c(1, 2))
  plot(xcoords[(abs(oldG-newG) > 50)], ycoords[(abs(oldG-newG) > 50)], pch = 16, xlab = "X", ylab = "Y", main = "entire array")
  points(col(TIFF)[TIFF < 400], row(TIFF)[TIFF < 400], col = "red", pch = 16)
  plot(xcoords[(abs(oldG-newG) > 50)], ycoords[(abs(oldG-newG) > 50)], pch = 16, xlim = c(1145,1180), ylim = c(15500,15580), xlab = "X", ylab = "Y", main = "zoomed in")
  points(col(TIFF)[TIFF < 400], row(TIFF)[TIFF < 400], col = "red", pch = 16)
  garbage <- dev.off()

  # To make image plots to identify spatial artifacts on the array surface that can occur from mis-handling or scanning problems
  pdf(file = paste(sapply(strsplit(tiff.Files[i], ".tif"), "[", 1), "imageplot.pdf", sep = "_"))
  par(mfrow = c(6, 2))
  par(mar = c(1, 1, 1, 1))
  print(imageplot(data, array = i, low = "lightgreen", high = "darkgreen", zlim = c(4,10), main = sectionNames(data)[i]))
  garbage <- dev.off()

  # To plot the location of outliers on the arrays with the most obvious spatial artifacts and plot their location
  pdf(file = paste(sapply(strsplit(tiff.Files[i], ".tif"), "[", 1), "outliers.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(outlierplot(data, array = i, main = paste(sectionNames(data)[i], "outliers")))
  garbage <- dev.off()

  # To store green channel intensities for all arrays
  if (dim(array.GreenInt)[1] == dim(as.matrix(logGreenChannelTransform(data, i)))[1]) {
    array.GreenInt[, i] = as.matrix(logGreenChannelTransform(data, i))
  } else {
    array.GreenInt[, i] = rbind(as.matrix(logGreenChannelTransform(data, i)), as.matrix(rep("NA", dim(array.GreenInt)[1]-dim(as.matrix(logGreenChannelTransform(data, i)))[1]), 1))
  }

  # To move up one directory and out of qc_details directory
  setwd("..")

}

# To change the directory to QC
setwd(raw_qc_dir)

# To create a box plot
pdf(file = paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "BoxPlot.pdf", sep = "_"), width = 11, height = 8.5)
boxplot(data, transFun = logGreenChannelTransform, col = "green", ylab = expression(log[2](intensity)), las = 2, outline = FALSE, main = sapply(strsplit(tiff.Files[1], "_"), "[", 1))
garbage <- dev.off()

# To store percentile values of log2 of number of beads, GreenChannel intensities percentiles, control probe detection function that returns the percentage of each control type that are significantly expressed above background level for all the array
array.qcvalues = matrix(, nrow = length(array.names), ncol = 62)

# BASH performs three types of artifact detection in the style of the affmetrix-oriented Harsh-light package: Compact analysis identifies large clusters of outliers, where each outlying bead must be an immediate neighbor of another outliers; Diffuse analysis finds regions that contain more outliers than would be anticipated by chance, and Extended analysis looks for chip-wide variation, such as a consistent gradient effect

# To loop all the arrays
for (j in 1:length(array.names)) {

  # To store percentile values of log2 of GreenChannel intensities
  array.qcvalues[j, 1] = log2(array.numBeads[j])
  array.qcvalues[j, 2] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.01, na.rm = TRUE)[[1]])
  array.qcvalues[j, 3] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.05, na.rm = TRUE)[[1]])
  array.qcvalues[j, 4] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.1, na.rm = TRUE)[[1]])
  array.qcvalues[j, 5] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.25, na.rm = TRUE)[[1]])
  array.qcvalues[j, 6] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.5, na.rm = TRUE)[[1]])
  array.qcvalues[j, 7] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.75, na.rm = TRUE)[[1]])
  array.qcvalues[j, 8] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.9, na.rm = TRUE)[[1]])
  array.qcvalues[j, 9] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.95, na.rm = TRUE)[[1]])
  array.qcvalues[j, 10] = suppressWarnings(quantile(as.numeric(array.GreenInt[, j]), 0.99, na.rm = TRUE)[[1]])
  array.qcvalues[j, 11] = suppressWarnings(mean(as.numeric(array.GreenInt[, j]), na.rm = TRUE)[[1]])
  array.qcvalues[j, 12] = suppressWarnings(sd(as.numeric(array.GreenInt[, j]), na.rm = TRUE)[[1]])

  # To run Bash, mask the affected beads and visualize the regions that have been exclude
  BASHoutput <- suppressMessages(BASH(data, array = j))
  data <- setWeights(data, wts = BASHoutput$wts, array = j)

  # To plot the location of regions that have been excluded by BASH
  pdf(file = paste(sapply(strsplit(tiff.Files[j], ".tif"), "[", 1), "bash_exclude.pdf", sep = "_"), width = 11, height = 8.5)
  showArrayMask(data, array = j, override = TRUE)
  garbage <- dev.off()

  # To store BASH QC information; the $wts with 0 indicating that a bead should be masked
  # Number of masked beads
  array.qcvalues[j, 13] = BASHoutput$QC[[1]]

  # Total number of beads
  array.qcvalues[j, 14] = table(getBeadData(data, what = "wts", array = j))[[2]]

  # Ratio of number of masked beads to total number of beads
  array.qcvalues[j, 15] = (BASHoutput$QC[[1]]/table(getBeadData(data, what = "wts", array = j))[[2]])

  # The Extended score returned by BASH is an indication of the level of variability across the entire surface of the chip. If this value is large it may indicate a significant gradient effect in the intensities
  array.qcvalues[j, 16] = BASHoutput$QC[[2]]

  # The Extended score returned by BASH in the previous use case gives an indication of the level of variability across the entire surface of the chip. If this value is large it may indicate a significant gradient effect in the intensities
  # The HULK function can be used to smooth out any gradients that are present
  # HULK uses information about neighboring beads, but rather than mask them out as in BASH, it adjusts the log-intensities by the weighted average of residual values within a local neighborhood
  HULKoutput <- suppressMessages(HULK(data, array = j, transFun = logGreenChannelTransform))
  data <- insertBeadData(data, array = j, data = HULKoutput,  what = "GrnHulk")

}

# To loop all the arrays
for (k in 1:length(array.names)) {

  # Two particular controls on expression arrays are housekeeping and biotin controls, which are expected to be highly expressed in any sample
  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "ControlProbes.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, main = paste(sectionNames(data)[k], "Positive Controls"), ylim = c(4, 15)))
  garbage <- dev.off()

  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "housekeeping.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, positiveControlTags = c("housekeeping"), colList = c("blue"), ylim = c(4, 15), controlProfile = makeControlProfile("Humanv4")))
  garbage <- dev.off()

  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "biotin.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, positiveControlTags = c("biotin"), colList = c("blue"), ylim = c(4, 15), controlProfile = makeControlProfile("Humanv4")))
  garbage <- dev.off()

  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "cy3_hyb.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, positiveControlTags = c("cy3_hyb"), colList = c("blue"), ylim = c(4, 15), controlProfile = makeControlProfile("Humanv4")))
  garbage <- dev.off()

  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "labeling.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, positiveControlTags = c("labeling"), colList = c("blue"), ylim = c(4, 15), controlProfile = makeControlProfile("Humanv4")))
  garbage <- dev.off()

  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "low_stringency_hyb.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, positiveControlTags = c("low_stringency_hyb"), colList = c("blue"), ylim = c(4, 15), controlProfile = makeControlProfile("Humanv4")))
  garbage <- dev.off()

  pdf(file = paste(sapply(strsplit(tiff.Files[k], ".tif"), "[", 1), "negative.pdf", sep = "_"), width = 11, height = 8.5)
  suppressMessages(poscontPlot(data, array = k, positiveControlTags = c("negative"), colList = c("blue"), ylim = c(4, 15), controlProfile = makeControlProfile("Humanv4")))
  garbage <- dev.off()

  # To store control probe detection function that returns the percentage of each control type that are significantly expressed above background level for all the array
  array.qcvalues[k, 17] = suppressMessages(as.matrix(controlProbeDetection(data, array = k, tagsToDetect = c("housekeeping"), negativeTag = "negative"))[1,1])
  array.qcvalues[k, 18] = suppressMessages(as.matrix(controlProbeDetection(data, array = k, tagsToDetect = c("biotin"), negativeTag = "negative"))[1,1])
  array.qcvalues[k, 19] = suppressMessages(as.matrix(controlProbeDetection(data, array = k, tagsToDetect = c("cy3_hyb"), negativeTag = "negative"))[1,1])
  array.qcvalues[k, 20] = suppressMessages(as.matrix(controlProbeDetection(data, array = k, tagsToDetect = c("labeling"), negativeTag = "negative"))[1,1])
  array.qcvalues[k, 21] = suppressMessages(as.matrix(controlProbeDetection(data, array = k, tagsToDetect = c("low_stringency_hyb"), negativeTag = "negative"))[1,1])

}


# Certain probes located on the Y chromosome are beneficial in discriminating the gender of samples in a population; these probes can be incorporated into a QC report
# Make sure that probes are the same between Humanv4 vs. Humanv3
cprof <- suppressMessages(makeControlProfile("Humanv4"))
sexprof <- data.frame("ArrayAddress" = c("5270068", "1400139", "6860102"), "Tag" = rep("Gender", 3))
cprof <- rbind(cprof, sexprof)

# To summarize the control intensities for the each array, then tabulate the mean and standard deviation of all control probes on every array
qcReport <- suppressMessages(makeQCTable(data, controlProfile = cprof))
data <- insertSectionData(data, what = "BeadLevelQC", data=qcReport)

# To store mean and stdev of several types of control probes
array.qcvalues[, 22] = as.matrix(qcReport)[, 1]
array.qcvalues[, 23] = as.matrix(qcReport)[, 8]
array.qcvalues[, 24] = as.matrix(qcReport)[, 2]
array.qcvalues[, 25] = as.matrix(qcReport)[, 9]
array.qcvalues[, 26] = as.matrix(qcReport)[, 3]
array.qcvalues[, 27] = as.matrix(qcReport)[, 10]
array.qcvalues[, 28] = as.matrix(qcReport)[, 4]
array.qcvalues[, 29] = as.matrix(qcReport)[, 11]
array.qcvalues[, 30] = as.matrix(qcReport)[, 5]
array.qcvalues[, 31] = as.matrix(qcReport)[, 12]
array.qcvalues[, 32] = as.matrix(qcReport)[, 6] 
array.qcvalues[, 33] = as.matrix(qcReport)[, 13]
array.qcvalues[, 34] = as.matrix(qcReport)[, 7] 
array.qcvalues[, 35] = as.matrix(qcReport)[, 14]


# All observations are extracted, transformed and then grouped together according to their ArrayAddressID. Outliers are removed and the mean and standard deviation of the remaining beads are calculated
# The default options of summarize apply a log2 transformation, remove outliers using the Illumina 3 M.A.D cut-off and report the mean and standard deviation for each bead type
datasumm <- summarize(BLData = data, useSampleFac = FALSE)

# The detection score, or detection p-value is a standard measure for Illumina expression experiments, and can be viewed as an empirical estimate of the p-value for the null hypothesis that a particular probe in not expressed
# These can be calculated for summarized data provided that the identity of the negative controls on the array is known using the function calculateDetection
array.qcvalues[, 36] = as.matrix(dim(datasumm))[1]
array.qcvalues[, 37] = as.matrix(dim(datasumm))[2]

# To create a box plot of summarized gene expression
pdf(file = paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "summarized_intensity_BoxPlot.pdf", sep = "_"), width = 11, height = 8.5)
boxplot(exprs(datasumm), ylab = expression(log[2](intensity)), las = 2, outline = FALSE)
garbage <- dev.off()

# To create a box plot of number of beads used for the summarization of gene expression
pdf(file = paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "numbers_of_beads_BoxPlot.pdf", sep = "_"), width = 11, height = 8.5)
boxplot(nObservations(datasumm),  ylab = "number of beads", las = 2, outline = FALSE)
garbage <- dev.off()

# The detection score, or detection p-value is a standard measure for Illumina expression experiments, and can be viewed as an empirical estimate of the p-value for the null hypothesis that a particular probe in not expressed
# These can be calculated for summarized data provided that the identity of the negative controls on the array is known using the function calculateDetection
det <- calculateDetection(datasumm)

# To loop all the arrays
for (l in 1:length(array.names)) {

  # To store the percentile of summarized gene expression intensities
  array.qcvalues[l, 38] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.01, na.rm = TRUE)[[1]])
  array.qcvalues[l, 39] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.05, na.rm = TRUE)[[1]])
  array.qcvalues[l, 40] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.1, na.rm = TRUE)[[1]])
  array.qcvalues[l, 41] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.25, na.rm = TRUE)[[1]])
  array.qcvalues[l, 42] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.50, na.rm = TRUE)[[1]])
  array.qcvalues[l, 43] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.75, na.rm = TRUE)[[1]])
  array.qcvalues[l, 44] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.90, na.rm = TRUE)[[1]])
  array.qcvalues[l, 45] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.95, na.rm = TRUE)[[1]])
  array.qcvalues[l, 46] = suppressWarnings(quantile(as.matrix(exprs(datasumm))[, l], 0.99, na.rm = TRUE)[[1]])
  array.qcvalues[l, 47] = suppressWarnings(mean(as.matrix(exprs(datasumm))[, l], na.rm = TRUE))
  array.qcvalues[l, 48] = suppressWarnings(sd(as.matrix(exprs(datasumm))[, l], na.rm = TRUE))

  # To store the percentile of the number of beads used for the summarization of gene expression
  array.qcvalues[l, 49] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.01, na.rm = TRUE)[[1]])
  array.qcvalues[l, 50] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.05, na.rm = TRUE)[[1]])
  array.qcvalues[l, 51] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.1, na.rm = TRUE)[[1]])
  array.qcvalues[l, 52] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.25, na.rm = TRUE)[[1]])
  array.qcvalues[l, 53] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.50, na.rm = TRUE)[[1]])
  array.qcvalues[l, 54] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.75, na.rm = TRUE)[[1]])
  array.qcvalues[l, 55] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.90, na.rm = TRUE)[[1]])
  array.qcvalues[l, 56] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.95, na.rm = TRUE)[[1]])
  array.qcvalues[l, 57] = suppressWarnings(quantile(as.matrix(nObservations(datasumm))[, l], 0.99, na.rm = TRUE)[[1]])
  array.qcvalues[l, 58] = suppressWarnings(mean(as.matrix(nObservations(datasumm))[, l], na.rm = TRUE))
  array.qcvalues[l, 59] = suppressWarnings(sd(as.matrix(nObservations(datasumm))[, l], na.rm = TRUE))

  # To store detection counts information with probes having P-value < 0.01
  # P-value < 0.01 was selected based on lumi package information
  n.count = 0

  for (m in 1:as.matrix(dim(datasumm))[1]) {

    if(!is.na(as.matrix(det)[m, l])) {

      if(as.matrix(det)[m, l] < 0.01) {
        n.count = n.count + 1
      }
    }
  }

  #To get a ratio of detection probes
  array.qcvalues[l, 60] = n.count
  array.qcvalues[l, 61] = as.matrix(dim(datasumm))[1]
  array.qcvalues[l, 62] = n.count/as.matrix(dim(datasumm))[1]

}

# To store detection counts information with probes having P-value < 0.05
Detection(datasumm) <- det

# To write an entire QC table
write.table(array.qcvalues, file = paste(sapply(strsplit(tiff.Files[1], "_"), "[", 1), "raw_qc_details.txt", sep = "_"), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = array.names, col.names = FALSE)

# To move up one directory and out of qc_details directory
setwd("..")
