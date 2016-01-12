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
}
