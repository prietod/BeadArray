# This script is written by Dr. Hemang Parikh as on February 12, 2016
# The Health Informatics Institute (HII) at the University of South Florida, Tampa, FL

# To compare several different normalization methods based on replicated data from donors
# To load libraries
library("psych")
library("corrplot")

#------------------------------------------------------------------------
# Get cli args and assign appropriate variables
#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

result_dir <- toString(args[1])

setwd(result_dir)

#------------------------------------------------------------------------
# Get cli args and assign appropriate variables
#------------------------------------------------------------------------

# To create a function to obtain various statistical values
pairwisecorrelationstat <- function(matrix.data) {

  # To create a matrix to save results
  matrix.data.results <- matrix(, nrow = (((dim(matrix.data)[2]) - 2) * ((dim(matrix.data)[2]) - 3))/2, ncol = 9)

  # To store column names
  matrix.data.cols = colnames(matrix.data)

  # To count row numbers
  n.row.loop = 1

  # To run a loop with a dimension of a data matrix
  for (x.loop in 3:dim(matrix.data)[2]) {

    # To check if it is the last column
    if(x.loop < dim(matrix.data)[2]) {

      # To run a loop with next column for correlation analysis
      for (y.loop in 4:dim(matrix.data)[2]) {

        # To compare only unique pairwise combinations by looking at second column to be larger than first column
        if (y.loop > x.loop) {

          # To save a row name by joining headers
          matrix.data.results[n.row.loop, 1] <- paste(colnames(matrix.data)[x.loop], colnames(matrix.data)[y.loop], sep = ".and.")

          # To perform Spearman correlation
          matrix.data.results[n.row.loop, 2] = cor(as.numeric(matrix.data[, x.loop]), as.numeric(matrix.data[, y.loop]), use = "pairwise.complete.obs", method = c("spearman"))

          # To obtain the Fisher's z' transformation for Spearman
          matrix.data.results[n.row.loop, 3] = fisherz(as.numeric(matrix.data.results[n.row.loop, 2]))

          # To obtain Spearman r square
          matrix.data.results[n.row.loop, 4] = (as.numeric(matrix.data.results[n.row.loop, 2])) ^ 2


          # To perform Pearson correlation
          matrix.data.results[n.row.loop, 5] = cor(as.numeric(matrix.data[, x.loop]), as.numeric(matrix.data[, y.loop]), use = "pairwise.complete.obs", method = c("pearson"))

          # To obtain the Fisher's z' transformation for Pearson
          matrix.data.results[n.row.loop, 6] = fisherz(as.numeric(matrix.data.results[n.row.loop, 5]))

          # To obtain Pearson r square
          matrix.data.results[n.row.loop, 7] = (as.numeric(matrix.data.results[n.row.loop, 5])) ^ 2


          # To perform Kolmogorov-Smirnov test
          matrix.data.results[n.row.loop, 8] = suppressWarnings(as.numeric(ks.test(as.numeric(matrix.data[, x.loop]), as.numeric(matrix.data[, y.loop]), alternative = "two.sided")$statistic))
          matrix.data.results[n.row.loop, 9] = suppressWarnings(as.numeric(ks.test(as.numeric(matrix.data[, x.loop]), as.numeric(matrix.data[, y.loop]), alternative = "two.sided")$p.value))

          # To increase a number for a loop
          n.row.loop = n.row.loop + 1

        }
      }
    }
  }

  # To return result matrix
  return(matrix.data.results)

}

# To read all the .txt files
file.names = list.files(pattern="*.txt")

# To store coefficient of variations
cv.values <- data.frame(NULL)

# To store correlation results
data.row.number = as.matrix(read.table(file = file.names[1], sep = "\t", header = TRUE))
donor.data.cor.results <- matrix(, nrow = (((dim(data.row.number)[2]) - 2) * ((dim(data.row.number)[2]) - 3))/2, ncol = 9*(length(file.names)))

# To define a column header
donor.col.header = character(0)

# To read each of the text file separately
for (i in 1:length(file.names)) {

  # To print the name of a file
  print(file.names[i])

  # To read each of the text file separately
  donor.data <- na.omit(as.matrix(read.table(file = file.names[i], sep = "\t", header = TRUE)))
  infinite.values <- is.infinite(donor.data)
  donor.data[infinite.values] <- 0

  # To make a correlation plots
  donor.data.cor <- matrix(, nrow = dim(donor.data)[1], ncol = (dim(donor.data)[2] - 2))
  donor.data.cor[, 1:dim(donor.data.cor)[2]] <- as.matrix(as.numeric(donor.data[, 3:dim(donor.data)[2]]))
  row.names(donor.data.cor) <- donor.data.cor[, 1]
  data.cor.plots <- cor(donor.data.cor)

  pdf(file = paste("method", letters[i], "correlation_plot.pdf", sep = "_"), width = 11, height = 8.5)
  corrplot.mixed(data.cor.plots, lower = "ellipse", upper = "circle")
  garbage <- dev.off()

  # To run for loop for each gene
  for (j.id in 1:dim(donor.data)[1]) {

    # To calculate the coefficient of variations for each gene from each file
    cv.values[j.id, i] = sd(as.numeric(donor.data[j.id, 3:dim(donor.data)[2]]), na.rm = TRUE) / mean(as.numeric(donor.data[j.id, 3:dim(donor.data)[2]]), na.rm = TRUE)

  }

  # To store correlation results
  donor.data.cor.results[, ((as.numeric(i)*9)-8):(as.numeric(i)*9)] = as.matrix(pairwisecorrelationstat(donor.data))[, 1:9]

  donor.col.header = c(donor.col.header, paste("method", letters[i], "ID", sep = "."))
  donor.col.header = c(donor.col.header, paste("method", letters[i], "r.Spearman", sep = "."), paste("method", letters[i], "Fisher.Z.Spearman", sep = "."), paste("method", letters[i], "rsquare.Spearman", sep = "."))
  donor.col.header = c(donor.col.header, paste("method", letters[i], "r.Pearson", sep = "."), paste("method", letters[i], "Fisher.Z.Pearson", sep = "."), paste("method", letters[i], "rsquare.Pearson", sep = "."))
  donor.col.header = c(donor.col.header, paste("method", letters[i], "D.K-S-test", sep = "."), paste("method", letters[i], "p.K-S-test", sep = "."))

}

# To make the coefficient of variations box plot
pdf(file = "coefficient_variations_boxplot.pdf", width = 11, height = 8.5)
boxplot(cv.values, ylab = "Coefficient of variations", names = letters[1:length(file.names)], xlab = "Different methods", main = "Coefficient of variations between different methods")
garbage <- dev.off()

# To make the coefficient of variations box plot
# DISABLED: pdf(file = "coefficient_variations_boxplot_a-j_m-p.pdf", width = 11, height = 8.5)
# DISABLED: boxplot(cv.values[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16)], ylab = "Coefficient of variations", names = letters[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16)], xlab = "Different methods", main = "Coefficient of variations between different methods")
# DISABLED: garbage <- dev.off()

# To give a header line
colnames(donor.data.cor.results) = donor.col.header

# To make a table with counts with different cut-offs of r
donor.r.cutoffs <- matrix(, nrow = 42, ncol = (length(file.names) + 1))

# To store first column of a r square cutoff matrix
donor.r.cutoffs.rows = character(0)
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Spearman.99.r.count", "Spearman.99.r.number", "Spearman.99.r.percentage", "Spearman.95.r.count", "Spearman.95.r.number", "Spearman.95.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Spearman.90.r.count", "Spearman.90.r.number", "Spearman.90.r.percentage", "Spearman.85.r.count", "Spearman.85.r.number", "Spearman.85.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Spearman.80.r.count", "Spearman.80.r.number", "Spearman.80.r.percentage", "Spearman.70.r.count", "Spearman.70.r.number", "Spearman.70.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Spearman.50.r.count", "Spearman.50.r.number", "Spearman.50.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Pearson.99.r.count", "Pearson.99.r.number", "Pearson.99.r.percentage", "Pearson.95.r.count", "Pearson.95.r.number", "Pearson.95.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Pearson.90.r.count", "Pearson.90.r.number", "Pearson.90.r.percentage", "Pearson.85.r.count", "Pearson.85.r.number", "Pearson.85.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Pearson.80.r.count", "Pearson.80.r.number", "Pearson.80.r.percentage", "Pearson.70.r.count", "Pearson.70.r.number", "Pearson.70.r.percentage")
donor.r.cutoffs.rows = c(donor.r.cutoffs.rows, "Pearson.50.r.count", "Pearson.50.r.number", "Pearson.50.r.percentage")

donor.r.cutoffs[, 1] <- donor.r.cutoffs.rows

donor.r.cutoffs.cols = c("Info")

# To for loop numbers of methods
for (k in 1:length(file.names)) {

  # Header line
  donor.r.cutoffs.cols = c(donor.r.cutoffs.cols, as.character(paste("method", letters[k], sep = ".")))

  # To store different cut-offs of r square results
  donor.r.cutoffs[1, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.99))
  donor.r.cutoffs[2, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[3, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.99))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[4, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.95))
  donor.r.cutoffs[5, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[6, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.95))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[7, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.90))
  donor.r.cutoffs[8, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[9, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.90))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[10, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.85))
  donor.r.cutoffs[11, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[12, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.85))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[13, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.80))
  donor.r.cutoffs[14, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[15, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.80))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[16, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.70))
  donor.r.cutoffs[17, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[18, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.70))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[19, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.50))
  donor.r.cutoffs[20, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[21, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Spearman", sep = "."))]) > 0.50))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[22, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.99))
  donor.r.cutoffs[23, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[24, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.99))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[25, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.95))
  donor.r.cutoffs[26, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[27, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.95))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[28, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.90))
  donor.r.cutoffs[29, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[30, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.90))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[31, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.85))
  donor.r.cutoffs[32, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[33, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.85))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[34, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.80))
  donor.r.cutoffs[35, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[36, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.80))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[37, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.70))
  donor.r.cutoffs[38, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[39, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.70))/dim(donor.data.cor.results)[1]
  donor.r.cutoffs[40, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.50))
  donor.r.cutoffs[41, k+1] = dim(donor.data.cor.results)[1]
  donor.r.cutoffs[42, k+1] = sum((as.numeric(donor.data.cor.results[, as.character(paste("method", letters[k], "r.Pearson", sep = "."))]) > 0.50))/dim(donor.data.cor.results)[1]

}

# To give a header line
colnames(donor.r.cutoffs) = donor.r.cutoffs.cols

# To write output files
write.table(cv.values, "coefficient_variations.cls", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(donor.data.cor.results, "correlation_statistics.cls", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(donor.r.cutoffs, "r_cutoff.cls", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
