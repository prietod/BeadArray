# This script is written by Dr. Hemang Parikh as on February 04, 2016
# The Health Informatics Institute (HII) at the University of South Florida, Tampa, FL

# To make a list of exclude samples based on cutoff from different QC variables

# To load library
# library("doBy")
# library("pracma")


#------------------------------------------------------------------------
# Get cli args and assign appropriate variables
#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# To create a directory for storing pdf files
data_dir <- toString(args[1])
raw_qc_file <- toString(args[2])

# To change the directory
setwd(data_dir)

# To filter out samples based on the different QC variables
# To read the qc data from a text file
array.qc.details <- read.table(file = raw_qc_file, sep = "\t", header = TRUE)

# To run a for loop to generate different plots
# To create a box plots for different QC variables for each column
pdf(file = "all_raw_qc_details_boxplot.pdf", width = 11, height = 8.5)

for (m in 2:dim(as.matrix(array.qc.details))[2]) {
  boxplot(array.qc.details[, m], main = "BoxPlot", ylab = names(array.qc.details)[m])
  stripchart(array.qc.details[, m], vertical = TRUE, add = TRUE, pch = 21, col = "blue")
}
garbage <- dev.off()

# To create a Q-Q plots for different QC variables for each column
pdf(file = "all_raw_qc_details_qqplot.pdf", width = 11, height = 8.5)

for (n in 2:dim(as.matrix(array.qc.details))[2]) {
  qqnorm(array.qc.details[, n], main = "Normal Q-Q Plot", ylab = names(array.qc.details)[n])
}
garbage <- dev.off()

# To create a histogram plots for different QC variables for each column
pdf(file = "all_raw_qc_details_histogram.pdf", width = 11, height = 8.5)

for (o in 2:dim(as.matrix(array.qc.details))[2]) {
  hist(array.qc.details[, o], main = "Histogram", breaks = 15, border = "black", col = "skyblue", ylab = names(array.qc.details)[o])
  lines(density(array.qc.details[, o]))
}
garbage <- dev.off()

# To create a density plots for different QC variables for each column
pdf(file = "all_raw_qc_details_density.pdf", width = 11, height = 8.5)

for (p in 2:dim(as.matrix(array.qc.details))[2]) {
  plot(density(array.qc.details[, p]), main = "Density", ylab = names(array.qc.details)[p])
}
garbage <- dev.off()

# To create a box plots for gender specific probes (three probes are not that reliable)
# Gender specific probes are not useful for QC as there are only three probes and variability is very high with each of the group
# pdf(file = "all_raw_qc_details_gender_boxPlot.pdf", width = 11, height = 8.5)

# boxplot(array.qc.details$Gender_log2_mean ~ array.qc.details$Sex, main = "BoxPlot", ylab = "Gender_log2_mean")
# stripchart(array.qc.details$Gender_log2_mean ~ array.qc.details$Sex, vertical = TRUE, add = TRUE, pch = 21, col = "blue")
# garbage <- dev.off()

# To obtain median values of male and female gender specific probes (three probes are not that reliable)
# gender.median = as.matrix(summaryBy(array.qc.details$Gender_log2_mean ~ array.qc.details$Sex, data = array.qc.details, FUN = list(quantile), na.rm=TRUE))

# if (gender.median[1, 1] == "Female") {
#   female.cutoff = gender.median[1, 5]
# }

# if (gender.median[2, 1] == "Male") {
#   male.cutoff = gender.median[2, 3]
# }

# To filter out samples based on following
array.qc.filter <- data.frame(NULL, NULL, NULL)
n_filter_count = 1
n_pass_count = 1

# To run a for loop for all the samples
for (q in 1:dim(as.matrix(array.qc.details))[1]) {

  # To filter out samples based on having percentage of detection probes values < 0.20
  if (as.numeric(array.qc.details$Percent_detect_probes[q]) < 0.20) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Low percentage of detected probes"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Percent_detect_probes[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on having percentage of Housekeeping genes expressed above the background level of the array < 80
  else if (as.numeric(array.qc.details$Housekeeping_percentage[q]) < 80) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Low percentage of Housekeeping genes expressed above the background level of the array"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Housekeeping_percentage[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on having percentage of Biotin genes expressed above the background level of the array < 80
  else if (as.numeric(array.qc.details$Biotin_percentage[q]) < 80) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Low percentage of Biotin genes expressed above the background level of the array"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Biotin_percentage[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on having high percentage of beads masked by the BASH > 0.25
  else if (as.numeric(array.qc.details$Bash_percentage_masked[q]) > 0.25) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "High percentage of beads are masked by the BASH method"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Bash_percentage_masked[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on having a significant gradient effect in the intensities. It gives an indication of the level of variability across the entire surface of the chip. Bash_extended_score > 0.35
  else if (as.numeric(array.qc.details$Bash_extended_score[q]) > 0.35) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Having a significant gradient effects in the intensities"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Bash_extended_score[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on median numbers of beads used to create the summary values < 16
  else if (as.numeric(array.qc.details$Nos_beads_50_percentile[q]) < 16) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Low median number of beads used to create the summary values for each probe on each array after outliers removal"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Nos_beads_50_percentile[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on mean numbers of beads used to create the summary values < 16
  else if (as.numeric(array.qc.details$Nos_beads_mean[q]) < 16) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Low mean number of beads used to create the summary values for each probe on each array after outliers removal"
    array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Nos_beads_mean[q])
    n_filter_count = n_filter_count + 1
  }

  # To filter out samples based on total numbers of beads < 100000
  else if (2^(as.numeric(array.qc.details$Nos_beads_log2[q])) < 100000) {
    array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
    array.qc.filter[n_filter_count, 2] = "Low total number of beads on the array"
    array.qc.filter[n_filter_count, 3] = as.numeric(2^(as.numeric(array.qc.details$Nos_beads_log2[q])))
    n_filter_count = n_filter_count + 1
  }

  # Gender specific probes are not useful for QC as there are only three probes and variability is very high with each of the group
  # To filter out samples based on female cutoff > female.cutoff
#   else if ((!is.na(as.character(array.qc.details$Sex[q]))) & (strcmp(as.character(array.qc.details$Sex[q]), "Female")) & (as.numeric(array.qc.details$Gender_log2_mean[q]) > female.cutoff)) {
#     array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
#     array.qc.filter[n_filter_count, 2] = "High expression of gender probes for a female"
#     array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Gender_log2_mean[q])
#     n_filter_count = n_filter_count + 1
#   }

  # To filter out samples based on male cutoff < male.cutoff
#   else if ((!is.na(as.character(array.qc.details$Sex[q]))) & (strcmp(as.character(array.qc.details$Sex[q]), "Male")) & (as.numeric(array.qc.details$Gender_log2_mean[q]) < male.cutoff)) {
#     array.qc.filter[n_filter_count, 1] = array.qc.details$Array_ID[q]
#     array.qc.filter[n_filter_count, 2] = "Low expression of gender probes for a male"
#     array.qc.filter[n_filter_count, 3] = as.numeric(array.qc.details$Gender_log2_mean[q])
#     n_filter_count = n_filter_count + 1
#   }

  # To pass to next sample; to calculate total numbers of passed samples
  else {
    n_pass_count = n_pass_count + 1
  }

}

# To write a table with samples IDs to be removed with quality information
write.table(array.qc.filter, file = "exclude_sample_list.txt", append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = c("Array_ID", "Quality_info", "Quality_value"))
