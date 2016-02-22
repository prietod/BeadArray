# Read from tab-delimited Donor_data.txt
donor.data <- as.matrix(read.table("Donor_1_Spearman.txt", header = TRUE, sep = "\t"))

# Compute the largest y value used in the data
max.y <- as.numeric(max(donor.data[, 3:5]))

# Define colors to be used for 0.90, 0.80 and 0.70 r cut-offs
plot.colors <- c("blue", "red", "forestgreen")

# To make the correlation coefficient plots
pdf(file = "Spearman_correlation_coefficient_donor_1.pdf", width = 11, height = 8.5)

# Graph autos using y axis that ranges from 0 to max.y
# Turn off axes and annotations (axis labels) so it can be specified
plot(as.numeric(donor.data[, "Spearman.90.r.percentage"]), type = "o", col = plot.colors[1], ylim = c(0.7, max.y), axes = FALSE, ann = FALSE)

# Make x axis using methods labels
axis(1, at = 1:16, lab = letters[1:16])

# Make y axis with horizontal labels that display ticks at every 4 marks
axis(2, las = 1, at = c(0.7, 0.8, 0.9, 1))

# Create box around plot
box()

# Graph values with red dashed line and square points
lines(as.numeric(donor.data[, "Spearman.85.r.percentage"]), type = "o", pch = 22, lty = 2, col = plot.colors[2])

# Graph values with green dotted line and diamond points
lines(as.numeric(donor.data[, "Spearman.80.r.percentage"]), type = "o", pch = 23, lty = 3, col = plot.colors[3])

# Label the x and y axes with dark green text
title(xlab = "Different methods")
title(ylab = "Proportions of pairwise correlation pairs")

# Create a legend at (1, max.y) that is slightly smaller (cex) and uses the same line colors and points used by the actual plots
legend(0.7, 0.75, c("Spearman: r cutoff = 0.90", "Spearman: r cutoff = 0.85", "Spearman: r cutoff = 0.80"), cex = 0.8, col = plot.colors, pch = 21:23, lty = 1:3)

# Turn off device driver (to flush output to pdf)
dev.off()