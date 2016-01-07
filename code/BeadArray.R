#

#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("illuminaHumanv4.db"))

suppressMessages(library(beadarray))
#suppressMessages(library(illuminaHumanv4.db))

# get cli args and assign appropriate variables

args <- commandArgs(trailingOnly = TRUE)

data_dir <- toString(args[1])
results_dir <- toString(args[2])
qc_dir <- toString(args[3])


#comment #dir.create(paste(getwd(),"/results/",sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
#comment #dir.create(paste(getwd(),"/qc_details/",sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
#comment
#comment
#comment #There is a huge issue with sampleSheet so just move everything into one folder and then read it
#comment #chipPath <- "/home/h/hparikh/IL_2_G_9224523070"
#comment #list.files(chipPath)
#comment #sampleSheetFile <- paste(chipPath, "/sampleSheet.csv",sep="")
#comment #readLines(sampleSheetFile)
#comment
#comment #Known issue; there is a problem of using .csv file instead it is better to import all the arrays from a directory
#comment #http://www.rdocumentation.org/packages/beadarray/functions/readIllumina
#comment #https://github.com/Bioconductor-mirror/beadarray/blob/master/R/readIllumina.R
#comment #https://support.bioconductor.org/p/61023/
#comment
#comment
#comment #To read the data using readIllumina(), which will extract the intensities from the .txt, .locs and also TIFF files
#comment #Annotation from https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html
#comment data <- readIllumina(dir=".", useImages=TRUE, illuminaAnnotation="Humanv4")
#comment
#comment
#comment #slotNames(data)
#comment #data@sectionData
#comment Array_Names=as.matrix(sectionNames(data))
#comment Array_NumBeads=as.matrix(numBeads(data))
#comment #head(data[[1]])
#comment #getBeadData(data, array=1, what="Grn")[1:5]
#comment #getBeadData(data, array=1, what="GrnX" )[1:5]
#comment #getBeadData(data, array=1, what="ProbeID")[1:5]
#comment
#comment
#comment #log2(data[[1]][1:10, "Grn"])
#comment #log2(getBeadData(data, array=1, what="Grn")[1:10])
#comment #logGreenChannelTransform(data, array=1)[1:10]
#comment #logGreenChannelTransform
#comment
#comment #To list all the TIFF files
#comment TIFF_Files = list.files(pattern="*.tif")
#comment
#comment #To store green intensities for all the array which is used for the box plot
#comment Array_GreenInt <- matrix(, nrow = max(Array_NumBeads), ncol = length(TIFF_Files))
#comment
#comment #To read each of the TIFF file separately
#comment for (i in 1:length(TIFF_Files))
#comment {
#comment    TIFF <- readTIFF(TIFF_Files[i])
#comment   #cbind (col(TIFF)[which(TIFF==0)], row(TIFF)[which(TIFF==0)])
#comment   xcoords <- getBeadData(data, array=i, what="GrnX")
#comment   ycoords <- getBeadData(data, array=i, what="GrnY")
#comment   #par(mfrow = c(1, 3))
#comment   #offset <- 1
#comment
#comment   #plotTIFF(TIFF+offset, c(1517, 1527), c(5507, 5517), values=T, textCol="yellow", main=expression(log[2](intensity+1)))
#comment   #points(xcoords[503155], ycoords[503155], pch=16, col="red")
#comment   #plotTIFF(TIFF+offset, c(1202, 1212), c(13576, 13586), values=T, textCol="yellow", main=expression(log [2](intensity+1)))
#comment   #points(xcoords[625712], ycoords[625712], pch=16, col="red")
#comment   #plotTIFF(TIFF+offset, c(1613, 1623), c(9219, 9229), values=T, textCol="yellow", main=expression(log[2](intensity+1)))
#comment   #points(xcoords[767154], ycoords[767154], pch=16, col="red")
#comment
#comment   #To calculate a robust measure of background for each array using median of the five lowest pixel values
#comment   Brob<-medianBackground(TIFF,cbind(xcoords,ycoords))
#comment   data<-insertBeadData(data,array=i,what="GrnRB",Brob)
#comment
#comment   #To calculate foreground values in the normal way
#comment   TIFF2<-illuminaSharpen(TIFF)
#comment
#comment   #To calculate foreground values
#comment   IllF<-illuminaForeground(TIFF2, cbind(xcoords,ycoords))
#comment   data<-insertBeadData(data,array=i,what="GrnF",IllF)
#comment
#comment   #To substract the median background values to get locally background corrected intensities
#comment   data<-backgroundCorrectSingleSection(data, array = i, fg="GrnF", bg="GrnRB", newName = "GrnR")
#comment
#comment   #To compare the Illumina intensity with the robust intensity, plot the locations of beads whose expressions change substantially and overlay the locations of the implausibly low-intensity pixels in red
#comment   oldG<-getBeadData(data,array=i,"Grn")
#comment   newG<-getBeadData(data,array=i,"GrnR")
#comment   #summary(oldG-newG)
#comment
#comment   #To change the directory to QC
#comment   setwd(paste(getwd(),"/qc_details/",sep=""))
#comment
#comment   pdf(file=paste(sapply(strsplit(TIFF_Files[i], ".tif"), "[", 1), "image_low_intensity.pdf",sep="_"), width=11, height=8.5)
#comment   par(mfrow=c(1,2))
#comment   plot(xcoords[(abs(oldG-newG)>50)],ycoords[(abs(oldG-newG)>50)],pch=16,xlab="X",ylab="Y",main="entire array")
#comment   points(col(TIFF)[TIFF<400],row(TIFF)[TIFF<400],col="red",pch=16)
#comment   plot(xcoords[(abs(oldG-newG)>50)],ycoords[(abs(oldG-newG)>50)],pch=16,xlim=c(1145,1180),ylim=c(15500,15580),xlab="X",ylab="Y",main="zoomed in")
#comment   points(col(TIFF)[TIFF<400],row(TIFF)[TIFF<400],col="red",pch=16)
#comment   garbage <- dev.off()
#comment
#comment   #To make image plots to identify spatial artifacts on the array surface that can occur from mis-handling or scanning problems
#comment   pdf(file=paste(sapply(strsplit(TIFF_Files[i], ".tif"), "[", 1), "imageplot.pdf",sep="_"))
#comment   par(mfrow=c(6,2))
#comment   par(mar=c(1,1,1,1))
#comment   print(imageplot(data, array=i, low="lightgreen", high="darkgreen", zlim=c(4,10), main=sectionNames(data)[i]))
#comment   garbage <- dev.off()
#comment
#comment
#comment   #To plot the location of outliers on the arrays with the most obvious spatial artifacts and plot their location
#comment   pdf(file=paste(sapply(strsplit(TIFF_Files[i], ".tif"), "[", 1), "outliers.pdf",sep="_"), width=11, height=8.5)
#comment   outlierplot(data, array=i, main=paste(sectionNames(data)[i], "outliers"))
#comment   garbage <- dev.off()
#comment
#comment   #To move up one directory and out of qc_details directory
#comment   setwd("..")
#comment
#comment   #To store green channel intensities for all arrays
#comment   if (dim(Array_GreenInt)[1]==dim(as.matrix(logGreenChannelTransform(data,i)))[1])
#comment   {
#comment     Array_GreenInt[,i] = as.matrix(logGreenChannelTransform(data,i))
#comment   }
#comment   else
#comment   {
#comment     Array_GreenInt[,i] = rbind(as.matrix(logGreenChannelTransform(data,i)),as.matrix(rep("NA", dim(Array_GreenInt)[1]-dim(as.matrix(logGreenChannelTransform(data,i)))[1]),1))
#comment   }
#comment }
#comment
#comment
#comment #To change the directory to QC
#comment setwd(paste(getwd(),"/qc_details/",sep=""))
#comment
#comment #To create a box plot
#comment pdf(file=paste(sapply(strsplit(TIFF_Files[1], "_"), "[", 1), "BoxPlot.pdf",sep="_"), width=11, height=8.5)
#comment
#comment boxplot(data, transFun = logGreenChannelTransform, col = "green", ylab=expression(log[2](intensity)), las = 2, outline = FALSE, main = sapply(strsplit(TIFF_Files[1], "_"), "[", 1))
#comment
#comment garbage <- dev.off()
#comment
#comment #To move up one directory and out of qc_details directory
#comment setwd("..")
#comment
#comment #-----------------------------------------------------------------------------
#comment # Notes
#comment #source("http://www.bioconductor.org/biocLite.R")
#comment #biocLite(c("beadarray", "limma", "GEOquery", "illuminaHumanv1.db", "illuminaHumanv2.db", "illuminaHumanv3.db", "illuminaHumanv4.db", "BeadArrayUseCases", "GOstats", "GenomicRanges", "Biostrings"))
#comment
#comment #To follow Use Cases
#comment #https://bioconductor.org/packages/release/data/experiment/html/BeadArrayUseCases.html
#comment
#comment
#comment #Quality assessment using scanner metrics
#comment #We do not have this data available
#comment #ht12metrics <- read.table(system.file("extdata/Chips/Metrics.txt", package="BeadArrayUseCases"), sep="\t", header=TRUE, as.is=TRUE)
#comment #ht12snr <- ht12metrics$P95Grn/ht12metrics$P05Grn
#comment #labs <- paste(ht12metrics[, 2], ht12metrics[, 3], sep = "_")
#comment #par(mai = c(1.5, 0.8, 0.3, 0.1))
#comment #plot (1:12, ht12snr, pch = 19, ylab = "P95/P05", xlab = "", main = "Signal-to-noise ratio for HT12 data", axes = FALSE, frame.plot=TRUE)
#comment #axis(2)
#comment #axis(1, 1:12, labs, las = 2)
