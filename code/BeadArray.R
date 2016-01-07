#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("beadarray", "limma", "GEOquery", "illuminaHumanv1.db", "illuminaHumanv2.db", "illuminaHumanv3.db", "illuminaHumanv4.db", "BeadArrayUseCases", "GOstats", "GenomicRanges", "Biostrings"))

#To follow Use Cases
#https://bioconductor.org/packages/release/data/experiment/html/BeadArrayUseCases.html


#Quality assessment using scanner metrics
#We do not have this data available
#ht12metrics <- read.table(system.file("extdata/Chips/Metrics.txt", package="BeadArrayUseCases"), sep="\t", header=TRUE, as.is=TRUE)
#ht12snr <- ht12metrics$P95Grn/ht12metrics$P05Grn
#labs <- paste(ht12metrics[, 2], ht12metrics[, 3], sep = "_")
#par(mai = c(1.5, 0.8, 0.3, 0.1))
#plot (1:12, ht12snr, pch = 19, ylab = "P95/P05", xlab = "", main = "Signal-to-noise ratio for HT12 data", axes = FALSE, frame.plot=TRUE)
#axis(2)
#axis(1, 1:12, labs, las = 2)

#Data import and storage
library(beadarray)
library(illuminaHumanv4.db)


#To create two directories within
dir.create(paste(getwd(),"/results/",sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create(paste(getwd(),"/qc_details/",sep=""), showWarnings = TRUE, recursive = FALSE, mode = "0777")


#There is a huge issue with sampleSheet so just move everything into one folder and then read it
#chipPath <- "/home/h/hparikh/IL_2_G_9224523070"
#list.files(chipPath)
#sampleSheetFile <- paste(chipPath, "/sampleSheet.csv",sep="")
#readLines(sampleSheetFile)

#Known issue; there is a problem of using .csv file instead it is better to import all the arrays from a directory
#http://www.rdocumentation.org/packages/beadarray/functions/readIllumina
#https://github.com/Bioconductor-mirror/beadarray/blob/master/R/readIllumina.R
#https://support.bioconductor.org/p/61023/


#To read the data using readIllumina(), which will extract the intensities from the .txt, .locs and also TIFF files
#Annotation from https://bioconductor.org/packages/release/data/annotation/html/illuminaHumanv4.db.html
data <- readIllumina(dir=".", useImages=TRUE, illuminaAnnotation="Humanv4")


#slotNames(data)
#data@sectionData
Array_Names=as.matrix(sectionNames(data))
Array_NumBeads=as.matrix(numBeads(data))
#head(data[[1]])
#getBeadData(data, array=1, what="Grn")[1:5]
#getBeadData(data, array=1, what="GrnX" )[1:5]
#getBeadData(data, array=1, what="ProbeID")[1:5]


#log2(data[[1]][1:10, "Grn"])
#log2(getBeadData(data, array=1, what="Grn")[1:10])
#logGreenChannelTransform(data, array=1)[1:10]
#logGreenChannelTransform

#To list all the TIFF files
TIFF_Files = list.files(pattern="*.tif")

#To store green intensities for all the array which is used for the box plot
Array_GreenInt <- matrix(, nrow = max(Array_NumBeads), ncol = length(TIFF_Files))

#To read each of the TIFF file separately
for (i in 1:length(TIFF_Files))
{
	#To read each of the TIFF file
	TIFF <- readTIFF(TIFF_Files[i])
	#cbind (col(TIFF)[which(TIFF==0)], row(TIFF)[which(TIFF==0)])
	xcoords <- getBeadData(data, array=i, what="GrnX")
	ycoords <- getBeadData(data, array=i, what="GrnY")
	#par(mfrow = c(1, 3))
	#offset <- 1

	#plotTIFF(TIFF+offset, c(1517, 1527), c(5507, 5517), values=T, textCol="yellow", main=expression(log[2](intensity+1)))
	#points(xcoords[503155], ycoords[503155], pch=16, col="red")
	#plotTIFF(TIFF+offset, c(1202, 1212), c(13576, 13586), values=T, textCol="yellow", main=expression(log [2](intensity+1)))
	#points(xcoords[625712], ycoords[625712], pch=16, col="red")
	#plotTIFF(TIFF+offset, c(1613, 1623), c(9219, 9229), values=T, textCol="yellow", main=expression(log[2](intensity+1)))
	#points(xcoords[767154], ycoords[767154], pch=16, col="red")

	#To calculate a robust measure of background for each array using median of the five lowest pixel values
	Brob<-medianBackground(TIFF,cbind(xcoords,ycoords))
	data<-insertBeadData(data,array=i,what="GrnRB",Brob)

	#To calculate foreground values in the normal way
	TIFF2<-illuminaSharpen(TIFF)
	
	#To calculate foreground values
	IllF<-illuminaForeground(TIFF2, cbind(xcoords,ycoords))
	data<-insertBeadData(data,array=i,what="GrnF",IllF)

	#To substract the median background values to get locally background corrected intensities
	data<-backgroundCorrectSingleSection(data, array = i, fg="GrnF", bg="GrnRB", newName = "GrnR")

	#To compare the Illumina intensity with the robust intensity, plot the locations of beads whose expressions change substantially and overlay the locations of the implausibly low-intensity pixels in red
	oldG<-getBeadData(data,array=i,"Grn")
	newG<-getBeadData(data,array=i,"GrnR")
	#summary(oldG-newG)

	#To change the directory to QC
	setwd(paste(getwd(),"/qc_details/",sep=""))

	pdf(file=paste(sapply(strsplit(TIFF_Files[i], ".tif"), "[", 1), "image_low_intensity.pdf",sep="_"), width=11, height=8.5)
	par(mfrow=c(1,2))
	plot(xcoords[(abs(oldG-newG)>50)],ycoords[(abs(oldG-newG)>50)],pch=16,xlab="X",ylab="Y",main="entire array")
	points(col(TIFF)[TIFF<400],row(TIFF)[TIFF<400],col="red",pch=16)
	plot(xcoords[(abs(oldG-newG)>50)],ycoords[(abs(oldG-newG)>50)],pch=16,xlim=c(1145,1180),ylim=c(15500,15580),xlab="X",ylab="Y",main="zoomed in")
	points(col(TIFF)[TIFF<400],row(TIFF)[TIFF<400],col="red",pch=16)
	garbage <- dev.off()

	#To make image plots to identify spatial artifacts on the array surface that can occur from mis-handling or scanning problems
	pdf(file=paste(sapply(strsplit(TIFF_Files[i], ".tif"), "[", 1), "imageplot.pdf",sep="_"))
	par(mfrow=c(6,2))
	par(mar=c(1,1,1,1))
	print(imageplot(data, array=i, low="lightgreen", high="darkgreen", zlim=c(4,10), main=sectionNames(data)[i]))
	garbage <- dev.off()


	#To plot the location of outliers on the arrays with the most obvious spatial artifacts and plot their location
	pdf(file=paste(sapply(strsplit(TIFF_Files[i], ".tif"), "[", 1), "outliers.pdf",sep="_"), width=11, height=8.5)
	outlierplot(data, array=i, main=paste(sectionNames(data)[i], "outliers"))
	garbage <- dev.off()

	#To move up one directory and out of qc_details directory
	setwd("..")

	#To store green channel intensities for all arrays
	if (dim(Array_GreenInt)[1]==dim(as.matrix(logGreenChannelTransform(data,i)))[1])
	{
		Array_GreenInt[,i] = as.matrix(logGreenChannelTransform(data,i))
	}
	else
	{
		Array_GreenInt[,i] = rbind(as.matrix(logGreenChannelTransform(data,i)),as.matrix(rep("NA", dim(Array_GreenInt)[1]-dim(as.matrix(logGreenChannelTransform(data,i)))[1]),1))
	}
}


#To change the directory to QC
setwd(paste(getwd(),"/qc_details/",sep=""))

#To create a box plot
pdf(file=paste(sapply(strsplit(TIFF_Files[1], "_"), "[", 1), "BoxPlot.pdf",sep="_"), width=11, height=8.5)

boxplot(data, transFun = logGreenChannelTransform, col = "green", ylab=expression(log[2](intensity)), las = 2, outline = FALSE, main = sapply(strsplit(TIFF_Files[1], "_"), "[", 1))

garbage <- dev.off()

#To move up one directory and out of qc_details directory
setwd("..")
