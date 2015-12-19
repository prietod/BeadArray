#!/usr/bin/env Rscript

library(beadarray)

#------------------------------------------------------------------------------------------------------------------
# Parse the arguments
#------------------------------------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE) # get command-line arguments
dirname <- toString(head(args, 1))       # first command-line argument as a string
section_names <- tail(args, 1)           # second, third, etc... command-line arguments as a vector

#------------------------------------------------------------------------------------------------------------------
# Load the array into the data variable 
#------------------------------------------------------------------------------------------------------------------

data <- readIllumina(dir=dirname, useImages=FALSE, illuminaAnnotation="Humanv4", SectionNames = section_names)


