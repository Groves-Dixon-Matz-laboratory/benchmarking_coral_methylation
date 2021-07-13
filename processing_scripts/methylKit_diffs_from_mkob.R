#!/usr/bin/env Rscript
#methylKit_diffs_from_mkob.R
##############################################
#### for installing the package on TACC: #####
##############################################
# #(only needs to be done once)
# source("https://bioconductor.org/biocLite.R")
# myTaccRLib="/home1/02260/grovesd/R/x86_64-pc-linux-gnu-library/3.4" #find lib by typing installed.packages() with R open
# install.packages('matrixStats')
# install.packages('futile.logger')
# install.packages('Rcurl')
# biocLite("GenomeInfoDb", lib = myTaccRLib, lib.loc = myTaccRLib)
# biocLite("DelayedArray", lib = myTaccRLib, lib.loc = myTaccRLib)
# biocLite("methylKit", lib = myTaccRLib, lib.loc = myTaccRLib) 
# biocLite("genomation", lib = myTaccRLib, lib.loc = myTaccRLib) 
##############################################

#load packages
suppressMessages(library(methylKit))
suppressMessages(library(optparse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
# library(plyr)

#PARSE ARGUMENTS

option_list = list(
  
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input R data file", metavar="input_file"),

  make_option(c("-m", "--metadata.file"), type="character", default=NULL, 
              help="Table with sample data for comparisons", metavar="metadata_file"),
  
  make_option(c("-N", "--ncore"), type="integer", default=1, 
              help="Number of cores to use", metavar="number_of_cores")
)

print("Parsing arugments...")

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

inputFile = opt$input
metadata.file = opt$metadata.file
N.CORES = opt$ncore
outputPrefix = sub('.Rdata', '', inputFile)
outputFile = paste(outputPrefix, 'methylKit.Rdata', sep='_')
mdat = read.table(metadata.file, header = T, stringsAsFactors=F)
treats0 = mdat$treat
treats = c()
for (t in treats0){
  treats = append(treats, c(t,t))
}

print(paste('loading input file ', inputFile, '...', sep=''))
ll=load(inputFile)
countsCols = colnames(my.reduced.counts)[grep('^num', colnames(my.reduced.counts))]
cdat = data.frame(my.reduced.counts)[countsCols]
cdat.t0 = cdat[,treats==0]
cdat.t1 = cdat[,treats==1]
nNA.t0 = apply(cdat.t0, 1, function(x) sum(is.na(x)))
nNA.t1 = apply(cdat.t1, 1, function(x) sum(is.na(x)))
table(nNA.t0)
table(nNA.t1)
bad.t0 = nNA.t0==(ncol(cdat.t0)-4)
bad.t1 = nNA.t1==(ncol(cdat.t1)-4)
bad = bad.t0 | bad.t1


nZero.t0 = apply(cdat.t0, 1, function(x) sum(x==0, na.rm=TRUE))
nZero.t1 = apply(cdat.t1, 1, function(x) sum(x==0, na.rm=TRUE))
table(nZero.t0)
table(nZero.t1)
bad.t0 = nZero.t0 > (ncol(cdat.t0)-4)
bad.t1 = nZero.t1 > (ncol(cdat.t1)-4)
bad = bad.t0 | bad.t1

good.counts = my.reduced.counts[!bad, ]



# rsums = apply(cdat, 1, function(x) sum(x, na.rm=TRUE))
# zero = rsums==0
# nonZero = my.reduced.counts[!zero,]


#CALCULATE DIFFERENCES
print('Running DiffMeth on windows...')
myDiff=calculateDiffMeth(good.counts, mc.cores=N.CORES)

print("Done.")
print("Results head:")
print(head(myDiff))


#SAVE THE RESULTS
windowOutName = paste(outfile.prefix, "_methylKit.Rdata", sep="")
print(paste("Saving window results as", windowOutName))
save(wbounds, totCounts, myDiff, file=windowOutName)