#!/usr/bin/env Rscript
#methylKit1.R
#This script gives summary statsistics of the 
#differential methylation analysis performed by methylKit0.R



#load packages
library(methylKit)
library(optparse)
# library(plyr)

#PARSE ARGUMENTS

option_list = list(
  
  make_option(c("-d", "--methylDiff.object"), type="character", default=NULL, 
              help="The methyl diff object output from methylKit0.R", metavar="meth_diff_object"),
	
	make_option(c("-q", "--q.value"), type="character", default="0.05", 
              help="Number of cores to use", metavar="q_value_cutoff"),

	make_option(c("-p", "--pct.diff"), type="character", default="25", 
              help="output file name [default= %default]", metavar="meth_difference_cutoff"),

	make_option(c("-s", "--outfile.name"), type="character", default="methylKitOutput.txt", 
              help="output file name [default= %default]", metavar="save_as")
)

print("Parsing arugments...")

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dat.in = opt$methylDiff.object
q.cut = as.numeric(opt$q.value)
pctd.cut = as.numeric(opt$pct.diff)
outfile.name = opt$outfile.name


print("--------------")
print("Run parameters:")
print(paste("Dataset =", dat.in))
print(paste("q-value cutoff =", q.cut))
print(paste("output file name =", outfile.name))


#SET UP INPUT
print("Loading data...")
ll=load(dat.in)








