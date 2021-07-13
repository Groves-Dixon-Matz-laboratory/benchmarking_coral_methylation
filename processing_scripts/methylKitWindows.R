#!/usr/bin/env Rscript
#methylKit.R
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
  
  make_option(c("-m", "--metadata.file"), type="character", default=NULL, 
              help="Tab-delimited table with sample data for comparisons", metavar="metadata_file"),
	
	make_option(c("-N", "--ncore"), type="integer", default=1, 
              help="Number of cores to use", metavar="number_of_cores"),
	
	make_option(c("-s", "--prefix"), type="character", default="methylKitOutput", 
              help="output file name [default= %default]", metavar="output_prefix"),
	
	make_option(c("-lowCount", "--minimum_depth"), type="integer", default=10, 
              help="Minimum fold coverage to keep a position (see ?filterByCoverage)", metavar="minimum_depth"),
	
	make_option(c("-hiPct", "--hi.pct"), type="integer", default=99.9, 
              help="output file name [default= %default]", metavar="highest_percentile"),
  
  make_option(c("-a", "--assembly"), type="character", default='assembly', 
              help="string handle for assembly", metavar="assembly_tag"),

  make_option(c("-mpg", "--min_per_group"), type="integer", default=2, 
              help="minimum samples per treatment group that need data to keep site
              (passed to min.per.group in unite() function)", metavar="min_per_group"),

  make_option(c("-w", "--windows_bed"), type="character",
              help="name of the bed file with windows (expect 4 column bed files with chr,start,end,", metavar="windows_bed")
)

print("Parsing arugments...")

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadata.file = opt$metadata.file
N.CORES = opt$ncore
outfile.prefix = opt$prefix
LOW.COUNT = opt$minimum_depth
HI.PCT = opt$hi.pct
assembly.tag = opt$assembly
windowBoundFile = opt$windows_bed
minPerGroup = opt$min_per_group
saveWindows = opt$save_windows



print("--------------")
print("Run parameters:")
print(paste("meta data =", metadata.file))
print(paste("number of cores to use =", N.CORES))
print(paste("output file prefix =", outfile.prefix))
print(paste("low count cutoff =", LOW.COUNT))
print(paste("minimum per group =", minPerGroup))
print(paste("high percentile cutoff =", HI.PCT))
print(paste("windows input bed file =", windowBoundFile))



#SET UP INPUT
print(paste("Reading in metadata file", metadata.file))
mdat = read.table(metadata.file, header = TRUE, stringsAsFactors=FALSE, sep='\t')
fl = mdat$run
sid = as.list(mdat$id)
file.list=as.list(fl)
treat=mdat$treat
print("Sample traits:")
print(mdat)


#READ IN WINDOW BOUNDARIES
print('Reading in window boundaries...')
wbounds = read.table(windowBoundFile, col.names=c('chr', 'start', 'end', 'name'))
print(head(wbounds))


#READ IN THE DATA
print("Reading in data...")
myobj=methRead(file.list,
           sample.id=sid,
           assembly=assembly.tag,
           treatment=treat,
           context="CpG",
           pipeline = "bismarkCoverage",
           header = FALSE
)


#FILTER CPGS BASED ON COVERAGE

filtered.myobj=filterByCoverage(myobj,lo.count=LOW.COUNT,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=HI.PCT)

#ASSEMBLE DATA TOGETHER INCLUDING CPGS WITH COVERAGE IN ALL SAMPLES

#from manual: "In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses). In addition, setting destrand=TRUE will only work when operating on base-pair resolution, otherwise setting this option TRUE will have no effect. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples."

print("Merging up the CpG sites...")
print(paste('minimum per group = ', minPerGroup))
methForLvl=methylKit::unite(filtered.myobj,
                      min.per.group = 0L)
meth=methylKit::unite(filtered.myobj,
                      min.per.group = as.integer(minPerGroup))

print("Results:")
print(head(meth))
print("Dimentions:")
print(dim(meth))

# print("Saving merged CpGs as mdf.tsv")
# meth.out = paste(outfile.prefix, "unit_object.Rdata", sep="_")
# print(paste("Saving the united cpg table as", meth.out))
# save(meth, file=meth.out)
# write.table(data.frame(meth), file="mdf.tsv", quote=FALSE, row.names=FALSE, sep="\t")


#SET UP WINDOWS

#assign coordinates based on window inputs
SEQNAMES = wbounds$chr
START=wbounds$start
END=wbounds$end
RANGES = IRanges(start=START, end=END)
NAMES=wbounds$name

#created the Granges object specifying the regions
my.win=GRanges(seqnames=SEQNAMES,
               ranges=RANGES)

#create counts object for regions

my.reg.countsForLvl = regionCounts(object = methForLvl,
                             regions = my.win)

my.reg.counts = regionCounts(object = meth,
                             regions = my.win)


#GET SUMS FROM REGION COUNTS

#make dataframe
rcdat = as_tibble(data.frame(my.reg.countsForLvl))
rcdat2 = rcdat %>% 
  select(-strand, -grep('coverage', colnames(rcdat)))

#gather ncs and nts by sample
ncs = rcdat2 %>% 
  select(chr, start, end, grep('numCs', colnames(rcdat2))) %>% 
  gather(key='sample', value='nC', -chr, -start, -end) %>% 
  mutate(sample=sub('numCs', '', sample))

nts = rcdat2 %>% 
  select(chr, start, end, grep('numTs', colnames(rcdat2))) %>% 
  gather(key='sample', value='nT', -chr, -start, -end) %>% 
  mutate(sample=sub('numTs', '', sample))

#merge with wbounds
pos = ncs %>% 
  full_join(nts, by=c('chr', 'start', 'end', 'sample'))
totCounts = merge(pos, wbounds, by = c('chr', 'start', 'end'))[,c('chr', 'start', 'end', 'name', 'sample', 'nC', 'nT')]
head(totCounts)



#calculate differences
print('Running DiffMeth on windows...')
myDiff=calculateDiffMeth(my.reg.counts, mc.cores=N.CORES)

print("Done.")
print("Results head:")
print(head(myDiff))


#SAVE THE RESULTS

#save the main results
windowOutName = paste(outfile.prefix, "_methylKit.Rdata", sep="")
print(paste("Saving window results as", windowOutName))
save(wbounds, totCounts, myDiff, file=windowOutName)

#also save the count split by window
wcOut = paste(windowOutName, '_methylKit_regionCounts.Rdata', sep='')
print(paste('Saving window counts as', wcOut))
save(my.reg.counts, file=wcOut)


#ALSO SAVE AS DATAFRAME
ddf = data.frame(myDiff)
ddf = ddf[order(ddf$pvalue),] %>%
  filter(qvalue < 0.1)
dataFrameOut = paste(outfile.prefix, "significant.tsv", sep="_")
print(head(ddf, n=20))
print(paste("Writing out significant results (qvalue < 0.1) to file:", dataFrameOut))
write.table(ddf, file=dataFrameOut, row.names=FALSE, quote=FALSE, sep='\t')





