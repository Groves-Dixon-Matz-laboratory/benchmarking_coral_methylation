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

  make_option(c("-ppg", "--min_prop_group"), type="double", default=0.75, 
              help="minimum proportion of samples per treatment group that need data to keep site
              (This argument uses the metadata.file to get the count to pass to the min.per.group in unite() function)", metavar="min_prop_group"),

  make_option(c("-w", "--windows_bed"), type="character",
              help="name of the bed file with windows (expect 4 column bed files with chr,start,end,", metavar="windows_bed")
)

print("Parsing arugments...")

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## #For debugging:
metadata.file = 'bisulfite_treat_table.txt'
N.CORES = 12
outfile.prefix = 'TEST'
LOW.COUNT = 3
HI.PCT = 99.9
assembly.tag = '/work/02260/grovesd/lonestar/bumblebee_genome/GCF_000214255.1_Bter_1.0_genomic.fasta'
windowBoundFile = '/work/02260/grovesd/lonestar/bumblebee_genome/cdsBoundaries.bed'
propPerGroup = 0.75




metadata.file = opt$metadata.file
N.CORES = opt$ncore
outfile.prefix = opt$prefix
LOW.COUNT = opt$minimum_depth
HI.PCT = opt$hi.pct
assembly.tag = opt$assembly
windowBoundFile = opt$windows_bed
propPerGroup = opt$min_prop_group
saveWindows = opt$save_windows



print("--------------")
print("Run parameters:")
print(paste("meta data =", metadata.file))
print(paste("number of cores to use =", N.CORES))
print(paste("output file prefix =", outfile.prefix))
print(paste("low count cutoff =", LOW.COUNT))
print(paste("minimum proportion of group =", propPerGroup))
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

#GET THE min_per_group integer
print("Calcualting the min_per_group integer from --min_prop_group and --metadata.file...")
mnGroupSize = round(mean(table(mdat$treatInfo)), digits=0)
minPerGroup = ceiling(propPerGroup*mnGroupSize)
print(paste('Mean group size (rounded to integer) =', mnGroupSize))
print(paste('min_per_group (ceiling integer of ', propPerGroup, '*', mnGroupSize,') = ', minPerGroup, sep=''))


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

filtered.myobj=filterByCoverage(myobj,
  lo.count=as.integer(LOW.COUNT))

#ASSEMBLE DATA TOGETHER INCLUDING CPGS WITH COVERAGE IN ALL SAMPLES

#from manual: "In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses). In addition, setting destrand=TRUE will only work when operating on base-pair resolution, otherwise setting this option TRUE will have no effect. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples."

print("Merging up the CpG sites...")
print(paste('minimum per group = ', minPerGroup))
meth=methylKit::unite(filtered.myobj,
                      min.per.group = as.integer(minPerGroup))
# methForLvl=methylKit::unite(filtered.myobj,
#                       min.per.group = 0L) #decided not to use this, but leaving for reference

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
my.reg.counts = regionCounts(object = meth,
                             regions = my.win)
# my.reg.countsForLvl = regionCounts(object = methForLvl,
#                              regions = my.win) #decided not to use this, but leaving for reference


#GET SUMS FROM REGION COUNTS

#make dataframe
rcdat = as_tibble(data.frame(my.reg.counts))
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


#CHECK REGION COUNTS FOR EXCESSIVE MISSING VALUES BEFORE RUNNING calculateDiffMeth()

print("Checking for windows with too many missing values...")
mdat$num = 1:nrow(mdat)

get_coverages = function(treatNum){
  grp = mdat %>%
    filter(treat==treatNum) %>%
    pull(num)
  keep = paste('coverage', grp, sep='')
  rcSub = rcdat %>%
    select(keep)
  notNaCount = apply(rcSub, 1, function(x) sum(!is.na(x)))
  goodRows = notNaCount >= minPerGroup
  return(goodRows)
}

good0 = get_coverages(0)
good1 = get_coverages(1)
good = good0 & good1

print(paste(nrow(my.reg.counts), 'total regions in dataset'))
print(paste(sum(good), 'total regions passed missing data test'))
filt.reg.counts = my.reg.counts[good,]


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





