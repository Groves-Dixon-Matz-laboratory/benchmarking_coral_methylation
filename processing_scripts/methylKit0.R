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
library(methylKit)
library(optparse)
# library(plyr)

#PARSE ARGUMENTS

option_list = list(
  
  make_option(c("-m", "--metadata.file"), type="character", default=NULL, 
              help="Table with sample data for comparisons", metavar="metadata_file"),
	
	make_option(c("-N", "--ncore"), type="integer", default=1, 
              help="Number of cores to use", metavar="number_of_cores"),
	
	make_option(c("-s", "--prefix"), type="character", default="methylKitOutput", 
              help="output file name [default= %default]", metavar="output_prefix"),
	
	make_option(c("-lowCount", "--low.count"), type="integer", default=10, 
              help="Minimum fold coverage to keep a position (see ?filterByCoverage)", metavar="minimum_depth"),
	
	make_option(c("-hiPct", "--hi.pct"), type="integer", default=99.9, 
              help="output file name [default= %default]", metavar="highest_percentile"),

	make_option(c("-u", "--summary"), type="logical", default=FALSE, 
              help="Boolean for whether to run summary stats", metavar="run_summary_stats"),

	make_option(c("-a", "--assembly"), type="character", default='assembly', 
              help="string handle for assembly", metavar="assembly_tag")
)

print("Parsing arugments...")

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

metadata.file = opt$metadata.file
N.CORES = opt$ncore
outfile.prefix = opt$prefix
LOW.COUNT = opt$low.count
HI.PCT = opt$hi.pct
do.summary = opt$summary
assembly.tag = opt$assembly


print("--------------")
print("Run parameters:")
print(paste("meta data =", metadata.file))
print(paste("number of cores to use =", N.CORES))
print(paste("output file prefix =", outfile.prefix))
print(paste("low count cutoff =", LOW.COUNT))
print(paste("high percentile cutoff =", HI.PCT))
print(paste("run summary =", do.summary))


#SET UP INPUT
print(paste("Reading in metadata file", metadata.file))
mdat = read.table(metadata.file, header = T, stringsAsFactors=F)
fl = mdat$run
sid = as.list(mdat$id)
file.list=as.list(fl)
treat=mdat$treat
print("Sample traits:")
print(mdat)




#READ IN THE DATA
print("Reading in data...")
myobj=methRead(file.list,
           sample.id=sid,
           assembly=assembly.tag,
           treatment=treat,
           context="CpG",
           pipeline = "bismarkCoverage",
           header = F
)


#GET SUMMARY STATISTICS
if (do.summary==TRUE){
	print("Looking at sample summary statistics...")
	for (i in 1:length(file.list)){
		file.name = file.list[[i]][1]
		print(paste(file.name, ":", sep=""))
		print(getMethylationStats(myobj[[i]],plot=FALSE,both.strands=FALSE))
		
		#plot cpg meth percentages
		image1.out = paste(file.name, ".cpgMethPctHistogram.pdf", sep='')
		pdf(file=image1.out)
		getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
		dev.off()
		
		#plot cpg coverages
		getCoverageStats(myobj[[i]],plot=FALSE,both.strands=F)
		image2.out = paste(file.name, ".cpgCoverageHistogram.pdf", sep='')
		pdf(file=image2.out)
		getCoverageStats(myobj[[i]],plot=TRUE,both.strands=F)
		dev.off()
	}
} else{
		print("Skipping prefiltered summary stats.")
} 



#FILTER CPGS BASED ON COVERAGE

filtered.myobj=filterByCoverage(myobj,lo.count=LOW.COUNT,lo.perc=NULL,
                                      hi.count=NULL,hi.perc=HI.PCT)

#GET POST-FILTERING SUMMARY STATISTICS
if (do.summary==TRUE){
	print("Looking at post-filtered sample summary statistics...")
	for (i in 1:length(file.list)){
		file.name = file.list[[i]][1]
		print(paste(file.name, ":", sep=""))
		print(getMethylationStats(filtered.myobj[[i]],plot=FALSE,both.strands=FALSE))
		
		#plot cpg meth percentages
		image1.out = paste(file.name, ".POST_FILT_cpgMethPctHistogram.pdf", sep='')
		pdf(file=image1.out)
		getMethylationStats(filtered.myobj[[i]],plot=TRUE,both.strands=FALSE)
		dev.off()
		
		#plot cpg coverages
		getCoverageStats(filtered.myobj[[i]],plot=FALSE,both.strands=F)
		image2.out = paste(file.name, ".POST_FILT_cpgCoverageHistogram.pdf", sep='')
		pdf(file=image2.out)
		getCoverageStats(filtered.myobj[[i]],plot=TRUE,both.strands=F)
		dev.off()
	}
}	else{
		print("Skipping post filtering summary stats.")
} 


#ASSEMBLE DATA TOGETHER INCLUDING CPGS WITH COVERAGE IN ALL SAMPLES

# from manual: "In order to do further analysis, we will need to get the bases covered in all samples. The following function will merge all samples to one object for base-pair locations that are covered in all samples. Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses). In addition, setting destrand=TRUE will only work when operating on base-pair resolution, otherwise setting this option TRUE will have no effect. The unite() function will return a methylBase object which will be our main object for all comparative analysis. The methylBase object contains methylation information for regions/bases that are covered in all samples."

print("Merging up the CpG sites with coverage from all samples...")
meth=unite(filtered.myobj, destrand=F)
print("Results:")
print(head(meth))
print("Dimentions:")
print(dim(meth))

print("Saving merged CpGs as mdf.tsv")
meth.out = paste(outfile.prefix, "unit_object.Rdata", sep="_")
print(paste("Saving the united cpg table as", meth.out))
save(meth, file=meth.out)
write.table(data.frame(meth), file="mdf.tsv", quote=FALSE, row.names=FALSE, sep="\t")


#RUN SUMMARY STATS

if (do.summary==TRUE){
	print("Running summary statstics...")
	#GET SAMPLE CORRELATIONS
	print("Calculating sample correlations...")
	print(getCorrelation(meth,plot=F))
	pdf(file="correlation_plot.pdf")
	getCorrelation(meth,plot=TRUE)
	dev.off()


	#CLUSTER SAMPLES
	print("Running sample clustering...")
	pdf("clustering_plot.pdf")
	clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
	dev.off()


	#PCA SAMPLES

	#save the pca for closer examination on pc
	print("Running PCA...")
	pdf('sample_pca.pdf')
	pca = PCASamples(meth, obj.return=TRUE)
	save(pca, file="sample_pca.Rdata")
	dev.off()
} else {
	print("Skipping running summary statistics.")
}


#DIFFERENTIAL METHYLATION ANALYSIS
print("Starting differential methylation tests...")
print(paste(c("Using", N.CORES, "parallel cores"), collapse=" "))
print("Groups:")
print(mdat)
myDiff=calculateDiffMeth(meth, mc.cores=N.CORES)
print("Done.")
print("Results head:")
print(head(myDiff))


#SAVE THE METH DIFF OBJECT
objOutName = paste(outfile.prefix, "Rdata", sep=".")
print(paste("Saving meth diff object as", objOutName))
save(myDiff, file=objOutName)


#ALSO SAVE AS DATAFRAME
ddf = data.frame(myDiff)
ddf = ddf[order(ddf$pvalue),]
dataFrameOut = paste(outfile.prefix, "tsv", sep=".")
print(head(ddf, n=20))
print(paste("Writing out results to file:", dataFrameOut))
write.table(ddf, file=dataFrameOut, row.names=F, quote=F)





