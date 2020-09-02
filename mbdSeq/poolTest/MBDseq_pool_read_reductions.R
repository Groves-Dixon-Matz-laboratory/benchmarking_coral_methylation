#MBDseq_read_reductions.R
rm(list=ls())
library(tidyverse)
library(cowplot)
source('mbdSeq/MBD_functions.R')

#read in the counts
infile = 'mbdSeq/poolTest/mbd_gene_multicov.tsv'
counts = read.table(infile, header = TRUE)
head(counts)
pdat = counts[,1:4]
lengths = pdat$stop - pdat$start
cdat = counts[,5:ncol(counts)]
rownames(cdat) = pdat$name
windowNames = rownames(cdat)


# BUILD REDUCED COUNT SETS ------------------------------------------------

#build a step-wise reduced counts matrices
reductionMultipliers0 = sapply(12:1, function(x) return(2^-x))
reductionMultipliers = append(reductionMultipliers0, c(0.75, 0.90, 1.0))

resampledList = list()
csums = apply(cdat, 2, sum)

for (redM in reductionMultipliers){
  pctLab = paste('pct', redM*100, sep='.')
  print('---------------')
  print(paste('reduction =', redM))
  propSub = csums*redM
  resampled=c()
  for (s in 1:ncol(cdat)) {
    probs= cdat[,s]/csums[s]
    subSample = sample(c(1:nrow(cdat)),
                       propSub[s],
                       prob=probs,
                       replace=TRUE)
    cts = hist(subSample,
               breaks=c(0:nrow(cdat)),
               plot=FALSE)$counts
    resampled = data.frame(cbind(resampled, cts))
  }
  colnames(resampled)=names(csums)
  rsums = apply(resampled, 2, sum)
  print('mean resampled proportion:')
  print(mean(rsums/csums))
  rownames(resampled) = rownames(cdat)
  resampledList[[pctLab]] = resampled
}
names(resampledList)
# resampledList[['pct.100']]=cdat


# SET UP COLDATA ----------------------------------------------------------
counts=cdat
sample = colnames(counts)
fraction = sapply(sample, function(x) strsplit(x, '_', fixed=TRUE)[[1]][1])
genotype = sapply(sample, function(x) strsplit(x, '_', fixed=TRUE)[[1]][2])
coldata = data.frame(Run=sample,
                     genotype=factor(genotype, levels = c('L5', 'N12')),
                     fraction=factor(fraction, levels=c('m', 'ub')))
rownames(coldata)=colnames(counts)
coldata


# GET DIFFERENCES FOR REDUCED COUNT SETS ----------------------------------
library(DESeq2)

#modified version of function to work with pools since they are oblivious to tissue
get_response_with_ub_pool = function(counts){
  
  print('Running DESeq for dataset:')
  print(head(counts))
  #set up input matrix for DESeq
  int.ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                       colData = coldata,
                                       design = formula(~ genotype + fraction + genotype:fraction))
  #run DESeq
  int.dds = DESeq(int.ddsHTSeq,
                  fitType = 'local')
  print('Results names:')
  print(resultsNames(int.dds))
  print(paste('Using', RESULTSNAME))
  INDEPDENDENT_FILTERING=FALSE
  res = results(int.dds, name=RESULTSNAME, independentFiltering=INDEPDENDENT_FILTERING)
  return(res)
}



#get genotype results
RESULTSNAME = 'genotypeN12.fractionub'
g.ub.resList = lapply(resampledList, function(x) get_response_with_ub_pool(x))
names(g.ub.resList)=names(resampledList)

#assemble into single dataframe
names(g.ub.resList)
change_names = function(n){
  dat = g.ub.resList[[n]]
  dat = data.frame(dat)[,c('log2FoldChange', 'padj')]
  colnames(dat)=paste(n, colnames(dat), sep='_')
  dat$name = rownames(dat)
  return(dat)
}
colnamedList = lapply(names(g.ub.resList), function(x) change_names(x))
mbd.dat = purrr::reduce(colnamedList, full_join, by='name')
head(mbd.dat)


#compare them to results from before
source('benchmarking_functions.R')
ll=load('mbdSeq/datasets/genotype_response_ub.Rdata')
ll
odat = g.ub.resListPos[['gene']]


lfcs = mbd.dat %>%
  dplyr::select('name', grep('log2FoldChange', colnames(mbd.dat))) %>%
  left_join(odat, by = 'name')
reductions = colnames(lfcs)[grep('pct', colnames(lfcs))]
pltList = list()
for (rd in reductions){
  print('--------')
  print(rd)
  plt=plot_scatter_r2_annotated(dat = lfcs,
                                      xcol='log2FoldChange',
                                      ycol=rd,
                                      xlab='original',
                                      ylab=rd) +
    labs(subtitle=sub('_log2FoldChange', '', rd)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  pltList[[rd]]=plt
}

length(pltList)
plot_grid(plotlist=pltList, nrow=2)


#save
save(mbd.dat, csums, file='mbdSeq/poolTest/countReducedDifferencesWub.Rdata')



