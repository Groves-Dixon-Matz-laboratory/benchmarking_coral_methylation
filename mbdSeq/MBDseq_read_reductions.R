#MBDseq_read_reductions.R
rm(list=ls())
library(tidyverse)
library(cowplot)
source('mbdSeq/MBD_functions.R')

#read in the counts
infile = 'mbdSeq/datasets/mbd_gene_multicov.tsv'
counts = read.table(infile, header = TRUE)
head(counts)
pdat = counts[,1:4]
lengths = pdat$stop - pdat$start
cdat = counts[,5:ncol(counts)]
rownames(cdat) = pdat$name
windowNames = rownames(cdat)


#read in the pipeline counts
plCounts = read_tsv('mbdSeq/pipelineCounts/pipelineCounts.tsv',
                    col_names=c('run', 'value', 'stat'))
dedupMapped = plCounts %>% 
  filter(stat=="dedupMapped")
csums = dedupMapped %>% 
  pull(value)
names(csums) = dedupMapped$run
csums
save(csums, file='mbdSeq/datasets/csums.Rdata')

# BUILD REDUCED COUNT SETS ------------------------------------------------

#build a step-wise reduced counts matrices
reductionMultipliers0 = sapply(12:1, function(x) return(2^-x))
reductionMultipliers = append(reductionMultipliers0, c(0.75, 0.90, 1.0))

resampledList = list()

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
sample = sapply(colnames(counts), function(x) strsplit(x, '_')[[1]][1])
fraction = sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1])
genotype = sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
tissue = substr(sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3]),
                1,
                1)
coldata = data.frame(Run=sample,
                     genotype=factor(genotype, levels = c('L5', 'N12')),
                     tissue=factor(tissue, levels=c('t', 's')),
                     fraction=factor(fraction, levels=c('m', 'ub')))
rownames(coldata)=colnames(counts)
coldata




# GET METHYLATION LEVEL FOR REDUCED SETS -----------------------------------------------
library(DESeq2)

#same function from process_MBDseq.R works here
lvlList = lapply(resampledList, function(x) get_meth_lvl_mbd(x))
names(lvlList)=names(resampledList)


#check bimodals
hists = list()
for (n in names(lvlList)){
  h=lvlList[[n]] %>% 
    data.frame() %>% 
    ggplot(aes(x=log2FoldChange)) +
    geom_histogram() +
    labs(subtitle=n)
  hists[[n]]=h
}
plot_grid(plotlist=hists)
  
  
#assemble together
rdat = pdat
for (n in names(lvlList)){
  print(n)
  rdat[,n]=data.frame(lvlList[[n]])$log2FoldChange
}
head(rdat)



#check the correlations
source('benchmarking_functions.R')
resampled = colnames(rdat)[grep('pct', colnames(rdat))]
pltList = list()
for (rs in resampled){
  pltList[[rs]]=plot_scatter_r2_annotated(dat = rdat,
                            xcol='pct.100',
                            ycol=rs,
                            xlab='original',
                            ylab=rs) + 
    labs(subtitle=rs) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

plot_grid(plotlist = pltList, nrow=2)

#save
save(rdat, csums, file='mbdSeq/datasets/countReducedLvls.Rdata')



# GET DIFFERENCES FOR REDUCED COUNT SETS ----------------------------------

#this also works with function from process_MBDseq.R

#get genotype results
RESULTSNAME = 'genotypeN12.fractionub'
g.ub.resList = lapply(resampledList, function(x) get_response_with_ub(x))
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
save(mbd.dat, csums, file='mbdSeq/datasets/countReducedDifferencesWub.Rdata')

# 
# # GET DIFFERENCES FOR REDUCED COUNT SET WITHOUT UNBOUND ----------------------------
# #left for reference
# 
# CONTRAST = c('genotype', 'L5', 'N12')
# g.b.resList = lapply(resampledList, function(x) get_response_without_ub(x))
# names(g.b.resList)=names(resampledList)
# 
# 
# 
# 
# 
# 
# #assemble into single dataframe
# names(g.b.resList)
# change_names = function(n){
#   dat = g.b.resList[[n]]
#   dat = data.frame(dat)[,c('log2FoldChange', 'padj')]
#   colnames(dat)=paste(n, colnames(dat), sep='_')
#   dat$name = rownames(dat)
#   return(dat)
# }
# colnamedList2 = lapply(names(g.b.resList), function(x) change_names(x))
# mbd.dat2 = purrr::reduce(colnamedList2, full_join, by='name')
# head(mbd.dat2)
# 
# 
# #compare them to results from before
# ll=load('mbdSeq/datasets/genotype_response_ub.Rdata')
# ll
# odat = g.ub.resListPos[['gene']]
# 
# 
# lfcs2 = mbd.dat2 %>%
#   dplyr::select('name', grep('log2FoldChange', colnames(mbd.dat2))) %>%
#   left_join(odat, by = 'name')
# reductions = colnames(lfcs)[grep('pct', colnames(lfcs))]
# pltList = list()
# for (rd in reductions){
#   print('--------')
#   print(rd)
#   plt=plot_scatter_r2_annotated(dat = lfcs2,
#                                 xcol='log2FoldChange',
#                                 ycol=rd,
#                                 xlab='original',
#                                 ylab=rd) +
#     labs(subtitle=sub('_log2FoldChange', '', rd)) +
#     theme(axis.title.x = element_blank(),
#           axis.title.y = element_blank())
#   pltList[[rd]]=plt
# }
# 
# length(pltList)
# plot_grid(plotlist=pltList, nrow=2)
# 
# 
# #save
# csums2 = csums[grep("^m", names(csums))]
# save(mbd.dat2, csums2, file='mbdSeq/datasets/countReducedDifferencesWoutUB.Rdata')
