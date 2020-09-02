#methylRAD_read_reductions.R
#artificially reduce read counts to assess effects of agreement with other assays

rm(list=ls())
library(tidyverse)
source('benchmarking_functions.R')

#----- SELECT DATASET TO RUN -----------#
#full datset
infile = 'methylRAD/datasets/mr_gene_multicov.tsv';SUFFIX='FULL' #original full dataset

#or

#8-sample subset
infile = 'methylRAD/datasets/eightSample_gene_multicov.tsv';SUFFIX='' #8-sample reduction made in process_methylRAD.R
#---------------------------------------#

#READ IN THE COUNTS
counts = read.table(infile, header = TRUE)
head(counts)
pdat = counts[,1:4]
lengths = pdat$end - pdat$start
cdat = counts[,5:ncol(counts)]
rownames(cdat) = pdat$name
windowNames = rownames(cdat)

#read in the pipeline counts
pcdat = read_tsv('methylRAD/pipelineCounts/pipelineCounts.tsv',
                 col_names=c('sample', 'value', 'stat')) %>% 
  filter(stat=='predupPropPaired')
csums = pcdat$value/2
mod.names = sub('_flagstats.txt', '' , pcdat$sample)
mod.names = sub('-', '.', mod.names, fixed = TRUE)
mod.names = sub('-', '.', mod.names, fixed = TRUE)
names(csums) = mod.names
sum(colnames(counts) %in% mod.names)
csums = csums[names(csums) %in% colnames(counts)]
length(csums)
csumOut=paste('methylRAD/datasets/csums', SUFFIX, '.Rdata', sep = '')
save(csums, file=csumOut)

# BUILD REDUCED COUNT SETS ------------------------------------------------

# #build a step-wise reduced counts matrices
reductionMultipliers0 = sapply(12:1, function(x) return(2^-x))
reductionMultipliers = append(reductionMultipliers0, c(0.75, 0.90, 1))

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

resampledList[['original']]=cdat #load the original as the set you started with (no resampling)

resampledListFile=paste('methylRAD/datasets/resampledList', SUFFIX, '.Rdata', sep='')
save(resampledList, file=resampledListFile)

# GET FPKM FOR REDUCED SETS -----------------------------------------------

# #optionally load the reduced datset from the full mdRAD set (instead of 8 sample reduction)
# ll=load('methylRAD/datasets/resampledListFULL.Rdata')

#get fpkm for methyl RAD counts (note this different from )
get_fpkm_reduction = function(rcounts){
  m = apply(rcounts, 2, function(x) sum(x)/1e6)
  fpkm=sweep(rcounts, 1, k, `/`) %>% 
    sweep(2, m, `/`)
  return(fpkm)
}


k=lengths/1e3
fpkmList = lapply(resampledList, function(x) get_fpkm_reduction(x))
names(fpkmList) = names(resampledList)

#average them
get_means = function(df){
  mns=apply(df, 1, mean)
  res=data.frame(tag=windowNames,
                 mn=log(mns, 2))
  return(res)
}
mnList = lapply(fpkmList, function(x) get_means(x))
names(mnList)=names(fpkmList)


#assemble together
rdat = pdat
for (n in names(mnList)){
  rdat[,n]=mnList[[n]]$mn
}
head(rdat)


#check the correlations
resampled = colnames(rdat)[grep('pct', colnames(rdat))]
pltList = list()
for (rs in resampled){
  pltList[[rs]]=plot_scatter_r2_annotated(dat = rdat,
                            xcol='original',
                            ycol=rs,
                            xlab='original',
                            ylab=rs) + 
    labs(subtitle=rs) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

plot_grid(plotlist = pltList, nrow=2)

#save
countReducedLvlsFile=paste('methylRAD/datasets/countReducedLvls',SUFFIX,'.Rdata',sep='')
save(rdat, csums, file=countReducedLvlsFile)


# GET DIFFERENCES FOR REDUCED COUNT SETS ----------------------------------
names(resampledList)

#SET UP COLDATA
counts=resampledList[[1]]
sample = sapply(colnames(counts), function(x) strsplit(x, '_')[[1]][1])
enzyme = sub('mr',
             '',
             sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]))
genotype = sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
tissue = substr(sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3]),
                1,
                1)


coldata = data.frame(Run=sample,
                     genotype=genotype,
                     tissue=tissue,
                     enzyme=enzyme)


#RUN DIFF COVERAGE
library(DESeq2)

#get responses for both (includes option to skip those with too few reads left to analyze)
get_response_reduction = function(counts){
  
  print('Running DEseq on data:')
  print(head(counts))
  
  #set up input matrix for DESeq
  ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                   colData = coldata,
                                   design = formula(~ tissue + genotype))
  # Nzeros = apply(counts, 1, function(x) return(0 %in% x))
  medCount =  apply(counts, 1, median)
  over1 = sum(medCount > 1)
  if (over1 > nrow(counts)/10){
    #run DESeq
    dds = DESeq(ddsHTSeq,
                fitType='local')
    #get DEseq results
    resultsNames(dds)
    res = results(dds, contrast = CONTRAST, independentFiltering=FALSE)

  } else {
    print('WARNING. LESS THAN 90% OF THE GENES HAVE MEDIAN COUNT LESS THAN 1')
    res = 'TOO FEW READS TO ANALYZE'
  }
  return(res)
}


#for genotype
CONTRAST=c('genotype', 'L5', 'N12')
g.resList = lapply(resampledList, function(x) get_response_reduction(x))
names(g.resList)=names(resampledList)
length(g.resList)
length(resampledList)
lapply(g.resList, head)
g.resList2=list()
for (n in names(g.resList)){
  r=g.resList[[n]]
  if (class(r)!='character'){
    g.resList2[[n]]=r
  }
}
length(g.resList2)

#assemble into single dataframe
names(g.resList2)
change_names = function(n){
  dat = g.resList[[n]]
  dat = data.frame(dat)[,c('log2FoldChange', 'padj')]
  colnames(dat)=paste(n, colnames(dat), sep='_')
  dat$name = rownames(dat)
  return(dat)
}
colnamedList = lapply(names(g.resList2), function(x) change_names(x))
mr.dat = purrr::reduce(colnamedList, full_join, by='name')
head(mr.dat)


#compare them to results from before
ll=load('methylRAD/datasets/genotype_responses.Rdata')
ll
odat = mr.g.list[['gene']] 
  

lfcs = mr.dat %>% 
  select('name', grep('log2FoldChange', colnames(mr.dat))) %>% 
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
    labs(subtitle=rd) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  pltList[[rd]]=plt
}

length(pltList)
plot_grid(plotlist=pltList, nrow=2)


#save
countReducedDifferencesFile=paste('methylRAD/datasets/countReducedDifferences', SUFFIX, '.Rdata', sep='')
save(mr.dat, csums, file=countReducedDifferencesFile)

