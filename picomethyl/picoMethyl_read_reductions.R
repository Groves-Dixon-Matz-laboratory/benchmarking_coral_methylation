#picoMethyl_read_reductions.R

#note this starts with read reductions already done on TACC
#originally did it with sections titled: 'RESAMPLE FOR GBM LEVEL' and 'RESAMPLE COV FILES FOR SENSITIVITY' in walkthrough
#then redid it by resampleing the fastq files under 'SIMULATE REDUCED SEQUENCING'
#the original sections were moved to an appendix. Data from the originals are saved under picomethyl/covResampled_readReductions


rm(list=ls())
library(tidyverse)
source('benchmarking_functions.R')
source('picomethyl/picomethyl_functions.R')

# ASSEMBLE GLM LEVELS FROM BASIC ------------------------------------------

#set up paths and file names for reductions
redDir = 'picomethyl/datasets/readReduced_basicGBM_lvl/'
fileList = list.files(path=redDir, pattern='pct*')
filePaths = list.files(path=redDir, pattern='pct*', full.names=TRUE)
pcts = sub('_gene_basicStatsBed.tsv', '', fileList)


#get csums
cdat = read_tsv('picomethyl/pipeline_counts/pipeline_counts.txt',
                col_names=c('run', 'count', 'step')) %>% 
  filter(step=='mappedCount')
csums = cdat %>% 
  pull(count)
names(csums)=cdat$run
csums
save(csums, file='picomethyl/datasets/csums.Rdata')

#upload the window files
redList = list()
for (i in 1:length(pcts)){
  filePath = filePaths[i]
  pct = pcts[i]
  fres = read_tsv(filePath)
  rdf = data.frame(fres)
  rdf[,pct]=log2_convesion_noZero(rdf$fracMeth)
  redList[[pct]]=rdf[c('chr', 'start', 'end', 'name', pct)]
}
names(redList)

#merge into dataframe
rdat = purrr::reduce(redList, full_join, by=c('chr', 'start', 'end', 'name'))
head(rdat)


#CHECK CORRELATIONS

#upload original
ll=load('comparisons/datasets/gbmLvl.Rdata')
fres = read_tsv('picomethyl/datasets/basic_stats/gene_basicStatsBed.tsv') %>% 
  mutate(l.fracMeth = log2_convesion_noZero(fracMeth)) %>% 
  dplyr::select(chr, start, end, name, l.fracMeth)
head(fres)

#setup columns
resampled0 = colnames(rdat)[grep('pct', colnames(rdat))]
props = as.numeric(sub('pct.', '', resampled0))
resampled = tibble(props, resampled0) %>% 
  arrange(props) %>% 
  pull(resampled0)


#plot
pltList = list()
for (rs in resampled){
  pdat = rdat %>% 
    full_join(fres, by = c('chr', 'start', 'end', 'name'))
  pltList[[rs]]=plot_scatter_r2_annotated(dat = pdat,
                                          xcol='l.fracMeth',
                                          ycol=rs,
                                          xlab='original',
                                          ylab=rs) + 
    labs(subtitle=rs) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}
# plot_grid(plotlist = pltList)

#save
save(rdat, csums, file='picomethyl/datasets/readReduced_basicGBM_lvl/countReducedLvls.Rdata')


# ASSEMBLE GBM RESPONSE REDUCTIONS ----------------------------------------

#set up paths and file names for reductions
redDir = 'picomethyl/datasets/readReduced_gbm_diff/'


fileList = list.files(path=redDir, pattern='pct*')
filePaths = list.files(path=redDir, pattern='pct*', full.names=TRUE)
pcts = sub('_genotype_genes_methylKit.Rdata', '', fileList)
diffList = list()
for (i in 1:length(pcts)){
  pct=pcts[i]
  filePath = filePaths[i]
  ll=load(filePath)
  wbounds$chr=as.character(wbounds$chr)
  d.df = data.frame(myDiff)
  d.df[,paste(pct, 'meth.diff', sep='_')]=d.df$meth.diff
  d.df[,paste(pct, 'qvalue', sep='_')]=d.df$qvalue
  d.df2 = d.df %>% 
    mutate(chr=as.character(chr)) %>% 
    dplyr::select(-strand, -meth.diff, -qvalue, -pvalue) %>% 
    left_join(wbounds, by = c('chr', 'start', 'end'))
  diffList[[pct]] = d.df2
}


#merge into dataframe
pm.dat = purrr::reduce(diffList, full_join, by=c('chr', 'start', 'end', 'name'))
head(pm.dat)

#check the number of significant genes
p.cols = colnames(pm.dat)[grep('qvalue', colnames(pm.dat))]
nsigs = c()
for (p in p.cols){
  nsigs = append(nsigs, sum(pm.dat[,p]<0.1, na.rm=TRUE) )
}
data.frame(p.cols, nsigs)




#compare them to results from before
ll=load('picomethyl/datasets/methylKit_results/genotype_response.Rdata')
ll
odat = g.diffListPos[['gene']] 


lfcs = pm.dat %>% 
  dplyr::select('chr', 'start', 'end', 'name', grep('_meth.diff', colnames(pm.dat))) %>% 
  left_join(odat, by = c('chr', 'start', 'end', 'name'))
reductions0 = colnames(lfcs)[grep('pct', colnames(lfcs))]
props0 = sub('_meth.diff', '', reductions0)
props = as.numeric(sub('pct.', '', props0))
reductions = tibble(reductions0,
                    props) %>% 
  arrange(props) %>% 
  pull(reductions0)
  
  
#double-check against original
pltList = list()
for (rd in reductions){
  print('--------')
  print(rd)
  plt=plot_scatter_r2_annotated(dat = lfcs,
                                xcol='meth.diff',
                                ycol=rd,
                                xlab='original',
                                ylab=rd) + 
    labs(subtitle=rd) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  pltList[[rd]]=plt
}

length(pltList)
# plot_grid(plotlist=pltList, nrow=2)

#replace the 100 with original (don't want to do this but leaving for refence)
m = pm.dat %>%
  left_join(odat, by = c('chr','start','end'))
pm.dat$pct.100_meth.diff = m$meth.diff
pm.dat$pct.100_qvalue = m$qvalue


#SAVE RESULTS
save(pm.dat, csums, file='picomethyl/datasets/readReduced_gbm_diff/countReducedDifferences.Rdata')


