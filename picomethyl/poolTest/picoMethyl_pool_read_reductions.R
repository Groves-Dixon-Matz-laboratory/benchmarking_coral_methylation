#picoMethyl_pool_read_reductions.R

#same as picoMethyl_read_reductions.R but only for the pool test
#(also skips the level reductions part)
rm(list=ls())
library(tidyverse)
source('benchmarking_functions.R')
source('picomethyl/picomethyl_functions.R')

# ASSEMBLE GBM RESPONSE REDUCTIONS ----------------------------------------

#set up paths and file names for reductions
redDir = 'picomethyl/poolTest/'


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
#(*note comparing to one with replicates instead of pools)
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
plot_grid(plotlist=pltList, nrow=2)


#get csums
cdat = read_tsv('picomethyl/poolTest/pipeline_counts.txt',
                col_names=c('run', 'count', 'step')) %>% 
  filter(step=='mappedCount')
csums = cdat %>% 
  pull(count)
names(csums)=cdat$run
csums


#SAVE RESULTS
save(pm.dat, csums, file='picomethyl/poolTest/countReducedDifferences.Rdata')


