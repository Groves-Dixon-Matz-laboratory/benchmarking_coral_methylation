#correlate_GBM_measures.R

library(tidyverse)
library(cowplot)
library(DESeq)
source('benchmarking_functions.R')
# library(boot)




# load RNAseq -------------------------------------------------------------

ll=load('rnaseq/datasets/tissue_results.Rdata')
rna.t = format_deseq_res(res)
ll=load('rnaseq/datasets/genotype_results.Rdata')
rna.g = format_deseq_res(res)

# load methylRAD data -----------------------------------------------------

format_mr_window = function(mr.res){
  colnames(mr.res) = paste('mr', colnames(mr.res), sep='.')
  data.frame(mr.res) %>% 
    dplyr::select(mr.log2FoldChange, mr.pvalue, mr.padj) %>% 
    mutate(gene=rownames(mr.res))
  
}

ll = load('methylRAD/host/window_results/methylRAD_10kb_window_counts_diffMeth.Rdata')
ll

#format tissue
mr.t = format_mr_window(t.res)
head(mr.t)

#format genotype
mr.g = format_mr_window(g.res)
head(mr.g)



# load picomethyl data ----------------------------------------------------

format_pm_window = function(pmdat){
  colnames(pmdat)=paste('pm', colnames(pmdat), sep='.')
  pmdat %>% 
    tidyr::unite(gene, pm.chr, pm.start, pm.end, sep='_') %>% 
    dplyr::select(gene, pm.pvalue, pm.qvalue, pm.meth.diff)
}

#tissue
ll=load('picomethyl/window_results/tipVside_windows_10000_10000_formatted_diff.Rdata')
ll
pm.t = format_pm_window(pm)
head(pm.t)

#genotype
ll=load('picomethyl/window_results/genotype_windows_10000_10000_formatted_diff.Rdata')
ll
pm.g = format_pm_window(pm)
head(pm.g)


# load mbdseq from recip meth ---------------------------------------------

format_mbd_window = function(mbd.dat){
  colnames(mbd.dat)=paste('mbd', colnames(mbd.dat), sep='.')
  mbd.dat %>% 
    tidyr::unite(gene, mbd.chr, mbd.start, mbd.end, sep='_') %>% 
    dplyr::select(gene, mbd.log2FoldChange, mbd.pvalue, mbd.padj)
}


ll=load('mbdSeq/host/window_results/mbd_10kb_window_diff.Rdata')
ll

mbd.ig = format_mbd_window(mbd.ig)
mbd.bg = format_mbd_window(mbd.bg)
mbd.it = format_mbd_window(mbd.it)
mbd.bt = format_mbd_window(mbd.bt)
head(mbd.bt)



# merge up the results ----------------------------------------------------

tList = list(mr.t,
             pm.t,
             mbd.it,
             mbd.bt)
gList = list(mr.g,
             pm.g,
             mbd.ig,
             mbd.bg)


mtdat = purrr::reduce(tList, full_join, by='gene')
mgdat = purrr::reduce(gList, full_join, by='gene')
head(mtdat)




# plot correlations -------------------------------------------------------


mtdat %>% 
  ggplot(aes(x=pm.meth.diff, y=mbd.log2FoldChange.x)) +
  geom_point() +
  geom_smooth()



mgdat %>% 
  ggplot(aes(x=pm.meth.diff, y=mbd.log2FoldChange.x)) +
  geom_point() +
  geom_smooth(method='lm')

mgdat %>% 
  ggplot(aes(x=pm.meth.diff, y=mbd.log2FoldChange.y)) +
  geom_point() +
  geom_smooth(method='lm')



mgdat %>% 
  ggplot(aes(x=pm.meth.diff, y=mr.log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm')



mgdat %>% 
  ggplot(aes(x=mbd.log2FoldChange.x, y=mr.log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm')


mgdat %>% 
  ggplot(aes(x=mbd.log2FoldChange.y, y=mr.log2FoldChange)) +
  geom_point() +
  geom_smooth(method='lm')












