#correlate_GBM_measures.R

library(tidyverse)
library(cowplot)
library(DESeq)
source('benchmarking_functions.R')
# library(boot)




# load RNAseq -------------------------------------------------------------

ll=load('rnaseq/host/datasets/')
rna.t = format_deseq_res(res)
ll=load('rnaseq/datasets/genotype_results.Rdata')
rna.g = format_deseq_res(res)

# load methylRAD data -----------------------------------------------------

ll = load('methylRAD/gbm_results/deseq_results.Rdata')
ll
mr.t = format_deseq_res(t.res)
mr.g = format_deseq_res(g.res)


# load picomethyl data ----------------------------------------------------

#funciton to load results output from treatEffect_by_glm_geneRobjectInputs.R (see picoMethyl_data_processing_pipeline.txt)
ll=load('metadata/mRNA_to_gene.Rdata')
load_pm_gbm_response = function(filePath, mtg){
  read.table(filePath, col.names=c('transcript', 'coeff', 'pvalue'), stringsAsFactors=FALSE) %>% 
    mutate(coeff = if_else(coeff=='Response_is_constant_0nM',
                           '0',
                           coeff),
           coeff = as.numeric(coeff)) %>% 
    filter(coeff > -10 & coeff < 10) %>% 
    left_join(mtg, by = 'transcript') %>% 
    as_tibble()
}

pm.t = load_pm_gbm_response('picomethyl/gbm_results/tipVside_all_gbm_treat_effects.tsv', mtg)
pm.g = load_pm_gbm_response('picomethyl/gbm_results/genotype_all_gbm_treat_effects.tsv', mtg)


# load mbdseq from recip meth ---------------------------------------------

ll=load('mbdSeq/gbm_results/gbm_responses.Rdata')
ll

imbd.t = format_deseq_res(it.res)
bmbd.t = format_deseq_res(bt.res)
imbd.g = format_deseq_res(ig.res)
bmbd.g = format_deseq_res(bg.res)



# correlate gbm change with ge change -------------------------------------

gdat = pm.g
rdat = rna.g
gbm.col = 'coeff'

ALPHA=0.05

#PLOT AGAINST GENE EXPRESSION

#for tissue
pm.ge.t = gbm_scatter(pm.g, rna.g, 'coeff', 'log2FoldChange', 'PM vs GE tissue effect', 'GBM difference', 'GE difference')
mbd.ge.t = gbm_scatter(imbd.t, rna.t, 'log2FoldChange.x', 'log2FoldChange.y', 'MBD vs GE tissue effect', 'GBM difference', 'GE difference')
mr.ge.t = gbm_scatter(mr.t, rna.t, 'log2FoldChange.x', 'log2FoldChange.y', 'MR vs GE tissue effect', 'GBM difference', 'GE difference')
plot_shared_x(list(pm.ge.t, mbd.ge.t, mr.ge.t), xlab=bquote(log[2]~GE))

#for genotype
pm.ge.g = gbm_scatter(pm.g, rna.g, 'coeff', 'log2FoldChange', 'PM vs GE tissue effect', 'GE difference', 'Picomethyl')
mbd.ge.g = gbm_scatter(bmbd.g %>% filter(gene %in% pm.g$gene), rna.g, 'log2FoldChange.x', 'log2FoldChange.y', 'MBD vs GE genotype effect', 'GE difference', 'MBD')
mr.ge.g = gbm_scatter(mr.g, rna.g, 'log2FoldChange.x', 'log2FoldChange.y', 'MR vs GE genotype effect', 'GE difference', 'methylRAD')
ge.gList=list(pm.ge.g,
             mbd.ge.g,
             mr.ge.g)
plot_shared_x(ge.gList, xlab=bquote(log[2]~GE))


#PLOT METH AGAINST METH

#for tissue
pm.mr.t = gbm_scatter(pm.t, mr.t, 'coeff', 'log2FoldChange', 'PM vs MR tissue effect', 'Picomethyl', 'MethylRAD')
pm.mbd.t = gbm_scatter(pm.t, bmbd.t, 'coeff', 'log2FoldChange', 'PM vs MBD tissue effect', 'Picomethyl', 'MBD-seq')
mbd.mr.t = gbm_scatter(bmbd.t, mr.t, 'log2FoldChange.x', 'log2FoldChange.y', 'MBD vs MR tissue effect', 'MBD', 'MethylRAD')
pltList = list(pm.mr.t,
               pm.mbd.t,
               mbd.mr.t)
plot_grid(pm.mr.t,
          pm.mbd.t,
          mbd.mr.t,
          nrow=1)

#for genotype
pm.mr.g = gbm_scatter(pm.g, mr.g, 'coeff', 'log2FoldChange', 'PM vs MR tissue effect', 'Picomethyl', 'MethylRAD')
pm.mbd.g = gbm_scatter(pm.g, imbd.g, 'coeff', 'log2FoldChange', 'PM vs MBD tissue effect', 'Picomethyl', 'MBD-seq')
mbd.mr.g1 = gbm_scatter(imbd.g, mr.g, 'log2FoldChange.x', 'log2FoldChange.y', 'MBD vs MR tissue effect', 'MBD w/ ub', 'MethylRAD')
mbd.mr.g2 = gbm_scatter(bmbd.g, mr.g, 'log2FoldChange.x', 'log2FoldChange.y', 'MBD vs MR tissue effect', 'MBD no ub', 'MethylRAD')
pltList = list(pm.mr.g,
               pm.mbd.g,
               mbd.mr.g2)
plot_grid(plotlist=pltList,
          nrow=1)

















# temp correlate methylRAD and picomethyl ---------------------------------

ALPHA=0.5
fm %>% 
  left_join(pm, by = 'gene') %>% 
  ggplot(aes(x=mnPM, y=mnMR)) +
  geom_point(alpha=ALPHA)

mm %>% 
  left_join(pm, by = 'gene') %>% 
  ggplot(aes(x=mnPM, y=mnMR)) +
  geom_point(alpha=ALPHA)

bm %>% 
  left_join(pm, by = 'gene') %>% 
  ggplot(aes(x=mnPM, y=mnMR)) +
  geom_point(alpha=ALPHA)




# plot bimodals -----------------------------------------------------------
pmbreaks = logit(c(0.005, 0.25, 0.98))
ybreaks = c(0,500,1500)
pmhist = pdat %>% 
  ggplot(aes(x=picoMethyl)) +
  geom_histogram() +
  scale_x_continuous(labels = inv_logit_labs, breaks = pmbreaks) +
  labs(subtitle='WGBS', x='Methylation %')

mbdhist = mbd %>% 
  ggplot(aes(x=mbdGbm)) +
  geom_histogram() +
  labs(subtitle='MBDseq', x=bquote(log[2]*"(capture vs flowthru)")) +
  theme(axis.title.y = element_blank())
  
mrhist = mmrdat %>% 
  ggplot(aes(x=mnMR)) +
  geom_histogram() +
  labs(subtitle='methylRAD', x=mrAxisLab) +
  theme(axis.title.y = element_blank()) +
  lims(x=c(-6,8)) +
  scale_y_continuous(breaks=c(0,1000,2000))
  

plot_grid(pmhist, mbdhist, mrhist, nrow=1, rel_widths=c(1,.9,.9))

ALPHA=0.015

# plot N genes with data --------------------------------------------------

Nmbd = length(unique(mbd$gene))
Nwgbs = length(unique(pdat$gene))
Nmr = length(unique(na.omit(mmrdat)$gene))

ndat = data.frame(N=c(Nmbd, Nwgbs, Nmr), method=c('MBD', 'WGBS', 'mRAD')) %>% 
  as_tibble()

#how many genes have data from each?
# ndat %>% 
#   ggplot(aes(x=method,y=N)) +
#   geom_bar(stat='identity') +
#   labs(y='N genes with any data')

# plot picomethyl and mbd -------------------------------------------------
pmbreaks = logit(c(0.01, 0.4, 0.98))
wVmbd = mbd %>% 
  full_join(pdat, by = 'gene') %>% 
  ggplot(aes(x=mbdGbm, y=picoMethyl)) +
  geom_point(alpha=ALPHA) +
  # geom_smooth(method='lm') +
  labs(y='Methylation %', x=bquote(log[2]*"(capture vs flowthru)"), subtitle='WGBS ~ MBDseq') +
  scale_y_continuous(labels = inv_logit_labs, breaks=pmbreaks) 
  



# plot methylRAD and mbd -----------------------------------------------
mrBreaks=c(-4,0,4)
mrLims = c(-6,6)
#averaged over all
mrVmbd = mmrdat %>% 
  full_join(mbd, by = 'gene') %>% 
  ggplot(aes(x=mbdGbm, y=mnMR)) +
  geom_point(alpha=ALPHA) +
  # geom_smooth(method='lm') +
  labs(y=mrAxisLab, x=bquote(log[2]*"(capture vs flowthru)"), subtitle='MethylRAD ~ MBDseq') +
  scale_y_continuous(breaks=mrBreaks, limits=mrLims)
  
  

# #split by library
# mrdat %>% 
#   full_join(mbd, by = 'gene') %>% 
#   ggplot(aes(x=mbdGbm, y=mnMR, color=library)) +
#   geom_point(alpha=0.1) +
#   # geom_smooth(method='lm') +
#   labs(y=mrAxisLab, x=bquote(log[2]*"(capture vs flowthru)"), color='enzyme_recognitionSiteSelection', subtitle='MethylRAD ~ MBDseq')
# 

# plot methylRAD and picoMethyl -----------------------------------------------
mrBreaks=c(-4,0,4)
mrLims = c(-6,6)
#averaged over all
mrVw = mmrdat %>% 
  full_join(pdat, by = 'gene') %>% 
  ggplot(aes(x=picoMethyl, y=mnMR)) +
  geom_point(alpha=ALPHA) +
  # geom_smooth(method='lm') +
  labs(y=mrAxisLab, x='Methylation %', subtitle='MethylRAD ~ WGBS') +
  scale_x_continuous(labels = inv_logit_labs, breaks=pmbreaks) +
  scale_y_continuous(breaks=mrBreaks, limits=mrLims) +
  theme(plot.margin=unit(c(5.5,15,5.5,5.5),"pt"))

# #split by library
# mrdat %>% 
#   full_join(pdat, by = 'gene') %>% 
#   ggplot(aes(x=picoMethyl, y=mnMR, color=library)) +
#   geom_point(alpha=0.1) +
#   geom_smooth(method='lm') +
#   labs(y=mrAxisLab, x='WGBS', color='enzyme_recognitionSiteSelection', subtitle='MethylRAD ~ WGBS') +
#   scale_x_continuous(labels = inv_logit_labs) 




#plot all three
plot_grid(wVmbd, mrVmbd, mrVw, nrow=1, rel_widths=c(1,1,1.07))


  