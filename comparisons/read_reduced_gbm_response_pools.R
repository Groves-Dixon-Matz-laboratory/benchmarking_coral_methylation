#read_reduced_gbm_response.R
rm(list=ls())
library(tidyverse)
library(cowplot)
source('benchmarking_functions.R')


# LOAD DATA ---------------------------------------------------------------

ll=load('comparisons/datasets/genotype_gbm_response.Rdata')
ll

#format and subset the gbm results
#note that methylKit doesn't actually use log2fold differences, 
#but changing to match DESeq output for coding
gbm.dat = g.gbm.dat
head(gbm.dat)
main.gbm = gbm.dat %>% 
  mutate(name=as.character(name)) %>% 
  dplyr::rename(pm.pvalue=pvalue,
                pm.padj=qvalue,
                pm.log2FoldChange = meth.diff) %>% 
  dplyr::select(chr, start, end, name, pm.pvalue, pm.padj, pm.log2FoldChange, mr.log2FoldChange, mr.pvalue, mr.padj, mbd.log2FoldChange, mbd.pvalue, mbd.padj)


#double-check correlations
cor(main.gbm[,c('pm.log2FoldChange', 'mr.log2FoldChange', 'mbd.log2FoldChange')],use='complete.obs')
mr = prep_for_cor(main.gbm$mr.log2FoldChange)
mbd = prep_for_cor(main.gbm$mbd.log2FoldChange)
pm = prep_for_cor(main.gbm$pm.log2FoldChange)
cor(mr,mbd, use='complete.obs')
cor(pm,mbd, use='complete.obs')


# SET SOME GLOBAL VARIABLES -----------------------------------------------

minToPlot = 300 #minimum number of complete observations going into a correlation to plot it
minSigToPlot = 20 #minimum number of significant genes to plot uncorroborated proportion
p.false.cut = 0.3 #raw p-value cutoff to infer a false positive if other two assay are greater
xCols = c('pm.log2FoldChange',
          'mbd.log2FoldChange',
          'mr.log2FoldChange')
nameSwaps = c('WGBS',
              'MBD',
              'MR')
names(nameSwaps)=xCols
assayAbrevs = c('pm',
                'mbd',
                'mr')
USE_MEDIAN=TRUE
#set up colors
pm.color = gg_color_hue(6)[1]
mbd.color = gg_color_hue(6)[2]
mr.color = gg_color_hue(6)[3]
alt2.color = gg_color_hue(6)[4]
plus2.color = gg_color_hue(6)[5]
all3.color = gg_color_hue(6)[6]
color.list = c(pm.color,
                mbd.color,
                mr.color,
                alt2.color,
                plus2.color,
                all3.color)
a.color.set=c(pm.color,
            mbd.color,
            mr.color)
a.shape.set=c(19,17,15)


# MBD REDUCED PLOTS CHANGE CORRELATION WITH UNBOUND ------------------------------------

#upload the methylRAD reduced read measures
ll=load('mbdSeq/poolTest/countReducedDifferencesWub.Rdata')
ll
mbd.csums = csums
mbd.dat = mbd.dat %>% 
  full_join(main.gbm, by = 'name')

#get correlation reduction
mbdCorRed = get_correlation_reduction(mbd.dat, xCols, mbd.csums, USE_MEDIAN) %>% 
  swap_names(nameSwaps)

#plot
mbd.diff.cors = plot_correlation_reduction(mbdCorRed, 'MBD-seq', '12 pools\n(6 captured, 6 unbound)', a.color.set)


#plot the total genes interrogated drop
mbd.only=mbd.dat %>% 
  dplyr::select(-mr.padj, -mbd.padj)
mbd.n.tested = get_n_tested(mbd.only, mbd.csums, useMedian=TRUE)

mbd.nt.plt = mbd.n.tested %>% 
  ggplot(aes(x=medReads/1e6,y=nTested)) +
  geom_point(color=mbd.color) +
  geom_line(color=mbd.color) +
  labs(title='MBD-seq')


# WGBS REDUCED PLOTS CHANGE CORRELATION ------------------------------------

#upload the methylRAD reduced read measures
ll=load('picomethyl/poolTest/countReducedDifferences.Rdata')
ll
pm.csums=csums
pm.dat = pm.dat %>% 
  mutate(name=as.character(name)) %>% 
  full_join(main.gbm, by = 'name')

#change the methylKit notation to resemble DESeq2
cn = colnames(pm.dat)
cnmod = sub('_meth.diff', '_log2FoldChange', cn)
cnmod = sub('_qvalue', '_padj', cnmod)
colnames(pm.dat) = cnmod

#get correlation reduction
pmCorRed = get_correlation_reduction(pm.dat, xCols, csums, USE_MEDIAN) %>% 
  swap_names(nameSwaps)


#plot the correlation drop
pm.diff.cors = plot_correlation_reduction(pmCorRed, 'WGBS', '6 pools\n  ', a.color.set)


#plot the total genes interrogated drop
pmOnly = pm.dat %>% 
  dplyr::select(-mbd.padj, -mr.padj)
pm.n.tested = get_n_tested(pmOnly, csums, useMedian=TRUE) %>% 
  arrange(prop)

pm.nt.plt = pm.n.tested %>% 
  ggplot(aes(x=medReads/1e6,y=nTested)) +
  geom_point(color=pm.color) +
  geom_line(color=pm.color) +
  labs(title='WGBS')


# PLOT CORR REDUCTION TOGETHER --------------------------------------------

#set up plots
fpltList0 = list(pm.diff.cors,
                 mbd.diff.cors)
fpltList = lapply(fpltList0, function(x) return(x + 
                                                  lims(y=c(0,1)) +
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        legend.position = 'none')))

#build legend
l.cor = cowplot::get_legend(mbd.diff.cors + labs(color='comparison', shape='comparison'))

#plot
corPlts = plot_grid(plotlist = fpltList, nrow=1, labels = LETTERS[1:3])
corPlts2 = plot_grid(corPlts, l.cor, nrow=1, rel_widths=c(5,1))
xlab = plot_grid(ggdraw() + draw_label('reads per pool (millions)'))
ylab = ggdraw() + draw_label('correlation', angle=90)
top = plot_grid(ylab, corPlts2, nrow=1, rel_widths=c(0.04, 1))
full = plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.1))
full



#PRINT OUT STATS
pm.n = length(pm.csums)
pm.select = 40
mbd.n = length(mbd.csums)
mbd.select = 8
pm.tot = pm.n*pm.select
mbd.tot = mbd.n*mbd.select
pm.tot
mbd.tot


#STUFF BELOW JUST FOR REFERENCE






# PLOT NUMBER OF GENES TESTED ---------------------------------------------

ntList0 = list(pm.nt.plt,
              mbd.nt.plt)
maxGenes = max(nrow(mbd.dat),
               nrow(pm.dat))
clean_nts = function(x){
  clean=x + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    lims(y=c(0,maxGenes)) +
    labs(title='')
}
ntList = lapply(ntList0, function(x) clean_nts(x))


ntPlts=plot_grid(plotlist=ntList, nrow=1)
ntPlts2=plot_grid(ntPlts, legend, nrow=1, rel_widths=c(5,1))
xlab = plot_grid(ggdraw() + draw_label('reads per sample (millions)'))
ylab = ggdraw() + draw_label('genes tested', angle=90)
top = plot_grid(ylab, ntPlts2, nrow=1, rel_widths=c(0.04, 1))
ntFull = plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.1))
ntFull



# SET UP SIGNIFICANT REGION SETS ------------------------------------------

CUT=0.1

#individual assays
sig.pm = main.gbm %>% 
  filter(pm.padj < CUT) %>% 
  pull(name)
sig.mbd = main.gbm %>% 
  filter(mbd.padj < CUT) %>% 
  pull(name)
sig.mr = main.gbm %>% 
  filter(mr.padj < CUT) %>% 
  pull(name)

#pairs
sig.pm.mbd = sig.pm[sig.pm %in% sig.mbd]
sig.pm.mr = sig.pm[sig.pm %in% sig.mr]
sig.mbd.mr = sig.mbd[sig.mbd %in% sig.mr]

#any two-way
sdf = data.frame(pm=main.gbm$pm.padj,
                 mbd=main.gbm$mbd.padj,
                 mr=main.gbm$mr.padj)
rownames(sdf)=main.gbm$name
nsig = apply(sdf, 1, function(x) return(sum(x<CUT, na.rm=TRUE)))
sig.2 = names(nsig)[nsig>1]

#threeway
sig.all = sig.pm.mbd[sig.pm.mbd %in% sig.mr]

#nSig table
sigNames = c('WGBS',
             'MBD',
             'MR',
             'WGBS & MBD',
             'WGBS & MR',
             'MBD & MR',
             'any 2',
             'all 3')
sigSets = list(sig.pm,
               sig.mbd,
               sig.mr,
               sig.pm.mbd,
               sig.pm.mr,
               sig.mbd.mr,
               sig.2,
               sig.all)
names(sigSets)=sigNames
sigNs = sapply(sigSets, length)
sigDf = data.frame(group = sigNames,
                   Nsig = sigNs)
alternate.pairings = c('WGBS & MBD',
                       'WGBS & MR',
                       'MBD & MR') #these go in as the two 'other' assays for each individual. Parsed in get_significance_correspondance 
alt.string='alt. 2'
color.names = c('WGBS',
                'MBD',
                'MR',
                alt.string,
                'any 2',
                'all 3')
names(color.list) = color.names
shape.list=c(19,17,15,3,7,8)
names(shape.list)=color.names
color.list
shape.list


#CALL LIKELY FALSE POSITIVES
pm.fps = infer_false_positives('pm', sig.pm, p.false.cut)
pm.dir.fps = pm.fps[['direction']]
pm.either.fps = pm.fps[['either']]
mbd.fps = infer_false_positives('mbd', sig.mbd, p.false.cut)
mbd.dir.fps = mbd.fps[['direction']]
mbd.either.fps = mbd.fps[['either']]
# mr.fps = infer_false_positives('mr', sig.mr, p.false.cut)
# mr.dir.fps = mr.fps[['direction']]
# mr.either.fps = mr.fps[['either']]


# MBD REDUCED SENSITIVITY -------------------------------------------------
head(mbd.dat)

#get sensitivity and precision estimates
res2=get_significance_correspondance(mbd.dat, sigSets, mbd.dir.fps, mbd.either.fps, mbd.csums, useMedian=TRUE)
res3 = res2 %>% 
  filter(!grepl('MBD', group)) %>% 
  na.omit()

#sensitivity
mbd.sens = plot_sensitivity(res3, 'MBD', color.list, shape.list)

#inferred false positives
mbd.fps = plot_fp_rate(res3, a.shape.set[2])

#2-way to false positive ratio
mbd.ratio = plot_tp_to_fp_ratio(res3, color.choice=mbd.color, shape.choice=a.shape.set[2])

# WGBS REDUCED SENSITIVITY -------------------------------------------------
head(pm.dat)

#get sensitivity and precision estimates
res2=get_significance_correspondance(pm.dat, sigSets, pm.dir.fps, pm.either.fps, pm.csums, useMedian=TRUE)
res3 = res2 %>% 
  filter(!grepl('WGBS', group)) %>% 
  na.omit()

#sensitivity
pm.sens = plot_sensitivity(res3, 'WGBS', color.list, shape.list)

#inferred false positives
pm.fps = plot_fp_rate(res3, a.shape.set[1])

#2-way to false positive ratio
pm.ratio = plot_tp_to_fp_ratio(res3, color.choice=pm.color, shape.choice = a.shape.set[1])


# PLOT SENSITIVITY PLOTS --------------------------------------------------

sens.list0 = list(pm.sens,
                mbd.sens)
clean_sens = function(x){
  z=x+theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position='none') +
    lims(y=c(0,1))
}
sense.list = lapply(sens.list0, function(x) clean_sens(x))

#build legend
r.shape.list=shape.list[c('all 3', 'any 2', 'alt. 2', 'WGBS', 'MBD', 'MR')]
l.sens = plot_sensitivity(res2, 'none', color.list, r.shape.list) + labs(color='comparison',
                                                                         shape='comparison')
l.sens = cowplot::get_legend(l.sens)


# PLOT FALSE POSITIVES ----------------------------------------------------

fp.list0 = list(pm.fps,
               mbd.fps)
clean_fp = function(x){
  z=x+theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position='none') +
    lims(y=c(0,0.17))
  return(z)
}
fp.list = lapply(fp.list0, function(x) clean_fp(x))

#build legend
r.shape.list=shape.list[c('all 3', 'any 2', 'alt. 2', 'WGBS', 'MBD', 'MR')]
l.fp= plot_fp_rate(res3) +
  scale_color_manual(values=c('red', 'firebrick'),
                     labels=c('direction',
                              'p-value')) +
  labs(color='false\npositive')
l.fp = cowplot::get_legend(l.fp)



# PLOT FALSE POSITIVE RATIO -----------------------------------------------

ratio.list0 = list(pm.ratio,
                  mbd.ratio)
clean_ratio = function(x){
  z=x+theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position='none') +
    lims(y=c(0,7.5))
  return(z)
}
ratio.list = lapply(ratio.list0, function(x) clean_ratio(x))

#build legend
l.rat = cowplot::get_legend(mbd.diff.cors + labs(color='asssay', shape='asssay'))

# FINAL PLOT --------------------------------------------------------------

#check initial plots
corPlts0 = plot_grid(plotlist = fpltList, nrow=1, labels = LETTERS[1:3])
sensPlts0 = plot_grid(plotlist=sense.list, nrow=1, labels = LETTERS[4:6])
fpPlts0 = plot_grid(plotlist=fp.list, nrow=1, labels = LETTERS[7:9])
ratPlts0 = plot_grid(plotlist=ratio.list, nrow=1, labels = LETTERS[10:12])

#build without axes
fpltList2 = fpltList %>% 
  remove_y() %>% 
  remove_x()
sense.list2 = sense.list %>% 
  remove_y() %>% 
  remove_x()
fp.list2 = fp.list %>% 
  remove_y() %>% 
  remove_x()
ratio.list2 = remove_y(ratio.list)
ratio.list2[[1]] = ratio.list2[[1]] + scale_y_continuous(breaks=c(0,2,4,6), labels=c('0.00', '2.00', '4.00', '6.00'))

#make new panels
labelx = c(0,-0.05, -0.05)
rel.widths = c(1.2,1,1)
corPlts0 = plot_grid(plotlist = fpltList2, nrow=1, labels = LETTERS[1:3], label_x= labelx, rel_widths=rel.widths)
sensPlts0 = plot_grid(plotlist=sense.list2, nrow=1, labels = LETTERS[4:6], label_x= labelx, rel_widths=rel.widths)
fpPlts0 = plot_grid(plotlist=fp.list2, nrow=1, labels = LETTERS[7:9], label_x= labelx, rel_widths=rel.widths)
ratPlts0 = plot_grid(plotlist=ratio.list2, nrow=1, labels = LETTERS[10:12], label_x= labelx, rel_widths=rel.widths)

#add legends
corPlts1 = plot_grid(corPlts0, l.cor, nrow=1, rel_widths=c(5,1))
sensPlts1 = plot_grid(sensPlts0, l.sens, nrow=1, rel_widths=c(5,1))
fpPlts1 = plot_grid(fpPlts0, l.fp, nrow=1, rel_widths=c(5,1))
ratPlts1 = plot_grid(ratPlts0, l.rat, nrow=1, rel_widths=c(5,1))

#ADD Y AXIS TITLES
#correlation
corY = ggdraw() + draw_label('correlation', angle=90)
corPlts2 = plot_grid(corY, corPlts1, nrow=1, rel_widths=c(0.04, 1))
#sensitivity
sensY = ggdraw() + draw_label('sensitivity', angle=90)
sensPlts2 = plot_grid(sensY, sensPlts1, nrow=1, rel_widths=c(0.04, 1))
#false positive
fpY = ggdraw() + draw_label('inferred FP', angle=90)
fpPlts2 = plot_grid(fpY, fpPlts1, nrow=1, rel_widths=c(0.04, 1))
#ratio 
ratY = ggdraw() + draw_label('any 2 : inferred FP', angle=90)
ratPlts2 = plot_grid(ratY, ratPlts1, nrow=1, rel_widths=c(0.04, 1))

#ADD X AXIS
allX = ggdraw() + draw_label('median read count per sample (millions)')


#FINAL PLOT
top = plot_grid(corPlts2, sensPlts2, fpPlts2, ratPlts2, nrow=4, align='h', rel_heights = c(1.2, 1, 1, 1))
plot_grid(top, allX, nrow=2,  rel_heights=c(30, 1))

