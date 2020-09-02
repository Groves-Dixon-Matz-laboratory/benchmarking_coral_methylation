#read_reduced_gbm_response.R
rm(list=ls())
source('benchmarking_functions.R')
source('figure_plotting/read_reduction_functions.R')

#SELECT WHICH WAY TO PLOT:
#full mdRAD full stats
source('figure_plotting/linePltVersions/fullMdRAD_fullCompare.R')
#full mdRAD simple stats
source('figure_plotting/linePltVersions/fullMdRAD_simpleCompare.R')
#8 sample mdRAD full stats
source('figure_plotting/linePltVersions/8MdRAD_fullCompare.R')
#8 samplw mdRAD simple stats
source('figure_plotting/linePltVersions/8MdRAD_simpleCompare.R')

# LOAD DATA ---------------------------------------------------------------

ll=load('comparisons/datasets/genotype_gbm_response.Rdata')
ll

# #change the full mdRAD dataset results to the 8-sample reduction
# #optionally comment this to run compared to the full dataset
# g.gbm.dat$mr.log2FoldChange = g.gbm.dat$esmr.log2FoldChange
# g.gbm.dat$mr.pvalue = g.gbm.dat$esmr.pvalue
# g.gbm.dat$mr.padj = g.gbm.dat$esmr.padj



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

PCUT=0.1
RCUT = 1.96

minToPlot = 300 #minimum number of complete observations going into a correlation to plot it
minSigToPlot = 20 #minimum number of significant genes to plot uncorroborated proportion
p.false.cut = 0.3 #raw p-value cutoff to infer a false positive if other two assay are greater
xCols = c('pm.log2FoldChange',
          'mbd.log2FoldChange',
          'mr.log2FoldChange')
nameSwaps = c('WGBS',
              'MBD',
              'mdRAD')
names(nameSwaps)=xCols
assayAbrevs = c('pm',
                'mbd',
                'mr')
USE_MEDIAN=TRUE


#PCA TO SUMMARIZE METHYLATION --------------------------------------------

#SCALE DATASET AND GET PCA
log2folds = main.gbm %>%
  dplyr::select(pm.log2FoldChange, mr.log2FoldChange, mbd.log2FoldChange)
scaled = map_dfc(log2folds, ~ (.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)) %>%
  data.frame()
rownames(scaled)=main.gbm$name
scaled_nona = na.omit(scaled)
colnames(scaled_nona)=c('pm', 'mr', 'mbd')

#BUILD PCA
mpca = prcomp(scaled_nona)
percentVar <- mpca$sdev^2/sum(mpca$sdev^2)
mcoords = mpca$x %>%
  data.frame() %>%
  mutate(PC1=PC1*-1)
pdat = cbind(scaled_nona, mcoords)

#GET RESIDUALS
get_resid = function(var) {
  dat = dplyr::select(pdat, PC1, y = {{ var }})
  as.numeric(lm(y ~ PC1, data = dat)$resid)
}

rdat = pdat %>%
  mutate(pm.resid = get_resid(pm),
         mr.resid = get_resid(mr),
         mbd.resid = get_resid(mbd),
         name=rownames(scaled_nona)) %>%
  as_tibble()

#CALL FALSE POSITIVES
fpdat = rdat %>%
  mutate(pm.fp = abs(pm.resid) > RCUT,
         mr.fp = abs(mr.resid) > RCUT,
         mbd.fp = abs(mbd.resid) > RCUT)

pm.fps = fpdat %>%
  filter(pm.fp) %>%
  pull(name)

mbd.fps = fpdat %>%
  filter(mbd.fp) %>%
  pull(name)

mr.fps = fpdat %>%
  filter(mr.fp) %>%
  pull(name)


#PLOT CORRELATIONS WITH PC1 AND INFERRED FALSE POSITIVES
pct=round(percentVar[1], digits=2)*100
XLAB=paste('PC1 (', pct, '% of variance)', sep='')
pm.plt = plot_scatter_pearsonCor_annotated(fpdat, 'PC1', 'pm', 'PC1', 'WGBS', ALPHA=1) +
  geom_point(aes(color=pm.fp)) +
  scale_color_manual(values=c('black', 'red')) +
  labs(color=paste('|residual| > ', RCUT, sep=''),
       title='WGBS',
       subtitle=paste('N=', sum(fpdat$pm.fp))) +
  theme(plot.subtitle=element_text(color='red'),
        plot.title=element_text(hjust=0.5))
mbd.plt = plot_scatter_pearsonCor_annotated(fpdat, 'PC1', 'mbd', 'PC1', 'MBD-Seq', ALPHA=1) +
  geom_point(aes(color=mbd.fp)) +
  scale_color_manual(values=c('black', 'red')) +
  labs(color='|residual| > 1.96',
       title='MBD-seq',
       subtitle=paste('N=', sum(fpdat$mbd.fp))) +
  theme(plot.subtitle=element_text(color='red'),
        plot.title = element_text(hjust=0.5))
mr.plt = plot_scatter_pearsonCor_annotated(fpdat, 'PC1', 'mr', 'PC1', 'mdRAD', ALPHA=1) +
  geom_point(aes(color=mr.fp)) +
  scale_color_manual(values=c('black', 'red')) +
  labs(color='|residual| > 1.96',
       title='MethylRAD',
       subtitle=paste('N=', sum(fpdat$mr.fp))) +
  theme(plot.subtitle=element_text(color='red'),
        plot.title = element_text(hjust=0.5))
l=cowplot::get_legend(pm.plt)
pltList0 = list(pm.plt, mbd.plt, mr.plt)
pltList1 = lapply(pltList0, function(x) return(x + theme(axis.title.x = element_blank(), legend.position = 'none')))
plt=plot_grid(plotlist = pltList1, nrow=1)
xlab = ggdraw() + draw_label(paste('PC1 (', pct, '% of variance explained)', sep=''))
plt2=plot_grid(plt, xlab, nrow=2, rel_heights = c(12,1))
plot_grid(plt2, l, rel_widths=c(5,1))


#add PC1 to main.gbm
main.gbm = main.gbm %>%
  left_join(fpdat[,c('name', 'PC1', 'pm.fp', 'mr.fp', 'mbd.fp')], by = 'name')

# METHYLRAD REDUCED PLOTS CHANGE CORRELATION -------------------------------

#Handle selection of full or 8-sample subset
if (FULL_mdRAD){
  mrSamples='12 samples'
  mrLibs = '24 libraries'
  ll=load('methylRAD/datasets/countReducedDifferencesFULL.Rdata')
} else {
  mrSamples='8 samples'
  mrLibs = '8 libraries'
  ll=load('methylRAD/datasets/countReducedDifferences.Rdata')
}
#set csums to toal per sample
mr.csums = csums
names(mr.csums) = sub('mrF.', '', names(Fsums), fixed=TRUE)
length(mr.csums) #n samples
millReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
mdTitle=paste(mrSamples, mrLibs, millReads, sep='\n')

#merge in the data
mr.dat = mr.dat %>% 
  mutate(name=as.character(name)) %>% 
  full_join(main.gbm, by = 'name')

#get correlation reduction
mrCorRed = get_correlation_reduction(mr.dat, xCols, mr.csums, USE_MEDIAN) %>% 
  swap_names(nameSwaps)

#build plot

t=paste(mrSamples, mrLibs, millReads, sep='\n')
mr.diff.cors = plot_correlation_reduction(mrCorRed, 'mdRAD',t, a.color.set)

#plot the total genes interrogated drop
mr.only=mr.dat %>% 
  dplyr::select(-mr.padj, -mbd.padj)

# mr.n.tested = get_n_tested(mr.only, mr.csums, useMedian=TRUE)
# mr.nt.plt = mr.n.tested %>% 
#   ggplot(aes(x=medReads/1e6,y=nTested)) +
#   geom_point(color=mr.color) +
#   geom_line(color=mr.color) +
#   labs(title='MethylRAD')


# MBD REDUCED PLOTS CHANGE CORRELATION WITH UNBOUND ------------------------------------

#upload the methylRAD reduced read measures
ll=load('mbdSeq/datasets/countReducedDifferencesWub.Rdata')
ll
mbd.dat = mbd.dat %>% 
  full_join(main.gbm, by = 'name')
millReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
mbdTitle=paste('8 samples', '16 libraries', millReads, sep='\n')
mCsums = csums[grep('^m', names(csums))]
uCsums = csums[grep('^ub', names(csums))]
mbd.csums = mCsums+uCsums
names(mbd.csums)=sub('ub-', '', names(uCsums))


#get correlation reduction
mbdCorRed = get_correlation_reduction(mbd.dat, xCols, mbd.csums, USE_MEDIAN) %>% 
  swap_names(nameSwaps)

#plot
mbd.diff.cors = plot_correlation_reduction(mbdCorRed, 'MBD-seq', mbdTitle, a.color.set)


#plot the total genes interrogated drop
mbd.only=mbd.dat %>% 
  dplyr::select(-mr.padj, -mbd.padj)



# mbd.n.tested = get_n_tested(mbd.only, mbd.csums, useMedian=TRUE)
# 
# mbd.nt.plt = mbd.n.tested %>% 
#   ggplot(aes(x=medReads/1e6,y=nTested)) +
#   geom_point(color=mbd.color) +
#   geom_line(color=mbd.color) +
#   labs(title='MBD-seq')


# WGBS REDUCED PLOTS CHANGE CORRELATION ------------------------------------

#upload the methylRAD reduced read measures
ll=load('picomethyl/datasets/readReduced_gbm_diff/countReducedDifferences.Rdata')
ll
pm.csums=csums
pm.dat = pm.dat %>% 
  mutate(name=as.character(name)) %>% 
  full_join(main.gbm, by = 'name')
pmMillReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
pmTitle=paste('8 samples', '8 libraries', pmMillReads, sep='\n')


#change the methylKit notation to resemble DESeq2
cn = colnames(pm.dat)
cnmod = sub('_meth.diff', '_log2FoldChange', cn)
cnmod = sub('_qvalue', '_padj', cnmod)
colnames(pm.dat) = cnmod

#get correlation reduction
pmCorRed = get_correlation_reduction(pm.dat, xCols, csums, USE_MEDIAN) %>% 
  swap_names(nameSwaps)


#plot the correlation drop
pm.diff.cors = plot_correlation_reduction(pmCorRed, 'WGBS', pmTitle, a.color.set)


#plot the total genes interrogated drop
pmOnly = pm.dat %>% 
  dplyr::select(-mbd.padj, -mr.padj)



# pm.n.tested = get_n_tested(pmOnly, csums, useMedian=TRUE) %>% 
#   arrange(prop)
# 
# pm.nt.plt = pm.n.tested %>% 
#   ggplot(aes(x=medReads/1e6,y=nTested)) +
#   geom_point(color=pm.color) +
#   geom_line(color=pm.color) +
#   labs(title='WGBS')


# PLOT CORR REDUCTION TOGETHER --------------------------------------------

#set up plots
fpltList0 = list(pm.diff.cors,
                 mbd.diff.cors,
                 mr.diff.cors)
fpltList = lapply(fpltList0, function(x) return(x + 
                                                  lims(y=c(0,1)) +
                                                  theme(axis.title.x = element_blank(),
                                                        axis.title.y = element_blank(),
                                                        plot.title = element_text(hjust=0.5),
                                                        legend.position = 'none')))

#build legend
l.cor = cowplot::get_legend(mbd.diff.cors + labs(color='comparison', shape='comparison'))

#plot
corPlts = plot_grid(plotlist = fpltList, nrow=1, labels = LETTERS[1:3])
corPlts2 = plot_grid(corPlts, l.cor, nrow=1, rel_widths=c(5,1))
xlab = plot_grid(ggdraw() + draw_label('reads per sample (millions)'))
ylab = ggdraw() + draw_label('correlation', angle=90)
top = plot_grid(ylab, corPlts2, nrow=1, rel_widths=c(0.04, 1))
full = plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.1))
full



# SET UP SIGNIFICANT REGION SETS ------------------------------------------

#individual assays
sig.pm = main.gbm %>% 
  filter(pm.padj < PCUT) %>% 
  pull(name)
sig.mbd = main.gbm %>% 
  filter(mbd.padj < PCUT) %>% 
  pull(name)
sig.mr = main.gbm %>% 
  filter(mr.padj < PCUT) %>% 
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
nsig = apply(sdf, 1, function(x) return(sum(x<PCUT, na.rm=TRUE)))
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


# CALL LIKELY FALSE POSITIVES ---------------------------------------------


#make new df with wrong calls
pm.ww = add_wrong_to_df(pm.dat)
mbd.ww = add_wrong_to_df(mbd.dat)
mr.ww = add_wrong_to_df(mr.dat)

  


# METHYLRAD REDUCED SENSITIVITY AND FPS -------------------------------------------

#get sensitivity and precision estimates
sdf0=get_significance_correspondance(mr.ww,
                                     sigSets,
                                     mr.csums,
                                     useMedian=TRUE,
                                     alternate.pairings=alternate.pairings)
sdf = sdf0 %>% 
  # filter(!grepl('MR', group)) %>% 
  na.omit()

#sensitivity
mr.sens = plot_sensitivity(sdf, 'MR', color.list, shape.list)

#precision
mr.prec = plot_precision(sdf)

#plot multi sensitivity
mr.sens = plot_sig_agreement(sdf, 'sensitivity', "WGBS & MBD", color.list, shape.list, groups_to_exclude = EXCLUDE_GROUPS)

#plot multi precision
mr.prec = plot_sig_agreement(sdf, 'precision', "WGBS & MBD", color.list, shape.list, groups_to_exclude = EXCLUDE_GROUPS)



#false pos estimate
mr.fdf = gather_fp_data(mr.ww, mr.csums)
mr.fp = plot_false_pos(mr.fdf)





#number tested
mr.ntest = plot_n_tested(mr.fdf, color.list[['MR']], 'MethylRAD')


# MBD REDUCED SENSITIVITY -------------------------------------------------

#get sensitivity and precision estimates
res2=get_significance_correspondance(mbd.ww, sigSets, mbd.csums, useMedian=TRUE)
res3 = res2 %>% 
  # filter(!grepl('MBD', group)) %>%
  na.omit()

#sensitivity
mbd.sens = plot_sensitivity(res3, 'MBD', color.list, shape.list)

#precision
mbd.prec = plot_precision(res3)

#plot multi sensitivity
mbd.sens = plot_sig_agreement(res3,
                              'sensitivity',
                              "WGBS & MR",
                              color.list,
                              shape.list,
                              groups_to_exclude = EXCLUDE_GROUPS)

#plot multi precision
mbd.prec = plot_sig_agreement(res3,
                              'precision',
                              "WGBS & MR",
                              color.list,
                              shape.list,
                              groups_to_exclude = EXCLUDE_GROUPS)

#false pos estimate
mbd.fdf = gather_fp_data(mbd.ww, mbd.csums)
mbd.fp = plot_false_pos(mbd.fdf)

#number tested
mbd.ntest = plot_n_tested(mbd.fdf, color.list[['MBD']], 'MBD-seq')


# WGBS REDUCED SENSITIVITY -------------------------------------------------

#get sensitivity and precision estimates
res2=get_significance_correspondance(pm.ww, sigSets, pm.csums, useMedian=TRUE)
res3 = res2 %>% 
  # filter(!grepl('WGBS', group)) %>% 
  na.omit()

#sensitivity
pm.sens = plot_sensitivity(res3, 'WGBS', color.list, shape.list)

#precision
pm.prec = plot_precision(res3)

#plot multi sensitivity
alternate.pairings
pm.sens = plot_sig_agreement(res3, 'sensitivity', "MBD & MR", color.list, shape.list, groups_to_exclude = EXCLUDE_GROUPS)

#plot multi precision
pm.prec = plot_sig_agreement(res3, 'precision', "MBD & MR", color.list, shape.list, groups_to_exclude = EXCLUDE_GROUPS)

#false pos estimate
pm.fdf = gather_fp_data(pm.ww, pm.csums)
pm.fp = plot_false_pos(pm.fdf)

#number tested
pm.ntest = plot_n_tested(pm.fdf, color.list[['WGBS']], 'WGBS')

# PLOT SENSITIVITY PLOTS --------------------------------------------------

sens.list0 = list(pm.sens,
                  mbd.sens,
                  mr.sens)
clean_sens = function(x){
  z=x+theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position='none') +
    lims(y=c(0,1))
}
sense.list = lapply(sens.list0, function(x) clean_sens(x))

#build legend
l.sens = cowplot::get_legend(pm.prec+labs(color='comparison',
                                          shape='comparison'))


# PLOT PRECISION PLOTS ----------------------------------------------------

precYlims = c(0, 1)
prec.list0 = list(pm.prec,
                   mbd.prec,
                   mr.prec)
clean_prec = function(x){
  z=x+theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position='none') +
    lims(y=precYlims)
  return(z)
}
prec.list = lapply(prec.list0, function(x) clean_prec(x))

#build legend ()
# p.shape.list=shape.list[c('any 2')]
# p.color.list = list('firebrick')
# l.sens = plot_sensitivity(res2, 'none', color.list, r.shape.list) + labs(color='comparison',
#                                                                          shape='comparison')
# l.sens = cowplot::get_legend(l.sens)






# PLOT FALSE POSITIVE RATE -----------------------------------------------
ratioYlims = c(0, 0.03)
ratio.list0 = list(pm.fp,
                   mbd.fp,
                   mr.fp)
clean_ratio = function(x){
  z=x+theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position='none') +
    lims(y=ratioYlims)
  return(z)
}
ratio.list = lapply(ratio.list0, function(x) clean_ratio(x))

#build legend
l.rat = cowplot::get_legend(ratio.list0[[1]])

# FINAL PLOT --------------------------------------------------------------

#build without axes
fpltList2 = fpltList %>% 
  remove_y() %>% 
  remove_x()
sense.list2 = sense.list %>% 
  remove_y() %>% 
  remove_x()
sense.list2[[1]] = sense.list2[[1]] + scale_y_continuous(breaks = seq(0,1,0.25),
                                                         labels=c('0.00', '0.25', '0.50', '0.75', ''),
                                                         limits=c(0,1))

prec.list1.5 = remove_y(prec.list)
prec.list2 = lapply(prec.list1.5, function(x) return(x + scale_y_continuous(breaks=c(0, 0.3, 0.6,0.9),
                                                                            labels = c('0.00', '0.30', '0.60','0.90'),
                                                                            limits=c(0,1))))

prec.list2[[1]] = prec.list2[[1]] + 
  scale_y_continuous(breaks=c(0, 0.3, 0.6,0.9),
                     labels = c('0.00', '0.30', '0.60','0.90'),
                     limits=c(0,1)) +
  scale_x_continuous(breaks=c(0, 15, 30,45))
prec.list2[[2]] = prec.list2[[2]] + 
  scale_x_continuous(breaks=c(0, 15, 30,45))
  
  
# ratio.list2[[2]] = ratio.list2[[2]] + scale_x_continuous(breaks=c(0, 3, 6, 9))
# 


#make new panels
labelx = c(0,-0.05, -0.05)
rel.widths = c(1.2,1,1)
corPlts0 = plot_grid(plotlist = fpltList2, nrow=1, labels = LETTERS[1:3], label_x= labelx, rel_widths=rel.widths)
sensPlts0 = plot_grid(plotlist=sense.list2, nrow=1, labels = LETTERS[4:6], label_x= labelx, rel_widths=rel.widths)
precPlts0 = plot_grid(plotlist=prec.list2, nrow=1, labels = LETTERS[7:9], label_x= labelx, rel_widths=rel.widths)
# ratPlts0 = plot_grid(plotlist=ratio.list2, nrow=1, labels = LETTERS[7:9], label_x= labelx, rel_widths=rel.widths)

#add legends
corPlts1 = plot_grid(corPlts0, l.cor, nrow=1, rel_widths=c(5,1))
sensPlts1 = plot_grid(sensPlts0, l.sens, nrow=1, rel_widths=c(5,1))
precPlts1 = plot_grid(precPlts0, l.sens, ncol=2, rel_widths=c(5,1))
# ratPlts1 = plot_grid(ratPlts0, l.rat, ncol=2, rel_widths=c(5,1))




#ADD Y AXIS TITLES
YLAB_WIDTHS=c(0.08, 1)
#correlation
corY = ggdraw() + draw_label('\ncorrelation', angle=90)
corPlts2 = plot_grid(corY, corPlts1, nrow=1, rel_widths=YLAB_WIDTHS)
#sensitivity
sensY = ggdraw() + draw_label('comparative\nsensitivity', angle=90)
sensPlts2 = plot_grid(sensY, sensPlts1, nrow=1, rel_widths=YLAB_WIDTHS)
#precision 
precY = ggdraw() + draw_label('comparative\nprecision', angle=90)
precPlts2 = plot_grid(precY, precPlts1, nrow=1, rel_widths=YLAB_WIDTHS)
# #ratio 
# ratY = ggdraw() + draw_label('estimated FP rate', angle=90)
# ratPlts2 = plot_grid(ratY, ratPlts1, nrow=1, rel_widths=c(0.04, 1))

#ADD X AXIS
allX = ggdraw() + draw_label('median read count per sample (millions)')


#FINAL PLOT
top = plot_grid(corPlts2, sensPlts2, precPlts2, nrow=3, align='h', rel_heights = c(1.2, 1, 1))
plot_grid(top, allX, nrow=2,  rel_heights=c(20, 1))



# PLOT N TESTED PLOT ------------------------------------------------------

#note only the reductions are plotted
n.testedList = list(pm.ntest, mbd.ntest, mr.ntest+scale_x_continuous(breaks=c(0,2,4,6), limits=c(0,6)))
noy = remove_y(n.testedList)
nox = lapply(noy, function(x) return(x+theme(axis.title.x = element_blank())))
plts = plot_grid(plotlist = nox, nrow=1, rel_widths = c(1,.75, .75))
xlab =  ggdraw() + draw_label('median mapped counts')
full = plot_grid(plts, xlab, nrow=2, rel_heights = c(10,1))
full


