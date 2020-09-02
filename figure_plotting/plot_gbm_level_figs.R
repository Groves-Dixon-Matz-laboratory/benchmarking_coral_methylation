#plot_gbm_level_figs.R
rm(list=ls())
source('benchmarking_functions.R')


# LOAD GBM LEVEL DATA -----------------------------------------------------

#load data from compare_level.R
ll=load('comparisons/datasets/gbmLvl.Rdata')
ll
head(gbm.dat)
aldat = gbm.dat %>% 
  dplyr::select(chr, start, end, name, l.fracMeth, pm.glmLvl, l.mCpG_per_CpG, l.mCpG_per_bp, mrB, mrF, mrM, mrF.s, mrM.s, mrB.s, mbd.score, mbd.FPKM, mbd.FPCpGM,cpgoe)
dim(aldat)

# BUILD GBM HISTOGRAMS ----------------------------------------------------

#FOR MAIN FIGURE

#wgbs
Npm = sum(!is.na(aldat$l.fracMeth))
pmHist = aldat %>% 
  ggplot(aes(x=l.fracMeth)) +
  geom_histogram() +
  scale_x_continuous(labels = log2_to_percent) +
  labs(title='WGBS', subtitle=paste(Npm, 'genes'), x='% methylation')

#mbdseq
Nmbd = sum(!is.na(aldat$mbd.score))
mbdHist = aldat %>% 
  ggplot(aes(x=mbd.score)) +
  geom_histogram() +
  labs(title='MBD-seq', subtitle=paste(Nmbd, 'genes'),x='MBD-score')

#methylRAD
Nmr = sum(!is.na(aldat$mrB))
mrHist = aldat %>% 
  ggplot(aes(x=mrB)) +
  geom_histogram() +
  labs(title='mdRAD', subtitle=paste(Nmr, 'genes'), x=bquote(log[2]~FPKM))


# BUILD LEVEL CORRELATIONS ------------------------------------------------

#FOR MAIN FIGURE

ALPHA=0.01

pm.mbd.lvl = plot_scatter_pearsonCor_annotated(aldat,
                                               xcol='l.fracMeth',
                                               ycol='mbd.score',
                                               xlab='WGBS',
                                               ylab='MBD-seq',
                                               ALPHA=ALPHA) + scale_x_continuous(labels = log2_to_percent)
pm.mr.lvl = plot_scatter_pearsonCor_annotated(aldat,
                                               xcol='l.fracMeth',
                                               ycol='mrB',
                                               xlab='WGBS',
                                               ylab='mdRAD',
                                               ALPHA=ALPHA) + scale_x_continuous(labels = log2_to_percent)
mbd.mr.lvl = plot_scatter_pearsonCor_annotated(aldat,
                                               xcol='mbd.score',
                                               ycol='mrB',
                                               xlab='MBD-seq',
                                               ylab='mdRAD',
                                               ALPHA=ALPHA)



# plot figure 1 -----------------------------------------------------------
pmAxis = c(0.01, 0.10, 0.9)
pmLimits = c(log(0.001,2), log(1, 2))
mbdAxis = c(-2, 0, 2, 4)
mrAxis = c(-6, 0, 6)
plot_grid(pmHist + scale_x_continuous(breaks=log(pmAxis, 2),
                                      labels=pmAxis*100,
                                      limits = pmLimits) +
            labs(x='% methylation'),
          mbdHist + scale_x_continuous(breaks=mbdAxis),
          mrHist + scale_x_continuous(breaks=mrAxis),
          pm.mbd.lvl + 
            scale_x_continuous(breaks=log(pmAxis, 2), labels=pmAxis*100,limits = pmLimits) +
            scale_y_continuous(breaks=mbdAxis),
          pm.mr.lvl +
            scale_x_continuous(breaks=log(pmAxis, 2), labels=pmAxis*100,limits = pmLimits) +
            scale_y_continuous(breaks=mrAxis),
          mbd.mr.lvl +
            scale_x_continuous(breaks=mbdAxis) +
            scale_y_continuous(breaks=mrAxis),
          nrow=2,
          labels = LETTERS[1:6])



# PLOT CORRELATIONS WITH CPGOE --------------------------------------------

mcols = c('l.fracMeth',
          'mbd.score',
          'mrB')
mnames = c('WGBS',
           'MBD-seq',
           'mdRAD')
ylabs = list('% methylation',
          'MBD-score',
          bquote(log[2]*FPKM))
cpgplots = list()
for (i in 1:length(mcols)){
  ycol=mcols[i]
  ylab = ylabs[[i]]
  name=mnames[i]
  print(paste(name, '..', sep='.'))
  plt=plot_scatter_pearsonCor_annotated(aldat,
                                     xcol='cpgoe',
                                     ycol=ycol,
                                     xlab='CpGo/e',
                                     ylab=ylab,
                                     ALPHA=ALPHA) + 
    theme(axis.title.x = element_blank()) +
    labs(title=name)
  cpgplots[[ycol]]=plt
  
}
cpgplots[['l.fracMeth']] = cpgplots[['l.fracMeth']] + scale_y_continuous(labels=log2_to_percent)
plot_shared_x(cpgplots, bquote(CpG[o/e]))


# MBD-SEQ MEASURES --------------------------------------------------------

#choose the columns, titles, and subtitles
mbd.cols = c('l.fracMeth',
          'mrB',
          'mbd.score',
          'mbd.FPKM',
          'mbd.FPCpGM')
mbd.titles = c('WGBS',
           'mdRAD',
           'MBD-seq',
           'MBD-seq',
           'MBD-seq')
mbd.subtitles = list('% methylation',
                 bquote(log[2]*'FPKM'),
                 'MBD-score',
                 bquote(log[2]*'FPKM'),
                 bquote(log[2]*'FPCpGM'))

#build plots
mbd.pltList = plot_pairs_hist_diagonal(mbd.cols, mbd.titles, mbd.subtitles)
toChange = names(mbd.pltList)[grep('^l.fracMeth', names(mbd.pltList))]
for (n in toChange){
  mbd.pltList[[n]]=mbd.pltList[[n]] + scale_x_continuous(labels = log2_to_percent)
}
plot_grid(plotlist = mbd.pltList, nrow = length(mbd.cols))



# METHYLRAD MEASURES ------------------------------------------------------

mr.cols = c('l.fracMeth',
          'mbd.score',
          'mrF',
          'mrM',
          'mrF.s',
          'mrM.s')
mr.titles = c('WGBS',
           'MBD-seq',
           'mdRAD',
           'mdRAD',
           'mdRAD',
           'mdRAD')
mr.subtitles = list('% methylation',
              'MBD-score',
              'FspEI',
              'MspJI',
              'FspEI/site',
              'MspJI/site')

mr.pltList = plot_pairs_hist_diagonal(mr.cols, mr.titles, mr.subtitles)
toChange = names(mr.pltList)[grep('l.fracMeth', names(mr.pltList))]
for (n in toChange){
  mr.pltList[[n]]=mr.pltList[[n]] + scale_x_continuous(labels = log2_to_percent)
}
plot_grid(plotlist = mr.pltList, nrow = length(mr.cols))


# PICOMETHYL MEASURES -----------------------------------------------------

pm.cols = c('l.fracMeth',
          'pm.glmLvl',
          'l.mCpG_per_CpG',
          'l.mCpG_per_bp',
          'mbd.score',
          'mrB')
pm.titles = c('WGBS',
           'WGBS',
           'WGBS',
           'WGBS',
           'MBD-seq',
           'mdRAD')
pm.subtitles = list('% methylation',
           'glm coef. (logit)',
           bquote(log[2]*'mCpG/CpG'),
           bquote(log[2]*'mCpG/bp'),
           'MBD-score',
           bquote(log[2]*'FPKM'))

pm.pltList = plot_pairs_hist_diagonal(pm.cols, pm.titles, pm.subtitles)
toChange = names(pm.pltList)[grep('l.fracMeth', names(pm.pltList))]
for (n in toChange){
  pm.pltList[[n]]=pm.pltList[[n]] + scale_x_continuous(labels = log2_to_percent)
}
plot_grid(plotlist = pm.pltList, nrow = length(pm.cols))


pltList=list()
for (i in 1:length(mcols)){
  for (j in 1:length(mcols)){
    xcol=mcols[i]
    ycol=mcols[j]
    xname = mnames[i]
    yname = mnames[j]
    plt = plot_scatter_pearsonCor_annotated(aldat,
                                            xcol=xcol,
                                            ycol=ycol,
                                            xlab=xname,
                                            ylab=yname,
                                            ALPHA=ALPHA) + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
    pltName = paste(xcol,ycol,sep='_')
    pltList[[pltName]]=plt
  }
}

#replace diagonal with histograms
histList = list()
for (i in 1:length(mcols)){
  col=mcols[i]
  name=mnames[i]
  subt = subtitles[[i]]
  Ngenes = sum(!is.na(aldat[,col]))
  hist = aldat %>% 
    ggplot(aes_string(x=col)) +
    geom_histogram() +
    labs(title=name, x=name, subtitle=subt) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  histList[[col]]=hist
}

for (col in mcols){
  selfName = paste(col,col,sep='_')
  pltList[[selfName]] = histList[[col]]
}

plot_grid(plotlist = pltList, nrow = length(mcols))

