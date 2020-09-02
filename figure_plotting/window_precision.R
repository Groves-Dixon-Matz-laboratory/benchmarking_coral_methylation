#window_precision.R
rm(list=ls())
library(tidyverse)
source('benchmarking_functions.R')

############################################################
#################### SET UP LEVEL PLOTS ####################
############################################################


# LOAD PICOMETHYL DATA ----------------------------------------------------

ll=load('picomethyl/datasets/windowPrecision/methLevelList.Rdata')
ll
pm.lvlList = lvlList
names(pm.lvlList)


# LOAD MBDSEQ DATA --------------------------------------------------------

ll=load('mbdSeq/datasets/windowPrecision/methLevelList.Rdata')
ll
mbd.lvlList = lvlList
names(mbd.lvlList)


# LOAD METHYLRAD DATA -----------------------------------------------------

ll=load('methylrad/datasets/windowPrecision/methLevelList.Rdata')
ll
mr.lvlList = lvlList
names(mr.lvlList)



# MERGE UP ----------------------------------------------------------------

wnames = names(pm.lvlList)

all.lvlList=list()
for (wn in wnames){
  pm=pm.lvlList[[wn]] %>% 
    dplyr::rename(pm.lvl = fracMeth) %>% 
    dplyr::select(chr:name, pm.lvl) %>% 
    mutate(name=as.character(name))
  mbd=mbd.lvlList[[wn]] %>% 
    dplyr::rename(mbd.lvl = log2FoldChange) %>% 
    dplyr::select(chr:name, mbd.lvl) %>% 
    mutate(name=as.character(name)) 
  mr=mr.lvlList[[wn]] %>% 
    dplyr::rename(mr.lvl = mrB) %>% 
    dplyr::select(chr:name, mr.lvl) %>% 
    mutate(name=as.character(name)) 
  mdat = pm %>% 
    full_join(mbd, by='name') %>% 
    full_join(mr, by='name') %>% 
    dplyr::select(name, 
                  pm.lvl, mbd.lvl, mr.lvl)
  all.lvlList[[wn]]=mdat
}


# GET CORRELATIONS --------------------------------------------------------

wn.to.bp = list('100bp'=100,
                '500bp'=500,
                '1kb'=1000,
                '5kb'=5000,
                '10kb'=10000)

#Set up correlations for plotting by window size
l.cdat=data.frame()
for (wn in wnames){
  wdat = all.lvlList[[wn]]
  bp=wn.to.bp[[wn]]
  pm.mbd = get_pearson_cor(wdat, 'pm.lvl', 'mbd.lvl')
  pm.mr = get_pearson_cor(wdat, 'pm.lvl', 'mr.lvl')
  mbd.mr = get_pearson_cor(wdat, 'mbd.lvl', 'mr.lvl')
  res=data.frame(bp,
                 pm.mbd,
                 pm.mr,
                 mbd.mr)
  l.cdat=rbind(l.cdat, res)
}

l.lcdat1 = l.cdat %>% 
  pivot_longer(cols=c('pm.mbd', 'pm.mr', 'mbd.mr'),
               names_to = 'pair',
               values_to = 'cor') %>% 
  separate(pair, into=c('a1','a2')) %>% 
  mutate(a1=factor(a1),
         a2=factor(a2))
l.lcdat2=l.lcdat1
l.lcdat2$a1=l.lcdat1$a2
l.lcdat2$a2=l.lcdat1$a1
l.lcdat = rbind(l.lcdat1,l.lcdat2) %>% 
  mutate(a1=factor(a1, levels=c('pm','mbd','mr')),
         a2=factor(a2, levels=c('pm', 'mbd', 'mr')))

colors=gg_color_hue(3)
pm.col=colors[1]
mbd.col<-colors[2]
mr.col<-colors[3]
color.sets = list('pm'=c(mbd.col,
                         mr.col),
                  'mbd'=c(pm.col,
                          mr.col),
                  'mr'=c(pm.col,
                         mbd.col))
shape.sets = list('pm'=c(17,
                         15),
                  'mbd'=c(19,
                          15),
                  'mr'=c(19,
                         17))
title.sets=list('pm'='WGBS',
                'mbd'='MBD-seq',
                'mr'='mdRAD')


#BUILD PLOTS
assays=unique(l.lcdat$a1)
pltList=list()
for (a in assays){
  colors = color.sets[[a]]
  shapes=shape.sets[[a]]
  title=title.sets[[a]]
  plt=l.lcdat %>%
    filter(a1==a) %>%
    mutate(kb=bp/1e3) %>%
    ggplot(aes(x=kb, y=cor, color=a2, shape=a2)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values=colors) +
    scale_shape_manual(values=shapes) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none') 
  pltList[[a]]=plt
}
pltList[[3]]


###########################################################
################## SET UP RESPONSE PLOTS ##################
###########################################################


# LOAD PICOMETHYL DATA ----------------------------------------------------

ll=load('picomethyl/datasets/windowPrecision/genotype_response.Rdata')
ll
pm.gList = g.diffListPos
names(pm.gList)


# LOAD MBDSEQ DATA --------------------------------------------------------

ll=load('mbdSeq/datasets/windowPrecision/genotype_response_ub.Rdata')
ll
mbd.gList = g.ub.resListPos
names(mbd.gList)


# LOAD METHYLRAD DATA -----------------------------------------------------

ll=load('methylrad/datasets/windowPrecision/genotype_responses.Rdata')
ll
mr.gList = mr.g.list
names(mr.gList)



# MERGE UP ----------------------------------------------------------------

wnames = names(pm.gList)

all.gList=list()
for (wn in wnames){
  pm=pm.gList[[wn]] %>% 
    dplyr::rename(pm.log2FoldChange = meth.diff,
                  pm.pvalue = pvalue,
                  pm.padj = qvalue) %>% 
    mutate(name=as.character(name))
  mbd=mbd.gList[[wn]] %>% 
    dplyr::rename(mbd.log2FoldChange = log2FoldChange,
                  mbd.pvalue = pvalue,
                  mbd.padj = padj) %>% 
    mutate(name=as.character(name)) 
  mr=mr.gList[[wn]] %>% 
    dplyr::rename(mr.log2FoldChange = log2FoldChange,
                  mr.pvalue = pvalue,
                  mr.padj = padj) %>% 
    mutate(name=as.character(name)) 
  mdat = pm %>% 
    mutate(name=as.character(name)) %>% 
    full_join(mbd, by='name') %>% 
    full_join(mr, by='name') %>% 
    dplyr::select(name, 
           pm.log2FoldChange, pm.pvalue, pm.padj,
           mbd.log2FoldChange, mbd.pvalue, mbd.padj,
           mr.log2FoldChange, mr.pvalue, mr.padj)
  all.gList[[wn]]=mdat
}


# GET CORRELATIONS --------------------------------------------------------

w.pltList=list()
for (wn in wnames){
  wdat = all.gList[[wn]]
  pm.mbd = plot_scatter_pearsonCor_annotated(wdat, 'pm.log2FoldChange', 'mbd.log2FoldChange', 'WGBS', 'MBD-seq', ALPHA=0.1) + labs(title=wn)
  pm.mr = plot_scatter_pearsonCor_annotated(wdat, 'pm.log2FoldChange', 'mr.log2FoldChange', 'WGBS', 'MethylRAD', ALPHA=0.1) + labs(title=wn)
  plt=plot_grid(pm.mbd, pm.mr)
  w.pltList[[wn]]=plt
}

#DECIDED THESE WERE A LITTLE EXCESSIVE, BUT LEAVING FOR REFERENCE
# plot_grid(w.pltList[[1]],
#           w.pltList[[2]],
#           w.pltList[[3]],
#           w.pltList[[4]],
#           w.pltList[[5]],
#           nrow=5)
# 
# 
wn.to.bp = list('100bp'=100,
                '500bp'=500,
                '1kb'=1000,
                '5kb'=5000,
                '10kb'=10000)

#Set up correlations for plotting by window size
cdat=data.frame()
for (wn in wnames){
  wdat = all.gList[[wn]]
  bp=wn.to.bp[[wn]]
  pm.mbd = get_pearson_cor(wdat, 'pm.log2FoldChange', 'mbd.log2FoldChange')
  pm.mr = get_pearson_cor(wdat, 'pm.log2FoldChange', 'mr.log2FoldChange')
  mbd.mr = get_pearson_cor(wdat, 'mbd.log2FoldChange', 'mr.log2FoldChange')
  res=data.frame(bp,
                 pm.mbd,
                 pm.mr,
                 mbd.mr)
  cdat=rbind(cdat, res)
}

lcdat1 = cdat %>% 
  pivot_longer(cols=c('pm.mbd', 'pm.mr', 'mbd.mr'),
                     names_to = 'pair',
                     values_to = 'cor') %>% 
  separate(pair, into=c('a1','a2')) %>% 
  mutate(a1=factor(a1),
         a2=factor(a2))
lcdat2=lcdat1
lcdat2$a1=lcdat1$a2
lcdat2$a2=lcdat1$a1
lcdat = rbind(lcdat1,lcdat2) %>% 
  mutate(a1=factor(a1, levels=c('pm','mbd','mr')),
         a2=factor(a2, levels=c('pm', 'mbd', 'mr')))

colors=gg_color_hue(3)
pm.col=colors[1]
mbd.col<-colors[2]
mr.col<-colors[3]
color.sets = list('pm'=c(mbd.col,
                         mr.col),
                  'mbd'=c(pm.col,
                          mr.col),
                  'mr'=c(pm.col,
                         mbd.col))
shape.sets = list('pm'=c(17,
                         15),
                  'mbd'=c(19,
                          15),
                  'mr'=c(19,
                         17))
title.sets=list('pm'='WGBS',
                'mbd'='MBD-seq',
                'mr'='mdRAD')


#BUILD PLOTS
assays=unique(lcdat$a1)
pltList=list()
for (a in assays){
  colors = color.sets[[a]]
  shapes=shape.sets[[a]]
  title=title.sets[[a]]
  plt=lcdat %>%
    filter(a1==a) %>%
    mutate(kb=bp/1e3) %>%
    ggplot(aes(x=kb, y=cor, color=a2, shape=a2)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values=colors) +
    scale_shape_manual(values=shapes) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position='none') +
    scale_y_continuous(breaks=c(0, 0.1,0.2, 0.3, 0.4),
                       labels=c('0', '0.1','0.2', '0.3', ''),
                       limits = c(0,0.4))
  pltList[[a]]=plt
}
pltList[[1]]


# BUILD FINAL PLOT --------------------------------------------------------


#legend
lplt = lcdat %>% 
  mutate(a2=factor(a2, levels=c('pm', 'mbd', 'mr'))) %>% 
  ggplot(aes(x=bp, y=cor, color=a2)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values=c(pm.col,mbd.col,mr.col),
                     labels=c('WGBS', 'MBD', 'MR')) +
  labs(color='') +
  theme(legend.position='bottom',
        legend.justification = 'center')
l=cowplot::get_legend(lplt)

#PLOT PIECES
labelx = c(-0.05, -0.05, -0.05)
#level
l.plts=plot_grid(plotlist=l.pltList, nrow=1, labels=LETTERS[1:3])
l.ylab = ggdraw() + draw_label('level\ncorrelation', angle=90)
l.pltsY = plot_grid(l.ylab, l.plts, nrow=1, rel_widths=c(1,12))
#response
r.plts=plot_grid(plotlist=pltList, nrow=1, labels=LETTERS[4:6])
r.ylab = ggdraw() + draw_label('difference\ncorrelation', angle=90)
r.pltsY = plot_grid(r.ylab, r.plts, nrow=1, rel_widths=c(1,12))
#together
plts=plot_grid(l.pltsY, r.pltsY, nrow=2, rel_heights=c(1,0.8))
# outerY = ggdraw() + draw_label('correlation', angle=90)
# plts.y = plot_grid(outerY, plts, nrow=1, rel_widths=c(1,30))
xlab = plot_grid(ggdraw() + draw_label('window size (Kb)'))
plts.x = plot_grid(plts, xlab, nrow=2, rel_heights=c(15,1))
plts.f = plot_grid(plts.x, l, rel_heights=c(15,1), nrow=2)
plts.f


  

#BUILD PLOTS
assays=unique(l.lcdat$a1)
f.pltList=list()
for (a in assays){
  colors = color.sets[[a]]
  title=title.sets[[a]]
  shapes=shape.sets[[a]]
  lvlDat = l.lcdat %>% 
    filter(a1==a) %>% 
    mutate(kb=bp/1e3)
  difDat = lcdat %>% 
    filter(a1==a) %>% 
    mutate(kb=bp/1e3)
  plt = ggplot() +
    geom_point(data=lvlDat, aes(x=kb, y=cor, color=a2, shape=a2)) +
    geom_line(data=lvlDat, aes(x=kb, y=cor, color=a2)) +
    geom_point(data=difDat, aes(x=kb, y=cor, color=a2, shape=a2)) +
    geom_line(data=difDat, aes(x=kb, y=cor, color=a2), lty=2) +
    labs(title=title) +
    scale_color_manual(values=colors) +
    scale_shape_manual(values=shapes) +
    scale_x_continuous(breaks=seq(0,10, 2)) +
    lims(y=c(0,0.95)) +
    theme(axis.title.x = element_blank(),
          # axis.title.y = element_blank(),
          plot.title = element_text(hjust=0.5),
          legend.position='none')
  f.pltList[[a]]=plt
}

#legend
lcdat$test='difference'
l.lcdat$test = 'level'
lplt=rbind(lcdat, l.lcdat) %>% 
  mutate(a2=factor(a2, levels=c('pm', 'mbd', 'mr')),
         test=factor(test, levels=c('level', 'difference'))) %>% 
  ggplot(aes(x=bp, y=cor, color=a2)) +
  geom_point() +
  geom_line(aes(lty=test)) +
  scale_color_manual(values=c(pm.col,mbd.col,mr.col),
                     labels=c('WGBS', 'MBD', 'mdRAD')) +
  labs(color='', lty='') +
  theme(legend.position='bottom',
        legend.justification = 'center')
l=cowplot::get_legend(lplt)

#assemble
names(f.pltList)
noyList = remove_y(f.pltList)
noyList[[1]] = noyList[[1]] + labs(y='correlation')
plts = plot_grid(plotlist = noyList, nrow=1,  rel_widths = c(1.45,1,1))
xlab = plot_grid(ggdraw() + draw_label('window size (Kb)'))
withx = plot_grid(plts, xlab, nrow=2, rel_heights = c(12,1))
withl = plot_grid(withx, l, rel_heights=c(8,1), nrow=2)
withl



