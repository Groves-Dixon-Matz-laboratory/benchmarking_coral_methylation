#read_reduced_gbm_level.R
rm(list=ls())

library(tidyverse)
library(cowplot)


# SOME FUNCTIONS ----------------------------------------------------------

#funciton to clear away things that mess up cor
prep_for_cor = function(vector){
  bad=c(-Inf, Inf, NaN)
  vector[vector %in% bad]<-NA
  return(vector)
}

#function to get correlations for set of y columns
get_correlations_by_column = function(dat, xcol, ycols){
  corVec = c()
  xVec =  prep_for_cor(dat[,xcol])
  for (yc in ycols){
    yVec = prep_for_cor(dat[,yc])
    xyCor = cor(x=xVec, y=yVec, use='complete.obs')
    corVec = append(corVec, xyCor)
  }
  res = data.frame(y=ycols,
                   x=xcol,
                   yxCor = corVec)
  return(res)
}



# LOAD DATA ---------------------------------------------------------------

#the full dataset level measures
ll=load('comparisons/datasets/gbmLvl.Rdata')
ll
head(gbm.dat)
main.gbm = gbm.dat %>% 
  select(name, l.fracMeth, mrB, mbd.score)


# CHECK CORELATIONS FOR GLM N GENE TEST -----------------------------------

# IT DOESN'T NOT MATTER. see picomethyl/datasets/glmLvlNgeneTest/glmLvl_Ngene_test.R

# METHYLRAD REDUCED PLOTS -------------------------------------------------

#upload the methylRAD reduced read measures
ll=load('methylRAD/datasets/countReducedLvlsFULL.Rdata')
ll
head(rdat)


#merge with main gbm measures
mr.dat = rdat %>% 
  full_join(main.gbm, by='name')

#get correlations for each percentage
pcts = colnames(rdat)[grep('pct', colnames(rdat))]
# pcts = pcts[pcts != 'pct.100']
length(pcts)


#get WGBS correlations
pmCors = get_correlations_by_column(dat=mr.dat,
                           xcol='l.fracMeth',
                           ycols=pcts)

#get MBD correlations
mbdCors = get_correlations_by_column(dat=mr.dat,
                                    xcol='mbd.score',
                                    ycols=pcts)

#get MR correlations
mrCors = get_correlations_by_column(dat=mr.dat,
                                     xcol='mrB',
                                     ycols=pcts)

#assemble
medReadCount = sum(csums)
mrRes = rbind(pmCors,
               mbdCors,
               mrCors) %>% 
  mutate(pct = as.numeric(sub('pct.', '', y)),
         prop=pct/100,
         medReads = prop*medReadCount,
         logProp = log(prop, 2))


#plot
convert_to_counts = function(x){
  round((2^x)*medReadCount, digits=0)
}
mr.lvl = mrRes %>% 
  ggplot(aes(x=medReads, y=yxCor,color=x)) +
  geom_point() +
  geom_line() +
  labs(x='total reads',
       y='GBM level correlation',
       subtitle='mdRAD')



# MBD-SEQ REDUCED PLOTS ---------------------------------------------------

#upload the methylRAD reduced read measures
ll=load('mbdSeq/datasets/countReducedLvls.Rdata')
ll
head(rdat)


#merge with main gbm measures
mbd.dat = rdat %>% 
  full_join(main.gbm, by='name')

#get correlations for each percentage
pcts = colnames(rdat)[grep('pct', colnames(rdat))]
# pcts = pcts[pcts != 'pct.100']
length(pcts)


#get WGBS correlations
pmCors = get_correlations_by_column(dat=mbd.dat,
                                    xcol='l.fracMeth',
                                    ycols=pcts)

#get MBD correlations
mbdCors = get_correlations_by_column(dat=mbd.dat,
                                     xcol='mbd.score',
                                     ycols=pcts)

#get MR correlations
mrCors = get_correlations_by_column(dat=mbd.dat,
                                    xcol='mrB',
                                    ycols=pcts)

#assemble
medReadCount = sum(csums)
mbdRes = rbind(pmCors,
               mbdCors,
               mrCors) %>% 
  mutate(pct = as.numeric(sub('pct.', '', y)),
         prop=pct/100,
         medReads = prop*medReadCount,
         logProp = log(prop, 2))


#plot
mbd.lvl = mbdRes %>% 
  ggplot(aes(x=medReads, y=yxCor,color=x)) +
  geom_point() +
  geom_line() +
  labs(x='total reads',
       y='GBM level correlation',
       subtitle='MBD-seq')


# PICOMETHYL REDUCED PLOTS ------------------------------------------------

#upload the methylRAD reduced read measures
# ll=load('picomethyl/datasets/readReduced_glmGBM_lvl/countReducedLvls.Rdata') #based on GLM
ll=load('picomethyl/datasets/readReduced_basicGBM_lvl/countReducedLvls.Rdata') #based on basic 
ll
head(rdat)


#merge with main gbm measures
pm.dat = rdat %>% 
  select(-chr, -start, -end) %>% 
  full_join(main.gbm, by='name')

#get correlations for each percentage
pcts = colnames(rdat)[grep('pct', colnames(rdat))]
# pcts = pcts[pcts != 'pct.100']
length(pcts)


#get WGBS correlations
pmCors = get_correlations_by_column(dat=pm.dat,
                                    xcol='l.fracMeth',
                                    ycols=pcts)

#get MBD correlations
mbdCors = get_correlations_by_column(dat=pm.dat,
                                     xcol='mbd.score',
                                     ycols=pcts)

#get MR correlations
mrCors = get_correlations_by_column(dat=pm.dat,
                                    xcol='mrB',
                                    ycols=pcts)

#assemble
medReadCount = sum(csums)
pmRes = rbind(pmCors,
               mbdCors,
               mrCors) %>% 
  mutate(pct = as.numeric(sub('pct.', '', y)),
         prop=pct/100,
         medReads = prop*medReadCount,
         logProp = log(prop, 2))


#plot
pm.lvl = pmRes %>% 
  ggplot(aes(x=medReads, y=yxCor,color=x)) +
  geom_point() +
  geom_line() +
  labs(x='total reads',
       y='GBM level correlation',
       subtitle='WGBS')



# FINAL PLOT --------------------------------------------------------------
#set up legend
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
lplt=pm.lvl +
  scale_color_manual(labels = c('WGBS', 'MBD-seq', 'mdRAD'), values=gg_color_hue(3)) +
  labs(color='') +
  theme(legend.position = 'bottom')
lgd = cowplot::get_legend(lplt)

#set up panels
YLIM=c(0.3, 1.1)
pltList = list(pm.lvl,
               mbd.lvl,
               mr.lvl)

pltList = list(pm.lvl + scale_x_continuous(breaks=c(0,1e8,2e8,3e8,4e8), limits=c(0,4.2e8)),
               mbd.lvl + scale_x_continuous(breaks=c(0,1e8,2e8,3e8,4e8), limits=c(0,4.2e8)),
               mr.lvl + scale_x_continuous(breaks=c(0,4e7,8e7,1.2e8), limits=c(0,1.3e8)))

clean_plts = function(x){
  x + theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'none') +
    lims(y=YLIM)
}

pl2 = lapply(pltList, function(x) clean_plts(x))

# #alter axis ticks
# pl2[[3]] = pl2[[3]] + scale_x_continuous(breaks = c(0, 6e+07, 1.2e+08, 1.8e+08))

#plot
xlab = plot_grid(ggdraw() + draw_label('total reads'), lgd, nrow=2)
ylab = ggdraw() + draw_label('GBM correlation', angle=90)
plts = plot_grid(plotlist = pl2, nrow=3)
top = plot_grid(ylab, plts, nrow=1, rel_widths=c(0.06, 1))
full = plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.08))
full



# IDENTIFY INFLECTION POINTS ----------------------------------------------
inflections = c()
tops = c()
#pm
point=5
prop = pmRes %>% 
  arrange(desc(prop)) %>% 
  pull(medReads) %>% 
  unique()
inflect = prop[point]
pmRes %>% 
  ggplot(aes(x=medReads, y=yxCor,color=x)) +
  geom_point() +
  geom_line() +
  labs(x='median reads per sample',
       y='GBM level correlation',
       subtitle='WGBS') +
  geom_vline(xintercept = inflect, lty=2)
inflections = append(inflections, inflect)
top = pmRes %>% 
  filter(x=='l.fracMeth',
         medReads==inflect) %>% 
  pull(yxCor)
tops = append(tops, top)

#mbd
point=4
prop = mbdRes %>% 
  arrange(desc(prop)) %>% 
  pull(medReads) %>% 
  unique()
inflect = prop[point]
mbdRes %>% 
  ggplot(aes(x=medReads, y=yxCor,color=x)) +
  geom_point() +
  geom_line() +
  labs(x='total reads',
       y='GBM level correlation',
       subtitle='MBD-seq') +
  geom_vline(xintercept = inflect, lty=2)
inflections = append(inflections, inflect)
top = mbdRes %>% 
  filter(x=='mbd.score',
         medReads==inflect) %>% 
  pull(yxCor)
tops = append(tops, top)


#mr
point=6
prop = mrRes %>% 
  arrange(desc(prop)) %>% 
  pull(medReads) %>% 
  unique()
inflect = prop[point]
mrRes %>% 
  ggplot(aes(x=medReads, y=yxCor,color=x)) +
  geom_point() +
  geom_line() +
  labs(x='total reads',
       y='GBM level correlation',
       subtitle='MBD-seq') +
  geom_vline(xintercept = inflect, lty=2)
inflections = append(inflections, inflect)
top = mrRes %>% 
  filter(x=='mrB',
         medReads==inflect) %>% 
  pull(yxCor)
tops = append(tops, top)

#look at them
names(inflections) = c('pm', 'mbd', 'mr')
data.frame(inflections/1e6)


#PLOT WITH INFLECTION POINTS
pl3=list()
for (i in 1:length(inflections)){
  lineDf = data.frame(x=inflections[i],
                      xend=inflections[i],
                      y=0,
                      yend=tops[i])
  xcoor = inflections[i]
  ycoor = tops[i]
  pl3[[i]] = pl2[[i]] + geom_segment(aes(x = xcoor, y = YLIM[1], xend = xcoor, yend = ycoor), lwd=0.5, color='grey', lty=2)
}

xlab = plot_grid(ggdraw() + draw_label('total reads'), lgd, nrow=2)
ylab = ggdraw() + draw_label('GBM correlation', angle=90)
plts = plot_grid(plotlist = pl3, nrow=3)
top = plot_grid(ylab, plts, nrow=1, rel_widths=c(0.06, 1))
full = plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.08))
full