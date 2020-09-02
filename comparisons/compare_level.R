#compare_level.R

rm(list=ls())
source('benchmarking_functions.R')



# load picomethyl data ----------------------------------------------------
#generate this list of dataframes using picomethyl/process_picoMethyl.R
ll=load('picomethyl/datasets/methLevelList.Rdata')
ll
pmList=lvlList
names(pmList)
lapply(pmList, head)

# load mdRAD data -----------------------------------------------------

ll = load('methylRAD/datasets/methLevelList.Rdata')
ll
mrList = lvlList
remove(lvlList)
lapply(mrList, head)

# load mbdseq from recip meth ---------------------------------------------

ll=load('mbdSeq/datasets/methLevelList.Rdata')
ll
mbdList = lvlList
remove(lvlList)
lapply(mbdList, head)

#remove single outlier from exons
mbdList[['exon']] %>% 
  ggplot(aes(x=log2FoldChange)) +
  geom_histogram()
w = max(abs(mbdList[['exon']]$log2FoldChange), na.rm=TRUE)
mbdList[['exon']] = mbdList[['exon']] %>% 
  filter(abs(log2FoldChange) < w)
mbdList[['exon']] %>% 
  ggplot(aes(x=log2FoldChange)) +
  geom_histogram()

# merge them into single dataframes ---------------------------------------

#temp fix to match col names
change_stop_to_end = function(df){
  colnames(df)[colnames(df)=='stop']<-'end'
  return(df)
}


#build merged df list
mList = list()
names=names(pmList)
for (n in names){
  print(n)
  nList = list(
    pmList[[n]],
    mrList[[n]],
    mbdList[[n]]
  )
  nList2 = lapply(nList, function(x) change_stop_to_end(x))
  ndf = purrr::reduce(nList2, full_join, by=c('chr', 'start', 'end', 'name'))
  colnames(ndf)[colnames(ndf)=='log2FoldChange']<-'mbd.score'
  mList[[n]] = ndf
}


#new list of merged dataframes
lapply(mList, head)


# REMOVE DUPLICATE ROWS ---------------------------------------------------

#these somehow got included at some point
for (x in names(mList)){
  print('-------')
  print(x)
  dat = mList[[x]]
  print(nrow(dat))
  print(length(unique(dat$name)))
  reg = paste(paste(dat$chr, dat$start, sep='_'), dat$end, sep='_')
  print(length(unique(reg)))
  newdat = dat[!duplicated(reg),]
  mList[[x]]=newdat
}

lapply(mList, nrow) #these match the bed files


# SAVE DATASETS -------------------------------------------------

#pull gene bounds from list
gbm.dat = mList[['gene']]
save(gbm.dat, file='comparisons/datasets/gbmLvl.Rdata')



#pull genes based only on exons
edat = mList[['exon']]
eg.dat = edat %>% 
  mutate(name=sub('-RA:cds', '', name)) %>% 
  group_by(name) %>% 
  summarize(fracMeth = mean(fracMeth, na.rm=TRUE),
            mrB = mean(mrB, na.rm=TRUE),
            mbd.score=mean(mbd.score, na.rm=TRUE)) %>% 
  
  mutate(lfracMeth = log2_convesion_noZero(fracMeth)) %>% 
  dplyr::select(-fracMeth)

eg.dat %>% 
  pivot_longer(mrB:lfracMeth, 'stat', 'value') %>% 
  ggplot(aes(x=value)) +
  facet_grid(~stat) +
  geom_histogram()
save(eg.dat, file='comparisons/datasets/exonGbmLvl.Rdata')



# BUILD FIGURE 1 FOR NON GENES --------------------------------------------
lvlCompList = list()
names(mList)


for (region in names(mList)){
  print(paste(region, '...',sep=''))

  aldat = mList[[region]]
  
  #BUILD HISTOGRAMS
  #wgbs
  Npm = sum(!is.na(aldat$l.fracMeth))
  pmHist = aldat %>% 
    ggplot(aes(x=l.fracMeth)) +
    geom_histogram() +
    scale_x_continuous(labels = log2_to_percent) +
    labs(title='WGBS', subtitle=paste(Npm, 'windows'), x='% methylation')
  
  #mbdseq
  Nmbd = sum(!is.na(aldat$mbd.score))
  mbdHist = aldat %>% 
    ggplot(aes(x=mbd.score)) +
    geom_histogram() +
    labs(title='MBD-seq', subtitle=paste(Nmbd, 'windows'),x='MBD-score')
  #mdRAD
  Nmr = sum(!is.na(aldat$mrB))
  mrHist = aldat %>% 
    ggplot(aes(x=mrB)) +
    geom_histogram() +
    labs(title='mdRAD', subtitle=paste(Nmr, 'windows'), x=bquote(log[2]~FPKM))
  
  #BUILD SCATTERPLOTS
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
  
  #ASSEMBLE MULTIPANEL
  rplt = plot_grid(pmHist,
            mbdHist,
            mrHist,
            pm.mbd.lvl,
            pm.mr.lvl,
            mbd.mr.lvl,
            nrow=2,
            labels = LETTERS[1:6])
  lvlCompList[[region]]=rplt
}

#PLOT THEM

#or manually
lvlCompList[['exon']]
lvlCompList[['promoter']]
lvlCompList[['1kb']]


# #or in a loop
# for (reg in c('exon', 'promoter', '1kb')){
#   fn = paste('/Users/grovesdixon/lab_files/projects/benchmark_meth/writing/bm_figure_bits/supplemental/', reg, '_fig1.png', sep='')
#   png(filename=fn,
#       width=1548,
#       height=980)
#   plot(lvlCompList[[reg]])
#   dev.off()
# }


# compare level accross measures ---------------------------------


df=gene
xcol='fracMeth'
ycol='mbd.score'
ALPHA=0.5 

plot_scatter = function(df, xcol, ycol){
  bad = c(NA, NaN, Inf, -Inf)
  badx = df[,xcol] %in% bad
  bady = df[,ycol] %in% bad
  badr = badx | bady
  subdf = df[,c(ycol, xcol)] %>% 
    filter(!badr)
  lm1 = lm(subdf[,ycol] ~ subdf[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  df %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha=ALPHA) +
    labs(subtitle=bquote(R^2*' = '*.(r2)))
}
plot_scatter(gene, xcol, ycol)

my_pairs = function(df){
  print(colnames(df))
  pltList = list()
  for (xcol in colnames(df)){
    for (ycol in colnames(df)){
      pair=paste(xcol, ycol, sep='-')
      print(pair)
      if (xcol==ycol){
        plt= ggdraw() + draw_label(xcol)
      } else {
        plt = plot_scatter(df, xcol, ycol)
      }
      pltList[[pair]]=plt
    }
  }
  return(pltList)
}

sub_and_plot = function(df){
  df$lFracMeth=logit(df$fracMeth)
  pltList = my_pairs(df[,c('lFracMeth', 'mrB', 'mbd.score')])
  return(pltList)
    
}


npltList = list()
for (n in names){
  print('------')
  print(n)
  npltList[[n]]=sub_and_plot(mList[[n]])
}

##### PLOT SCATTERPLOTS

# #genes
# plot_grid(plotlist=npltList[['gene']], nrow=3)
# 
# #exons
# plot_grid(plotlist=npltList[['exon']], nrow=3)
# 
# #promoters
# plot_grid(plotlist=npltList[['promoter']], nrow=3)
# 
# #tss
# plot_grid(plotlist=npltList[['tss']], nrow=3)


# CHECK METHYL RAD LEVELS -------------------------------------------------

mrCols = c('mrF', 'mrM', 'mrB', 'mrF.s', 'mrM.s', 'mrB.s')

ALPHA=0.01
pltList=list()
for (mcol in mrCols){
  print(mcol)
  plt=pm.mbd.lvl = plot_scatter_pearsonCor_annotated(gbm.dat,
                                                 xcol='pm.glmLvl',
                                                 ycol=mcol,
                                                 xlab='WGBS',
                                                 ylab=mcol,
                                                 ALPHA=ALPHA) + scale_x_continuous(labels = inv_logit_labs)
  pltList[[mcol]]=plt
}

plot_grid(plotlist=pltList)



# COMPARE WITH EXPRESSION LEVEL -------------------------------------------
ll=load('rnaseq/datasets/rld.Rdata')
ll
head(rld.df)
ge = data.frame(gene=rownames(rld.df),
                gelvl=apply(rld.df, 1, function(x) mean(x, na.rm=TRUE)),
                stringsAsFactors=FALSE)







# COMPARE WITH OTHER METH LVL MEASURES ------------------------------------


#load new gbm levl estimates
ll=load('picomethyl/datasets/new_genes_glmLvls.Rdata')
ll
colnames(fres)[colnames(fres)=='name']<-'gene'
head(fres)
zed = fres %>% 
  filter(mnLvl > -10) %>% 
  pull(mnLvl) %>% 
  min()
fres = fres %>% 
  mutate(modLvl = if_else(mnLvl < -10,
                          zed,
                          mnLvl)) %>% 
  filter(modLvl < 10)
fres %>% 
  ggplot(aes(x=modLvl)) +
  geom_density()


#load previous
ll=load('picomethyl/old/gbm_results/mean_gbm_level.Rdata')
head(pm)

pm %>% 
  ggplot(aes(x=mnPM)) +
  geom_density()

glms = fres %>% 
  dplyr::select(gene, modLvl) %>% 
  left_join(pm[,c('gene', 'mnPM')], by = 'gene')

glms %>% 
  ggplot(aes(x=mnPM, y=modLvl)) +
  geom_point(alpha=0.1)


#load alternative measures
bdat = read.table('picomethyl/old/gbm_results/basicLvl/all_basic_gbmLvl_by_gene.tsv', header = TRUE) %>% 
  mutate(gene=sub('-RA', '', parentGene),
         mCpG_per_bp = ifelse(mCpG_per_bp > 0.1,
                               NA,
                               mCpG_per_bp)) %>% 
  dplyr::select(gene,
         mCpG_per_CpG,
         mCpG_per_bp,
         mnMeth,
         medMeth) %>% 
  as_tibble()

omeths = glms %>% 
  left_join(bdat, by = 'gene')


gdat = mList[['gene']]
head(gdat)


#merge all together
gdat2 = gdat %>% 
  mutate(lFracMeth=logit(fracMeth)) %>% 
  dplyr::select(name, lFracMeth, mrB, mrF, mrM, mbd.score, pm.glmLvl) %>% 
  dplyr::rename(gene = name) %>% 
  left_join(ge, by = 'gene') %>% 
  left_join(omeths, by = 'gene') %>% 
  as_tibble()
head(gdat2)


lm1=lm(gdat2$mnPM~gdat2$pm.glmLvl)
summary(lm1)

ll=load('picomethyl/datasets/quickGlm_fromCov.Rdata')
ll
head(fres)
ndat = fres %>% 
  select(name, glmLvl) %>% 
  dplyr::rename(gene=name,
                nglmLvl=glmLvl)
gdat3=gdat2 %>% 
  left_join(ndat, by = 'gene')
lm1=lm(gdat3$mnPM~gdat3$nglmLvl)
head(gdat3)
x=gdat3$mnPM - gdat3$nglmLvl
summary(lm1)
gdat3 %>% 
  ggplot(aes(x=mnPM, y=pm.glmLvl)) +
  geom_point(alpha=0.1)


#see how new methylation lvl estimate not as good
bad=c(NA, Inf, -Inf)
corMeth = gdat2 %>% 
  mutate(mrB = ifelse(mrB %in% bad,
                      NA,
                      mrB)) %>% 
  dplyr::select(mrB, mbd.score, mnPM, modLvl) %>% 
  dplyr::rename(oldWGBS = mnPM,
                newWGBS = modLvl)

ALPHA=0.3
pltList = my_pairs(data.frame(corMeth))
plot_grid(plotlist = pltList, nrow=4)


#get Z-scores for summarry meth estimate
zdat = gdat2 %>% 
  mutate(mrB = ifelse(mrB %in% bad,
                       NA,
                       mrB)) %>% 
  dplyr::select(mrB, mbd.score, mnPM) %>% 
  scale() %>% 
  data.frame()

# pltList = my_pairs(zdat)
# plot_grid(plotlist=pltList, nrow=3)

sumry = apply(zdat, 1, mean)
sdat = data.frame(gene=gdat2$gene,
                   smryMeth = sumry)


lgdat = gdat2 %>% 
  left_join(sdat, by = 'gene') %>% 
  pivot_longer(-c('gene', 'gelvl'),
               names_to = 'assay',
               values_to = 'methLvl')

#look at densities
lgdat %>% 
  ggplot() +
  geom_density(aes(x=methLvl)) +
  facet_wrap(~assay, scales='free')


#relationship to expression level
lgdat %>% 
  ggplot(aes(x=methLvl,y=gelvl)) +
  geom_point(alpha=0.01) +
  geom_smooth(method='lm') +
  facet_wrap(~assay, scales='free')

plot_scatter(data.frame(lgdat), 'methLvl', 'gelvl') +
  facet_wrap(~assay, scales='free')

assays = unique(lgdat$assay)
rs = c()
bad=c(NA, Inf, -Inf)
scatList = list()
ALPHA=0.01
for (a in assays){
  adat = lgdat %>% 
    filter(assay==a,
           !gelvl %in% bad,
           !methLvl %in% bad) 
  lma=lm(gelvl~methLvl, data=adat)
  r2=summary(lma)$r.squared
  rs = append(rs, r2)
  plt=plot_scatter(data.frame(adat), 'methLvl', 'gelvl') + labs(title=a)
  scatList[[a]]=plt
}
res=data.frame(assays, rs)
res
plot_grid(plotlist=scatList, nrow=4)

lgdat %>% 
  filter(assay=='mCpG_per_bp') %>% 
  ggplot(aes(x=methLvl, y=gelvl)) +
  geom_point(alpha=0.01) +
  lims(x=c(0,0.08))

# COMPARE WITH EXPRESSION DIFFERENCES -----------------------------------------

ll=load('rnaseq/datasets/tissue_results.Rdata')
rna.t = format_deseq_res(res)
ll=load('rnaseq/datasets/genotype_results.Rdata')
rna.g = format_deseq_res(res)
ll=load('rnaseq/datasets/rld.Rdata')



#COEFFICIENT OF VARIATION

#calculate coefficient of variance for the genes
cv=data.frame(cv=apply(rld.df, 1, function(x) sd(x) / mean(x))) %>% 
  mutate(gene=rownames(rld.df)) %>% 
  as_tibble() %>% 
  full_join(lgdat, by = 'gene')
head(cv)

rs = c()
bad=c(NA, Inf, -Inf)
cvList = list()
ALPHA=0.01
for (a in assays){
  adat = cv %>% 
    filter(assay==a,
           !gelvl %in% bad,
           !methLvl %in% bad) 
  lma=lm(gelvl~methLvl, data=adat)
  r2=summary(lma)$r.squared
  rs = append(rs, r2)
  plt=plot_scatter(data.frame(adat), 'methLvl', 'cv') + labs(title=a)
  cvList[[a]]=plt
}
res=data.frame(assays, rs)
res
plot_grid(plotlist=cvList, nrow=4)


  



pm.cv.ge = gbm_scatter(pm, cv, 'mnPM', 'cv', 'PM vs GE tissue', 'picomethyl', 'coef. var.')
mbd.cv.ge = gbm_scatter(mbd, cv, 'mbd.score', 'cv', 'MBD vs GE tissue', 'MBD-score', 'coef. var.')
mr.cv.ge = gbm_scatter(bm, cv, 'mnMR', 'cv', 'MR vs GE tissue', 'methylRAD', 'coef. var.')

cvList = list(pm.cv.ge,
             mbd.cv.ge,
             mr.cv.ge)
plot_grid(plotlist=cvList, nrow=1)


#BETWEEN TISSUES

pm.t.ge = gbm_scatter(pm, rna.t, 'mnPM', 'absDiff', 'PM vs GE tissue', 'picomethyl', bquote(log[2]~GE))
mbd.t.ge = gbm_scatter(mbd, rna.t, 'mbd.score', 'absDiff', 'MBD vs GE tissue', 'MBD-score', bquote(log[2]~GE))
mr.t.ge = gbm_scatter(bm, rna.t, 'mnMR', 'absDiff', 'MR vs GE tissue', 'methylRAD', bquote(log[2]~GE))

tList = list(pm.t.ge,
              mbd.t.ge,
              mr.t.ge)
plot_grid(plotlist=tList, nrow=1)


#BETWEEN COLONIES

pm.g.ge = gbm_scatter(pm, rna.g, 'mnPM', 'absDiff', 'PM vs GE genotype', 'picomethyl', bquote(log[2]~GE))
mbd.g.ge = gbm_scatter(mbd, rna.g, 'mbd.score', 'absDiff', 'MBD vs GE genotype', 'MBD-score', bquote(log[2]~GE))
mr.g.ge = gbm_scatter(bm, rna.g, 'mnMR', 'absDiff', 'MR vs GE genotype', 'methylRAD', bquote(log[2]~GE))

gList = list(pm.g.ge,
             mbd.g.ge,
             mr.g.ge)
plot_grid(plotlist=gList, nrow=1)


#DEVELOPMENTAL STAGE
#(from GE meta)
ll=load('~/gitreps/Acropora_gene_expression_meta/deseqResults/adultVlarva_deseqResults.Rdata')
ll
rna.dev = format_deseq_res(res)

pm.dev.ge = gbm_scatter(pm, rna.dev, 'mnPM', 'absDiff', 'PM vs GE age', 'picomethyl', bquote(log[2]~GE))
mbd.dev.ge = gbm_scatter(mbd, rna.dev, 'mbd.score', 'absDiff', 'MBD vs GE age', 'MBD-score', bquote(log[2]~GE))
mr.dev.ge = gbm_scatter(bm, rna.dev, 'mnMR', 'absDiff', 'MR vs GE age', 'methylRAD', bquote(log[2]~GE))

gList = list(pm.dev.ge,
             mbd.dev.ge,
             mr.dev.ge)
plot_grid(plotlist=gList, nrow=1)


#GRAND COEFFICIENT OF VAR

ll=load('~/gitreps/Acropora_gene_expression_meta/largeIgnored/fullDataset_vsd.Rdata')
ll
rld.df=t(datExpr)
rld.df[1:10,1:10]


#calculate coefficient of variance for the genes
cv=data.frame(cv=apply(rld.df, 1, function(x) sd(x) / mean(x))) %>% 
  mutate(gene=rownames(rld.df)) %>% 
  as_tibble()

pm.cv.ge = gbm_scatter(pm, cv, 'mnPM', 'cv', 'PM vs GE tissue', 'picomethyl', 'coef. var.')
mbd.cv.ge = gbm_scatter(mbd, cv, 'mbd.score', 'cv', 'MBD vs GE tissue', 'MBD-score', 'coef. var.')
mr.cv.ge = gbm_scatter(bm, cv, 'mnMR', 'cv', 'MR vs GE tissue', 'methylRAD', 'coef. var.')

cvList = list(pm.cv.ge,
              mbd.cv.ge,
              mr.cv.ge)
plot_grid(plotlist=cvList, nrow=1)



  