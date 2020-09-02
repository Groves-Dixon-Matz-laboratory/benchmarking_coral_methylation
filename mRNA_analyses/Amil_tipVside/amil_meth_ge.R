#amil_meth_ge.R
rm(list=ls())

source('meth_ge_functions.R')

# LOAD GBM DATA -----------------------------------------------------------

ll=load('bioprojects/amillepora_PRJNA601565/level_results/gbmLvl.Rdata')
ll
gbm.dat = gbm.dat %>% 
  dplyr::select(name, l.fracMeth, mrB, mbd.score) %>% 
  as_tibble()

ll=load('bioprojects/amillepora_PRJNA601565/level_results/exonGbmLvl.Rdata')
ll
eg.dat=eg.dat %>% 
  dplyr::rename(l.fracMeth=lfracMeth)


# LOAD EXPRESSION DATA --------------------------------------------------

ll=load('bioprojects/amillepora_PRJNA601565/rna_results/rld.Rdata')
ll

ll=load('bioprojects/amillepora_PRJNA601565/rna_results/previous/genotype_results.Rdata')
ll
g.res=res %>% 
  data.frame() %>% 
  mutate(name=rownames(res)) %>% 
  dplyr::select(name, log2FoldChange, pvalue, padj) %>% 
  as_tibble()

ll=load('bioprojects/amillepora_PRJNA601565/rna_results/previous/tissue_results.Rdata')
ll
t.res=res %>% 
  data.frame() %>% 
  mutate(name=rownames(res)) %>% 
  dplyr::select(name, log2FoldChange, pvalue, padj) %>% 
  as_tibble()



# PLOT HISTOGRAMS ---------------------------------------------------------

#wgbs
Npm = sum(!is.na(eg.dat$l.fracMeth))
PM.BREAKS = log(c(0.002, .015, .12,1), 2)
pmHist = eg.dat %>% 
  ggplot(aes(x=l.fracMeth)) +
  geom_histogram() +
  scale_x_continuous(breaks = PM.BREAKS,
                     labels = log2_to_percent) +
  labs(title='WGBS', subtitle=paste(Npm, 'genes'), x='% methylation')

#mbdseq
Nmbd = sum(!is.na(gbm.dat$mbd.score))
mbdHist = gbm.dat %>% 
  ggplot(aes(x=mbd.score)) +
  geom_histogram() +
  labs(title='MBD-seq', subtitle=paste(Nmbd, 'genes'),x='MBD-score')

#mdRAD
Nmr = sum(!is.na(gbm.dat$mrB))
mrHist = gbm.dat %>% 
  ggplot(aes(x=mrB)) +
  geom_histogram() +
  labs(title='mdRAD', subtitle=paste(Nmr, 'genes'), x=bquote(log[2]~FPKM))

histList = list(pmHist, mbdHist, mrHist)

# COMPARE WITH EXPRESION LEVEL --------------------------------------------

#for genes
rld.df = data.frame(assay(rld))
genes = rownames(rld.df)
ge.lvl = apply(rld.df, 1, function(x) mean(x, na.rm=TRUE)) %>% 
  data.frame() %>% 
  mutate(name=genes) %>% 
  as_tibble() %>% 
  set_names('mnGe', 'name')


gbm.ge = gbm.dat %>% 
  left_join(ge.lvl, by = 'name') 

#WGBS for genes
plot_scatter_pearsonCor_annotated(gbm.ge, 'l.fracMeth', 'mnGe', 'GBM', 'expression level', ALPHA=0.1) +
  geom_smooth(method='lm', se=FALSE)

#mbd-seq genes
mbdLvl = plot_scatter_pearsonCor_annotated(gbm.ge, 'mbd.score', 'mnGe', 'MBD-score GBM', 'expression level', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

#methylRAD for genes
radLvl = plot_scatter_pearsonCor_annotated(gbm.ge, 'mrB', 'mnGe', 'mdRAD GBM', 'expression level', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='mdRAD') +
  theme(plot.title = element_text(hjust=0.5))


#for exons
eg.dat
e.ge = eg.dat %>% 
  left_join(ge.lvl, by = 'name')

wgbsLvl = plot_scatter_pearsonCor_annotated(e.ge, 'l.fracMeth', 'mnGe', 'GBM', 'expression level', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))


#mbd-exons
eMbdLvl=plot_scatter_pearsonCor_annotated(e.ge, 'mbd.score', 'mnGe', 'MBD-seq GBM', 'expression level', ALPHA=0.1) +
  geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) 

#mbd-exons
eMdRADLvl=plot_scatter_pearsonCor_annotated(e.ge, 'mrB', 'mnGe', 'mdRAD GBM', 'expression level', ALPHA=0.1) +
  geom_smooth(method='lm', se=FALSE) 


#PLOT BEST TOGETHER
lvlList = list(wgbsLvl + scale_x_continuous(breaks = PM.BREAKS,
                                            labels = log2_to_percent),
               mbdLvl + lims(x=c(-2.2,3)),
               radLvl)

plot_shared_x_y(pltList,
                xlab='gbM level',
                ylab='mRNA level',
                relXlab = 1/15,
                relYlab = 1/20)


# COMPARE WITH DIFFERENTIAL EXPRESSION ------------------------------------

gTissue = t.res %>% 
  left_join(gbm.dat, by='name') %>% 
  mutate(absl = abs(log2FoldChange))

eTissue = t.res %>% 
  left_join(eg.dat, by='name') %>% 
  mutate(absl = abs(log2FoldChange))
  

#plot for wgbs
wgbs.t = plot_scatter_pearsonCor_annotated(gTissue, 'l.fracMeth', 'absl', 'WGBS', 'tissue difference', ALPHA=0.1, ylim=c(0,2)) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))

#plot for mbd-seq
mbd.t = plot_scatter_pearsonCor_annotated(gTissue, 'mbd.score', 'absl', 'MBD-seq GBM', 'tissue difference', ALPHA=0.1, ylim=c(0,2)) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

#plot for methylRAD
mr.t = plot_scatter_pearsonCor_annotated(gTissue, 'mrB', 'absl', 'MethylRAD GBM', 'tissue difference', ALPHA=0.1, ylim=c(0,2)) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='mdRAD') +
  theme(plot.title = element_text(hjust=0.5))

#plot for wgbs exons
wgbs.t=plot_scatter_pearsonCor_annotated(eTissue, 'l.fracMeth', 'absl', 'WGBS GBM', 'tissue difference', ALPHA=0.1, ylim=c(0,2)) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))

#plot for mbeseq exons
eMbdTip=plot_scatter_pearsonCor_annotated(eTissue, 'mbd.score', 'absl', 'MBD-score GBM', 'tissue difference', ALPHA=0.1, ylim=c(0,2)) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2)

#plot for mdRAD
eMdTip = plot_scatter_pearsonCor_annotated(eTissue, 'mrB', 'absl', 'mdRAD GBM', 'tissue difference', ALPHA=0.1, ylim=c(0,2)) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2)


#plot together 
tipList = list(wgbs.t, mbd.t, mr.t)
plot_grid(plotlist = tipList, nrow=1)



# ASSEMBLE TOGETHER -------------------------------------------------------

PM.BREAKS = log(c(0.002, .015, .12,1), 2)
MBD.BREAKS = c(-2, 0, 2, 4)
MR.BREAKS = c(-6,0,6)

histList2 = lapply(histList, function(x) return(x+theme(axis.title = element_blank(),
                                                        plot.title = element_text(hjust=0.5),
                                                        plot.subtitle = element_text(hjust=0.5))))
lvlList2 = lapply(lvlList, function(x) return(x+theme(axis.title = element_blank(),
                                                      plot.title = element_blank())))
tipList2 = lapply(tipList, function(x) return(x+theme(axis.title = element_blank(),
                                                      plot.title = element_blank())))

hists = plot_grid(plotlist = histList2, nrow=1)
lvls = plot_grid(plotlist = lvlList2, nrow=1)
tips = plot_grid(plotlist = tipList2, nrow=1)

hy =  ggdraw() + draw_label('count', angle=90)
ly =  ggdraw() + draw_label('mRNA level', angle=90)
ty =  ggdraw() + draw_label('mRNA difference', angle=90)
ywidth=1/30
yhists = plot_grid(hy, hists, nrow = 1, rel_widths = c(ywidth,1))
ylvls = plot_grid(ly, lvls, nrow = 1, rel_widths = c(ywidth,1))
ytips = plot_grid(ty, tips, nrow = 1, rel_widths = c(ywidth,1))

plts = plot_grid(yhists,
          ylvls,
          ytips,
          nrow=3)
xlab = ggdraw() + draw_label('gbM level')
plot_grid(plts, xlab, nrow=2, rel_heights = c(1,1/25))











# DIFFERENTIAL METH AND GE FOR -----------------------------------

#FOR GBM
ll=load('comparisons/datasets/tissue_gbm_response.Rdata')
ll  
ll=load('comparisons/datasets/genotype_gbm_response.Rdata')
ll  


merge_meth_and_ge = function(mdat, gedat){
  mdat$name = as.character(mdat$name)
  gedat$name = as.character(gedat$name)
  rdat = mdat %>% 
    left_join(gedat, by='name') %>% 
    as_tibble()
}

t.gbm = merge_meth_and_ge(t.gbm.dat, t.res)


wgbs.tissue = plot_scatter_pearsonCor_annotated(t.gbm, 'meth.diff', 'log2FoldChange', 'WGBS', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))

mbd.tissue = plot_scatter_pearsonCor_annotated(t.gbm, 'mbd.log2FoldChange', 'log2FoldChange', 'MBD-seq', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

mr.tissue = plot_scatter_pearsonCor_annotated(t.gbm, 'mr.log2FoldChange', 'log2FoldChange', 'MethylRAD', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE) +
  labs(title='MethylRAD') +
  theme(plot.title = element_text(hjust=0.5))


plot_grid(wgbs.tissue, mbd.tissue, mr.tissue, nrow=1)



# FOR PROMOTERS -----------------------------------------------------------

ll=load('comparisons/datasets/tissue_promoter_response.Rdata')
ll  
ll=load('comparisons/datasets/genotype_promoter_response.Rdata')
ll  

t.prom = merge_meth_and_ge(t.dat, t.res)


wgbs.tissue = plot_scatter_pearsonCor_annotated(t.prom, 'meth.diff', 'log2FoldChange', 'WGBS', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))

mbd.tissue = plot_scatter_pearsonCor_annotated(t.prom, 'mbd.log2FoldChange', 'log2FoldChange', 'MBD-seq', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

mr.tissue = plot_scatter_pearsonCor_annotated(t.prom, 'mr.log2FoldChange', 'log2FoldChange', 'MethylRAD', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='MethylRAD') +
  theme(plot.title = element_text(hjust=0.5))


plot_grid(wgbs.tissue, mbd.tissue, mr.tissue, nrow=1)




# FOR 500 BP WINDOWDS -----------------------------------------------------

ll=load('comparisons/datasets/tissue_500bp_response.Rdata')
ll  
ll=load('comparisons/datasets/genotype_500bp_response.Rdata')
ll  


#GET CLOSEST GENE


gdat = read_tsv('metadata/geneBounds.tsv',
                col_names = c('chr', 'start', 'end', 'description'))

CUT=0.2
bdat = t.dat %>% 
  dplyr::select(chr:end, meth.diff, pvalue, mr.log2FoldChange, mr.pvalue, mbd.log2FoldChange, mbd.pvalue) %>% 
  as_tibble() %>% 
  filter(pvalue < CUT,
         mr.pvalue < CUT,
         mbd.pvalue < CUT)
bdat

#get unique chrs
gchrs = unique(gdat$chr)
bchrs = unique(bdat$chr)
uchrs = gchrs[gchrs %in% bchrs]

for (chrom in uchrs){
  print(paste(chrom, '..', sep='.'))
  gsub = gdat %>% 
    filter(chr == chrom)
  bsub = bdat %>% 
    filter(chr == chrom)
  rdf = data.frame()
  
  for (i in 1:nrow(bsub)){
    s = as.numeric(bsub[i,'start'])
    e = as.numeric(bsub[i,'end'])
    leftGenes = gsub %>% 
      filter(start <= s) %>% 
      mutate(dist=start - s)
    rightGenes = gsub %>% 
      filter(start >= e) %>% 
      mutate(dist = start + e)
    closeLeft = leftGenes %>% 
      filter(dist == min(dist))
    closeRight = rightGenes %>% 
      filter(dist == min(dist))
    rres =bsub[i,]
    if (nrow(closeLeft)>0){
      rres$closeLeft = closeLeft$description
      rres$leftDist = closeLeft$dist 
    } else{
      rres$closeLeft = 'none'
      rres$leftDist = 'none'
    }
    if (nrow(closeRight)>0){
      rres$closeRight = closeRight$description
      rres$rightDist = closeRight$dist
    } else{
      rres$closeRight = 'none'
      rres$rightDist = 'none'
    }
    rdf = rbind(rdf, rres)
  }
}

split_out_name = function(sVector){
  s1 = sapply(sVector, function(x) strsplit(x, ';')[[1]][1])
  r = sub('ID=', '', s1)
  return(r)
}

rdf$name = split_out_name(rdf$closeRight)  
# rdf$name = split_out_name(rdf$closeLeft)

rdat = rdf %>% 
  left_join(t.res, by= 'name') %>% 
  dplyr::select(chr:mbd.pvalue, name:padj)


absdat = rdat %>% 
  mutate(ge = abs(log2FoldChange),
         pm = abs(meth.diff),
         mbd = abs(mbd.log2FoldChange),
         mr = abs(mr.log2FoldChange)) %>% 
  dplyr::select(ge:mr)

absdat %>% 
  ggplot(aes(x=pm, y=ge)) +
  geom_point(alpha=0.2)

wgbs.tissue = plot_scatter_pearsonCor_annotated(absdat, 'pm', 'ge', 'WGBS', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='WGBS') +
  theme(plot.title = element_text(hjust=0.5))

mbd.tissue = plot_scatter_pearsonCor_annotated(t.prom, 'mbd.log2FoldChange', 'log2FoldChange', 'MBD-seq', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='MBD-seq') +
  theme(plot.title = element_text(hjust=0.5))

mr.tissue = plot_scatter_pearsonCor_annotated(t.prom, 'mr.log2FoldChange', 'log2FoldChange', 'MethylRAD', 'GE', ALPHA=0.1) +
  # geom_smooth(se=FALSE) +
  geom_smooth(method='lm', se=FALSE, color='red', lty=2) +
  labs(title='MethylRAD') +
  theme(plot.title = element_text(hjust=0.5))




