#plot_gbm_response_figs.R
library(limma)
library(ggforce)
rm(list=ls())
source('benchmarking_functions.R')


# SELECT COMPARISON TO PLOT --------------------------------------------------

#TISSUE COMPARISON
ll=load('comparisons/datasets/tissue_gbm_response.Rdata')
ll
PM.YLIM=c(0,6)
rdat = t.gbm.dat
head(rdat)

#OR

#GENOTYPE COMPARISON
ll=load('comparisons/datasets/genotype_gbm_response.Rdata')
ll
PM.YLIM=c(0,400)
rdat = g.gbm.dat


# SELECT mdRAD SET TO PLOT ------------------------------------------------

#8-SAMPLE DATASET
mdRAD_lfc = 'esmr.log2FoldChange'
mdRAD_pval = 'esmr.pvalue'
mdRAD_padj = 'esmr.padj'
mdRAD_sample = '8 samples'
ll=load('methylRAD/datasets/csums.Rdata')
millReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
mdTitle=paste('8 samples', '8 libraries', millReads, sep='\n')

#OR

#FULL DATASET (24 samples)
mdRAD_lfc = 'mr.log2FoldChange'
mdRAD_pval = 'mr.pvalue'
mdRAD_padj = 'mr.padj'
mdRAD_sample = '24 samples'
ll=load('methylRAD/datasets/csumsFULL.Rdata')
millReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
mdTitle=paste('12 samples', '24 libraries', millReads, sep='\n')



# BUILD VOLCANO PLOTS ------------------------------------------------------

#SET UP TITLES
#pm
ll=load('picomethyl/datasets/csums.Rdata')
pmMillReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
pmTitle=paste('8 samples', '8 libraries', pmMillReads, sep='\n')
#mbd
ll=load('mbdSeq/datasets/csums.Rdata')
millReads=paste(round(sum(csums)/1e6, digits = 0),' million total reads', sep='')
mbdTitle=paste('8 samples', '16 libraries', millReads, sep='\n')


#handle any zero pvalues
rdat$pvalue = convesion_noZero(rdat$pvalue)
rdat$mbd.pvalue = convesion_noZero(rdat$mbd.pvalue)
rdat$mr.pvalue = convesion_noZero(rdat$mr.pvalue)

HJUST=0.5
#WGBS
v.pm = plot_volcano_general(rdat,
                            xcol='meth.diff',
                            ycol='pvalue',
                            sigcol='qvalue',
                            sigcut=0.1,
                            title='WGBS',
                            subtitle=pmTitle,
                            xlab='log odds') +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))

#MBD-seq
v.mbd = plot_volcano_general(rdat,
                             xcol='mbd.log2FoldChange',
                             ycol='mbd.pvalue',
                             sigcol='mbd.padj',
                             sigcut=0.1,
                             title='MBD-seq',
                             subtitle=mbdTitle,
                             xlab=bquote(log[2]*'fold difference'))+
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))

#mdRAD
v.mr = plot_volcano_general(rdat,
                            xcol=mdRAD_lfc,
                            ycol=mdRAD_pval,
                            sigcol=mdRAD_padj,
                            sigcut=0.1,
                            title=mdRAD_label,
                            subtitle = mdTitle,
                            xlab=bquote(log[2]*'fold difference')) +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))



# BUILD RESPONSE CORRELATIONS ------------------------------------------------

#FOR MAIN FIGURE

ALPHA=0.4

pm.mbd.resp = plot_scatter_pearsonCor_annotated(rdat,
                                                xcol='meth.diff',
                                                ycol='mbd.log2FoldChange',
                                                xlab='WGBS',
                                                ylab='MBD-seq',
                                                ALPHA=ALPHA)
pm.mr.resp = plot_scatter_pearsonCor_annotated(rdat,
                                               xcol='meth.diff',
                                               ycol=mdRAD_lfc,
                                               xlab='WGBS',
                                               ylab=mdRAD_label,
                                               ALPHA=ALPHA)
mbd.mr.resp = plot_scatter_pearsonCor_annotated(rdat,
                                                xcol='mbd.log2FoldChange',
                                                ycol=mdRAD_lfc,
                                                xlab='MBD-seq',
                                                ylab=mdRAD_label,
                                                ALPHA=ALPHA)



# PLOT FIGURE 2 -----------------------------------------------------------
top = plot_grid(v.pm,
                v.mbd,
                v.mr,
                nrow=1,
                labels = LETTERS[1:3])
bottom = plot_grid(pm.mbd.resp,
                   pm.mr.resp,
                   mbd.mr.resp,
                   nrow=1,
                   labels = LETTERS[4:6])
plot_grid(top,bottom, nrow=2, rel_heights = c(1,0.95))



plot_grid(v.pm,
          v.mbd,
          v.mr,
          pm.mbd.resp,
          pm.mr.resp,
          mbd.mr.resp,
          nrow=2,
          labels = LETTERS[1:6])



# GET SOME STATS ----------------------------------------------------------

srdat = rdat %>% 
  mutate(sig.mbd=mbd.padj < 0.1,
         sig.mr = mr.padj < 0.1,
         sig.pm = qvalue < 0.1) %>% 
  dplyr::select(sig.mbd, sig.mr, sig.pm)
nsig = apply(srdat, 1, function(x) sum(x, na.rm=TRUE))
table(nsig)




# VENN DIAGRAM ------------------------------------------------------------

#isolate sig genes
PCUT=0.1
sigdat=rdat %>% 
  mutate(pm.sig=as.numeric(qvalue<PCUT & !is.na(qvalue)),
         mbd.sig=as.numeric(mbd.padj<PCUT & !is.na(mbd.padj)),
         mr.sig=as.numeric(mr.padj<PCUT & !is.na(mr.padj))) %>% 
  as_tibble() %>% 
  dplyr::select(mr.sig,
                mbd.sig,
                pm.sig)


FILL=c(pm.color,mbd.color,mr.color)
LABELS=c('WGBS', 'MBD', 'mdRAD')

#df for drawing circles
df.venn <- data.frame(x = c(0, 0.866, -0.866),
                      y = c(1, -0.5, -0.5),
                      labels = factor(LABELS, levels=LABELS))

#use limma to get venn counts
vdc <- vennCounts(sigdat)
class(vdc) <- 'matrix'
df.vdc <- as.data.frame(vdc)[-1,] %>%
  mutate(x = c(0, 1.2, 0.8, -1.2, -0.8, 0, 0),
         y = c(1.2, -0.6, 0.5, -0.6, 0.5, -1, 0))

#get solo calls
cdf = df.vdc[,1:3]
solo = df.vdc[apply(cdf,1,sum)==1,]


#build venn
ggplot(df.venn) +
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels), alpha = .3, size = 0.5, colour = 'black') +
  coord_fixed() +
  theme_void() +
  theme(legend.position = 'bottom') +
  scale_fill_manual(values = FILL) +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)

#get some stats
in.pm = as.logical(df.vdc$pm.sig)
in.mbd = as.logical(df.vdc$mbd.sig)
in.mr = as.logical(df.vdc$mr.sig)
nsig = df.vdc$Counts
tpm = sum(nsig[in.pm])
tmbd = sum(nsig[in.mbd])
tmr = sum(nsig[in.mr])


#pm corroboration
sum(nsig[in.pm & in.mbd])/tmbd
sum(nsig[in.pm & in.mr])/tmr
pm.solo = round(sum(nsig[in.pm & !in.mr & !in.mbd]) / tpm, digits=2)*100

#mbd corroboration
sum(nsig[in.pm & in.mbd])/tpm
sum(nsig[in.mr & in.mbd])/tmr
mbd.solo = round(sum(nsig[in.mbd & !in.mr & !in.pm]) / tmbd, digits=2)*100

#mr corroboration
sum(nsig[in.pm & in.mr])/tpm
sum(nsig[in.mr & in.mbd])/tmbd
mr.solo = round(sum(nsig[in.mr & !in.mbd & !in.pm]) / tmr, digits=2)*100

#STATS FOR OVERLAP
sigdat
fisher.test(sigdat$pm.sig, factor(sigdat$mbd.sig))
fisher.test(sigdat$pm.sig, factor(sigdat$mr.sig))
fisher.test(sigdat$mbd.sig, factor(sigdat$mr.sig))




# CHECK FOR CORRELATION WITH SMALL EFFECTS --------------------------------
#point here is to show that even for low level differeneces, (like those seen in tissue),
#we see much stronger correlations for genotype.

head(rdat)
lrdat = rdat %>% 
  filter(abs(meth.diff) < 15)
v.pm = plot_volcano_general(lrdat,
                            xcol='meth.diff',
                            ycol='pvalue',
                            sigcol='qvalue',
                            sigcut=0.1,
                            title='WGBS',
                            subtitle=pmTitle,
                            xlab='log odds') +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))
v.mbd = plot_volcano_general(lrdat,
                             xcol='mbd.log2FoldChange',
                             ycol='mbd.pvalue',
                             sigcol='mbd.padj',
                             sigcut=0.1,
                             title='MBD-seq',
                             subtitle=mbdTitle,
                             xlab=bquote(log[2]*'fold difference'))+
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))
v.mr = plot_volcano_general(lrdat,
                            xcol=mdRAD_lfc,
                            ycol=mdRAD_pval,
                            sigcol=mdRAD_padj,
                            sigcut=0.1,
                            title=mdRAD_label,
                            subtitle = mdTitle,
                            xlab=bquote(log[2]*'fold difference')) +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))
pm.mbd.resp = plot_scatter_pearsonCor_annotated(lrdat,
                                                xcol='meth.diff',
                                                ycol='mbd.log2FoldChange',
                                                xlab='WGBS',
                                                ylab='MBD-seq',
                                                ALPHA=ALPHA)
pm.mr.resp = plot_scatter_pearsonCor_annotated(lrdat,
                                               xcol='meth.diff',
                                               ycol=mdRAD_lfc,
                                               xlab='WGBS',
                                               ylab=mdRAD_label,
                                               ALPHA=ALPHA)
mbd.mr.resp = plot_scatter_pearsonCor_annotated(lrdat,
                                                xcol='mbd.log2FoldChange',
                                                ycol=mdRAD_lfc,
                                                xlab='MBD-seq',
                                                ylab=mdRAD_label,
                                                ALPHA=ALPHA)
plot_grid(v.pm,
          v.mbd,
          v.mr,
          pm.mbd.resp,
          pm.mr.resp,
          mbd.mr.resp,
          nrow=2,
          labels = LETTERS[1:6])





# -------------------------------------------------------------------------
#------------------- REPEAT WITH ONLY METHYLATED GENES -------------------#
# -------------------------------------------------------------------------

#USE PICOMETHYL TO IDENTIFY METHYLATED GENES
ll=load('comparisons/datasets/gbmLvl.Rdata')
gbm.dat = gbm.dat %>% 
  dplyr::select(name, l.fracMeth)
gbm.dat %>% 
  ggplot(aes(x=l.fracMeth)) +
  geom_histogram() +
  geom_vline(xintercept = -5)
mgenes = gbm.dat %>% 
  filter(l.fracMeth>-5) %>% 
  pull(name)

#SUBSET RDAT 
mrdat = rdat %>% 
  filter(name %in% mgenes)

HJUST=0.5
#WGBS
v.pm = plot_volcano_general(mrdat,
                            xcol='meth.diff',
                            ycol='pvalue',
                            sigcol='qvalue',
                            sigcut=0.1,
                            title='WGBS',
                            xlab='log odds')+
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))

#MBD-seq
v.mbd = plot_volcano_general(mrdat,
                             xcol='mbd.log2FoldChange',
                             ycol='mbd.pvalue',
                             sigcol='mbd.padj',
                             sigcut=0.1,
                             title='MBD-seq',
                             xlab=bquote(log[2]*'fold difference'))+
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))

#MethylRAD
v.mr = plot_volcano_general(mrdat,
                            xcol='mr.log2FoldChange',
                            ycol='mr.pvalue',
                            sigcol='mr.padj',
                            sigcut=0.1,
                            title='MethylRAD',
                            xlab=bquote(log[2]*'fold difference')) +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))



# BUILD RESPONSE CORRELATIONS ------------------------------------------------

#FOR MAIN FIGURE

ALPHA=0.4

pm.mbd.resp = plot_scatter_pearsonCor_annotated(mrdat,
                                                xcol='meth.diff',
                                                ycol='mbd.log2FoldChange',
                                                xlab='WGBS',
                                                ylab='MBD-seq',
                                                ALPHA=ALPHA)
pm.mr.resp = plot_scatter_pearsonCor_annotated(mrdat,
                                               xcol='meth.diff',
                                               ycol='mr.log2FoldChange',
                                               xlab='WGBS',
                                               ylab='MethylRAD',
                                               ALPHA=ALPHA)
mbd.mr.resp = plot_scatter_pearsonCor_annotated(mrdat,
                                                xcol='mbd.log2FoldChange',
                                                ycol='mr.log2FoldChange',
                                                xlab='MBD-seq',
                                                ylab='MethylRAD',
                                                ALPHA=ALPHA)



# PLOT FIGURE 1 -----------------------------------------------------------

plot_grid(v.pm + lims(y=PM.YLIM),
          v.mbd,
          v.mr,
          pm.mbd.resp,
          pm.mr.resp,
          mbd.mr.resp,
          nrow=2,
          labels = LETTERS[1:6])



# GET SOME STATS ----------------------------------------------------------

smrdat = mrdat %>% 
  mutate(sig.mbd=mbd.padj < 0.1,
         sig.mr = mr.padj < 0.1,
         sig.pm = qvalue < 0.1) %>% 
  dplyr::select(sig.mbd, sig.mr, sig.pm)
nsig = apply(smrdat, 1, function(x) sum(x, na.rm=TRUE))
table(nsig)






############################################################################
######################## REPEAT FOR ALL WINDOW TYPES########################
############################################################################

ll=load('comparisons/datasets/all_genotype_response.Rdata')
ll

region='exon'
region='promoter'
region='1kb'

rdat = gList[[region]]


# BUILD VOLCANO PLOTS ------------------------------------------------------

#handle any zero pvalues
rdat$pvalue = convesion_noZero(rdat$pvalue)
rdat$mbd.pvalue = convesion_noZero(rdat$mbd.pvalue)
rdat$mr.pvalue = convesion_noZero(rdat$mr.pvalue)

#WGBS
v.pm = plot_volcano_general(rdat,
                            xcol='meth.diff',
                            ycol='pvalue',
                            sigcol='qvalue',
                            sigcut=0.1,
                            title='WGBS',
                            subtitle=pmTitle,
                            xlab='log odds') +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))
#MBD-seq
v.mbd = plot_volcano_general(rdat,
                             xcol='mbd.log2FoldChange',
                             ycol='mbd.pvalue',
                             sigcol='mbd.padj',
                             sigcut=0.1,
                             title='MBD-seq',
                             subtitle=mbdTitle,
                             xlab=bquote(log[2]*'fold difference'))+
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))

#MethylRAD
v.mr = plot_volcano_general(rdat,
                            xcol='mr.log2FoldChange',
                            ycol='mr.pvalue',
                            sigcol='mr.padj',
                            sigcut=0.1,
                            title=mdRAD_label,
                            subtitle = mdTitle,
                            xlab=bquote(log[2]*'fold difference')) +
  theme(plot.title = element_text(hjust = HJUST),
        plot.subtitle = element_text(hjust=HJUST))


# BUILD RESPONSE CORRELATIONS ------------------------------------------------

#FOR MAIN FIGURE

ALPHA=0.05

pm.mbd.resp = plot_scatter_pearsonCor_annotated(rdat,
                                                xcol='meth.diff',
                                                ycol='mbd.log2FoldChange',
                                                xlab='WGBS',
                                                ylab='MBD-seq',
                                                ALPHA=ALPHA)
pm.mr.resp = plot_scatter_pearsonCor_annotated(rdat,
                                               xcol='meth.diff',
                                               ycol='mr.log2FoldChange',
                                               xlab='WGBS',
                                               ylab='mdRAD',
                                               ALPHA=ALPHA)
mbd.mr.resp = plot_scatter_pearsonCor_annotated(rdat,
                                                xcol='mbd.log2FoldChange',
                                                ycol='mr.log2FoldChange',
                                                xlab='MBD-seq',
                                                ylab='mdRAD',
                                                ALPHA=ALPHA)



# PLOT FIGURE 1 -----------------------------------------------------------

plot_grid(v.pm,
          v.mbd,
          v.mr,
          pm.mbd.resp,
          pm.mr.resp,
          mbd.mr.resp,
          nrow=2,
          labels = LETTERS[1:6])



# GET SOME STATS ----------------------------------------------------------

srdat = rdat %>% 
  mutate(sig.mbd=mbd.padj < 0.1,
         sig.mr = mr.padj < 0.1,
         sig.pm = qvalue < 0.1) %>% 
  dplyr::select(sig.mbd, sig.mr, sig.pm)
nsig = apply(srdat, 1, function(x) sum(x, na.rm=TRUE))
table(nsig)







