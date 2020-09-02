#compare_responses.R
#assemble data from each assay for comparing how 
#they measure differences in methylation between 
#treatment groups.

rm(list=ls())
source('benchmarking_functions.R')




# LOAD PICOMETHYL DATA ----------------------------------------------------


ll=load('picomethyl/datasets/methylKit_results/tissue_response.Rdata')
ll
ll=load('picomethyl/datasets/methylKit_results/genotype_response.Rdata')
ll
pm.t = t.diffListPos
pm.g = g.diffListPos
remove(t.diffListPos)
remove(g.diffListPos)
lapply(pm.t, head)
lapply(pm.g, head)
names(pm.t)
names(pm.g)


# LOAD mdRAD DATA -----------------------------------------------------

change_deseq_cols = function(rdf){
  colnames(rdf)[colnames(rdf)=='log2FoldChange']<-paste(CATSTRING, 'log2FoldChange', sep='')
  colnames(rdf)[colnames(rdf)=='pvalue']<-paste(CATSTRING, 'pvalue', sep='')
  colnames(rdf)[colnames(rdf)=='padj']<-paste(CATSTRING, 'padj', sep='')
  return(rdf)
}


ll = load('methylRAD/datasets/tissue_responses.Rdata')
ll
ll = load('methylRAD/datasets/genotype_responses.Rdata')
ll

CATSTRING='mr.'
mr.t = lapply(mr.t.list, function(x) change_deseq_cols(x))
mr.g = lapply(mr.g.list, function(x) change_deseq_cols(x))
remove(mr.t.list)
remove(mr.g.list)
lapply(mr.t, head)
lapply(mr.g, head)


#ALSO LOAD THE 8-SAMPLE SUBSET FOR mdRAD
ll = load('methylRAD/datasets/eightSample_tissue_responses.Rdata')
ll
ll = load('methylRAD/datasets/eightSample_genotype_responses.Rdata')
ll
CATSTRING='esmr.'
es.mr.t = lapply(es.mr.t.list, function(x) change_deseq_cols(x))
es.mr.g = lapply(es.mr.g.list, function(x) change_deseq_cols(x))
remove(es.mr.t.list)
remove(es.mr.g.list)
lapply(es.mr.t, head)
lapply(es.mr.g, head)

# 
# # # #temp check on volcano
gdat = data.frame(es.mr.g[['gene']])
es=plot_volcano_general(gdat,
                       xcol = 'esmr.log2FoldChange',
                       ycol='esmr.pvalue',
                       sigcol = 'esmr.padj')
gdat = data.frame(mr.g[['gene']])
x=plot_volcano_general(gdat,
                       xcol = 'mr.log2FoldChange',
                       ycol='mr.pvalue',
                       sigcol = 'mr.padj')
plot_grid(es,x)

# # # #temp check on volcano
gdat = data.frame(es.mr.t[['gene']])
es=plot_volcano_general(gdat,
                        xcol = 'esmr.log2FoldChange',
                        ycol='esmr.pvalue',
                        sigcol = 'esmr.padj')
gdat = data.frame(mr.t[['gene']])
x=plot_volcano_general(gdat,
                       xcol = 'mr.log2FoldChange',
                       ycol='mr.pvalue',
                       sigcol = 'mr.padj')
plot_grid(es,x)



# LOAD MBDSEQ DATA --------------------------------------------------------


ll=load('mbdSeq/datasets/tissue_response_ub.Rdata')
ll
ll=load('mbdSeq/datasets/genotype_response_ub.Rdata')
ll

CATSTRING='mbd.'
mbd.t = lapply(t.ub.resListPos, function(x) change_deseq_cols(x))
mbd.g = lapply(g.ub.resListPos, function(x) change_deseq_cols(x))
remove(t.ub.resListPos)
remove(g.ub.resListPos)
lapply(mbd.t, head)
lapply(mbd.g, head)


# MERGE METH DATA INTO SINGLE DFS ---------------------------------------------------

change_stop_to_end = function(df){
  colnames(df)[colnames(df)=='stop']<-'end'
  return(df)
}


#tissue
tAssayList = list(pm.t, mr.t, mbd.t, es.mr.t)
tnames = names=names(pm.t)
tList = merge_from_assay_lists(tnames, tAssayList)

#genotype
gAssayList = list(pm.g, mr.g, mbd.g, es.mr.g)
gnames = names=names(pm.g)
gList = merge_from_assay_lists(gnames, gAssayList)



# SAVE SELECT DATASETS ----------------------------------------------------

#GBM RESULTS
t.gbm.dat = tList[['gene']]
g.gbm.dat = gList[['gene']]
save(g.gbm.dat, file='comparisons/datasets/genotype_gbm_response.Rdata')
save(t.gbm.dat, file='comparisons/datasets/tissue_gbm_response.Rdata')

#PROMOTER RESULTS
t.dat = tList[['promoter']]
g.dat = gList[['promoter']]
save(g.dat, file='comparisons/datasets/genotype_promoter_response.Rdata')
save(t.dat, file='comparisons/datasets/tissue_promoter_response.Rdata')

# #500bp window RESULTS decided this was too much of a pain
# t.dat = tList[['500bp']]
# g.dat = gList[['500bp']]
# save(t.dat, file='comparisons/datasets/tissue_500bp_response.Rdata')
# save(g.dat, file='comparisons/datasets/genotype_500bp_response.Rdata')


#SAVE all genotype
save(gList, file='comparisons/datasets/all_genotype_response.Rdata')

#OUTPUT SELECTION FOR GO_MWU
#write out just WGBS data for tissue gbm
head(t.gbm.dat)

gene = t.gbm.dat$name
logp = -log(t.gbm.dat$pvalue, 10)
godat = data.frame(gene,logp)
godat %>% 
  write_csv(path='mRNA_analyses/Amil_tipVside/GO_MWU/wgbs_tissue_diff.csv')



# write out missing genes for Fisher GO -----------------------------------

#idea here is to assess if there is functional enrichment among genes that are NOT measured with each assay
head(g.gbm.dat)
length(unique(g.gbm.dat$name)) #total number of genes
nrow(g.gbm.dat)                #one row for each
sum(is.na(g.gbm.dat$meth.diff))
sum(is.na(g.gbm.dat$mbd.log2FoldChange))
sum(is.na(g.gbm.dat$mr.log2FoldChange))

#write out for go enrichment
output_go = function(col){
  go_df = data.frame('gene' = g.gbm.dat$name,
                     'missing' = as.numeric(is.na(g.gbm.dat[,col])))
  print('total missing genes:')
  print(sum(go_df$missing))
  outname = paste('./go_mwu/', col, '.csv', sep='')
  go_df %>% 
    write_csv(path = outname)
}
output_go('meth.diff')
output_go('mbd.log2FoldChange')
output_go('mr.log2FoldChange')


# LOAD TAGSEQ DATA --------------------------------------------------------

format_tagseq = function(res){
  res2 = change_deseq_cols(data.frame(res))
  res2$name = as.character(rownames(res2))
  return(res2)
}

CATSTRING='rna.'

#tissue
ll=load('rnaseq/datasets/tissue_results.Rdata')
ll
rna.t = format_tagseq(res)

#genotype
ll=load('rnaseq/datasets/genotype_results.Rdata')
ll
rna.g = format_tagseq(res)



# ADD RNA DATA TO SELECT SETS ------------------------------------------------------
names=names(tList)
select = c('gene','promoter', 'tss')

add_rna_responses = function(responseList, rnaRes, snames){
  newList = list()
  for (s in snames){
    print(s)
    newDf = responseList[[s]] %>% 
      mutate(name=as.character(name)) %>% 
      full_join(rnaRes, by = 'name')
    newList[[s]]=newDf
  }
  return(newList)
}

trnaList = add_rna_responses(tList, rna.t, select)
grnaList = add_rna_responses(gList, rna.g, select)



# COMPARE RESONSES --------------------------------------------------------

sub_and_plot = function(df){
  pltList = my_pairs(df[,c('meth.diff', 'mr.log2FoldChange', 'mbd.log2FoldChange')])
  return(pltList)
  
}


#build lists of plots
ALPHA=0.3
t.npltList = list()
g.npltList = list()
for (n in names){
  print('------')
  print(n)
  t.npltList[[n]]=sub_and_plot(tList[[n]])
  g.npltList[[n]]=sub_and_plot(gList[[n]])
}


##### PLOT
NROW=3

#GENES
#tissue
png('comparisons/figures/tissue_genes_response_cors.png')
plot_grid(plotlist=t.npltList[['gene']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/genotype_genes_response_cors.png')
plot_grid(plotlist=g.npltList[['gene']], nrow=NROW)
dev.off()

#EXONS
#tissue
png('comparisons/figures/tissue_exon_resonse_cors.png')
plot_grid(plotlist=t.npltList[['exon']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/genotype_exon_resonse_cors.png')
plot_grid(plotlist=g.npltList[['exon']], nrow=NROW)
dev.off()

#PROMOTERS
#tissue
png('comparisons/figures/tissue_promoter_response_cors.png')
plot_grid(plotlist=t.npltList[['promoter']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/genotype_promoter_response_cors.png')
plot_grid(plotlist=g.npltList[['promoter']], nrow=NROW)
dev.off()

#TSS
#tissue
png('comparisons/figures/tissue_tss_response_cors.png')
plot_grid(plotlist=t.npltList[['tss']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/genotype_tss_response_cors.png')
plot_grid(plotlist=g.npltList[['tss']], nrow=NROW)
dev.off()

#1KB
#tissue
png('comparisons/figures/tissue_1kb_response_cors.png')
plot_grid(plotlist=t.npltList[['1kb']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/genotype_1kb_response_cors.png')
plot_grid(plotlist=g.npltList[['1kb']], nrow=NROW)
dev.off()




# COMPARE RESPONSES WITH GE -----------------------------------------------



sub_and_plot_ge = function(df){
  pltList = my_pairs(df[,c('meth.diff', 'mr.log2FoldChange', 'mbd.log2FoldChange',  'rna.log2FoldChange')])
  return(pltList)
  
}


#build lists of plots
ALPHA=0.3
tge.npltList = list()
gge.npltList = list()
for (n in select){
  print('------')
  print(n)
  tge.npltList[[n]]=sub_and_plot_ge(trnaList[[n]])
  gge.npltList[[n]]=sub_and_plot_ge(grnaList[[n]])
}


##### PLOT
NROW=4

#GENES
#tissue
png('comparisons/figures/GEtissue_genes_response_cors.png')
plot_grid(plotlist=tge.npltList[['gene']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/GEgenotype_genes_response_cors.png')
plot_grid(plotlist=gge.npltList[['gene']], nrow=NROW)
dev.off()

#PROMOTERS
#tissue
png('comparisons/figures/GEtissue_promoter_response_cors.png')
plot_grid(plotlist=tge.npltList[['promoter']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/GEgenotype_promoter_response_cors.png')
plot_grid(plotlist=gge.npltList[['promoter']], nrow=NROW)
dev.off()

#TSS
#tissue
png('comparisons/figures/GEtissue_tss_response_cors.png')
plot_grid(plotlist=tge.npltList[['tss']], nrow=NROW)
dev.off()
#genotype
png('comparisons/figures/GEgenotype_tss_response_cors.png')
plot_grid(plotlist=gge.npltList[['tss']], nrow=NROW)
dev.off()












# SIG ONLY CORS -----------------------------------------------------------

genes = grnaList[['promoter']]
head(genes)


CUT=0.4
sig = genes %>% 
  filter(qvalue < CUT,
         mbd.padj < CUT)


nrow(sig)


lmpm = lm(sig$rna.log2FoldChange~sig$meth.diff)
r2 = round(summary(lmpm)$r.squared, digits=2)
pm = sig %>% 
  ggplot(aes(x=meth.diff, y=rna.log2FoldChange)) +
  geom_point() +
  labs(subtitle=bquote(R^2~"="~.(r2)))
lmmr = lm(sig$rna.log2FoldChange~sig$mr.log2FoldChange)
r2 = round(summary(lmmr)$r.squared, digits=2)
mr = sig %>% 
  ggplot(aes(x=mr.log2FoldChange, y=rna.log2FoldChange)) +
  geom_point() +
  labs(subtitle=bquote(R^2~"="~.(r2)))
lmmbd = lm(sig$rna.log2FoldChange~sig$mbd.log2FoldChange)
r2 = round(summary(lmmbd)$r.squared, digits=2)
mbd = sig %>% 
  ggplot(aes(x=mbd.log2FoldChange, y=rna.log2FoldChange)) +
  geom_point() +
  labs(subtitle=bquote(R^2~"="~.(r2)))

plot_grid(pm, mbd, mr, nrow=1)




