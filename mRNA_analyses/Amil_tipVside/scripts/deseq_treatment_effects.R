#deseq_treatment_effects.R 

#SET UP THE DATA TO RUN DESEQ
library(DESeq2)
library(cowplot)
library(tidyverse)
rm(list=ls())

source('rnaseq/Amil_tipVside/scripts/rnaseq_functions.R')


#load the data
ll = load('rnaseq/Amil_tipVside/datasets/deseqInput.Rdata')
ll
coldata
head(counts)
sum(colnames(counts) == coldata$sample)==ncol(counts)


#get variance stabilized counts for timepoint 2
dds<-DESeqDataSetFromMatrix(counts,
	colData = coldata, 
	design = formula(~ tissue +
	                   genotype)
	)

INDEPDENDENT_FILTERING = FALSE
dds <- DESeq(dds)
resultsNames(dds)
t.res = results(dds, contrast = c('tissue', 't', 's'), independentFiltering=INDEPDENDENT_FILTERING)
g.res = results(dds, contrast = c('genotype', 'L5', 'N12'), independentFiltering=INDEPDENDENT_FILTERING)


#volcano for tissue
tvol = data.frame(t.res) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
	geom_point(alpha=0.1) +
	scale_color_manual(values=c('black', 'red')) + 
	labs(subtitle='Tip vs side',
	     y=bquote(-log[10]~pvalue),
	     x=bquote(log[2]~fold~difference))  +
  theme(legend.position='none')

#volcano for genotype
gvol=data.frame(g.res) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) + 
  labs(subtitle='Genotype',
       y=bquote(-log[10]~pvalue),
       x=bquote(log[2]~fold~difference)) +
  theme(legend.position='none')

plot_grid(tvol, gvol)


#save results
save(dds, file='rnaseq/Amil_tipVside/datasets/dds.Rdata')
res=t.res
save(res, file='rnaseq/Amil_tipVside/datasets/tissue_results.Rdata')
res=g.res
save(res, file='rnaseq/Amil_tipVside/datasets/genotype_results.Rdata')


#write out for GO_MWU
head(t.res)
goOut = data.frame(gene=rownames(t.res), logp=-log(t.res$pvalue, 10)) %>% 
  mutate(logp=if_else(t.res$log2FoldChange<0,
                      logp*-1,
                      logp))
hist(goOut$logp)
write.csv(goOut, file='rnaseq/Amil_tipVside/GO_MWU/tipVsideGo.csv', row.names=F, quote=F)

