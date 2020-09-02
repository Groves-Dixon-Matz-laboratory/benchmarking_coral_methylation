#initialize_counts.R
#prepare reads for DESeq and get gist of data

library('DESeq2')
library('cowplot')
theme_set(theme_cowplot())
library('tidyverse')
rm(list=ls())


#upload read counts from htseq
counts = read.table("mRNA_analyses/Amil_tipVside/datasets/rnaseq_gene_counts_DEDUPPED.txt", header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]
revisedNames = sapply(colnames(counts), function(x) strsplit(x,'_')[[1]][1])
colnames(counts) = revisedNames
head(counts)
length(unique(colnames(counts)))


#remove non-gene objects from counts data
dim(counts)
tots = apply(counts, 2, sum)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))


#get count totals
tots = apply(counts, 2, sum)
hist(tots)
write.table(data.frame(tots), file="total_counted_on_genes.tsv", sep="\t", quote=F) 
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))


#remove genes with low coverage
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]



#SET UP COLDATA
#first revise the miss-sexed individuals
sample = colnames(counts)
coldata=tibble(
  sample=sample,
  tissue=factor(if_else(grepl('T', sample),
                   't',
                   's')),
  genotype=factor(if_else(grepl('L5', sample),
                   'L5',
                   'N12'))
  ) %>% 
  data.frame()
coldata


#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
#set up input matrix for DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~1))

#run DESeq
dds = DESeq(ddsHTSeq)

#get DEseq results
res = results(dds)

#get variance stabilized counts and save them
rld = rlog(dds)
rld.df=assay(rld)
colnames(rld.df) = colnames(counts)

#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================

datExpr0 = as.data.frame(t(rld.df));

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have geneswith too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
rld.df=t(datExpr0)
rld=rld[rownames(rld.df),]
dim(rld.df)
dim(rld)
nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#build sample heatmaps 
library(pheatmap)
phm=pheatmap(cor(rld.df))

#plot pca
NTOP=nrow(counts)
group='tissue'
pca <- plotPCA(rld, intgroup=c(group), returnData=TRUE, ntop=2000)
addx=3
addy=2
pca %>% 
  mutate(genotype=substr(name, 1,1),
         tissue = if_else(tissue=='s',
                          'side',
                          'tip')) %>% 
  ggplot(aes(x=PC1, y=PC2, shape=tissue, color=genotype)) +
  geom_point(size=5) +
  lims(x=c(min(pca$PC1)-addx, max(pca$PC1+addx)),
       y=c(min(pca$PC2)-addy, max(pca$PC2+addy))) +
  labs(title = 'RNA')
  



#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
# 
#=====================================================================================

#Remove outliers by setting a branch cut threshold
# Plot a line to show the cut
cut.height = 110
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = cut.height, col = "red", lty = 2);

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.height, minSize = 4)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
keepSampleNames = rownames(datExpr0)[keepSamples]
outlierNames = rownames(datExpr0)[clust==0]
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr) #number of samples left after outlier removal
print(paste(length(outlierNames), "samples were flagged as outliers and removed:"))
outlierNames
print(paste(nSamples, "samples were kept"))


#replot heatmap without outlier
rld.df = rld.df[, !colnames(rld.df) %in% outlierNames]
coldata = coldata[!rownames(coldata) %in% outlierNames, ]
pheatmap(cor(rld.df), labels_col= substr(coldata$treatment, start=1,stop=1))



#save the outlier names so you can optionally remove them in other analyses
# save(outlierNames, file = 'datasets/outliers.Rdata')
counts=counts[,!colnames(counts) %in% outlierNames]
coldata=coldata[!coldata$sample %in% outlierNames,]
save(counts, coldata, file="rnaseq/Amil_tipVside/datasets/deseqInput.Rdata")
save(rld.df, rld, file="rnaseq/Amil_tipVside/datasets/rld.Rdata")

