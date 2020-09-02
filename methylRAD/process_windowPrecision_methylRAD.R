#process_windowPrecision_methylRAD.R
rm(list=ls())
library('tidyverse')
library('DESeq2')
library('cowplot')
source('benchmarking_functions.R')
source('methylRAD/methylRAD_functions.R')

# READ IN BEDTOOLS COUNTS -------------------------------------------------

bedtoolsFileList=c('methylRAD/datasets/windowPrecision/100bp_chr1_windowRes.tsv',
                   'methylRAD/datasets/windowPrecision/500bp_chr1_windowRes.tsv',
                   'methylRAD/datasets/windowPrecision/1kb_chr1_windowRes.tsv',
                   'methylRAD/datasets/windowPrecision/5kb_chr1_windowRes.tsv',
                   'methylRAD/datasets/windowPrecision/10kb_chr1_windowRes.tsv')
names=c('100bp',
        '500bp',
        '1kb',
        '5kb',
        '10kb')


datList = lapply(bedtoolsFileList, function(x) read_bedtools(x, useName = TRUE))
posList = lapply(datList, function(x) return(x[['pos']]))
lengthList = lapply(posList, function(x) return(x$end-x$start))
countsList = lapply(datList, function(x) return(x[['counts']]))
names(countsList)=names
names(posList)=names
names(lengthList)=names

#split out the two enzymes
fCountsList = lapply(countsList, function(x) return(x[,grep('mrF', colnames(x))]))
mCountsList = lapply(countsList, function(x) return(x[,grep('mrM', colnames(x))]))
names(fCountsList)=names
names(mCountsList)=names


# GET FPKM -----------------------------------------------------

#run for entire set, and each enzyme individually
fpkmList = lapply(names, function(x) get_fpkm(x))
f.fpkmList = lapply(fpkmList, function(x) return(x[,grep('mrF', colnames(x))]))
m.fpkmList = lapply(fpkmList, function(x) return(x[,grep('mrM', colnames(x))]))

#get mean fkpm accross each dataset
b.lvlList = lapply(fpkmList, function(x) get_means(x))
f.lvlList = lapply(f.fpkmList, function(x) get_means(x))
m.lvlList = lapply(m.fpkmList, function(x) get_means(x))
names(b.lvlList)=names
names(f.lvlList)=names
names(m.lvlList)=names


#merge into single dataframes
lvlList = list()
for (n in names){
  blvl = b.lvlList[[n]]
  flvl = f.lvlList[[n]]
  mlvl = m.lvlList[[n]]
  colnames(blvl)=c('tag', 'mrB')
  colnames(flvl)=c('tag', 'mrF')
  colnames(mlvl)=c('tag', 'mrM')
  rList = list(blvl, flvl, mlvl)
  res = purrr::reduce(rList, full_join, by='tag')
  lvlList[[n]]=res
}
names(lvlList)=names


#re-append positions
lapply(lvlList, head)
for (n in names){
  print(n)
  lvlList[[n]] = cbind(posList[[n]], lvlList[[n]])
  # lvlList[[n]]$tag<-NULL
}
lapply(lvlList, head)


#SAVE
save(lvlList, file='methylRAD/datasets/windowPrecision/methLevelList.Rdata')




# SET UP COLDATA ----------------------------------------------------------


counts=countsList[[1]]
sample = sapply(colnames(counts), function(x) strsplit(x, '_')[[1]][1])
enzyme = sub('mr',
             '',
             sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]))
genotype = sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
tissue = substr(sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3]),
                1,
                1)


coldata = data.frame(Run=sample,
                     genotype=genotype,
                     tissue=tissue,
                     enzyme=enzyme)



# GET RESPONSES -----------------------------------------------------------

#for tissue
t.resList = lapply(countsList, function(x) get_response(x, coldata, selectContrast = c('tissue', 't', 's')))
names(t.resList)=names


#for genotype
g.resList = lapply(countsList, function(x) get_response(x, coldata, selectContrast = c('genotype', 'L5', 'N12')))
names(g.resList)=names


#CHECK VOLCANOS
#for tissue
tpltList=list()
for (n in names){
  print(n)
  rdf = t.resList[[n]]
  plt=plot_volcano_general(data.frame(rdf), 'log2FoldChange', 'pvalue', 'padj')
  tpltList[[n]]=plt
}
plot_grid(plotlist=tpltList, nrow=3)

#for genotypes
gpltList=list()
for (n in names){
  print(n)
  rdf = g.resList[[n]]
  plt=plot_volcano_general(data.frame(rdf), 'log2FoldChange', 'pvalue', 'padj')
  gpltList[[n]]=plt
}
plot_grid(plotlist=gpltList, nrow=3)



#RE-APPEND POSITIONS

#cehck names match up
lapply(t.resList, head)
lapply(t.resList, nrow)
lapply(posList, nrow)


#RE-APPEND POSITIONS
#re-append
mr.t.list = reappend_positions_to_res2(t.resList, names)
mr.g.list = reappend_positions_to_res2(g.resList, names)
lapply(mr.t.list, head)
lapply(mr.g.list, head)



#SAVE
save(mr.t.list, file='methylRAD/datasets/windowPrecision/tissue_responses.Rdata')
save(mr.g.list, file='methylRAD/datasets/windowPrecision/genotype_responses.Rdata')
