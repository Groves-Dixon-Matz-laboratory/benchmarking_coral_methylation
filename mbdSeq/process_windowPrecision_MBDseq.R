#process_MBDseq.R

rm(list=ls())
library('DESeq2')
library('tidyverse')
library('cowplot')
source('benchmarking_functions.R')
source('mbdSeq/MBD_functions.R')


# READ IN BEDTOOLS COUNTS -------------------------------------------------

bedtoolsFileList=c('mbdSeq/datasets/windowPrecision/100bp_chr1_windowRes.tsv',
                   'mbdSeq/datasets/windowPrecision/500bp_chr1_windowRes.tsv',
                   'mbdSeq/datasets/windowPrecision/1kb_chr1_windowRes.tsv',
                   'mbdSeq/datasets/windowPrecision/5kb_chr1_windowRes.tsv',
                   'mbdSeq/datasets/windowPrecision/10kb_chr1_windowRes.tsv')
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
lapply(countsList, head)



# SET UP COLDATA ----------------------------------------------------------
counts=countsList[[1]]
sample = sapply(colnames(counts), function(x) strsplit(x, '_')[[1]][1])
fraction = sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1])
genotype = sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
tissue = substr(sapply(sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3]),
                1,
                1)
coldata = data.frame(Run=sample,
                      genotype=factor(genotype, levels = c('L5', 'N12')),
                      tissue=factor(tissue, levels=c('t', 's')),
                      fraction=factor(fraction, levels=c('m', 'ub')))
rownames(coldata)=colnames(counts)
coldata



# GET METHYLATION LEVEL ---------------------------------------------------
#here use DESeq2 to calculate 'mbd-score' the log2 fold difference between
#methyated 'm' and unbound 'ub' fractions from the library preps


#Run
methLvlList = lapply(countsList, function(x) get_meth_lvl_mbd(x))
names(methLvlList)=names


#plot distributions
pltList = lapply(names(methLvlList), function(x) plot_dist(x))
plot_grid(plotlist=pltList, nrow=3)



#re-append positions
lapply(methLvlList, nrow)
lapply(posList, nrow)
lvlList = list()
#re-append
for (n in names){
  print(n)
  lvlList[[n]] = cbind(posList[[n]],
                       data.frame(methLvlList[[n]]))
}
lapply(lvlList, head)
lapply(posList, head)


#SAVE
save(lvlList, file='mbdSeq/datasets/windowPrecision/methLevelList.Rdata')



# GET RESPONSES USING UNBOUND ---------------------------------------------

#get tip vs side results
RESULTSNAME = 'tissues.fractionub'
t.ub.resList = lapply(countsList, function(x) get_response_with_ub(x))
names(t.ub.resList)=names


#get genotype results
RESULTSNAME = 'genotypeN12.fractionub'
g.ub.resList = lapply(countsList, function(x) get_response_with_ub(x))
names(g.ub.resList)=names

#tissue
resList=t.ub.resList
tvList = lapply(names, function(x) plot_volcano(x))
plot_grid(plotlist=tvList, nrow=3)

#genotype
resList=g.ub.resList
gvList = lapply(names, function(x) plot_volcano(x))
plot_grid(plotlist=gvList, nrow=3)



#RE-APPEND POSITIONS
#re-append
t.ub.resListPos = reappend_positions_to_res2(t.ub.resList, names)
g.ub.resListPos = reappend_positions_to_res2(g.ub.resList, names)
lapply(t.ub.resListPos, head)
lapply(g.ub.resListPos, head)



#SAVE
save(t.ub.resListPos, file='mbdSeq/datasets/windowPrecision/tissue_response_ub.Rdata')
save(g.ub.resListPos, file='mbdSeq/datasets/windowPrecision/genotype_response_ub.Rdata')

#HAVE NOT REPEATED THE WINDOW PRECISION ANALYSES FOR THE BASIC WAY (COMMENTED OUT BELOW), SINCE IT DIDN'T WORK WELL ANYWAY

# 
# #------- GET GBM DIFFERENCES USING IN BASIC WAY -----------#
# 
# #tissue
# CONTRAST = c('tissue', 't', 's')
# t.b.resList = lapply(countsList, function(x) get_response_without_ub(x))
# 
# #genotype
# CONTRAST = c('genotype', 'L5', 'N12')
# g.b.resList = lapply(countsList, function(x) get_response_without_ub(x))
# 
# 
# #VOLCANOS
# #tissue
# resList=t.b.resList
# tvList = lapply(names, function(x) plot_volcano(x))
# plot_grid(plotlist=tvList, nrow=3)
# 
# #genotype
# resList=g.b.resList
# gvList = lapply(names, function(x) plot_volcano(x))
# plot_grid(plotlist=gvList, nrow=3)
# 
# 
# #RE-APPEND POSITIONS
# #re-append
# t.b.resListPos = reappend_positions_to_res(t.b.resList, names)
# g.b.resListPos = reappend_positions_to_res(g.b.resList, names)
# lapply(t.ub.resListPos, head)
# lapply(g.ub.resListPos, head)
# 
# 
# #SAVE
# save(t.b.resListPos, file='mbdSeq/datasets/windowPrecision/tissue_response_basic.Rdata')
# save(g.b.resListPos, file='mbdSeq/datasets/genotype_response_basic.Rdata')
# 
# # COMPARE TWO METHODS -----------------------------------------------------
# 
# #tissue
# tissueList = list()
# for (n in names){
#   ubres = t.ub.resList[[n]]
#   bres = t.b.resList[[n]]
#   plt=compare_res(ubres, bres, xlab='Interaction (using UB)', ylab='Basic method (meth only)', subtitle=n)
#   tissueList[[n]]=plt
# }
# png(file='mbdSeq/figures/tissue_basic_v_ub.png')
# plot_grid(plotlist=tissueList, nrow=3)
# dev.off()
# 
# #genotype
# genoList = list()
# for (n in names){
#   ubres = g.ub.resList[[n]]
#   bres = g.b.resList[[n]]
#   plt=compare_res(ubres, bres, xlab='Interaction (using UB)', ylab='Basic method (meth only)', subtitle=n)
#   genoList[[n]]=plt
# }
# png(file='mbdSeq/figures/genotype_basic_v_ub.png')
# plot_grid(plotlist=genoList, nrow=3)
# dev.off()



