#process_MBDseq.R

rm(list=ls())
library('DESeq2')
library('tidyverse')
library('cowplot')
source('benchmarking_functions.R')
source('mbdSeq/MBD_functions.R')


# READ IN BEDTOOLS COUNTS -------------------------------------------------

#set up file names
file_dir = './mbdSeq/datasets'
bedtoolsFileList=list.files(file_dir,
                            pattern = 'bed_multicov.tsv',
                            full.names = TRUE)
names = sub(paste(file_dir, '/mbd_', sep=''), '', bedtoolsFileList)
names = sub('.bed_multicov.tsv', '', names)
names = sub('_Boundaries', '', names)
names = sub('Boundaries', '', names)
#fix some names for compatability
names[names=="cds"] = 'exon'
names[names=="window_1kb"] = '1kb'
print('Counts matrices to analyze:')
print(names)

#upload
datList = lapply(bedtoolsFileList, function(x) read_bedtools(x))
countsList = lapply(datList, function(x) return(x[['counts']]))
posList = lapply(datList, function(x) return(x[['pos']]))
names(countsList)=names
names(posList)=names
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



# GET RLD DIFFERENCES -----------------------------------------------------

#use rld to get estimate of absolute methylation level at individual
#level for a set of windows by taking difference between captured and unboudn fractions.
subCounts = list(countsList[['gene']],
                 countsList[['1kb']])
vsdList = lapply(subCounts, function(x) get_rld_diff(x))
mcoldata = coldata %>%
  filter(fraction=='m')
save(vsdList, mcoldata, file='mbdSeq/datasets/vsd_results.Rdata')



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


# ADD ADDITIONAL LEVEL MEASURES -------------------------------------------

#choose regions you want to add measures for
addRegions = c('gene')


#upload window stats for those regions
wStatList = list()
wStatList[['gene']]=read_window_stats('windowStats/geneBoundaries_nucStats.tsv')


#get fpkm
fpkmList = lapply(addRegions, function(x) get_mbd_fpkm(x))
names(fpkmList)=addRegions

#get fragments per CpG per million (fpcpgm)
fpcpgList = lapply(addRegions, function(x) get_fpCpGm(x))
names(fpcpgList)=addRegions


#merge these into the other level measures
for (n in addRegions){
  ldat = lvlList[[n]]
  fpkm = fpkmList[[n]]
  fpcpg = fpcpgList[[n]]
  new = merge(ldat, fpkm, by = 0)
  rownames(new)=new$Row.names
  new$Row.names<-NULL
  new2 = merge(new, fpcpg, by=0)
  rownames(new2)=new2$Row.names
  new2$Row.names<-NULL
  lvlList[[n]]=new2
}


#SAVE
save(lvlList, file='mbdSeq/datasets/methLevelList.Rdata')



# GET RESPONSES USING UNBOUND ---------------------------------------------

#get tip vs side results
RESULTSNAME = 'tissues.fractionub'
t.ub.resList = lapply(countsList, function(x) get_response_with_ub(x))
names(t.ub.resList)=names


x=countsList[[1]]
res=get_response_with_ub(x)
head(res)

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
t.ub.resListPos = reappend_positions_to_res(t.ub.resList, names)
g.ub.resListPos = reappend_positions_to_res(g.ub.resList, names)
lapply(t.ub.resListPos, head)
lapply(g.ub.resListPos, head)

#ADD ON ANY ADDITION DATASETS THAT WERE DONE ON TACC
#this is for window sizes too small to run directly
#these were run on TACC.
#See 'split by chromosome for parallelizing small windows' subsection in walkthrough

t500 = upload_response_to_add('mbdSeq/datasets/mbd_500bp_window_tissue_allResponse.tsv')
g500 = upload_response_to_add('mbdSeq/datasets/mbd_500bp_window_genotype_allResponse.tsv')
t.ub.resListPos[['500bp']] = t500
g.ub.resListPos[['500bp']] = g500
lapply(t.ub.resListPos, head)
lapply(g.ub.resListPos, head)
names(t.ub.resListPos)
names(g.ub.resListPos)


#SAVE
save(t.ub.resListPos, file='mbdSeq/datasets/tissue_response_ub.Rdata')
save(g.ub.resListPos, file='mbdSeq/datasets/genotype_response_ub.Rdata')

#------- GET GBM DIFFERENCES USING IN BASIC WAY -----------#

#tissue
CONTRAST = c('tissue', 't', 's')
t.b.resList = lapply(countsList, function(x) get_response_without_ub(x))

#genotype
CONTRAST = c('genotype', 'L5', 'N12')
g.b.resList = lapply(countsList, function(x) get_response_without_ub(x))


#VOLCANOS
#tissue
resList=t.b.resList
tvList = lapply(names, function(x) plot_volcano(x))
plot_grid(plotlist=tvList, nrow=3)

#genotype
resList=g.b.resList
gvList = lapply(names, function(x) plot_volcano(x))
plot_grid(plotlist=gvList, nrow=3)


#RE-APPEND POSITIONS
#re-append
t.b.resListPos = reappend_positions_to_res(t.b.resList, names)
g.b.resListPos = reappend_positions_to_res(g.b.resList, names)
lapply(t.ub.resListPos, head)
lapply(g.ub.resListPos, head)


#SAVE
save(t.b.resListPos, file='mbdSeq/datasets/tissue_response_basic.Rdata')
save(g.b.resListPos, file='mbdSeq/datasets/genotype_response_basic.Rdata')

# COMPARE TWO METHODS -----------------------------------------------------

#tissue
tissueList = list()
for (n in names){
  ubres = t.ub.resList[[n]]
  bres = t.b.resList[[n]]
  plt=compare_res(ubres, bres, xlab='Interaction (using UB)', ylab='Basic method (meth only)', subtitle=n)
  tissueList[[n]]=plt
}
plot_grid(plotlist=tissueList, nrow=3)


#genotype
genoList = list()
for (n in names){
  ubres = g.ub.resList[[n]]
  bres = g.b.resList[[n]]
  plt=compare_res(ubres, bres, xlab='Interaction (using UB)', ylab='Basic method (meth only)', subtitle=n)
  genoList[[n]]=plt
}
plot_grid(plotlist=genoList, nrow=3)


