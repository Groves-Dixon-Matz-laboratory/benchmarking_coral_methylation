#process_picoMethyl.R
#organize the WGBS data from TACC

rm(list=ls())
library('DESeq2')
library('tidyverse')
library('cowplot')
source('benchmarking_functions.R')
source('picomethyl/picomethyl_functions.R')


# PROCESS WINDOW METHYLATION COUNTS BY SAMPLE -----------------------------

#these are for the multivatiate analyses PCA and adonis
tt = read_tsv('picomethyl/datasets/methylKit_treatment_tables/genotype_methylKit0_inputTable.txt')
samples = tt$id
rgCountFiles = c('picomethyl/datasets/regionCounts/genotype_genes_methylKit_regionCounts.Rdata',
                 'picomethyl/datasets/regionCounts/genotype_1KbWindows_methylKit_regionCounts.Rdata')
rgNames = c('gene',
            '1kb')

get_meth_proportion_from_reg_counts = function(rgFile){
  ll=load(rgFile)
  rc = data.frame(my.reg.counts) %>% 
    as_tibble()
  pos = rc %>% 
    dplyr::select(c('chr', 'start', 'end'))
  ccols = colnames(rc)[grep('numCs', colnames(rc))]
  snums = as.numeric(sub('numCs', '', ccols))
  meth = rc %>% 
    dplyr::select(grep('numCs', colnames(rc))) %>% 
    data.frame()
  unmeth = rc %>% 
    dplyr::select(grep('numTs', colnames(rc))) %>% 
    data.frame()
  mdat = pos
  samples = ccols
  for (i in snums){
    sample = samples[i]
    methCol = paste('numCs', i, sep='')
    unmethCol = paste('numTs', i, sep='')
    mVec = meth[,methCol]
    umVec = unmeth[,unmethCol]
    rVec = mVec / (mVec + umVec)
    mdat[,sample] = rVec
  }
  return(mdat)
}
  
mList = lapply(rgCountFiles, function(x) get_meth_proportion_from_reg_counts(x))
names(mList) = rgNames
#also get samples
save(mList, samples, file='picomethyl/datasets/regionCounts/mprops_by_sample.Rdata')


# GET METHYLATION LEVEL ---------------------------------------------------

#set up file names
file_dir = './picomethyl/datasets/basic_stats'
basicFiles=list.files(file_dir,
                      pattern = '_basicStatsBed.tsv',
                      full.names = TRUE)
regionNames = sub(paste(file_dir, '/', sep=''), '', basicFiles)
regionNames = sub('_basicStatsBed.tsv', '', regionNames)
regionNames = sub('_Boundaries', '', regionNames)
regionNames = sub('Boundaries', '', regionNames)
#fix some names for compatability
regionNames[regionNames=="cds"] = 'exon'
regionNames[regionNames=="window_1kb"] = '1kb'
print('Counts matrices to analyze:')
print(regionNames)
names(basicFiles) = regionNames

#build set of window files from the set of basicStatsBed.tsv files
#build these files with nucleotide_stats_from_bed.py (see picoMethyl_data_prcocessing_pipeline.txt)
windowFiles = sub('./picomethyl/datasets/basic_stats/', '', basicFiles)
windowFiles = sub('_basicStatsBed.tsv', '', windowFiles)
windowFiles = paste('windowStats/', windowFiles, sep='')
windowFiles = paste(windowFiles, 'nucStats.tsv', sep='_')
names(windowFiles)=regionNames

#also add single glm file for gbm (this was just one other way to measure gbm)
glmFiles = c('picomethyl/datasets/glm_lvls/genes_glmLvls.Rdata') #(output from meth_lvl_by_glm.R)
names(glmFiles)=c('gene')

#check it all lines up
data.frame('basic' = basicFiles,
           'window' = windowFiles,
           'regions' = regionNames)


#UPLOAD THE GBM LEVEL RESULTS

#list to keep the level estimates in
lvlList = list()

#loop through list and upload
for (i in 1:length(regionNames)){
  #grab the file names
  rName = regionNames[i]
  glmFile = glmFiles[rName]
  windowFile = windowFiles[rName]
  basicFile = basicFiles[rName]
  print(paste(rName, glmFile, sep=' : '))
  #read in the datasets
  new.basic.dat = read_tsv(basicFile) %>% 
    mutate(mCpG_per_bp = Nmcpg / (end-start),
           l.fracMeth = log2_convesion_noZero(fracMeth))
  #if we have window stats, use them to get more detailed methylation measures
  if (!is.na(windowFile)){
    new.w.dat = read_tsv(windowFile) %>% 
      mutate(nTot = nA+nG+nT+nC,
             cpgoe = (nCpG/nTot) / ( (nC/nTot)*(nG/nTot) ) )
    bw.dat = new.basic.dat %>% 
      left_join(new.w.dat, by = c('chr', 'start', 'end', 'name')) %>% 
      mutate(mCpG_per_bp = Nmcpg / (end-start),
             l.mCpG_per_bp = log2_convesion_noZero(mCpG_per_bp),
             l.mCpG_per_CpG = log2_convesion_noZero(mCpG_per_CpG),
             l.fracMeth = log2_convesion_noZero(fracMeth))
  } else{
    bw.dat = new.basic.dat
  }
  #if we have glm estimates of methylation level, merge with those
  if (!is.na(glmFile)){
    print(paste('merging basic with glm data for:', rName))
    new.glm.dat = upload_glm(glmFile)
    new.dat = new.glm.dat %>% 
      full_join(bw.dat, by = c('chr', 'start', 'end', 'name'))
  } else {
    print(paste('no glm data found for:', rName))
    new.dat = bw.dat
  }
  lvlList[[rName]]=new.dat
}

#check results
lapply(lvlList, head)
lapply(lvlList, nrow)

#CHECK DISTRIBUTIONS

pltList = lapply(names(lvlList), function(x) plot_frac_dist(x))
plot_grid(plotlist=pltList, nrow=3)


#SAVE
save(lvlList, file='picomethyl/datasets/methLevelList.Rdata')



# READ IN THE METHYLKIT RESULTS -------------------------------------------------

#for tissue differences
tMethylKitFiles=c('picomethyl/datasets/methylKit_results/tissue_genes_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/tissue_exons_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/tissue_promoters_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/tissue_tss_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/tissue_1KbWindows_methylKit.Rdata')

#for genotype differences
gMethylKitFiles=c('picomethyl/datasets/methylKit_results/genotype_genes_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/genotype_exons_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/genotype_promoters_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/genotype_tss_methylKit.Rdata',
                  'picomethyl/datasets/methylKit_results/genotype_1KbWindows_methylKit.Rdata')

#names for datasets in lists
names=c('gene',
        'exon',
        'promoter',
        'tss',
        '1kb')

#load all the data
t.allList = lapply(tMethylKitFiles, function(x) load_methylKit_res(x))
g.allList = lapply(gMethylKitFiles, function(x) load_methylKit_res(x))
names(t.allList)=names
names(g.allList)=names



#separate out the position data (also fixes underscore in exon names back to dash)
get_pos = function(x){
  res=x[['bounds']]
  res$name = sub('_RA:cds', '-RA:cds', res$name, fixed = TRUE)
  return(res)
}
t.posList = lapply(t.allList, function(x) get_pos(x))
g.posList = lapply(g.allList, function(x) get_pos(x))
names(t.posList)=names
names(g.posList)=names




# GET RESPONSES --------------------------------------------------------

#APPEND POSITION DATA
add_names = function(diffList, posList){
  newList = list()
  for (n in names){
    pdat = posList[[n]]
    ddat = data.frame(diffList[[n]]$diff,
                      stringsAsFactors=FALSE)
    ddat$chr = as.character(ddat$chr)
    pdat$chr = as.character(pdat$chr)
    mdat = ddat %>% 
      left_join(pdat, by = c('chr', 'start', 'end'))
    newList[[n]]=mdat
  }
  return(newList)
}

t.diffListPos = add_names(t.allList, t.posList)
g.diffListPos = add_names(g.allList, g.posList)


#plot tissue
tpltList = list()
for (n in names){
  mk=t.diffListPos[[n]]
  plt=plot_volcano_general(mk, xcol = 'meth.diff', ycol='pvalue', sigcol = 'qvalue')
  tpltList[[n]]=plt
}
tpltList[['gene']]
# plot_grid(plotlist=tpltList, nrow=3)


#plot genotype
gpltList = list()
for (n in names){
  mk=g.diffListPos[[n]]
  plt=plot_volcano_general(mk, xcol = 'meth.diff', ycol='pvalue', sigcol = 'qvalue')
  gpltList[[n]]=plt
}
gpltList[['gene']]
# plot_grid(plotlist=gpltList, nrow=3)


#SAVE
save(t.diffListPos, file='picomethyl/datasets/methylKit_results/tissue_response.Rdata')
save(g.diffListPos, file='picomethyl/datasets/methylKit_results/genotype_response.Rdata')

