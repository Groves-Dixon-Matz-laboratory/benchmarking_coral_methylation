#process_picoMethyl.R
#organize the WGBS data from TACC

rm(list=ls())
library('DESeq2')
library('tidyverse')
library('cowplot')
source('benchmarking_functions.R')
source('picomethyl/picomethyl_functions.R')


# GET METHYLATION LEVEL ---------------------------------------------------

#list to keep the level estimates in
lvlList = list()
regionNames=c('100bp',
              '500bp',
              '1kb',
              '5kb',
              '10kb')
basicFiles = c('picomethyl/datasets/windowPrecision/100bp_chr1_chr1.basicStats',
               'picomethyl/datasets/windowPrecision/500bp_chr1_chr1.basicStats',
               'picomethyl/datasets/windowPrecision/1kb_chr1_chr1.basicStats',
               'picomethyl/datasets/windowPrecision/5kb_chr1_chr1.basicStats',
               'picomethyl/datasets/windowPrecision/10kb_chr1_chr1.basicStats') #output from basic_methylation_from_bed.R
names(basicFiles)=regionNames

#make empty place holders for other types not included here
glmFiles=c('fake')
windowFiles=c('fake')


#UPLOAD DATA FOR METHYLATION LEVEL

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
save(lvlList, file='picomethyl/datasets/windowPrecision/methLevelList.Rdata')




# READ IN THE METHYLKIT RESULTS -------------------------------------------------

#for genotype differences
gMethylKitFiles=c('picomethyl/datasets/windowPrecision/windowBoundaries_100bp_chr1_methylKit.Rdata',
                  'picomethyl/datasets/windowPrecision/windowBoundaries_500bp_chr1_methylKit.Rdata',
                  'picomethyl/datasets/windowPrecision/windowBoundaries_1kb_chr1_methylKit.Rdata',
                  'picomethyl/datasets/windowPrecision/windowBoundaries_5kb_chr1_methylKit.Rdata',
                  'picomethyl/datasets/windowPrecision/windowBoundaries_10kb_chr1_methylKit.Rdata')

#names for datasets in lists
names=c('100bp',
        '500bp',
        '1kb',
        '5kb',
        '10kb')

#load all the data
# t.allList = lapply(tMethylKitFiles, function(x) load_methylKit_res(x))
# names(t.allList)=names
g.allList = lapply(gMethylKitFiles, function(x) load_methylKit_res(x))
names(g.allList)=names

#separate out the position data
get_pos = function(x){
  res=x[['bounds']]
  return(res)
}
# t.posList = lapply(t.allList, function(x) get_pos(x))
# names(t.posList)=names
g.posList = lapply(g.allList, function(x) get_pos(x))
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

# t.diffListPos = add_names(t.diffList, t.posList)
g.diffListPos = add_names(g.allList, g.posList)


#plot genotype
gpltList = list()
for (n in names){
  mk=g.diffListPos[[n]]
  plt=plot_volcano_general(mk, xcol = 'meth.diff', ycol='pvalue', sigcol = 'qvalue')
  gpltList[[n]]=plt
}

plot_grid(plotlist=gpltList, nrow=3)


#SAVE
# save(t.diffListPos, file='picomethyl/datasets/windowPrecision/tissue_response.Rdata')
save(g.diffListPos, file='picomethyl/datasets/windowPrecision/genotype_response.Rdata')

