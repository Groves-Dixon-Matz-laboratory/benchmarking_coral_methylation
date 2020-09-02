#process_methylRAD.R

rm(list=ls())
library('tidyverse')
library('DESeq2')
library('cowplot')
source('benchmarking_functions.R')
source('methylRAD/methylRAD_functions.R')


# PREPARE AN 8-SAMPLE DATASET  --------------------------------------------

#this is for comparison to MBD-seq and WGBS with the same sample sizes
#will be incorperated only during the
esdat = read_tsv('methylRAD/datasets/mr_geneBoundaries.bed_multicov.tsv',
                 col_names = TRUE)
cnames = colnames(esdat)
keep1 = cnames[1:4]
keep2 = cnames[grepl('^mrF', cnames)]
keep3 = keep2[!grepl('-t3', keep2) & !grepl('-s3', keep2)]
length(keep3) #here are the 8 samples to match with WGBS 2 replicates for each genotype-tissue, (one FspE1 and one MspJ1)
keep = append(keep1, keep3)
esdat2 = esdat %>%
  dplyr::select(keep)
dim(esdat2)

#write out the reduced set
esdat2 %>%
  write_tsv(path='methylRAD/datasets/eightSample_gene_multicov.tsv')


# READ IN BEDTOOLS COUNTS -------------------------------------------------

add_tag_for_pos = function(x){
  x['pos'] %>% 
    unite(tag=name,chr,start,stop)
}

#set up file names
file_dir = './methylRAD/datasets'
bedtoolsFileList=list.files(file_dir,
                            pattern = 'bed_multicov.tsv',
                            full.names = TRUE)
names = sub(paste(file_dir, '/mr_', sep=''), '', bedtoolsFileList)
names = sub('.bed_multicov.tsv', '', names)
names = sub('_Boundaries', '', names)
names = sub('Boundaries', '', names)
#fix some names for compatability
names[names=="cds"] = 'exon'
names[names=="1kbWindow"] = '1kb'
print('Counts matrices to analyze:')
print(names)

#load the files
datList = lapply(bedtoolsFileList, function(x) read_bedtools(x))
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

#---- uplaod the 8-sample set separately ----#
#This will be used downstream for the reduced sample size response analysis
es.FileList=c('methylRAD/datasets/eightSample_gene_multicov.tsv')
es.names=c('gene')
es.datList = lapply(es.FileList, function(x) read_bedtools(x))
es.posList = lapply(es.datList, function(x) return(x[['pos']]))
es.lengthList = lapply(es.posList, function(x) return(x$end-x$start))
es.countsList = lapply(es.datList, function(x) return(x[['counts']]))
names(es.countsList)=es.names
names(es.posList)=es.names
names(es.lengthList)=es.names


# GET FPKM -----------------------------------------------------

#run for entire set, and each enzyme individually
fpkmList = lapply(names, function(x) get_fpkm(x))
f.fpkmList = lapply(fpkmList, function(x) return(x[,grep('mrF', colnames(x))]))
m.fpkmList = lapply(fpkmList, function(x) return(x[,grep('mrM', colnames(x))]))
names(fpkmList)=names
names(f.fpkmList)=names
names(m.fpkmList)=names

#output the GBM fpkms for pca_and_adonis.R
geneFpkm = fpkmList[['gene']]
save(geneFpkm, file='methylRAD/datasets/mdRAD_geneFPKMs.Rdata')


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



# GET READS PER SITE ------------------------------------------------------

wdat = read_window_stats('windowStats/geneBoundaries_nucStats.tsv')


#assign region types to run on
fptmRegions = c('gene')

#get mean fptm for FSPE1
SITES=c('nCCG', 'nCGG') #cut sites for FSPe1
f.fptm0 = lapply(fptmRegions, function(x) get_fptm(x))
f.fptm = lapply(f.fptm0, function(x) return(x[,grep('mrF', colnames(x))]))
f.mn.fptm = lapply(f.fptm, function(x) get_means(x))
names(f.mn.fptm)=fptmRegions

#get means for mspj1
SITES=c('nCGNR', 'nYNCG') #cut sites for Mspj1
m.fptm0 = lapply(fptmRegions, function(x) get_fptm(x))
m.fptm = lapply(m.fptm0, function(x) return(x[,grep('mrM', colnames(x))]))
m.mn.fptm = lapply(m.fptm, function(x) get_means(x))
names(m.mn.fptm)=fptmRegions

#get means accross both
SITES=c('nCCG', 'nCGG', 'nCGNR', 'nYNCG') #cut sites for both
b.fptm = lapply(fptmRegions, function(x) get_fptm(x))
b.mn.fptm = lapply(b.fptm, function(x) get_means(x))
names(b.mn.fptm)=fptmRegions

#add these onto the lvlList
for (n in fptmRegions){
  ldat = lvlList[[n]]
  fdat = f.mn.fptm[[n]] %>% 
    dplyr::rename(mrF.s = mn)
  mdat = m.mn.fptm[[n]] %>% 
    dplyr::rename(mrM.s = mn)
  bdat = b.mn.fptm[[n]] %>% 
    dplyr::rename(mrB.s = mn)
  ndatList = list(fdat, mdat, bdat)
  ndat = purrr::reduce(ndatList, full_join, by='tag')
  lvlList[[n]] = ldat %>% 
    left_join(ndat, by = 'tag')
}


#SAVE
save(lvlList, file='methylRAD/datasets/methLevelList.Rdata')


#CHECK DISTRIBUTIONS
pltList = list()
for (n in names){
  ldf = lvlList[[n]]
  plt=ldf %>% 
    ggplot(aes(x=mrB)) +
    geom_density() +
    labs(subtitle=n)
  pltList[[n]]=plt
}
plot_grid(plotlist=pltList, nrow=3)


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
CONTRAST=c('tissue', 't', 's')
t.resList = lapply(countsList, function(x) get_response(x, coldata, CONTRAST))
names(t.resList)=names


#for genotype
CONTRAST=c('genotype', 'L5', 'N12')
g.resList = lapply(countsList, function(x) get_response(x, coldata, CONTRAST))
names(g.resList)=names


#CHECK VOLCANOS

#for tissue
tpltList=list()
for (n in names){
  print(n)
  rdf = data.frame(t.resList[[n]])
  plt=plot_volcano_general(rdf)
  tpltList[[n]]=plt
}
plot_grid(plotlist=tpltList, nrow=3)

#for genotypes
gpltList=list()
for (n in names){
  print(n)
  rdf = data.frame(g.resList[[n]])
  plt=plot_volcano_general(rdf)
  gpltList[[n]]=plt
}
plot_grid(plotlist=gpltList, nrow=3)



#RE-APPEND POSITIONS

#cehck names match up
lapply(t.resList, head)
lapply(t.resList, nrow)
lapply(posList, nrow)

for (n in names){
  pos=posList[[n]]
  tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
  match=sum(rownames(t.resList[[n]])==tag)==nrow(t.resList[[n]])
  print(paste(n,match,sep=' = '))
}


#RE-APPEND POSITIONS
#re-append
mr.t.list = reappend_positions_to_res(t.resList, names)
mr.g.list = reappend_positions_to_res(g.resList, names)
lapply(mr.t.list, head)
lapply(mr.g.list, head)


# REPEAT RESPONSE ANALYSIS FOR 8-SAMPLE SET -------------------------------

#SET UP NEW COLDATA
es.counts=es.countsList[[1]]
es.sample = sapply(colnames(es.counts), function(x) strsplit(x, '_')[[1]][1])
es.enzyme = sub('mr',
             '',
             sapply(es.sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]))
es.genotype = sapply(es.sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
es.tissue = substr(sapply(es.sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3]),
                1,
                1)


es.coldata = data.frame(Run=es.sample,
                     genotype=es.genotype,
                     tissue=es.tissue,
                     enzyme=es.enzyme)
es.coldata #coldata for the 8-sample reduced set
nrow(es.coldata)

#GET RESPONSES

#for tissue
CONTRAST=c('tissue', 't', 's')
es.t.resList = lapply(es.countsList, function(x) get_response_singleEnzyme(x, es.coldata, CONTRAST))
names(es.t.resList)=es.names


#for genotype
CONTRAST=c('genotype', 'L5', 'N12')
es.g.resList = lapply(es.countsList, function(x) get_response_singleEnzyme(x, es.coldata, CONTRAST))
names(es.g.resList)=es.names

#RE-APPEND POSITIONS

#cehck names match up
lapply(es.t.resList, head)
lapply(es.t.resList, nrow)
lapply(es.posList, nrow)

for (n in es.names){
  pos=es.posList[[n]]
  tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
  match=sum(rownames(es.t.resList[[n]])==tag)==nrow(es.t.resList[[n]])
  print(paste(n,match,sep=' = '))
}


#RE-APPEND POSITIONS
#re-append
es.mr.t.list = reappend_positions_to_res(es.t.resList, es.names)
es.mr.g.list = reappend_positions_to_res(es.g.resList, es.names)
lapply(es.mr.t.list, head)
lapply(es.mr.g.list, head)

#check volcano
gdat = data.frame(es.mr.g.list[['gene']])
x=plot_volcano_general(gdat)
x


#SAVE THE 8-SAMPLE REDUCED RESULTS
save(es.mr.t.list, file='methylRAD/datasets/eightSample_tissue_responses.Rdata')
save(es.mr.g.list, file='methylRAD/datasets/eightSample_genotype_responses.Rdata')



# ADD RESPONSE RESULTS FOR SMALL WINDOWS ----------------------------------
#PROBABLY NOT GOING TO USE THESE, SINCE THEY ARE SO CUMBERSOME BUT KEEPING FOR REFERENCE
filePath='methylRAD/datasets/mr_500bp_window_genotype_allResponse.tsv'
t500 = upload_response_to_add('methylRAD/datasets/mr_500bp_window_tissue_allResponse.tsv')
g500 = upload_response_to_add('methylRAD/datasets/mr_500bp_window_genotype_allResponse.tsv')
mr.t.list[['500bp']] = t500
mr.g.list[['500bp']] = g500



#SAVE
save(mr.t.list, file='methylRAD/datasets/tissue_responses.Rdata')
save(mr.g.list, file='methylRAD/datasets/genotype_responses.Rdata')



# GET VSD RESULTS ---------------------------------------------------------

#Didn't end up needing this but keeping for refernce

# #subset for just gene and 1kb
# subCounts = list(countsList[['gene']],
#                  countsList[['1kb']])
# 
# #run vst on each set of counts
# vsdList0 = lapply(subCounts, function(x) get_vsd(x))
# names(vsdList0)=c('gene', '1kb')
# vsdList=vsdList0
# 
# #save
# save(vsdList, coldata, file='methylRAD/datasets/vst_results.Rdata')


