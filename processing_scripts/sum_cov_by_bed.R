#!/usr/bin/env Rscript
#gbm_by_glm_onSumsV2_subIterated.R

#PARSE ARUGMENTS
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  
  make_option(c("--cov"), type="character", default=NULL, 
              help="infile"),

  make_option(c("--bed"), type="character", default=NULL, 
              help="bed file to use"),

  make_option(c("--o"), type="character", default=NULL, 
              help="Name for output prefix")
)



print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
covIn = opt$cov
bedIn = opt$bed
outName = opt$o


#READ IN
print('Reading in data...')
cdat = read_tsv(covIn, col_names=c('chr', 'pos1', 'pos2', 'pctMeth', 'nM', 'nU'))
bdat = read_tsv(bedIn, col_names=c('chr', 'start', 'end', 'name'))



#LOOP THROUGH BED AND GET SUMS
print('Getting sums...')
nBeds = nrow(bdat)
chrs = unique(cdat$chr)
mSums = c()
uSums = c()
print(paste(length(chrs), 'scaffolds found in cov file'))
rdat = data.frame()

for (uchr in chrs){
  chrBed = bdat[bdat$chr==uchr,]
  chrCov = cdat[cdat$chr==uchr,]
  nBeds = nrow(chrBed)
  print('================')
  print(paste('chromosome', uchr, '...'))
  print(paste(nBeds, 'windows'))
  if (nBeds==0){
    next
  }
  for (i in 1:nBeds){
    if(i %% 100 == 0){
      print(paste(i,'of',nBeds))
    }
    brow = chrBed[i,]
    chr=as.character(brow['chr'])
    start = as.numeric(brow['start'])
    end = as.numeric(brow['end'])
    inside = (chrCov$chr==chr) & (chrCov$pos1 >= start) & (chrCov$pos1 <= end)
    sub = chrCov[inside,]
    brow$nM = sum(sub$nM)
    brow$nU = sum(sub$nU)
    rdat = rbind(rdat, brow)
  }
}



#WRITE OUT
if (nrow(rdat) > 0){
  write.table(rdat, file=outName, row.names=FALSE, quote=FALSE, sep='\t')
} else {
  print('no genes in this scaffold. Done.')
}











  