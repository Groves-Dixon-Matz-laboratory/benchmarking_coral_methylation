#!/usr/bin/env Rscript
#mbdseq_differences_TACC.R
library(tidyverse)
library(DESeq2)

# PARSE ARGUMENTS ---------------------------------------------------------

suppressMessages(library(optparse))
option_list = list(
  
  make_option(c("--i"), type="character", default=NULL, 
              help="input counts"),
  
  make_option(c("--pCut"), type="character", default=NULL, 
              help="cutoff for raw pvalue (only those less will be output"),
  
  make_option(c("--o"), type="character", default=NULL, 
              help="output prefix"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
inFile = opt$i
rawPvalCut = opt$pCut
outPrefix = opt$o


# READ IN THE DATA --------------------------------------------------------


#read in bedtools results
read_bedtools = function(bedtoolsFile, useName=FALSE){
  print(paste('Reading in file', bedtoolsFile))
  counts = read.table(bedtoolsFile, header = TRUE)
  positions = counts[,1:4]
  #make a tag of the name_chr_start_end
  name_chr = paste(counts$name, counts$chr, sep='_')
  start_end = paste(counts$start, counts$end, sep='_')
  tag = paste(name_chr, start_end, sep='_')
  rownames(counts)=tag
  if (useName==TRUE){
    rownames(counts)=counts$name
  }
  counts=counts[,5:ncol(counts)]
  return(list('counts'=counts,
              'pos'=positions))
}


dat = read_bedtools(inFile)
counts = dat[['counts']]
pos = dat[['pos']]


# SET UP COLDATA ----------------------------------------------------------

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

get_response = function(counts){
  
  print('Running DEseq on data:')
  print(head(counts))
  
  #set up input matrix for DESeq
  ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                   colData = coldata,
                                   design = formula(~ tissue + genotype + enzyme))
  #run DESeq
  dds = DESeq(ddsHTSeq,
              fitType = 'local')
  #get DEseq results
  resultsNames(dds)
  res = results(dds, contrast = CONTRAST, independentFiltering=FALSE)
  
}


#for tissue
CONTRAST=c('tissue', 't', 's')
t.res = get_response(counts)


#for genotype
CONTRAST=c('genotype', 'L5', 'N12')
g.res = get_response(counts)



# RE-APPEND POSITIONS -----------------------------------------------------
reappend = function(res.df, pos){
  tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
  match=sum(rownames(res.df)==tag)==nrow(res.df)
  print('all matched up?')
  print(match)
  newRes = cbind(pos, data.frame(res.df))
  return(newRes)
}

t.resPos = reappend(t.res, pos)
g.resPos = reappend(g.res, pos)


# WRITE OUT ---------------------------------------------------------------

#write out as tsvs with only those passing the cutoff
writeout = function(df, outname){
  df %>% 
    filter(pvalue < rawPvalCut) %>% 
    write_tsv(path=outname)
}

writeout(t.resPos, paste(outPrefix, 'tissue_response.tsv', sep='_'))
writeout(g.resPos, paste(outPrefix, 'genotype_response.tsv', sep='_'))