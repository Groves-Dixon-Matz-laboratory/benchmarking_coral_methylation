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

# GET RESPONSES USING UNBOUND ---------------------------------------------

#function to get response of a region using unbound and captured fractions
#requires global variable:RESULTSNAME
#eg: RESULTSNAME = 'tissues.fractionub'
get_response_with_ub = function(counts){
  
  print('Running DESeq for dataset:')
  print(head(counts))
  #set up input matrix for DESeq
  int.ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                       colData = coldata,
                                       design = formula(~ tissue + genotype + fraction + tissue:fraction + genotype:fraction))
  #run DESeq
  int.dds = DESeq(int.ddsHTSeq,
                  fitType = 'local')
  print('Results names:')
  print(resultsNames(int.dds))
  print(paste('Using', RESULTSNAME))
  INDEPDENDENT_FILTERING=FALSE
  res = results(int.dds, name=RESULTSNAME, independentFiltering=INDEPDENDENT_FILTERING)
  return(res)
}



#get tip vs side results
RESULTSNAME = 'tissues.fractionub'
t.ub.res = get_response_with_ub(counts)

#get genotype results
RESULTSNAME = 'genotypeN12.fractionub'
g.ub.res= get_response_with_ub(counts)



#RE-APPEND POSITIONS


###reappend positions to deseq res (assumes posList)
reappend_positions_to_resTACC = function(res, pos){
  tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
  match=sum(rownames(res)==tag)==nrow(res)
  print('Everything matching?')
  print(match)
  newDat = cbind(pos, data.frame(res))
  return(newDat)
}

#re-append
t.ub.resPos = reappend_positions_to_resTACC(t.ub.res, pos)
g.ub.resPos = reappend_positions_to_resTACC(g.ub.res, pos)


#write out as tsvs with only those passing the cutoff
writeout = function(df, outname){
  df %>% 
    filter(pvalue < rawPvalCut) %>% 
    write_tsv(path=outname)
}

writeout(t.ub.resPos, paste(outPrefix, 'tissue_response.tsv', sep='_'))
writeout(g.ub.resPos, paste(outPrefix, 'genotype_response.tsv', sep='_'))

