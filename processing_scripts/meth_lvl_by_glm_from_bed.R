#!/usr/bin/env Rscript
#gbm_by_glm_onSumsV2_subIterated.R

#PARSE ARUGMENTS
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  
  make_option(c("--i"), type="character", default=NULL, 
              help="infile"),

  make_option(c("--n"), type="integer", default=3000, 
              help="size of chunks to use"),

  make_option(c("--iter"), type="double", default=FALSE, 
              help="N iterations to run (number of times to calculate coefficients from chunks to average over)"),

  make_option(c("--o"), type="character", default=NULL, 
              help="Name for output prefix")
)



print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile = opt$i
outputPrefix = opt$o
chunkSize = opt$n
nIter = opt$iter



# LOAD DATA ---------------------------------------------------------------

dat = read.table(infile, header = TRUE, stringsAsFactors=FALSE)



#set up chucks
genes = unique(dat$name)
nGenes = length(genes)
geneChunks = list()
leftOver = genes
nChunk = ceiling(nGenes / chunkSize)
for (i in 1:(nChunk-1)){
  chunk = sample(leftOver, chunkSize)
  geneChunks[[i]]=chunk
  leftOver = leftOver[!leftOver %in% chunk]
}
geneChunks[[nChunk]]=leftOver


#Iterate through chunks
res=data.frame()
for (iteration in 1:nIter){
  print('-----------')
  print(paste('iteration', iteration))
  for (i in 1:length(geneChunks)){
    print('==')
    print(paste(c('chunk', i, 'of', length(geneChunks)), collapse=' '))
    subKeep = geneChunks[[i]]
    df = dat %>% 
      filter(name %in% subKeep)
    mylogit <- glm(cbind(nM,nU) ~ name, 
                   data=df, 
                   family = "binomial")
    coefDf = data.frame(mylogit$coefficients)
    int = as.numeric(coefDf['(Intercept)','mylogit.coefficients'])
    coefDf$norm = coefDf$mylogit.coefficients + int
    coefDf = coefDf[!grepl('Intercept', rownames(coefDf)),]
    coefDf$name=sub('name', '', rownames(coefDf))
    coefDf=coefDf[,c('name', 'mylogit.coefficients', 'norm')]
    coefDf$iteration = i
    res=rbind(res, coefDf)
  }
}


#get final meth levl estimates
fres = res %>% 
  group_by(name) %>% 
  summarize(glmLvl = mean(norm, na.rm=TRUE)) %>%
  left_join(dat, by = 'name') %>%
  select(chr, start, end, name, nM, nU, glmLvl)

#Save results
save(fres, file=paste(outputPrefix, 'glmLvls.Rdata', sep='_'))
write.table(fres,
  file=paste(outputPrefix, 'glmLvls.tsv', sep='_'),
  sep='\t',
  quote=FALSE,
  col.names=TRUE,
  row.names=FALSE)
