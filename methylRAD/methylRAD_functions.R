#methylRAD_functions.R

#get fpkm for methyl RAD counts
get_fpkm = function(n){
  print(n)
  lengths = lengthList[[n]]
  counts = countsList[[n]]
  m = apply(counts, 2, function(x) sum(x)/1e6)
  k=lengths/1e3
  fpkm=sweep(counts, 1, k, `/`) %>% 
    sweep(2, m, `/`)
  return(fpkm)
}





#Similar function to fpkm, but uses total sites for the enzyme instead
#See 'GET STATS FOR WINDOWS' in walkthrough for getting wdat
#Note global variable SITES needs to be assigned properly
get_fptm = function(n){
  x=countsList[[n]]
  mdat = merge(x, wdat, by=0)
  siteSub = mdat[,SITES]
  totSites = apply(siteSub, 1, sum)
  m = apply(x, 2, function(x) sum(x)/1e6)
  totSites = apply(siteSub, 1, sum)
  fptm=sweep(x, 1, totSites, `/`) %>% 
    sweep(2, m, `/`)
  return(fptm)
}


#funciton to get means for fpkm
get_means = function(df){
  mns=apply(df, 1, mean)
  res=data.frame(tag=names(mns),
                 mn=log(mns, 2))
  return(res)
}


#get responses using DESeq
#note requires global variable CONTRAST
#eg:
#CONTRAST=c('tissue', 't', 's')

get_response = function(counts, coldata, selectContrast){
  
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
  res = results(dds, contrast = selectContrast, independentFiltering=FALSE)
  
}

get_response_singleEnzyme = function(counts, coldata, selectContrast){
  
  print('Running DEseq on data:')
  print(head(counts))
  
  #set up input matrix for DESeq
  ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                   colData = coldata,
                                   design = formula(~ tissue + genotype))
  #run DESeq
  dds = DESeq(ddsHTSeq,
              fitType = 'local')
  #get DEseq results
  resultsNames(dds)
  res = results(dds, contrast = selectContrast, independentFiltering=FALSE)
  
}


get_vsd = function(counts){
  
  print('Running DEseq on data:')
  print(head(counts))
  
  #set up input matrix for DESeq
  ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                   colData = coldata,
                                   design = formula(~ tissue + genotype + enzyme))
  vsd = rlog(ddsHTSeq,
             fitType='local')
  vsdDf = assay(vsd)
  return(vsdDf)
}
