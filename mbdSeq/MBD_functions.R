#MBD_functions.R

#function to run deseq to get methylation level by comparing unbound vs captured fractions
get_meth_lvl_mbd = function(counts){
  
  print('Getting meth levl for dataset:')
  print(head(counts))
  
  #set up input matrix for DESeq
  ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                   colData = coldata,
                                   design = formula(~ fraction + genotype + tissue))
  #run DESeq
  dds = DESeq(ddsHTSeq,
              fitType = 'local')
  #get DEseq results
  resultsNames(dds)
  res = results(dds, contrast = c('fraction', 'm', 'ub'), independentFiltering=FALSE)
  return(res)
}


#function to get fpkm from captured fraction counts
get_mbd_fpkm = function(n){
  print(paste(n,'...',sep=''))
  dat = countsList[[n]]
  dat = dat %>% 
    dplyr::select(grep('^m', colnames(dat)))
  wdat = wStatList[[n]]
  mdat = merge(dat, wdat, by=0)
  m = apply(dat, 2, function(x) sum(x)/1e6)
  lengths = mdat$end - mdat$start
  k=lengths/1e3
  fpkm=sweep(dat, 1, k, `/`) %>% 
    sweep(2, m, `/`)
  mns = apply(fpkm, 1, function(x) mean(x, na.rm=TRUE))
  nonZero = min(mns[mns>0])
  mns[mns==0]<-nonZero
  res = data.frame(mbd.FPKM=log(mns, 2))
  return(res)
}

#similar function for fragments per CpG per million reads
get_fpCpGm = function(n){
  dat = countsList[[n]]
  dat = dat %>% 
    dplyr::select(grep('^m', colnames(dat)))
  wdat = wStatList[[n]]
  mdat = merge(dat, wdat, by=0)
  nCpG = mdat$nCpG
  m = apply(dat, 2, function(x) sum(x)/1e6)
  fpcpg=sweep(dat, 1, nCpG, `/`) %>% 
    sweep(2, m, `/`)
  mns = apply(fpcpg, 1, function(x) mean(x, na.rm=TRUE))
  nonZero = min(mns[mns>0])
  mns[mns==0]<-nonZero
  res = data.frame(mbd.FPCpGM=log(mns, 2))
  return(res)
}


#check distributions
plot_dist = function(name){
  res=methLvlList[[name]]
  data.frame(res) %>% 
    ggplot(aes(x=log2FoldChange)) +
    geom_density() +
    labs(subtitle=name)
}

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


#function to get response using only captured fraction
#requires global variable: CONTRAST
#eg: CONTRAST = c('tissue', 't', 's')
get_response_without_ub = function(counts){
  mcounts = counts %>% 
    dplyr::select(grep('^m', colnames(counts)))
  mcoldata = coldata[colnames(mcounts),]
  print('Running Deseq on methylated counts only:')
  print(head(mcounts))
  #set up input matrix for DESeq
  m.ddsHTSeq<-DESeqDataSetFromMatrix(mcounts,
                                     colData = mcoldata,
                                     design = formula(~ genotype + tissue))
  
  #run DESeq
  m.dds = DESeq(m.ddsHTSeq,
                fitType = 'local')
  
  #get DEseq results
  resultsNames(m.dds)
  bg.res = results(m.dds, contrast = CONTRAST, independentFiltering=FALSE)
  
}



#function to plot volcano plot
plot_volcano = function(name){
  res=resList[[name]]
  data.frame(res) %>% 
    mutate(sig=padj<0.1) %>% 
    ggplot(aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
    geom_point(alpha=0.5) +
    labs(subtitle=name) +
    scale_color_manual(values=c('black', 'red'))
}


#function to compare results files from DESEQ
compare_res = function(res1, res2, xlab, ylab, subtitle=''){
  df1 = data.frame(res1)
  df2 = data.frame(res2)
  mdf = merge(df1, df2, by = 0)
  lm1=lm(mdf$log2FoldChange.y ~ mdf$log2FoldChange.x)
  print(summary(lm1))
  mdf %>% 
    ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
    geom_point() +
    geom_smooth(method='lm') +
    labs(x=xlab, y=ylab)
  
}

#get rld differences unbound vs flowthru
#this provides an estimate of absolute methylation level
#for each individual
get_rld_diff = function(counts){
  print('Running DEseq on data:')
  print(head(counts))
  
  #set up input matrix for DESeq
  ddsHTSeq<-DESeqDataSetFromMatrix(counts,
                                   colData = coldata,
                                   design = formula(~ tissue + genotype + fraction))
  #run rld on all the samples
  vsd = rlog(ddsHTSeq,
             fitType='local')
  vsdf = data.frame(assay(vsd))
  
  #split into captured and unbound matrices
  mvsd = vsdf %>% 
    dplyr::select(grep('^m', colnames(vsdf))) %>% 
    as.matrix()
  uvsd = vsdf %>% 
    dplyr::select(grep('^ub', colnames(vsdf))) %>% 
    as.matrix()
  #take their difference
  ddat = mvsd - uvsd
  s0 = colnames(ddat)
  s1 = sub('m.', '', s0, fixed=TRUE)
  samples = sapply(s1, function(x) strsplit(x, split='_')[[1]][1])
  colnames(ddat) = samples
  print('differnce worked?')
  print(sum(mvsd[,'m.N12.s1_S5'] - uvsd[,'ub.N12.s1_S13'] == ddat[,'N12.s1'])==nrow(ddat))
  return(ddat)
}
