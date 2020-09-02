#read_reduction_functions.R

#function to get the correlation reduction data
get_correlation_reduction = function(dat, xCols, sampleSums, useMedian=FALSE){
  pct.cols = colnames(dat)[grep('pct.', colnames(dat))]
  lfc.cols = pct.cols[grep('log2FoldChange', pct.cols)]
  pval.cols = pct.cols[grep('padj', pct.cols)]
  
  rList = list()
  for (xcol in xCols){
    cors = get_correlations_by_column(dat=dat,
                                      xcol=xcol,
                                      ycols=lfc.cols)
    rList[[xcol]]=cors
  }
  res = purrr::reduce(rList, rbind)
  #assemble
  if (useMedian==TRUE){
    print('using median reads per sample')
    medReadCount = median(sampleSums)
  } else {
    print('using total reads')
    medReadCount = sum(sampleSums)
  }
  res2 = res %>% 
    mutate(pct0 = sub('pct.', '', y),
           pct = as.numeric(sub('_log2FoldChange', '', pct0)),
           prop=pct/100,
           medReads = (prop*medReadCount)/1e6,
           logProp = log(prop, 2),
           x=factor(x, levels=xCols))
}


# #get N tested
# get_n_tested = function(dat, sampleSums, useMedian=FALSE){
#   p.cols = colnames(dat)[grep('padj', colnames(dat))]
#   nres = c()
#   for(pc in p.cols){
#     nTested = sum(!is.na(dat[,pc]))
#     nres=append(nres, nTested)
#   }
#   res0 = data.frame(pct=p.cols,
#                     nTested = nres)
#   #assemble
#   if (useMedian==TRUE){
#     print('using median reads per sample')
#     medReadCount = median(sampleSums)
#   } else {
#     print('using total reads')
#     medReadCount = sum(sampleSums)
#   }
#   
#   res = res0 %>% 
#     mutate(pct0 = sub('pct.', '', pct),
#            pct = as.numeric(sub('_padj', '', pct0)),
#            prop=pct/100,
#            medReads = prop*medReadCount,
#            logProp = log(prop, 2))
#   return(res)
# }

#get sensitivity, precision, and uncorroborated proporitons
get_significance_correspondance = function(dat,
                                           sigSets,
                                           sampleSums,
                                           useMedian=FALSE,
                                           alternate.pairings){
  pval.cols = colnames(dat)[grep('padj', colnames(dat))]
  pctString = sub('_padj', '', pval.cols[grep('pct', pval.cols)])
  res = data.frame()
  for (pcol in pval.cols){
    rsig = dat[,pcol] < PCUT
    rsig[is.na(rsig)]<-FALSE
    rgenes = dat[rsig,'name']
    #sensitivity (proportion of 1s in sigSet identified as 1s in reduced set)
    sensitivity = sapply(sigSets, function(x) return(sum(rgenes %in% x) / length(x) ))
    #precision (proportion of predicted 1s that are actually 1s)
    precision = sapply(sigSets, function(x) return(sum(rgenes %in% x) / length(rgenes) ))
    nTested = sum(!is.na(dat[,pcol]))
    nSig = sum(dat[,pcol]<0.1, na.rm=TRUE)
    sub.res = data.frame(pct=pcol,
                         group=names(sensitivity),
                         sensitivity=sensitivity,
                         precision=precision,
                         nSig=nSig,
                         nTested=nTested)
    res=rbind(res,sub.res)
  }
  #refromat
  if (useMedian==TRUE){
    print('using median reads per sample')
    medReadCount = median(sampleSums)
  } else {
    print('using total reads')
    medReadCount = sum(sampleSums)
  }
  res2 = res %>% 
    mutate(prop0 = sub('pct.', '', pct),
           prop = as.numeric(sub('_padj', '', prop0))/100,
           lprop = log(prop, 2),
           nReads = (prop*medReadCount)/1e6,
    )
  return(res2)
}

#plot wrong, fps, fp rate
gather_fp_data = function(dat, sampleSums, useMedian=TRUE){
  ps0 = dat %>% 
    dplyr::select(contains('wrong')) %>% 
    colnames()
  pctStrings = sub('_wrong', '', ps0)
  res=data.frame()
  if (useMedian==TRUE){
    # print('using median reads per sample')
    medReadCount = median(sampleSums)
  } else {
    # print('using total reads')
    medReadCount = sum(sampleSums)
  }
  for (ps in pctStrings){
    pct = as.numeric(sub('pct.', '', ps))
    pcol = paste(ps, 'padj', sep='_')
    wcol = paste(ps, 'wrong', sep='_')
    sig = dat[,pcol] < PCUT
    wrong = dat[,wcol]
    fp = sig & wrong
    nRegions = length(sig)
    nNoNa = sum(!is.na(sig))
    nSig = sum(sig, na.rm=TRUE)
    nWrong = sum(wrong, na.rm=TRUE)
    nFp = sum(fp, na.rm=TRUE)
    sub.res = data.frame(pct=ps,
                         prop=pct/100,
                         nRegions = nRegions,
                         nTested = nNoNa,
                         nSig=nSig,
                         nWrong=nWrong,
                         nFp=nFp,
                         nReads=pct/100*medReadCount/1e6)
    res=rbind(res, sub.res)
  }
  return(res)
}



# USED FOR FALSE POSITIVE INFERENCE ---------------------------------------


# get_resid_df = function(pctLfc) {
#   pctSub = pcdat %>% 
#     filter(pct==pctLfc) %>% 
#     na.omit()
#   resid = as.numeric(lm(lfc ~ PC1, data = pctSub)$resid)
#   pctSub$resid = resid
#   pctSub$wrong = abs(pctSub$resid) > RCUT
#   return(pctSub)
# }
# 
# get_wrong_names = function(pctString){
#   pctLfc = paste(pctString, 'log2FoldChange', sep='_')
#   rdat = get_resid_df(pctLfc)
#   rdat %>% 
#     filter(wrong) %>% 
#     pull(name)
# }
# 
# 
# plot_pct_pc_scatter = function(pctString){
#   pctLfc = paste(pctString, 'log2FoldChange', sep='_')
#   rdat = get_resid_df(pctLfc)
#   nWrong = sum(rdat$wrong)
#   subtitle = paste('n=', nWrong, sep='')
#   rdat %>% 
#     ggplot(aes(x=PC1, y=lfc, color=wrong)) +
#     geom_point() +
#     scale_color_manual(values=c('black', 'red')) +
#     labs(title=pctString, subtitle=subtitle) +
#     theme(legend.position='none')
# }
# 
# 
# get_pcdat = function(dat){
#   lfcDat = dat %>% 
#     dplyr::select(name, contains('log2FoldChange')) %>% 
#     dplyr::select(name, contains('pct'))
#   lfcScaled = lfcDat %>% 
#     dplyr::select(-name) %>% 
#     map_dfc( ~ (.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)) %>% 
#     mutate(name=lfcDat$name)
#   pcdat = dat %>% 
#     dplyr::select(name, PC1) %>% 
#     left_join(lfcScaled, by = 'name') %>% 
#     pivot_longer(contains('pct'), names_to = 'pct', values_to = 'lfc')
#   return(pcdat)
# }

get_pctStrings = function(dat){
  ps0=dat %>% 
    dplyr::select(contains('pct')) %>% 
    dplyr::select(contains('log2Fold')) %>% 
    colnames()
  pctStrings = sub('_log2FoldChange', '', ps0)
  return(pctStrings)
}

add_wrong_to_df = function(dat){
  pctStrings = get_pctStrings(dat)
  #revise the original input and return
  for (pctString in pctStrings){
    wrongCol = paste(pctString, 'wrong', sep='_')
    lfc = dat[,paste(pctString, 'log2FoldChange', sep='_')]
    scaled = (lfc - mean(lfc, na.rm=TRUE) ) / sd(lfc, na.rm=TRUE)
    lmdat = data.frame(name=dat$name,
                       scaled=scaled,
                       PC1=dat$PC1) %>% 
      na.omit()
    #some reductions from picomethyl may not have not have any genes, skip these
    if (nrow(lmdat)>10){
      resid = as.numeric(lm(scaled ~ PC1, data = lmdat)$resid)
      wrongSet = lmdat$name[abs(resid) > RCUT]
      dat[,wrongCol] = dat$name %in% wrongSet
    } else{
      next
    }
  }
  return(dat)
}




# FOR PLOTTING ------------------------------------------------------------


#funciton to plot correlation reduction
plot_correlation_reduction=function(cor.df, TITLE, SUBTITLE, a.color.set, exclude=c()){
  cor.df %>% 
    filter(NcompleteObs > minToPlot,
           !x %in% exclude) %>% 
    ggplot(aes(x=medReads, y=yxCor,color=x, shape=x)) +
    geom_point() +
    geom_line() +
    labs(x='median reads per sample',
         y='difference correlation',
         title=TITLE,
         subtitle = SUBTITLE) +
    theme(plot.subtitle = element_text(hjust = 0.5)) +
    scale_color_manual(values=a.color.set)
}



#function to plot sensitivity line plot
plot_sensitivity = function(rdf, assay.to.rm, color.list, shape.list){
  mod.color.list = color.list
  mod.color = color.list[names(color.list) != assay.to.rm]
  mod.shape = shape.list[names(shape.list) != assay.to.rm]
  mod.levels=names(mod.shape)
  rdf %>% 
    mutate(group=as.character(group),
           group = if_else(group %in% alternate.pairings,
                           alt.string,
                           group),
           group=factor(group,
                        levels=mod.levels)) %>% 
    ggplot(aes(x=nReads, y=sensitivity, color=group, shape=group)) +
    geom_point() +
    geom_line() + 
    labs(x='median read count',
         y='sensitivity',
         color='',
         shape='') +
    scale_color_manual(values=mod.color) +
    scale_shape_manual(values=mod.shape)
}

#function to plot sensitivity line plot 
plot_sig_agreement = function(rdf, col.to.plot, alt.pair, color.list, shape.list, groups_to_exclude=c()){
  NCUT=10
  mod.color = color.list[!names(color.list) %in% groups_to_exclude]
  mod.shape = shape.list[!names(shape.list) %in% groups_to_exclude]
  mod.levels=names(mod.shape)
  mod.levels=names(shape.list)
  print('Using these levels:')
  print(mod.levels)
  rdf %>% 
    filter(nTested > NCUT) %>% 
    #first swap the alternative pair group to 'alt-pair'
    mutate(group=as.character(group),
           group = if_else(group == alt.pair,
                           alt.string,
                           group)) %>% 
    #now remove the other alternate.pairings from data
    filter(!group %in% alternate.pairings,
           ! group %in% groups_to_exclude) %>% 
    mutate(group=factor(group,
                        levels=mod.levels)) %>% 
    ggplot(aes_string(x='nReads', y=col.to.plot, color='group', shape='group')) +
    geom_point() +
    geom_line() + 
    labs(x='median read count',
         y='precision',
         color='',
         shape='') +
    scale_color_manual(values=color.list) +
    scale_shape_manual(values=shape.list)
}


#function to plot the precision of an assay
#input for this is created by get_significance_correspondance
#precision here is the proportion of significant genes that were
#significant in at least one other assay
plot_precision = function(rdf, COL='firebrick', SHAPE=7, NCUT=10){
  rdf %>% 
    filter(group %in% c('any 2'),
           nTested > NCUT) %>% 
    mutate(uncorroborated = 1-precision) %>% 
    ggplot(aes(x=nReads, y=precision)) +
    geom_point(color=COL, shape=SHAPE) +
    geom_line(color=COL) + 
    labs(x='median read count',
         y='precision',
         color='',
         shape='')
}


#function to plot the precision of an assay against each other
#input for this is created by get_significance_correspondance
#precision here is the proportion of significant genes that were
#significant in at least one other assay
plot_precision_multi = function(rdf, alt.pair, color.list, shape.list, exclude=c()){
  NCUT=10
  mod.levels=names(shape.list)
  print('Using these levels:')
  print(mod.levels)
  rdf %>% 
    filter(nTested > NCUT) %>% 
    #first swap the alternative pair group to 'alt-pair'
    mutate(group=as.character(group),
           group = if_else(group == alt.pair,
                           alt.string,
                           group)) %>% 
    #now remove the other alternate.pairings from data
    filter(!group %in% alternate.pairings) %>% 
    mutate(group=factor(group,
                        levels=mod.levels)) %>% 
    mutate(uncorroborated = 1-precision) %>% 
    ggplot(aes(x=nReads, y=precision, color=group, shape=group)) +
    geom_point() +
    geom_line() + 
    labs(x='median read count',
         y='precision',
         color='',
         shape='') +
    scale_color_manual(values=color.list) +
    scale_shape_manual(values=shape.list)
}





#old function to plot false positive rate estimate
old_plot_false_pos = function(fdf){
  nSigCut = 0 #there need to be at least this many to plot the point
  legendTitle = 'inferred\nerror'
  fdf %>% 
    mutate(`deviant` = ifelse(nSig > nSigCut,
                              nWrong / nRegions,
                              NA),
           `false pos.` = ifelse(nSig > nSigCut,
                                     nFp / nRegions,
                                     NA)) %>% 
    dplyr::select(nReads, `deviant`, `false pos.`) %>% 
    pivot_longer(`deviant`:`false pos.`, names_to = 'fpType', values_to = 'value') %>% 
    mutate(fpType=factor(fpType, levels=c('deviant',
                                          'false pos.'))) %>% 
    ggplot(aes(x=nReads, y=value, color=fpType)) +
    geom_point(aes(shape=fpType)) +
    geom_line(aes(lty=fpType)) +
    labs(x='median read count',
         y='false positive rate',
         lty=legendTitle,
         shape=legendTitle,
         color=legendTitle) +
    scale_color_manual(values = c('firebrick', 'firebrick1'))
}

#fuction to get false positives 
plot_false_pos = function(fdf){
  nSigCut = 0 #there need to be at least this many to plot the point
  legendTitle = 'inferred\nerror'
  fdf %>% 
    mutate(`deviant` = ifelse(nSig > nSigCut,
                              nWrong / nRegions,
                              NA),
           `false pos.` = ifelse(nSig > nSigCut,
                                 nFp / nRegions,
                                 NA)) %>% 
    dplyr::select(nReads, `deviant`, `false pos.`) %>% 
    pivot_longer(`deviant`:`false pos.`, names_to = 'fpType', values_to = 'value') %>% 
    mutate(fpType=factor(fpType, levels=c('deviant',
                                          'false pos.'))) %>% 
    ggplot(aes(x=nReads, y=value, color=fpType)) +
    geom_point(aes(shape=fpType)) +
    geom_line(aes(lty=fpType)) +
    labs(x='median read count',
         y='false positive rate',
         lty=legendTitle,
         shape=legendTitle,
         color=legendTitle) +
    scale_color_manual(values = c('firebrick', 'firebrick1'))
}








plot_n_tested = function(fdf, col, main){
  fdf %>% 
    filter(pct!='pct.100') %>% 
    ggplot(aes(x=nReads, y=nTested)) +
    geom_point(color=col) +
    geom_line(color=col) +
    labs(x='median read count',
         y='regions tested',
         title=main) +
    theme(plot.title = element_text(hjust = 0.5))
}


# 
# #plot inferred false positive rate
# plot_fp_rate = function(rdf, shape.choice=19){
#   rdf %>% 
#     dplyr::select(nReads, dirFps, eitherFps) %>% 
#     filter(!duplicated(nReads)) %>% 
#     pivot_longer(cols=c('dirFps', 'eitherFps'),
#                  names_to = 'falsePosType',
#                  values_to = 'prop') %>% 
#     ggplot(aes(x=nReads, y=prop, color=falsePosType)) +
#     geom_point(shape=shape.choice) +
#     geom_line() + 
#     labs(x='median read count',
#          y='inferred false positive rate') +
#     scale_color_manual(values=c('red', 'firebrick'))
# }
# 
# 
# #plot ratio of corroborated positives to false positives
# plot_tp_to_fp_ratio = function(rdf,
#                                group.filter='any 2',
#                                color.choice='black',
#                                shape.choice=19 ){
#   res3 %>% 
#     filter(group==group.filter,
#            nSig > minSigToPlot) %>% 
#     mutate(ratio = sensitivity / eitherFps,
#            lratio = log(sensitivity / eitherFps), 2) %>% 
#     ggplot(aes(x=nReads, y=ratio)) +
#     geom_point(color=color.choice, shape=shape.choice) +
#     geom_line(color=color.choice) + 
#     labs(x='median read count',
#          y=bquote('corroborated : false positive'))
# }