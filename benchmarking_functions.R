#benchmarking_funcitons.R
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(boot)

mdRAD_label='mdRAD'


#function to set up same color scheme as ggplot2 default
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#COLOR SCHEME
pm.color = gg_color_hue(6)[1]
mbd.color = gg_color_hue(6)[2]
mr.color = gg_color_hue(6)[3]
alt2.color = gg_color_hue(6)[4]
plus2.color = gg_color_hue(6)[5]
all3.color = gg_color_hue(6)[6]
color.list = c(pm.color,
               mbd.color,
               mr.color,
               alt2.color,
               plus2.color,
               all3.color)
a.color.set=c(pm.color,
              mbd.color,
              mr.color)
a.shape.set=c(19,17,15)


#get inverse logit labels for axes
inv_logit_labs = function(x) {
  round(inv.logit(x), 2)
}

#get inverse log2 labels for axes
log2_to_percent = function(x) {
  round(2^x, 2)*100
}

#function to make a log 2 conversion using smallest non-zero value for zeros
log2_convesion_noZero = function(x){
  nonZero=min(x[x>0], na.rm=TRUE)
  x[x==0]<-nonZero
  print(paste('replacing zeros with next smallest value =', nonZero))
  return(log(x,2))
}

#convert zero to next lowest
convesion_noZero = function(x){
  nonZero=min(x[x>0], na.rm=TRUE)
  x[x==0]<-nonZero
  print(paste('replacing zeros with next smallest value =', nonZero))
  return(x)
}


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
  rownames(counts)=
  return(list('counts'=counts,
              'pos'=positions))
}


#read in window stats
read_window_stats = function(fileName){
  wdat = read_tsv(fileName) %>% 
    data.frame()
  rownames(wdat)=paste(paste(paste(wdat$name, wdat$chr, sep='_'), wdat$start, sep='_'), wdat$end, sep='_')
  return(wdat)
}


#function for general volcano
plot_volcano_general = function(df,
                                xcol='log2FoldChange',
                                ycol='pvalue',
                                sigcol='padj',
                                sigcut=0.1,
                                ALPHA=0.5,
                                title='',
                                subtitle='',
                                xlab='difference',
                                ylab=bquote(log[10]*pvalue)){
  naForX = is.na(df[,xcol])
  df=df[!naForX,]
  sigVec = df[,sigcol]
  logP = -log(df[,ycol], 10)
  df=df %>% 
    mutate(sig=sigVec<sigcut,
           sig=if_else(is.na(sig),
                       FALSE,
                       sig),
           sig=factor(sig, levels=c(TRUE,FALSE)),
           logp = logP)
  nTrue = sum(df$sig==TRUE)
  nFalse = sum(df$sig==FALSE)
  sigString = paste('N', nTrue, sep='=')
  nsString = paste('N', nFalse, sep='=')
  df=df %>% 
    mutate(sigCall=if_else(sig==TRUE,
                           sigString,
                           nsString),
           sigCall=factor(sigCall, levels=c(sigString, nsString)))
  plt=df %>% 
    ggplot(aes_string(x=xcol, y='logp', color='sigCall')) +
    geom_point(alpha=ALPHA) +
    theme(legend.position='top',
          legend.title = element_blank(),
          plot.subtitle = element_text(hjust = 0.5),
          plot.title = element_text(hjust = 0.5)) +
    labs(title=title,
         # subtitle=paste('N', sum(!is.na(df[,xcol])), sep='='), #old way to plot the nyumber of genes
         subtitle = subtitle,
         x=xlab,
         y=ylab)
  if (sum(as.logical(df$sig))>0){
    plt=plt+scale_color_manual(values=c('red', 'black'))
  } else {
    plt=plt+scale_color_manual(values=c('black', 'black'))
  }
}


#plot a pca from vsd-like dataframe + coldata
plot_vsd_pca <- function (df, coldat, ntop = 25000, color.col='genotype', shape.col='tissue', returnData = F, pcs = 2, pc1 = 1, pc2 = 2, main = "\n", SIZE = 2, legendTitle=NULL, xInvert=1) 
{
  rv <- rowVars(as.matrix(df))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(pca$x[,1:pcs])
  d$col=coldat[,color.col]
  d$shape=coldat[,shape.col]
  
  attr(d, "percentVar") <- percentVar[1:2]
  g=d %>% 
    ggplot(aes(x=PC1, y=PC2, color=col, shape=shape)) +
    geom_point(size=5) +
    theme(axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank())
  
  # g = g + theme_bw()
  print(g)
  if (returnData == T){
    return(d)
  }
  else{
    return(g)
  }
}


plot_adonis_from_vsd = function(df, coldat, color.col='genotype', shape.col='tissue'){
  color.var = coldat[,color.col]
  shape.var = coldat[,shape.col]
  ad=adonis(t(df)~color.var+shape.var)
  labs=c("genotype","tissue","residuals")
  aovTab = ad$aov.tab
  r2s = aovTab$R2[1:3]
  ps = aovTab$`Pr(>F)`[1:2]
  stars = stars.pval(ps)
  labs[1] = paste(labs[1], stars[1], sep='')
  labs[2] = paste(labs[2], stars[2], sep='')
  cols=append(gg_color_hue(2), 'grey')
  # pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="WGBS")
  pidat = data.frame(r2 = r2s/sum(r2s)*100,
                     p = append(ps, NA),
                     lab=factor(labs, levels=labs))
  #plot with ggplot
  pidat %>% 
    ggplot() +
    geom_bar(aes(fill=lab, y=r2, x=''), stat='identity', color='black') +
    coord_polar("y", start=.5) +
    scale_fill_manual(values=cols) +
    theme(axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
}




#function for scatterplot with R2 annotation
plot_scatter_r2_annotated = function(dat, xcol, ycol, xlab, ylab, ALPHA=0.1){
  dat=data.frame(dat)
  bad = c(NA, NaN, Inf, -Inf)
  x=dat[,xcol]
  y=dat[,ycol]
  rem = (x %in% bad) | (y %in% bad)
  totBad = sum(rem)
  if (totBad>0){
    print(paste('removing', totBad, 'rows with NA/NaN/Inf'))
  }
  dat=dat[!rem,]
  lm1=lm(dat[,ycol]~dat[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  print(summary(lm1))
  plt = dat %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha = ALPHA) +
    labs(x=xlab,
         y=ylab)
  pbuild = ggplot_build(plt)
  yrange = pbuild$layout$panel_params[[1]]$y.range
  xrange = pbuild$layout$panel_params[[1]]$x.range
  plt +
    annotate("text", x = xrange[1], y = yrange[2],
             label = paste('italic(R) ^ 2 ==', r2), parse=TRUE, color='black',
             hjust=0)
}


get_pearson_cor = function(dat, xcol, ycol){
  dat=data.frame(dat)
  bad = c(NA, NaN, Inf, -Inf)
  x=dat[,xcol]
  y=dat[,ycol]
  rem = (x %in% bad) | (y %in% bad)
  totBad = sum(rem)
  if (totBad>0){
    print(paste('removing', totBad, 'rows with NA/NaN/Inf'))
  }
  dat=dat[!rem,]
  lm1=lm(dat[,ycol]~dat[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  pearsonCor=cor(x=dat[,xcol],
                 y=dat[,ycol])
  r=round(pearsonCor, digits=2)
  return(r)
}


#same as r2 version but for pearson correlation 
plot_scatter_pearsonCor_annotated = function(dat, xcol, ycol, xlab, ylab, ALPHA=0.1){
  dat=data.frame(dat)
  bad = c(NA, NaN, Inf, -Inf)
  x=dat[,xcol]
  y=dat[,ycol]
  rem = (x %in% bad) | (y %in% bad)
  totBad = sum(rem)
  if (totBad>0){
    print(paste('removing', totBad, 'rows with NA/NaN/Inf'))
  }
  dat=dat[!rem,]
  lm1=lm(dat[,ycol]~dat[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  pearsonCor=cor(x=dat[,xcol],
                 y=dat[,ycol])
  r=round(pearsonCor, digits=2)
  print(summary(lm1))
  plt = dat %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha = ALPHA) +
    labs(x=xlab,
         y=ylab)
  pbuild = ggplot_build(plt)
  yrange = pbuild$layout$panel_params[[1]]$y.range
  xrange = pbuild$layout$panel_params[[1]]$x.range
  plt +
    annotate("text", x = xrange[1], y = yrange[2],
             label = paste('italic(r) ==', r), parse=TRUE, color='black',
             hjust=0)
}


#Funciton to plot pairs of scatterplots with histograms on diagonal
plot_pairs_hist_diagonal = function(cols, titles, subtitles, incLower=FALSE){
  pltList=list()
  for (i in 1:length(cols)){
    for (j in 1:length(cols)){
      xcol=cols[i]
      ycol=cols[j]
      xname = titles[i]
      yname = titles[j]
      print(xcol)
      print(ycol)
      plt = plot_scatter_pearsonCor_annotated(aldat,
                                              xcol=xcol,
                                              ycol=ycol,
                                              xlab=xname,
                                              ylab=yname,
                                              ALPHA=ALPHA) + 
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())
      pltName = paste(xcol,ycol,sep='_')
      pltList[[pltName]]=plt
    }
  }
  
  #replace diagonal with histograms
  histList = list()
  for (i in 1:length(cols)){
    col=cols[i]
    title=titles[i]
    subtitle = subtitles[[i]]
    Ngenes = sum(!is.na(aldat[,col]))
    hist = aldat %>% 
      ggplot(aes_string(x=col)) +
      geom_histogram() +
      labs(title=title, x=title, subtitle=subtitle) + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.subtitle = element_text(hjust = 0.5))
    histList[[col]]=hist
  }
  
  for (col in cols){
    selfName = paste(col,col,sep='_')
    pltList[[selfName]] = histList[[col]]
  }
  #optionally remove plots below diagal
  if (incLower == FALSE){
    print('Removing plots from below diagonal...')
    done=c()
    keepPairs = c()
    for (i in cols){
      remaining = cols[!cols %in% done]
      for(j in remaining){
        keepPairs = append(keepPairs, paste(i,j,sep='_'))
      }
      done=append(done, i)
    }
    badPairs = names(pltList)[!names(pltList) %in% keepPairs]
    print('pairs to remove:')
    for (bp in badPairs){
      print(bp)
      pltList[[bp]]=ggdraw()
    }
  }
  return(pltList)
}





###reappend positions to deseq res (assumes posList)
reappend_positions_to_res = function(resList, names){
  checks = c()
  print('checking matchup...')
  for (n in names){
    pos=posList[[n]]
    tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
    match=sum(rownames(resList[[n]])==tag)==nrow(resList[[n]])
    print(paste(n,match,sep=' = '))
    checks=append(checks, match)
  }
  if (sum(checks)!=length(checks)){
    print('ERROR! Check agreements between positions!')
    break
  } else {
    print('Good. All matched up.')
  }

  newList = list()
  for (n in names){
    print(n)
    newList[[n]] = cbind(posList[[n]],
                         data.frame(resList[[n]]))
  }
  return(newList)
}

#reappend positions to deseq res (assumes posList)
reappend_positions_to_res2 = function(resList, names){
  checks = c()
  print('checking matchup...')
  for (n in names){
    pos=posList[[n]]
    tag =pos$name
    match=sum(rownames(resList[[n]])==tag)==nrow(resList[[n]])
    print(paste(n,match,sep=' = '))
    checks=append(checks, match)
  }
  if (sum(checks)!=length(checks)){
    print('ERROR! Check agreements between positions!')
    break
  } else {
    print('Good. All matched up.')
  }

  newList = list()
  for (n in names){
    print(n)
    newList[[n]] = cbind(posList[[n]],
                         data.frame(resList[[n]]))
  }
  return(newList)
}


#For adding smaller window sizes that were too big to run directly
#see mbdseq_differences_TACC.R in MBD-seq and methylRAD walkthroughs
upload_response_to_add = function(filePath){
  df = read.table(filePath, sep='\t', header = TRUE)
  df$name = sub('-', '_', df$name)
  df$name = sub('-', '_', df$name)
  tag0 = paste(df$name, df$chr, sep='_')
  coords = paste(df$start, df$end, sep='_')
  tag = paste(tag0, coords, sep='_')
  rownames(df) = tag
  return(df)
}


#used in compare_resonses.R to merge together dataframes from different assays
merge_from_assay_lists = function(names, assayList){
  newList = list()
  for (n in names){
    print(n)
    nList = lapply(assayList, function(x) return(x[[n]]))
    for (i in 1:length(nList)){
      adat = nList[[i]]
      if (is.null(adat)){
        print(paste('Warning. One of the assays (index =', i, ') lacked data for region type =', n))
        nList[[i]]<-NULL
      }
    }
    # nList2 = lapply(nList, function(x) change_stop_to_end(x)) #this was a temporary fix no longer needed
    ndf = purrr::reduce(nList, full_join, by=c('chr', 'start', 'end', 'name'))
    # colnames(ndf)[colnames(ndf)=='log2FoldChange']<-'mbd.score'
    newList[[n]] = ndf
  }
  return(newList)
}



plot_scatter = function(df, xcol, ycol){
  bad = c(NA, NaN, Inf, -Inf)
  badx = df[,xcol] %in% bad
  bady = df[,ycol] %in% bad
  badr = badx | bady
  subdf = df[,c(ycol, xcol)] %>% 
    filter(!badr)
  lm1 = lm(subdf[,ycol] ~ subdf[,xcol])
  r2 = round(summary(lm1)$r.squared, digits=2)
  df %>% 
    ggplot(aes_string(x=xcol, y=ycol)) +
    geom_point(alpha=ALPHA) +
    labs(subtitle=bquote(R^2*' = '*.(r2)))
}


my_pairs = function(df){
  print(colnames(df))
  pltList = list()
  for (xcol in colnames(df)){
    for (ycol in colnames(df)){
      pair=paste(xcol, ycol, sep='-')
      print(pair)
      if (xcol==ycol){
        plt= ggdraw() + draw_label(xcol)
      } else {
        plt = plot_scatter(df, xcol, ycol)
      }
      pltList[[pair]]=plt
    }
  }
  return(pltList)
}





















#----------- old
compare_with_picomethyl = function(pmdat, ydat, pmCol, ycol, xlab='Fractional Methylation', ylab='', ALPHA=0.1){
  mdf = pmdat %>% 
    left_join(ydat, by = 'gene') %>% 
    data.frame(stringsAsFactors=FALSE)
  x=mdf[,pmCol]
  y=mdf[,ycol]
  lm1=lm(y~x)
  print(summary(lm1))
  r2=summary(lm1)$r.square
  minpm = min(x)
  miny = min(y)
  mdf %>% 
    ggplot(aes_string(x=pmCol, y=ycol)) +
    geom_point(alpha=ALPHA) +
    labs(y=ylab, x=xlab) +
    scale_x_continuous(labels = inv_logit_labs)
}

plot_shared_x = function(plotList, xlab){
  revList = lapply(plotList, function(x) return(x + theme(axis.title.x = element_blank())))
  top = plot_grid(plotlist=revList, nrow=1)
  xlab = ggdraw() + draw_label(xlab)
  plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.1))
}


gbm_scatter = function(dat1, dat2, col1, col2, subtitle, xlab, ylab){
  d=dat1 %>% 
    inner_join(dat2, by = 'gene') %>% 
    data.frame()
  x=d[,col1]
  y=d[,col2]
  d=d[!is.na(x),]
  d=d[!is.na(y),]
  d=d[x != Inf & y != Inf,]
  d=d[x != -Inf & y != -Inf,]
  x=d[,col1]
  y=d[,col2]
  lm1=lm(y~x)
  print(summary(lm1))
  r2=round(summary(lm1)$r.square, digits=2)
  xadj = (max(x, na.rm=TRUE) - min(x, na.rm=TRUE) ) / 8
  d %>% 
    ggplot(aes_string(x=col1, y=col2)) +
    # geom_smooth(method='lm', lwd=0.5) +
    geom_point(alpha=ALPHA) +
    labs(subtitle=subtitle, x=xlab, y=ylab) +
    annotate("text", x = min(x, na.rm=TRUE)+xadj, y = max(y, na.rm=TRUE),
             label = paste('italic(R) ^ 2 ==', r2), parse=TRUE, color='black')
}



format_deseq_res = function(res){
  data.frame(res) %>% 
    mutate(gene=row.names(res),
           absDiff=abs(log2FoldChange)) %>% 
    as_tibble()
}


#add the gene column to a dataset with only transcript annotations
#From picoMethyl_data_processing_pipeline.txt:
#get a gff transcript-to-gene table
# parent_table_from_gff.R --gff $GFF --feature $PARENT_FEATURE --idTag $PARENT_ID --parentTag $BLOCK_ID --o mRNA_to_gene.txt
transcript_to_gene = function(df, transcriptColumn){
  load('metadata/mRNA_to_gene.Rdata')
  colnames(df)[colnames(df)==transcriptColumn]<-'transcript'
  newdf=df %>% 
    left_join(mtg, by = 'transcript')
  order=c('gene', colnames(newdf)[colnames(newdf)!='gene'])
  return(newdf[,order])
}




# FUNCTIONS FOR READ REDUCTION COMPARISONS --------------------------------

#funciton to clear away things that mess up cor
prep_for_cor = function(vector){
  bad=c(-Inf, Inf, NaN)
  vector[vector %in% bad]<-NA
  return(vector)
}

#function to get correlations for set of y columns
get_correlations_by_column = function(dat, xcol, ycols){
  corVec = c()
  nCompleteObsVec = c()
  xVec =  prep_for_cor(dat[,xcol])
  for (yc in ycols){
    yVec = prep_for_cor(dat[,yc])
    xyCor = cor(x=xVec, y=yVec, use='complete.obs')
    indat = data.frame(xVec, yVec)
    NcompleteObs = nrow(na.omit(indat))
    nCompleteObsVec = append(nCompleteObsVec, NcompleteObs)
    corVec = append(corVec, xyCor)
  }
  res = data.frame(y=ycols,
                   x=xcol,
                   yxCor = corVec,
                   NcompleteObs=nCompleteObsVec)
  return(res)
}

swap_names= function(df, new.names){
  df$x=as.character(df$x)
  xs = unique(df$x)
  for (x in xs){
    df$x[df$x==x]<-new.names[x]
  }
  df$x=factor(df$x, levels=new.names)
  return(df)
}


#funciton to remove X axes from plots in a list
remove_x = function(pltList){
  npltList=pltList
  for (i in 1:length(pltList)){
    npltList[[i]]=pltList[[i]] + theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.ticks.x = element_blank())
  }
  return(npltList)
}


#funciton to remove Y axes from all but first of plots in a list
remove_y = function(pltList){
  npltList=pltList
  for (i in 2:length(pltList)){
    npltList[[i]]=pltList[[i]] + theme(axis.title.y = element_blank(),
                                       axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank())
  }
  return(npltList)
}


