
library(cowplot)
theme_set(theme_cowplot())

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
plot_scatter_pearsonCor_annotated = function(dat, xcol, ycol, xlab, ylab, ALPHA=0.1, ylim=FALSE){
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
  if (length(ylim)>1){
    print('setting y lmiits')
    plt=plt + lims(y=ylim) +
      annotate("text", x = xrange[1], y = ylim[2],
               label = paste('italic(r) ==', r), parse=TRUE, color='black',
               hjust=0)
  } else {
    plt=plt +
      annotate("text", x = xrange[1], y = yrange[2],
               label = paste('italic(r) ==', r), parse=TRUE, color='black',
               hjust=0)
  }
  return(plt)
}


#
plot_shared_x = function(plotList, xlab){
revList = lapply(plotList, function(x) return(x + theme(axis.title.x = element_blank())))
top = plot_grid(plotlist=revList, nrow=1)
xlab = ggdraw() + draw_label(xlab)
plot_grid(top, xlab, nrow=2, rel_heights=c(1, 0.1))
}









#
plot_shared_x_y = function(plotList, xlab, ylab, relXlab=1/15, relYlab=1/15){
  revList = lapply(plotList, function(x) return(x + theme(axis.title = element_blank())))
  top = plot_grid(plotlist=revList, nrow=1)
  pxlab = ggdraw() + draw_label(xlab)
  pylab = ggdraw() + draw_label(ylab, angle=90)
  wX = plot_grid(top, pxlab, nrow=2, rel_heights=c(1, relXlab))
  wXY = plot_grid(pylab, wX, nrow=1, rel_widths=c(relYlab,1))
  return(wXY)
}

#get inverse log2 labels for axes
log2_to_percent = function(x) {
  round(2^x, 2)*100
}

#get inverse log2 labels for axes
log10_to_percent = function(x) {
  round(10^x, 6)*100
}

