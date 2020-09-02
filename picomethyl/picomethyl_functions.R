
#funciton to load results from methylKit difference script
load_methylKit_res = function(methylKitResFile){
  print(paste('Loading file', methylKitResFile))
  ll=load(methylKitResFile)
  wbounds$name = sub('-', '_', wbounds$name)
  wbounds$name = sub('-', '_', wbounds$name)
  return(list('bounds'=wbounds,
              'totCounts'=totCounts,
              'diff'=myDiff))
}


#function to upload and replace extreme values
upload_glm = function(glmFileName, CUT=10){
  ll=load(glmFileName)
  fres$start = as.numeric(fres$start)
  fres$end = as.numeric(fres$end)
  noExtreme = fres %>% 
    filter(glmLvl > -CUT & glmLvl < CUT)
  minNoExtreme = min(noExtreme$glmLvl)
  maxNoExtreme = max(noExtreme$glmLvl)
  fres$glmLvl[fres$glmLvl < -CUT]<-minNoExtreme
  fres$glmLvl[fres$glmLvl > CUT]<-maxNoExtreme
  colnames(fres)[colnames(fres)=='glmLvl']<-'pm.glmLvl'
  return(fres)
}


#for straight fractional methylation
plot_frac_dist = function(name){
  res=lvlList[[name]]
  res %>% 
    ggplot(aes(x=l.fracMeth)) +
    geom_histogram() +
    labs(subtitle=name)
}
