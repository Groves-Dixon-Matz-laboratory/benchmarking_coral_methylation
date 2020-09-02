#glmLvl_Ngene_test.R
rm(list=ls())

library(tidyverse)
library(cowplot)


# SOME FUNCTIONS ----------------------------------------------------------

#funciton to clear away things that mess up cor
prep_for_cor = function(vector){
  bad=c(-Inf, Inf, NaN)
  vector[vector %in% bad]<-NA
  return(vector)
}

#function to get correlations for set of y columns
get_correlations_by_column = function(dat, xcol, ycols){
  corVec = c()
  xVec =  prep_for_cor(dat[,xcol])
  for (yc in ycols){
    yVec = prep_for_cor(dat[,yc])
    xyCor = cor(x=xVec, y=yVec, use='complete.obs')
    corVec = append(corVec, xyCor)
  }
  res = data.frame(y=ycols,
                   x=xcol,
                   yxCor = corVec)
  return(res)
}


#count NAs
count_NAs=function(df){
  sum(is.na(df[,2]))
}


# LOAD DATA ---------------------------------------------------------------

#the full dataset level measures
ll=load('comparisons/datasets/gbmLvl.Rdata')
ll
head(gbm.dat)
main.gbm = gbm.dat %>% 
  select(name, pm.glmLvl, mrB, mbd.score)




# CHECK CORELATIONS FOR GLM N GENE TEST -----------------------------------
source('picomethyl/picomethyl_functions.R')
pathList = list.files('picomethyl/datasets/glmLvlNgeneTest', '.Rdata', full.names = TRUE)
fileList = list.files('picomethyl/datasets/glmLvlNgeneTest', '.Rdata', full.names = FALSE)
nameList0 = sub('gene_', '', fileList)
nameList = sub('_glmLvls.Rdata', '', nameList0)
dfList = list()
for (i in 1:length(fileList)){
  fp = pathList[i]
  n=nameList[i]
  fres = upload_glm(fp, CUT=10)
  x=data.frame(fres)
  x[,n]=x$pm.glmLvl
  x2=x[,c('name', n)]
  dfList[[n]]=x2
}

#count NAs
naCounts=c()
for (n in nameList){
  nas=count_NAs(dfList[[n]])
  naCounts = append(naCounts, nas)
}
data.frame(nameList,
           naCounts)


#assemble
ndat = purrr::reduce(dfList, full_join, by='name') %>% 
  full_join(main.gbm, by='name')
head(ndat)


#get WGBS correlations
pmCors = get_correlations_by_column(dat=ndat,
                                    xcol='pm.glmLvl',
                                    ycols=nameList)

#get MBD correlations
mbdCors = get_correlations_by_column(dat=ndat,
                                     xcol='mbd.score',
                                     ycols=nameList)

#get MR correlations
mrCors = get_correlations_by_column(dat=ndat,
                                    xcol='mrB',
                                    ycols=nameList)

#SO IT DOES NOT MATTER!!!!