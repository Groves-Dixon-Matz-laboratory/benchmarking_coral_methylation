

library(tidyverse)

#read in rep files
infiles = read.table('glmRepFiles.txt', stringsAsFactors=FALSE)$V1
fileList = list()
for (infile in infiles){
  ll=load(infile)
  fileList[[infile]]=fres
}

#rbind them and get overall means
dat = purrr::reduce(fileList, rbind)
fres = dat %>% 
  group_by(chr, start, end, name) %>% 
  summarize(meanLvl = mean(mnLvl))

#save
save(fres, file='glmRepLvlMeans.Rdata')