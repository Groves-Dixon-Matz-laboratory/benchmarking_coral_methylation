#conversion_efficiency.R
#check conversion efficiency based on lambda and mitochondrial alignments


library(tidyverse)
library(cowplot)

l=read_tsv('picomethyl/conversion_efficiency/allLambda.cov', col_names=c('chr', 'pos', 'm', 'um', 'fileName'))
m=read_tsv('picomethyl/conversion_efficiency/allMito.cov', col_names=c('chr', 'pos', 'm', 'um', 'fileName'))


#identify the samples with lamda in them
have.lambda = l %>% 
  mutate(total=m+um) %>% 
  group_by(fileName, chr) %>% 
  summarize(tot=sum(total)) %>% 
  filter(tot > 200) %>% 
  pull(fileName)

#subset for those
l2 = l %>% 
  filter(fileName %in% have.lambda)


#overall
1-sum(l2$m) / sum(l2$m + l2$um)
1-sum(m$m) / sum(m$m + m$um)


#boxplot by same
#note that one job got so little lambda DNA that we didn't get good estimates, those are the three low efficiency ones

do_barplot = function(d, SUBTITLE){
  sample = sapply(d$fileName, function(x) strsplit(x, '.trim_bismark_bt2.bismark.cov', fixed=TRUE)[[1]][1])
  sample0 = sub('pm-', '', sample)
  sample = sapply(sample0, function(x) strsplit(x, '_')[[1]][1])
  d$sample = sample
  g=d %>% 
    group_by(sample) %>% 
    summarize(tm = sum(m),
              t = sum(m+um)) %>% 
    mutate(prop = tm / t,
           eff = 1-prop)
  print('---------')
  g %>% 
    summarize(mnEff=mean(eff)*100,
              stdErr = sd(eff)/sqrt(nrow(g))*100) %>% 
    print()
  print('----------')
  
  g %>% 
    ggplot(aes(x=sample, y=eff)) +
    geom_bar(stat='identity') +
    lims(y=c(0,1)) +
    labs(y='Conversion efficiency', subtitle=SUBTITLE)
}




lplt = do_barplot(l2, 'lambda')
mplt = do_barplot(m, 'mitochondrial')
plot_grid(lplt, mplt, nrow=1)


l %>% 
  mutate(total=m+um) %>% 
  group_by(fileName, chr) %>% 
  summarize(tot=sum(total))
