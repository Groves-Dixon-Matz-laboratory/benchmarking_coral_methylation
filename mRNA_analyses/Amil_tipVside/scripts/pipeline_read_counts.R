#pipeline_read_counts.R
#plot read counts through data processing steps

#libs
library(tidyverse)
library(cowplot)


#read in
input = 'rnaseq/datasets/all_pipeline_counts.txt'
statLevels = c("rawCounts", "trimmedCounts", "mapped", "geneCounted")
sdat0 = read.table(input, header=F, col.names=c('Run', 'value', 'stat'), stringsAsFactors=F) %>% 
  mutate(stat=factor(stat, levels=statLevels)) %>% 
  as_tibble()

#fix periods inserted by feature counts
length(unique(sdat0$Run))
sdat = sdat0 %>% 
  mutate(Run = sub('.', '-', Run, fixed=T))
length(unique(sdat$Run))


#order the stats in pipline order
order = sdat %>% 
  filter(stat=='rawCounts') %>% 
  arrange(by=value) %>% 
  pull(Run) %>% 
  rev() 
sdat$Run = factor(sdat$Run, levels=order)



#plot full barplot
bp<-sdat %>% 
  ggplot(aes(x=stat, y=value, color=Run, fill=Run)) +
    geom_bar(stat='identity', position='dodge') +
    labs(y='Read count', x='Pipeline step') +
    theme(axis.text.x=element_text(angle=20, vjust=0.75))
bp


#plot gene counted only
gc = sdat %>% 
  filter(stat=='geneCounted') %>% 
  arrange(by=value)
gc$Run = factor(gc$Run, levels = rev(pull(gc, Run)))

histMeans = gc %>% 
  ggplot(aes(x=value)) +
  geom_histogram(bins=6) +
  labs(x='Mean gene counted reads accross projects')

rankedMeans = gc %>% 
  ggplot(aes(x=stat, y=value, color=Run, fill=Run, group=Run)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y='Read count', x='Pipeline step', title='Gene counted')
rankedMeans

plot_grid(histMeans, rankedMeans, nrow=1, rel_widths=c(.5,1))

#plot scatter abs
quartz()
lp<-sdat %>% 
  ggplot(aes(x=stat, y=value, color=Run)) +
  geom_point() +
  geom_line(aes(group=Run)) +
  theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75)) +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')
  

#plot scatter prop raw
pp<-sdat %>% 
  group_by(Run) %>% 
  mutate(prop=value/max(value)) %>% 
  ggplot(aes(x=stat, y=prop, color=Run)) +
    geom_point() +
    geom_line(aes(group=Run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none', axis.text.x=element_text(angle=20, vjust=0.75))

plot_grid(lp, pp)
