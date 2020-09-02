#repeat_bed_from_gff.R
#one-off script to build a bed file of the repeats from Amil.repeat.gff, part of the Amil.v2.00.chrs reference genome by Zach Fuller

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
dat = read_tsv('Amil.repeat.gff',
               col_names = c('chr', 'source', 'x', 'start', 'end', 'xx', 'dir', 'xxx', 'description'))
head(dat)

#split out the repeat types
description = dat$description
genus1 = sapply(split1, function(x) strsplit(x, 'genus:', fixed=TRUE)[[1]][2])
unique(genus1) #lots of variants of just a few types, split again on space
genus2 = sapply(genus1, function(x) strsplit(x, ' ', fixed=TRUE)[[1]][1])
unique(genus2) #LINE, DNA, LRT, have some further descriptions split by % sign
genus3 = sapply(genus2, function(x) strsplit(x, '%', fixed=TRUE)[[1]][1])
unique(genus3)
genus = sapply(genus3, function(x) strsplit(x, ';', fixed=TRUE)[[1]][1])
table(genus) #RC is Rolling cirlce (helitron)

#plot frequencies of repeat types
data.frame(table(genus)) %>% 
  arrange(desc(Freq)) %>% 
  mutate(genus = sub('_', ' ', genus),
  mutate(genus = factor(genus, levels = genus)) %>% 
  ggplot(aes(x=genus, y=Freq)) +
  geom_bar(stat='identity') +
  labs(x='repeat type',
       y='Frequency')



#build the bed file
bdat = data.frame('chr' = dat$chr,
                  'start' = dat$start,
                  'end' = dat$end,
                  'type' = as.character(genus))
head(bdat)


#write out a bed file for each type
unique(bdat$type)
types = na.omit(unique(bdat$type))
sub_rep = function(t){
  sub = bdat %>% 
    filter(type == t) %>% 
    arrange(chr, start)
  sub$type = paste(sub$type, seq(1,nrow(sub),1), sep='_')
  sub %>% 
    write_tsv(path = paste(sub(' ', '_', t), 'repeats.bed', sep='_'))
}
for (t in types){
  sub_rep(t)
}

