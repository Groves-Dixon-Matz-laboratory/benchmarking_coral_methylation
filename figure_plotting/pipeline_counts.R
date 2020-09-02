#pipeline_counts.R
#get stats for pipeline counts from each of the three assays

library(tidyverse)


# LOOK AT WGBS ------------------------------------------------------------

pc = read_tsv('picomethyl/pipeline_counts/pipeline_counts.txt',
              col_names=c('sample', 'value', 'stat'))

#look at total and median per sample for raw counts, trimmed, and mapped
pc %>% 
  group_by(stat) %>% 
  summarize(total = sum(value)/1e6,
            med = median(value)/1e6,
            N=n())

#look at the mean mapping efficiency
pc %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(mappingEff=mappedCount / trimmedCounts) %>% 
  pull(mappingEff) %>% 
  mean()


# LOOK AT MBDSEQ ----------------------------------------------------------

pc = read_tsv('mbdSeq/pipelineCounts/pipelineCounts.tsv',
              col_names=c('sample', 'value', 'stat'))

#totals
pc %>% 
  group_by(stat) %>% 
  summarize(total = sum(value)/1e6,
            med = median(value)/1e6,
            summed = sum(value),
            N=n())
  
#totals by library type
ldat = pc %>% 
  mutate(libType = substr(sample, start=1, stop=1)) %>% 
  group_by(stat, libType) %>% 
  summarize(total = sum(value)/1e6,
            med = median(value)/1e6,
            summed = sum(value))
ldat

#final percentages
mraw = ldat$summed[ldat$stat=='rawCounts' & ldat$libType=='m']
uraw = ldat$summed[ldat$stat=='rawCounts' & ldat$libType=='u']
mfinal = ldat$summed[ldat$stat=='dedupMapped' & ldat$libType=='m']
ufinal = ldat$summed[ldat$stat=='dedupMapped' & ldat$libType=='u']
mfinal/mraw
ufinal/uraw


# LOOK AT METHYLRAD -------------------------------------------------------

pc = read_tsv('methylRAD/pipelineCounts/pipelineCounts.tsv',
              col_names=c('sample', 'value', 'stat')) %>% 
  mutate(value=if_else(stat %in% c('predupMapped', 'predupPropPaired'),
                       value/2,
                       value))

#grand total
pc %>% 
  filter(!grepl('R2', sample)) %>% 
  group_by(stat) %>% 
  summarize(total = sum(value)/1e6,
            med = median(value)/1e6,
            summed = sum(value))

#percentage with correct starts
csPcts = read.table('methylRAD/pipelineCounts/correct_start_rates.txt')$V1
mean(csPcts)

#percentage PCR duplication
r1DupRates = read.table('methylRAD/pipelineCounts/R1_duplication_rates.txt')$V1
mean(r1DupRates)

#percentage startChecked
passRates = read.table('methylRAD/pipelineCounts/passing_filter_methylGBS_rates.txt')$V1
mean(passRates)

#total reads after filter_methylGBS
pc %>% 
  filter(stat=='startChecked',
         !grepl('R2', sample)) %>% 
  pull(value) %>% 
  sum() /1e6

#totals
pc %>% 
  filter(!grepl('R2', sample)) %>% 
  group_by(stat) %>% 
  summarize(total = sum(value)/1e6,
            med = median(value)/1e6,
            summed = sum(value),
            N=n())


#by enzyme
edat = pc %>% 
  filter(!grepl('R2', sample)) %>% 
  mutate(libType = substr(sample, start=3, stop=3)) %>% 
  group_by(stat, libType) %>% 
  summarize(total = sum(value)/1e6,
            med = median(value)/1e6,
            N=n()) %>% 
  data.frame()
edat



#properly paired mapping efficiency
tF = edat$total[edat$stat=='trimmedCounts' & edat$libType=='F']
tM = edat$total[edat$stat=='trimmedCounts' & edat$libType=='M']
ppF = edat$total[edat$stat=='predupPropPaired' & edat$libType=='F']
ppM = edat$total[edat$stat=='predupPropPaired' & edat$libType=='M']
ppF/tF
ppM/tM


#overall portion of raw reads that get mapped
rF = edat$total[edat$stat=='rawCounts' & edat$libType=='F']
rM = edat$total[edat$stat=='rawCounts' & edat$libType=='M']

ppF/rF
ppM/rM
