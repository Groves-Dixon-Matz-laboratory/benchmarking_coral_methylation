#region_type_comparisons.R
#look at additional regions upon a reviewer's request
rm(list=ls())
source('benchmarking_functions.R')

# load data ---------------------------------------------------------------

#picoimethyl (build with process_picoMethyl.R)
ll=load('picomethyl/datasets/methLevelList.Rdata')
ll
pm_list = lvlList
names(pm_list)

#mbd seq (build with process_MBDseq.R)
ll = load('mbdSeq/datasets/methLevelList.Rdata')
ll
mbd_list = lvlList
names(mbd_list)


#mdRAD (build with process_methylRAD.R)
ll = load('methylRAD/datasets/methLevelList.Rdata')
ll
mdr_list = lvlList
names(mdr_list)


# build up a single dataframe for selected regions -------------------------

#choose regions
select_regions = c('gene',
                   'promoter',
                   'exon',
                   'intron',
                   'five_prime_UTR',
                   'three_prime_UTR',
                   'intergenic',
                   "LINE_repeats",
                   "SINE_repeats",
                   "LTR_repeats",
                   "RC_repeats",
                   "Low_complexity_repeats",
                   "Simple_repeat_repeats")

#loop through and pull from each assay's dataframe list
dat = tibble()
for (r in select_regions){
  print(paste(r,'...',sep=''))
  #assemble levels for pm
  pm = pm_list[[r]]
  pml = pm %>% 
    dplyr::select(chr, start, end, name, l.fracMeth) %>% 
    mutate(region = r,
           assay = 'WGBS') %>% 
    set_names(c('chr', 'start', 'end', 'name', 'lvl', 'region', 'assay'))
  dat = rbind(dat, pml)
  #assemble levels for mbd
  mbd = mbd_list[[r]]
  mbdl = mbd %>% 
    dplyr::select(chr, start, end, name, log2FoldChange) %>% 
    mutate(region = r,
           assay = 'MBD-seq') %>% 
    set_names(c('chr', 'start', 'end', 'name', 'lvl', 'region', 'assay')) %>% 
    as_tibble()
  dat = rbind(dat, mbdl)
  #assemble levels for mdrad
  mdr = mdr_list[[r]]
  mdrl = mdr %>% 
    dplyr::select(chr, start, end, name, mrB) %>% 
    mutate(region = r,
           assay = 'mdRAD') %>% 
    set_names(c('chr', 'start', 'end', 'name', 'lvl', 'region', 'assay')) %>% 
    tibble()
  dat = rbind(dat, mdrl)
}
dat$region = factor(dat$region, levels = select_regions)
dat$assay = factor(dat$assay, levels = c('WGBS', 'MBD-seq', 'mdRAD'))
dat$win_name = paste(dat$chr, dat$start, sep='_')
dat$win_name = paste(dat$win_name, dat$end, sep='_')
dat$win_name = paste(dat$win_name, dat$name, sep='_')

# plot level for each type -------------------------------------------
#plot violin plot for level of eacch region type

#get the means
mns = dat %>% 
  group_by(region, assay) %>% 
  summarize(mn = mean(lvl, na.rm=TRUE))
mns

#plot with violins
clean_region_names = function(input_names){
  form_regions = sub('_repeats', '', input_names)
  form_regions = sub('_repeats', '', form_regions)
  form_regions = sub('_repeat', '', form_regions)
  form_regions = sub('_', ' ', form_regions)
  form_regions = sub('_', ' ', form_regions)
  form_regions[form_regions=='three prime UTR']<-"3' UTR"
  form_regions[form_regions=='five prime UTR']<-"5' UTR"
  return(form_regions)
}
form_regions = clean_region_names(select_regions)
dat$form_region = clean_region_names(dat$region)
dat$form_region = factor(dat$form_region, levels = form_regions)
dat %>% 
  filter( !(assay == 'MBD-seq' & lvl < -5) ) %>% 
  ggplot(aes(x=form_region, y=lvl)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.size=0.2) +
  labs(y='methylation level',
       x = 'region type') +
  facet_wrap(~assay, scales = 'free_y', nrow=3)



# compare within regions --------------------------------------------------
#build scatterplots showing correlation for different region types

ALPHA=0.01

#build scatterplots for each region type
region_scatterplots = list()
for (r in form_regions){
  print('----------------')
  print(paste(r,'...',sep=''))
  pm = dat %>% 
    filter(assay == 'WGBS' & form_region == r)
    
  mbd = dat %>% 
    filter(assay == 'MBD-seq' & form_region == r)
  
  mdr = dat %>% 
    filter(assay == 'mdRAD' & form_region == r)
  
  
  #plot pm and mbd
  pm.mbd = pm %>% 
    left_join(mbd, by = 'name')
  pm.mbd.lvl = plot_scatter_pearsonCor_annotated(pm.mbd,
                                                 xcol='lvl.x',
                                                 ycol='lvl.y',
                                                 xlab='WGBS',
                                                 ylab='MBD-seq',
                                                 ALPHA=ALPHA) + scale_x_continuous(labels = log2_to_percent)
  #plot pm and mdRAD
  pm.mdr = pm %>% 
    left_join(mdr, by = 'name') %>% 
    mutate(lvl.y = if_else(lvl.y < -10,
                         -10,
                         lvl.y),
           lvl.y = if_else(lvl.y > 10,
                         10,
                         lvl.y))
  pm.mdr.lvl = plot_scatter_pearsonCor_annotated(pm.mdr,
                                                 xcol='lvl.x',
                                                 ycol='lvl.y',
                                                 xlab='WGBS',
                                                 ylab='mdRAD',
                                                 ALPHA=ALPHA) + 
    scale_x_continuous(labels = log2_to_percent)
  avg_len = round(mean(pm$end - pm$start), digits=0)
  pans = plot_grid(pm.mbd.lvl, pm.mdr.lvl, nrow=1)
  r_form = sub('_', ' ', r)
  r_form = sub('_', ' ', r_form)
  title1 = ggdraw() + draw_label(r_form, fontface = 'bold')
  title2 = ggdraw() + draw_label(paste('mean length =', avg_len, 'bp'), size=11)
  title = plot_grid(title1, title2, nrow=2)
  plt = plot_grid(title, pans, rel_heights = c(1,10), nrow=2)
  region_scatterplots[[r]] = plt
}


#plot the regions we haven't looked at yet
repeat_regions = c("LINE",
                   "SINE",
                   "LTR",
                   "RC")

plt_list = map(repeat_regions, function(x) region_scatterplots[[x]])
plot_grid(plotlist = plt_list, nrow=2)





# coverage analysis -------------------------------------------------------


#load the bed files
bed_files = list.files(path = './windowStats', pattern = '.bed', full.names = TRUE)
bed_regions = list.files(path = './windowStats', pattern = '.bed', full.names = FALSE)
bed_regions = sub('.bed', '', bed_regions)
bed_regions = sub('_Boundaries', '', bed_regions)
bed_regions = sub('Boundaries', '', bed_regions)
bed_regions[bed_regions=='cds']<-'exon'
names(bed_files) = bed_regions


read_bed = function(x){
  d = read_tsv(x,
               col_names = c('chr', 'start', 'end', 'name'))
}
bed_list = map(bed_files, function(x) read_bed(x))
names(bed_list)

#loop through regions and get the proportion of annotated regions measured for each assay
mdat = tibble()
for (r in bed_regions){
  print(paste(r,'...',sep=''))
  rsub = dat %>% 
    filter(region == r)
  
  r_names = bed_list[[r]]
  u_name = paste(r_names$start, r_names$end, sep='_')
  u_name = paste(u_name, r_names$name, sep='_')
  tot_names = length(unique(u_name))
  
  measured = rsub %>% 
    filter(!is.na(lvl)) %>% 
    group_by(assay) %>% 
    summarize(tot_measured = length(unique(win_name))) %>% 
    mutate(tot_names = tot_names,
           prop = tot_measured / tot_names,
           region = r)
  noninf = rsub %>% 
    filter(!is.na(lvl),
           abs(lvl) < Inf) %>% 
    group_by(assay) %>% 
    summarize(tot_noninf = length(unique(win_name))) %>% 
    mutate(region = r)
  res = measured %>% 
    left_join(noninf, by = c('assay', 'region'))
  mdat = rbind(mdat, res)
}
mdat

#clean up the reigon names for plotting
mdat$form_region = clean_region_names(mdat$region)


#stack into tidy dataframe for total estimated and total measured
#(where we have the total with any number, and total with non-Inf number handing zero read regions for mdRAD)
tm = mdat %>% 
  mutate(prop = tot_measured / tot_names,
         prop_type = 'total estimated') %>% 
  select(assay, form_region, prop_type, prop)
tm
tni = mdat %>% 
  mutate(prop = tot_noninf / tot_names,
         prop_type = 'total measured') %>% 
  select(assay, form_region, prop_type, prop)

#build barplot
rbind(tm, tni) %>% 
  mutate(form_region = factor(form_region, levels = form_regions)) %>% 
  filter(!is.na(assay),
         prop_type == 'total measured') %>% 
  ggplot(aes(x=form_region, fill=assay, y=prop)) +
  geom_bar(position = 'dodge', stat='identity') +
  labs(y='proportion measured',
       x='region type')
