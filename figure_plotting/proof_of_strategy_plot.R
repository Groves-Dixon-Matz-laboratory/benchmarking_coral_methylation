#proof_of_strategy_plot.R
rm(list=ls())
source('benchmarking_functions.R')

# LOAD DATA ---------------------------------------------------------------
#load the pool count reduced data for WGBS and MBDseq, and the normal counts for MethylRAD
ll=load('picomethyl/poolTest/countReducedDifferences.Rdata')
pm.csums=csums
ll
ll=load('mbdSeq/poolTest/countReducedDifferencesWub.Rdata')
mbd.csums=csums
ll
ll=load('methylRAD/datasets/countReducedDifferences.Rdata')
mr.csums=csums
ll



# SELECT DATA -------------------------------------------------------------
#select particular reductions based on correlation plataeus figure 4

#WGBS
pm = pm.dat %>% 
  dplyr::select(name, pct.50_meth.diff) %>% 
  set_names(c('name', 'pm')) %>% 
  mutate(name=as.character(name)) %>% 
  as_tibble
pm.nreads = 0.5*sum(pm.csums) / 1e6

#MBD
mbd = mbd.dat %>% 
  dplyr::select(name, pct.25_log2FoldChange) %>% 
  set_names(c('name', 'mbd')) %>% 
  as_tibble
mbd.nreads = 0.25*sum(mbd.csums)/1e6

#MethylRAD
mr = mr.dat %>% 
  dplyr::select(name, pct.100_log2FoldChange) %>% 
  set_names(c('name', 'mr')) %>% 
  as_tibble
mr.nreads = 0.125*sum(mr.csums)

dList = list(pm, mbd, mr)
dat = dList %>% 
  purrr::reduce(full_join, by='name')


# GET PCA -----------------------------------------------------------

log2folds = dat %>% 
  select(pm, mbd, mr)
scaled = map_dfc(log2folds, ~ (.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)) %>% 
  data.frame()
rownames(scaled)=dat$name
scaled_nona = na.omit(scaled)


#BUILD PCA
mpca = prcomp(scaled_nona)
percentVar <- mpca$sdev^2/sum(mpca$sdev^2)
mcoords = mpca$x %>% 
  data.frame() %>% 
  mutate(PC1=PC1*-1)
pdat = cbind(scaled_nona, mcoords)

#GET RESIDUALS
get_resid = function(var) {
  dat = select(pdat, PC1, y = {{ var }})
  as.numeric(lm(y ~ PC1, data = dat)$resid)
}

rdat = pdat %>% 
  mutate(pm.resid = get_resid(pm),
         mr.resid = get_resid(mr),
         mbd.resid = get_resid(mbd),
         name=rownames(scaled_nona)) %>% 
  as_tibble()

#CALL FALSE POSITIVES
RCUT=1.96
fpdat = rdat %>% 
  mutate(pm.fp = abs(pm.resid) > RCUT,
         mr.fp = abs(mr.resid) > RCUT,
         mbd.fp = abs(mbd.resid) > RCUT)


fpdat %>% 
  ggplot(aes(x=PC1, y=mr, color=mr.fp)) +
  geom_point() +
  scale_color_manual(values=c('black', 'red'))


# PLOT --------------------------------------------------------------------

assay_resid = function(var){
  d2 = dplyr::select(dat, name, y=mr, x={{var}}) %>% 
    na.omit()
  d2$resid = as.numeric(lm(y ~ x, data = d2)$resid)
  return(d2)
}

add_resids = function(r, df, colName){
  r=r %>% 
    dplyr::select(name, resid)
  colnames(r)[2]=colName
  df %>% 
    left_join({{r}}, by='name')
}

pmr = assay_resid(pm)
mbdr = assay_resid(mbd)

dat2 = add_resids(pmr, dat, 'pmr')
dat2 = add_resids(pmr, dat2, 'mbdr')



plot_scatter_pearsonCor_annotated(dat2, 'pm', 'mr', 'WGBS', 'MethylRAD', ALPHA=0.3)
plot_scatter_pearsonCor_annotated(dat2, 'mbd', 'mr', 'MBD-seq', 'MethylRAD', ALPHA=0.3)
nred = sum(fpdat$mr.fp)


#FINAL PLOT
plot_scatter_pearsonCor_annotated(fpdat, 'PC1', 'mr', 'PC1', 'MethylRAD', ALPHA=0.3) + 
  geom_point(aes(color=mr.fp)) +
  scale_color_manual(values=c('black', 'red')) +
  labs(color='|residual| > 1.96',
       subtitle=paste('N =', nred),
       x=paste('PC1 (', round(percentVar[1], digits=2)*100, '% variance explained)', sep='')) +
  theme(plot.subtitle = element_text(color='red'))




#VOLCANO
nsig = sum(mr.dat$pct.12.5_padj < 0.1, na.rm=TRUE)
mr.dat %>% 
  select(lf=pct.12.5_log2FoldChange, p=pct.12.5_padj) %>% 
  filter(!is.na(p)) %>% 
  mutate(sig=factor(p<0.1, levels=c(TRUE, FALSE))) %>% 
  ggplot(aes(x=lf, y=-log(p, 10), color=sig)) +
  geom_point(alpha=0.3) +
  scale_color_manual(values=c('red', 'black')) +
  labs(x=bquote(log[2]*'fold difference'),
       y=bquote("-"*log[10]*'pvalue'),
       color='FDR>0.1',
       subtitle=paste('N =', nsig)) +
  theme(plot.subtitle)
  
