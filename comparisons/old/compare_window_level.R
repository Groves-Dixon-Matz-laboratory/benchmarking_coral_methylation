#compare_window_level.R
source('benchmarking_functions.R')
library(methylKit)
library(tidyverse)

# load methylRAD data -----------------------------------------------------

ll = load("methylRAD/host/window_results/methylRAD_10kb_window_counts_mnLvl.Rdata")
ll
head(bm)

# load picomethyl data ----------------------------------------------------

ll=load("picomethyl/window_results/tipVside_Windows_10000_10000_tiles_formatted.Rdata")
ll
head(pm)

# load mbdseq from recip meth ---------------------------------------------

ll=load("mbdSeq/host/window_results/mbd_10kb_window_lvls.Rdata")
ll
head(mbd)


# merge them up -----------------------------------------------------------

datList = list(bm, pm, mbd)
mdat = purrr::reduce(datList, full_join, by = 'gene') %>% 
  dplyr::select(-ts, -cs)
mdat

# correlate ---------------------------------

cols = colnames(mdat)[2:ncol(mdat)]
pltList = list()
for (c1 in cols){
  for (c2 in cols){
    plt=mdat %>% 
      ggplot(aes_string(x=c1,y=c2)) +
      geom_point(alpha=0.2) +
      # geom_smooth(method='lm') +
      labs(x=c1, y=c2)
    pltList[[paste(c1,c2)]]=plt
  }
}

length(pltList)
plot_grid(plotlist=pltList, nrow=3)


names(pltList)

pmx = list(pltList[['pm mbd']], pltList[['pm mnMR']])
plot_grid(plotlist=pmx, nrow=1)


