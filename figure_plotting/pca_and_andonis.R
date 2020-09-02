#pca_and_adonis.R
#plot PCAs and adonis partitioning of variance piecharts 
#for GBM measure from each assay

rm(list=ls())
library(vegan)
library(gtools)
library(DESeq2)
source('benchmarking_functions.R')

# PICOMETHYL --------------------------------------------------------------

ll=load('picomethyl/datasets/regionCounts/mprops_by_sample.Rdata')
ll

#prepare data for pca/adonis
mdat = mList[[1]]

# samples = colnames(mdat)[4:ncol(mdat)]
genotype =substr(samples, start=1, stop=1)
tissue = if_else(substr(samples, start=2, stop=2)=='s',
                 'side',
                 'tip')

#get the regions
dat = mdat %>% 
  tidyr::unite('reg', chr,start,end, sep='_')
reg = dat %>% 
  pull(reg)
#separate the data
dat = dat %>% 
  dplyr::select(-reg) %>% 
  data.frame()
rownames(dat)=reg
head(dat)

mns=apply(dat, 1, mean)
nna = sum(is.na(mns))
pct=paste('(', round(nna/nrow(dat),digits=2)*100, '%)', sep='')
print(paste('WARNING', nna, 'of', nrow(dat), pct, 'have missing data and will be removed.' ))
noNa = dat[!is.na(mns),]

#check the distribution
hist(log(mns),breaks=100) #doesn't look great, shows loss of data from using methylKit to get props


#plot pca
pm.coldata = data.frame(Run=samples,
                        genotype=substr(samples, start=1, stop=1),
                        tissue=substr(samples, start=2, stop=2))
pm.pca = plot_vsd_pca(noNa, pm.coldata)

#adonis
ad=adonis(t(noNa)~genotype+tissue)
labs=c("genotype","tissue","residuals")
aovTab = ad$aov.tab
r2s = aovTab$R2[1:3]
ps = aovTab$`Pr(>F)`[1:2]
stars = stars.pval(ps)
labs[1] = paste(labs[1], stars[1], sep='')
labs[2] = paste(labs[2], stars[2], sep='')
cols=append(gg_color_hue(2), 'grey')
# pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="WGBS")
pidat = data.frame(r2 = r2s/sum(r2s)*100,
                   p = append(ps, NA),
                   lab=factor(labs, levels=labs))
#plot with ggplot
pm.pie = pidat %>% 
  ggplot() +
  geom_bar(aes(fill=lab, y=r2, x=''), stat='identity', color='black') +
  coord_polar("y", start=.5) +
  scale_fill_manual(values=cols) +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())



# MBDSEQ ------------------------------------------------------------------

ll=load('mbdSeq/datasets/vsd_results.Rdata')
ll
mdat = vsdList[[1]]
mdat2 = mdat+abs(min(mdat))


#check distribition
mns = apply(mdat, 1, mean)
hist(mns, breaks=100)


#plot pca to show importance of enzyme
print('Sample names lined up?')
sum(paste('m',colnames(mdat), sep='.')==mcoldata$Run)==ncol(mdat)
mbd.pca = plot_vsd_pca(mdat, mcoldata)


#plot adonis
mbd.pie = plot_adonis_from_vsd(mdat2, mcoldata)


# METHYLRAD ---------------------------------------------------------------

#LOAD THE FPKM DATA AND SET UP COLDATA
ll=load('methylRAD/datasets/mdRAD_geneFPKMs.Rdata')
ll
mrSamples=colnames(geneFpkm)
enzyme = substr(mrSamples, 3,3)
genotype = sapply(mrSamples, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
tissue0 = sapply(mrSamples, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3])
tissue = substr(tissue0, 1,1)
mrColdata = data.frame(Run=mrSamples,
                       genotype,
                       tissue,
                       enzyme)

#SUBSET FOR THE TWO ENZYMES
fColdata = mrColdata %>% 
  filter(enzyme=='F')
fFpkm = geneFpkm[,fColdata$Run]
mColdata = mrColdata %>% 
  filter(enzyme=='M')
mFpkm = geneFpkm[,mColdata$Run]

#PLOT OVERALL RESULTS
#overall pca
mr.pca = plot_vsd_pca(geneFpkm, mrColdata, shape.col = 'enzyme') + 
  theme(legend.position = 'right') +
  labs(shape='enzyme',
       color='colony')

#overal adonis

ad=adonis(t(geneFpkm)~genotype+tissue+enzyme)
labs=c("genotype","tissue", "enzyme", "residuals")
aovTab = ad$aov.tab
r2s = aovTab$R2[1:4]
ps = aovTab$`Pr(>F)`[1:3]
stars = stars.pval(ps)
labs[1] = paste(labs[1], stars[1], sep='')
labs[2] = paste(labs[2], stars[2], sep='')
labs[3] = paste(labs[3], stars[3], sep='')
cols=append(gg_color_hue(3), 'grey')
# pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="WGBS")
pidat = data.frame(r2 = r2s/sum(r2s)*100,
                   p = append(ps, NA),
                   lab=factor(labs, levels=labs))
#plot with ggplot
mr.pie = pidat %>% 
  ggplot() +
  geom_bar(aes(fill=lab, y=r2, x=''), stat='identity', color='black') +
  coord_polar("y", start=.5) +
  scale_fill_manual(values=cols) +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

#--- PLOT INDIVIDUAL ENZYMES ---#

#FSPE1
#pca
f.pca = plot_vsd_pca(fFpkm, fColdata) + 
  theme(legend.position = 'right') +
  labs(shape='tissue',
       color='colony')

#adonis
f.pie = plot_adonis_from_vsd(fFpkm, fColdata)

#MSPJ1
#pca
m.pca = plot_vsd_pca(mFpkm, mColdata) + 
  theme(legend.position = 'right') +
  labs(shape='tissue',
       color='colony')

#adonis
m.pie = plot_adonis_from_vsd(mFpkm, mColdata)



# ASSEMBLE ALL RESULTS ----------------------------------------------------

#assemble pcas
pcaList = list(pm.pca + 
                 labs(title='WGBS') +
                 theme(legend.position = 'none'),
               mbd.pca + 
                 labs(title='MBD-seq') +
                 theme(legend.position = 'none'),
               f.pca + 
                 labs(title='FspE1') +
                 theme(legend.position = 'none'),
               m.pca +
                 labs(title='MspJ1'))
pcaList = lapply(pcaList, function(x) return(x+theme(plot.title = element_text(hjust=0.5))))


#assemble pie charts
pieList = list(pm.pie + labs(title='WGBS'),
               mbd.pie + labs(title='MBD-seq'),
               f.pie + labs(title='FspE1'),
               m.pie + labs(title='MspJ1'))
pieList = lapply(pieList, function(x) x+theme(plot.title=element_text(hjust=0.5)))


plot_grid(plotlist = pcaList, nrow=1, rel_widths = c(1,1,1,1.1))
plot_grid(plotlist = pieList, nrow=2)


#plot full mdRAD
plot_grid(mr.pca,
          mr.pie)

