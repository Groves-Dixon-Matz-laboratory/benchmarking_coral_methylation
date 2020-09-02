#gbm_and_ge.R

# LOAD GBM LEVEL DATA -----------------------------------------------------

#load data from compare_level.R
ll=load('comparisons/datasets/gbmLvl.Rdata')
ll
head(gbm.dat)
aldat = gbm.dat %>% 
  select(chr, start, end, name, l.fracMeth, pm.glmLvl, mrB, mrF, mrM, mrF.s, mrM.s, mrB.s, mbd.score, mbd.FPKM, mbd.FPCpGM)
dim(aldat)


# LOAD GENE EXPRESSION ----------------------------------------------------


