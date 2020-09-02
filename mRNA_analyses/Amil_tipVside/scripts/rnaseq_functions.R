#modified version of the function provided in the DESeq package that can be run from a dataframe
#instead of DESeqTransform object output from
mod.plotPCA.df <- function (df, coldat, intgroup = "condition", ntop = 25000, returnData = F, pcs = 10, pc1 = 1, pc2 = 2, main = "\n", SIZE = 5) 
{
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- as.data.frame(coldat[, intgroup, 
                                      drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    coldat[[intgroup]]
  }
  d <- data.frame(pca$x[,1:pcs], group = group, 
                  intgroup.df, name = colnames(df))
  attr(d, "percentVar") <- percentVar[1:2]
  g = ggplot(data = d, aes_string(x = paste('PC', pc1, sep = ''),
                                  y = paste('PC', pc2, sep = ''), color = "group")) + 
    geom_point(size = SIZE) + xlab(paste0(paste0(paste0("PC", pc1), ": "), round(percentVar[pc1] * 
                                                                                   100), "% variance")) + ylab(paste0(paste0(paste0("PC", pc2), ": "), round(percentVar[pc2] * 
                                                                                                                                                               100), "% variance")) + coord_fixed()
  g = g + ggtitle(main)
  g = g + theme_bw()
  print(g)
  if (returnData == T){
    return(d)
  }
}