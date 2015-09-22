## name: ttest_of_cij.r
## date: 09/15/2015

# Here I perform a t.test for each cij edge.

library(bio3d)

load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")

source("/Users/hyangl/Desktop/2015/functions/get.signif.R")

# t.test for the whole matrix and find signif different edges
p_ras_pearson <- get.signif(array_cij_ca_pearson_gtp[1:166,1:166,],
  array_cij_ca_pearson_gdp[1:166,1:166,])
p_ras_lmi <- get.signif(array_cij_ca_lmi_gtp[1:166,1:166,],
  array_cij_ca_lmi_gdp[1:166,1:166,])

length(which(p_ras_pearson<0.05))
#[1] 2199 * 2
length(which(p_ras_lmi<0.05))
#[1] 2051 * 2
length(which(p_ras_pearson<0.01))
#[1] 822 * 2
length(which(p_ras_lmi<0.01))
#[1] 748 * 2
## we can see the number of edges left (0.01) are similar when use cutoff.cij=0.4/0.55
## It seems the number of significant positions of lmi is small.
## In terms of this, it seems pearson can distinguish two states a little better than lmi.

p_gt_pearson <- get.signif(net.md$gtp$cijs0, net.md$gdp$cijs0)

length(which(p_gt_pearson<0.05))
# [1] 6240
length(which(p_gt_pearson<0.01))
# [1] 1360

mean_cij_ras_pearson <- list(gtp=apply(array_cij_ca_pearson_gtp[1:166,1:166,],c(1:2),mean),
  gdp=apply(array_cij_ca_pearson_gdp[1:166,1:166,],c(1:2),mean))
sd_cij_ras_pearson <- list(gtp=apply(array_cij_ca_pearson_gtp[1:166,1:166,],c(1:2),sd),
  gdp=apply(array_cij_ca_pearson_gdp[1:166,1:166,],c(1:2),sd))

mean_cij_ras_lmi <- list(gtp=apply(array_cij_ca_lmi_gtp[1:166,1:166,],c(1:2),mean),
  gdp=apply(array_cij_ca_lmi_gdp[1:166,1:166,],c(1:2),mean))
sd_cij_ras_lmi <- list(gtp=apply(array_cij_ca_lmi_gtp[1:166,1:166,],c(1:2),sd),
  gdp=apply(array_cij_ca_lmi_gdp[1:166,1:166,],c(1:2),sd))

mean_cij_gt_pearson <- list(gtp=apply(net.md$gtp$cijs0,c(1:2),mean),
  gdp=apply(net.md$gdp$cijs0,c(1:2),mean))
sd_cij_gt_pearson <- list(gtp=apply(net.md$gtp$cijs0,c(1:2),sd),
  gdp=apply(net.md$gdp$cijs0,c(1:2),sd))

save(p_ras_pearson, p_ras_lmi,
     p_gt_pearson,
     mean_cij_ras_pearson, sd_cij_ras_pearson,
     mean_cij_ras_lmi, sd_cij_ras_lmi,
     mean_cij_gt_pearson, sd_cij_gt_pearson,
     file="ttest_of_cij.RData")




