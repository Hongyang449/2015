## name: network_of_ttest_cij.r
## date: 09/15/2015

## Here I only use edges that pass t.test to build networks.

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/ttest_of_cij.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna3.R")
source("/Users/hyangl/Desktop/2015/functions/get.signif.R")
source("/Users/hyangl/Desktop/2015/functions/get.community.cij.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.nets.R")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")


p.cutoff <- 0.01
cutoff <- c(0.05,0.01,0.005)

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])

########### ras #########
cij_gtp <- filter.dccm(array_cij_ca_pearson_gtp[,,1], cutoff.cij=0.6)
cij_gdp <- filter.dccm(array_cij_ca_pearson_gdp[,,1], cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp[1:166,1:166], cutoff.cij=0)
cna_gdp <- cna(cij_gdp[1:166,1:166], cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

cij_pearson <- list(gtp=array_cij_ca_pearson_gtp[1:166,1:166,],gdp=array_cij_ca_pearson_gdp[1:166,1:166,])
cij_lmi <- list(gtp=array_cij_ca_lmi_gtp[1:166,1:166,],gdp=array_cij_ca_lmi_gdp[1:166,1:166,])
cmap <- list(gtp=cmap_ca_gtp_dist10[1:166,1:166],gdp=cmap_ca_gdp_dist10[1:166,1:166])

## pearson
community_cij_ttest_ras_pearson <- lapply(as.list(cutoff), function(x) {
  cij_filter <- matrix(0, nrow=dim(p_ras_pearson), ncol=dim(p_ras_pearson))
  cij_filter[p_ras_pearson < x] <- 1
  get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=0,
    membership=membership_ras, cij.filter=cij_filter)
  })
names(community_cij_ttest_ras_pearson) <- cutoff

p_community_cij_ttest_ras_pearson <- lapply(community_cij_ttest_ras_pearson, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_ttest_ras_pearson) <- cutoff

nets_ttest_ras_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  cij_filter <- matrix(0, nrow=dim(p_ras_pearson), ncol=dim(p_ras_pearson))
  cij_filter[p_ras_pearson < cutoff[x]] <- 1
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=0, extra.filter=cij_filter,
    nets_dummy=nets_dummy, membership=membership_ras, 
    signif=p_community_cij_ttest_ras_pearson[[x]], p.cutoff=p.cutoff)
  })
names(nets_ttest_ras_pearson_remodel) <- cutoff

## lmi
community_cij_ttest_ras_lmi <- lapply(as.list(cutoff), function(x) {
  cij_filter <- matrix(0, nrow=dim(p_ras_lmi), ncol=dim(p_ras_lmi))
  cij_filter[p_ras_lmi < x] <- 1
  get.community.cij(cij=cij_lmi, cmap=cmap, cutoff.cij=0,
    membership=membership_ras, cij.filter=cij_filter)
  })
names(community_cij_ttest_ras_lmi) <- cutoff

p_community_cij_ttest_ras_lmi <- lapply(community_cij_ttest_ras_lmi, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_ttest_ras_lmi) <- cutoff

nets_ttest_ras_lmi_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  cij_filter <- matrix(0, nrow=dim(p_ras_lmi), ncol=dim(p_ras_lmi))
  cij_filter[p_ras_lmi < cutoff[x]] <- 1
  remodel.nets(cij=cij_lmi, cmap=cmap, cutoff.cij=0, extra.filter=cij_filter,
    nets_dummy=nets_dummy, membership=membership_ras,
    signif=p_community_cij_ttest_ras_lmi[[x]], p.cutoff=p.cutoff)
  })
names(nets_ttest_ras_lmi_remodel) <- cutoff

########### gt #########
cij_gtp <- filter.dccm(net.md$gtp$cijs0, cmap=cm$GTP, cutoff.cij=0.6)
cij_gdp <- filter.dccm(net.md$gdp$cijs0, cmap=cm$GDP, cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp, cutoff.cij=0)
cna_gdp <- cna(cij_gdp, cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

cij_pearson <- list(gtp=net.md$gtp$cijs0, gdp=net.md$gdp$cijs0)
cmap <- list(gtp=cm$GTP, gdp=cm$GDP)

## pearson
community_cij_ttest_gt_pearson <- lapply(as.list(cutoff), function(x) {
  cij_filter <- matrix(0, nrow=dim(p_gt_pearson), ncol=dim(p_gt_pearson))
  cij_filter[p_gt_pearson < x] <- 1
  get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=0,
    membership=membership_gt, cij.filter=cij_filter)
  })
names(community_cij_ttest_gt_pearson) <- cutoff

p_community_cij_ttest_gt_pearson <- lapply(community_cij_ttest_gt_pearson, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_ttest_gt_pearson) <- cutoff

nets_ttest_gt_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  cij_filter <- matrix(0, nrow=dim(p_gt_pearson), ncol=dim(p_gt_pearson))
  cij_filter[p_gt_pearson < cutoff[x]] <- 1
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=0, extra.filter=cij_filter,
    nets_dummy=nets_dummy, membership=membership_gt,
    signif=p_community_cij_ttest_gt_pearson[[x]], p.cutoff=p.cutoff)
  })
names(nets_ttest_gt_pearson_remodel) <- cutoff

save(community_cij_ttest_ras_pearson, community_cij_ttest_ras_lmi,
     p_community_cij_ttest_ras_pearson, p_community_cij_ttest_ras_lmi,
     nets_ttest_ras_pearson_remodel, nets_ttest_ras_lmi_remodel,
     file="network_of_ttest_cij.RData")

# plot!
library(bio3d)
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_ttest_cij.RData")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")
cutoff <- c(0.05,0.01,0.005)

plot.nets(nets_ttest_ras_pearson_remodel[as.character(cutoff)], layout_2d=layout_2d[1:9,])
mtext("Ras_ttest_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/ras_ttest_pearson_signif.pdf")

plot.nets(nets_ttest_ras_lmi_remodel[as.character(cutoff)], layout_2d=layout_2d[1:9,])
mtext("Ras_ttest_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/ras_ttest_lmi_signif.pdf")

plot.nets(nets_ttest_gt_pearson_remodel[as.character(cutoff)], layout_2d=layout_2d)
mtext("Gt_ttest_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/gt_ttest_pearson_signif.pdf")


