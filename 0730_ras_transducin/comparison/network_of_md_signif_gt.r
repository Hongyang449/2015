## name: network_of_md_signif_gt.r
## date: 09/05/2015

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna3.R")
source("/Users/hyangl/Desktop/2015/functions/get.signif.R")
source("/Users/hyangl/Desktop/2015/functions/get.community.cij.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.nets.R")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")

p.cutoff <- 0.05
cutoff <- seq(0.0,0.7,by=0.05)
layout_2d <- layout_2d

membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])

# prepare the dummy nets!
cij_gtp <- filter.dccm(net.md$gtp$cijs0, cmap=cm$GTP, cutoff.cij=0.6)
cij_gdp <- filter.dccm(net.md$gdp$cijs0, cmap=cm$GDP, cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp, cutoff.cij=0)
cna_gdp <- cna(cij_gdp, cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

cij_pearson <- list(gtp=net.md$gtp$cijs0, gdp=net.md$gdp$cijs0)
cmap <- list(gtp=cm$GTP, gdp=cm$GDP)

########## pearson ###########

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_gt_pearson <- lapply(as.list(cutoff), function(x) {
  get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=x, membership=membership_gt)
  })
names(community_cij_md_gt_pearson) <- cutoff

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_gt_pearson <- lapply(community_cij_md_gt_pearson, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_md_gt_pearson) <- cutoff

# remodel nets
nets_md_gt_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_gt, signif=p_community_cij_md_gt_pearson[[x]],
    layout_2d=layout_2d, p.cutoff=p.cutoff)
  })
names(nets_md_gt_pearson_remodel) <- cutoff

save(community_cij_md_gt_pearson,
     p_community_cij_md_gt_pearson,
     nets_md_gt_pearson_remodel,
     file="network_of_md_signif_gt.RData")

# plot!

plot.nets(nets_md_gt_pearson_remodel[as.character(seq(0.4,0.6,by=0.05))], layout_2d=layout_2d)
mtext("Gt_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_gt_pearson_signif_0.4_0.6.pdf")



  
