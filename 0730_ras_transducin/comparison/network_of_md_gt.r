## name: network_of_md_gt.r
## date: 09/03/2015

## Here I try to explore parameters(cutoff.cij) for Ras networks.

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna2.R")

membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])
# x11() fullscreen

# remodel.cna() only use the cij to remodel; I don't need to calculate cna but pass in the cij
# What I did is provide dummy nets_pearson but change the cij (minus.log!)
# This is much faster when j is very small

# prepare the dummy nets!
cij_gtp <- filter.dccm(net.md$gtp$cijs0,
  cmap=cm$GTP, cutoff.cij=0.6)
cij_gdp <- filter.dccm(net.md$gdp$cijs0,
  cmap=cm$GDP, cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp, cutoff.cij=0)
cna_gdp <- cna(cij_gdp, cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

nets_md_gt_pearson_remodel <- NULL
layout(matrix(1:10,nrow=2))
for (j in seq(0.3,0.7,by=0.1)) {

  #### pearson !

  # calculate consensu cij
  cij_gtp <- filter.dccm(net.md$gtp$cijs0, cmap=cm$GTP, cutoff.cij=j)
  cij_gdp <- filter.dccm(net.md$gdp$cijs0, cmap=cm$GTP, cutoff.cij=j)

  # abs!
  cij_gtp <- abs(cij_gtp)
  cij_gdp <- abs(cij_gdp)

  # minus.log cij
  cij_gtp[cij_gtp>=1] <- 0.9999
  cij_gtp[cij_gtp>0] <- -log(cij_gtp[cij_gtp>0])
  cij_gdp[cij_gdp>=1] <- 0.9999
  cij_gdp[cij_gdp>0] <- -log(cij_gdp[cij_gdp>0])

  # change the cij of dummy networks
  nets_pearson <- nets_dummy
  nets_pearson$gtp$cij <- cij_gtp
  nets_pearson$gdp$cij <- cij_gdp

  # remodel networks
  nets_md_pearson_remodel <- remodel.cna(nets_pearson, member=membership_gt, method="sum",
    col.edge="feature", scut=4, normalize=FALSE)

  # calculate weights
  w1_pearson = (E(nets_md_pearson_remodel[[1]]$community.network)$weight)*1
  w2_pearson = (E(nets_md_pearson_remodel[[2]]$community.network)$weight)*1

  # plot!
  plot.cna(nets_md_pearson_remodel[[1]], layout=layout_2d, weights = w1_pearson,
    vertex.label=NA, main=paste0("gtp_cutoff.cij=",j))
  plot.cna(nets_md_pearson_remodel[[2]], layout=layout_2d, weights = w2_pearson,
    vertex.label=NA, main=paste0("gdp_cutoff.cij=",j))

  # list and save networks
  nets_md_gt_pearson_remodel <- c(nets_md_gt_pearson_remodel,
    list(nets_md_pearson_remodel))
}
mtext("Transducin_MD_networks(pearson)", line=-40, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_gt_pearson_0.3_to_0.7.pdf")

names(nets_md_gt_pearson_remodel) <- seq(0.3,0.7,by=0.1)
# if you print nets_md_gt_pearson_remodel, the description of all 7 members are the same
# but if you check community.cij, they are totally different!


save(nets_md_gt_pearson_remodel,
     file="network_of_md_gt.RData")



