## name: network_of_md_ras.r
## date: 09/01/2015

## Here I try to explore parameters(cutoff.cij) for Ras networks.

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna2.R")

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
# x11() fullscreen

# remodel.cna() only use the cij to remodel; I don't need to calculate cna but pass in the cij
# What I did is provide dummy nets_pearson but change the cij (minus.log!)
# This is much faster when j is very small

# prepare the dummy nets!
cij_gtp <- filter.dccm(array_cij_ca_pearson_gtp,
  cmap=cmap_ca_gtp_dist10, cutoff.cij=0.6)
cij_gdp <- filter.dccm(array_cij_ca_pearson_gdp,
  cmap=cmap_ca_gdp_dist10, cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp[1:166,1:166], cutoff.cij=0)
cna_gdp <- cna(cij_gdp[1:166,1:166], cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

nets_md_ras_pearson_remodel <- NULL
layout(matrix(1:14,nrow=2))
for (j in seq(0.0,0.6,by=0.1)) {

  #### pearson !

  # calculate consensu cij
  cij_gtp <- filter.dccm(array_cij_ca_pearson_gtp[1:166,1:166,],
    cmap=cmap_ca_gtp_dist10[1:166,1:166], cutoff.cij=j)
  cij_gdp <- filter.dccm(array_cij_ca_pearson_gdp[1:166,1:166,],
    cmap=cmap_ca_gdp_dist10[1:166,1:166], cutoff.cij=j)

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
  nets_md_pearson_remodel <- remodel.cna(nets_pearson, member=membership_ras, method="sum",
    col.edge="feature", scut=4, normalize=FALSE)

  # calculate weights
  w1_pearson = (E(nets_md_pearson_remodel[[1]]$community.network)$weight)*1
  w2_pearson = (E(nets_md_pearson_remodel[[2]]$community.network)$weight)*1

  # plot!
  plot.cna(nets_md_pearson_remodel[[1]], layout=layout_2d[1:9,], weights = w1_pearson,
    vertex.label=NA, main=paste0("gtp_cutoff.cij=",j))
  plot.cna(nets_md_pearson_remodel[[2]], layout=layout_2d[1:9,], weights = w2_pearson,
    vertex.label=NA, main=paste0("gdp_cutoff.cij=",j))

  # list and save networks
  nets_md_ras_pearson_remodel <- c(nets_md_ras_pearson_remodel,
    list(nets_md_pearson_remodel))
}
mtext("Ras_MD_networks(pearson)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson.pdf")

names(nets_md_ras_pearson_remodel) <- seq(0.0,0.6,by=0.1)
# if you print nets_md_ras_pearson_remodel, the description of all 7 members are the same
# but if you check community.cij, they are totally different!

nets_md_ras_lmi_remodel <- NULL
layout(matrix(1:14,nrow=2))
for (j in seq(0.0,0.6,by=0.1)) {

  #### lmi !

  # calculate consensu cij
  cij_gtp <- filter.dccm(array_cij_ca_lmi_gtp[1:166,1:166,],
    cmap=cmap_ca_gtp_dist10[1:166,1:166], cutoff.cij=j)
  cij_gdp <- filter.dccm(array_cij_ca_lmi_gdp[1:166,1:166,],
    cmap=cmap_ca_gdp_dist10[1:166,1:166], cutoff.cij=j)

  # abs!
  cij_gtp <- abs(cij_gtp)
  cij_gdp <- abs(cij_gdp)

  # minus.log cij
  cij_gtp[cij_gtp>=1] <- 0.9999
  cij_gtp[cij_gtp>0] <- -log(cij_gtp[cij_gtp>0])
  cij_gdp[cij_gdp>=1] <- 0.9999
  cij_gdp[cij_gdp>0] <- -log(cij_gdp[cij_gdp>0])

  # change the cij of dummy networks
  nets_lmi <- nets_dummy
  nets_lmi$gtp$cij <- cij_gtp
  nets_lmi$gdp$cij <- cij_gdp

  # remodel networks
  nets_md_lmi_remodel <- remodel.cna(nets_lmi, member=membership_ras, method="sum",
    col.edge="feature", scut=4, normalize=FALSE)

  # calculate weights
  w1_lmi = (E(nets_md_lmi_remodel[[1]]$community.network)$weight)*1
  w2_lmi = (E(nets_md_lmi_remodel[[2]]$community.network)$weight)*1

  # plot!
  plot.cna(nets_md_lmi_remodel[[1]], layout=layout_2d[1:9,], weights = w1_lmi,
    vertex.label=NA, main=paste0("gtp_cutoff.cij=",j))
  plot.cna(nets_md_lmi_remodel[[2]], layout=layout_2d[1:9,], weights = w2_lmi,
    vertex.label=NA, main=paste0("gdp_cutoff.cij=",j))

  # list and save networks
  nets_md_ras_lmi_remodel <- c(nets_md_ras_lmi_remodel,
    list(nets_md_lmi_remodel))
}
mtext("Ras_MD_networks(lmi)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_lmi.pdf")

names(nets_md_ras_lmi_remodel) <- seq(0.0,0.6,by=0.1)

save(nets_md_ras_pearson_remodel,
     nets_md_ras_lmi_remodel,
     file="network_of_md_ras.RData")



