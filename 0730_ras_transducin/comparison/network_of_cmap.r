## name: network_of_cmap.r
## date: 09/02/2015

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna2.R")

#### Ras ####

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

cna_cmap_gtp <- cna(cmap_ca_gtp_dist10[1:166,1:166], cutoff.cij=0)
cna_cmap_gdp <- cna(cmap_ca_gdp_dist10[1:166,1:166], cutoff.cij=0)

# remodel networks
nets_cmap <- list(gtp=cna_cmap_gtp,gdp=cna_cmap_gdp)
nets_cmap_remodel <- remodel.cna(nets_cmap, member=membership_ras, method="sum",
  col.edge="feature", scut=4, normalize=FALSE)

# calculate weights
w1_cmap = (E(nets_cmap_remodel[[1]]$community.network)$weight)*1
w2_cmap = (E(nets_cmap_remodel[[2]]$community.network)$weight)*1

# plot!
par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_cmap_remodel[[1]], layout=layout_2d[1:9,], weights = w1_cmap,
  vertex.label=NA, main="cna_cmap_gtp")
plot.cna(nets_cmap_remodel[[2]], layout=layout_2d[1:9,], weights = w2_cmap,
  vertex.label=NA, main="cna_cmap_gdp")

dev.copy2pdf(file="figures/cna_cmap_ras.pdf")

#### Gt ####

membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])

cna_cmap_gtp <- cna(cm$GTP, cutoff.cij=0)
cna_cmap_gdp <- cna(cm$GDP, cutoff.cij=0)


# remodel networks
nets_cmap <- list(gtp=cna_cmap_gtp,gdp=cna_cmap_gdp)
nets_cmap_remodel <- remodel.cna(nets_cmap, member=membership_gt, method="sum",
  col.edge="feature", scut=4, normalize=FALSE)

# calculate weights
w1_cmap = (E(nets_cmap_remodel[[1]]$community.network)$weight)*1
w2_cmap = (E(nets_cmap_remodel[[2]]$community.network)$weight)*1

# plot!
par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_cmap_remodel[[1]], layout=layout_2d, weights = w1_cmap,
  vertex.label=NA, main="cna_cmap_gtp")
plot.cna(nets_cmap_remodel[[2]], layout=layout_2d, weights = w2_cmap,
  vertex.label=NA, main="cna_cmap_gdp")

dev.copy2pdf(file="figures/cna_cmap_transducin.pdf")

