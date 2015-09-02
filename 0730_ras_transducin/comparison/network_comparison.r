## name: network_comparison.r
## date: 08/18/2015

# use sd instead of var to calculate edge colors between networks
# color the edge "separately"

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/gt_networks.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij_consensus.RData")

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/0105_manuscript/pca/sse_for_plot.RData")

source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/draw.box.R")
source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/label_ras.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.cna2.R")

pdb_ras <- read.pdb("/Users/hyangl/Desktop/0105_manuscript/dccm/pdb/gtp_plot.pdb")

# cna; 1:166
cna_ras_gtp <- cna(cij_ca_pearson_gtp_consensus[1:166,1:166], cutoff.cij=0)
cna_ras_gdp <- cna(cij_ca_pearson_gdp_consensus[1:166,1:166], cutoff.cij=0)

nets_ras <- list(gtp=cna_ras_gtp,gdp=cna_ras_gdp)
nets_gt <- nets_gt_3grps[1:2]

## import membership from ali
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
name_ras <- c("b1-b3,a1","p-loop","SI","SII","b4-b6","a3","a4","a5","loop8")
names(name_ras) <- 1:length(name_ras)
membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])
name_gt <- c("b1-b3,a1","p-loop","SI","SII","b4-b6", "SIII,a3","a4","a5","loop8","H1","H2")
names(name_gt) <- 1:length(name_gt)

nets_ras_remodel <- remodel.cna(nets_ras, member=membership_ras, method="sum", col.edge=
"feature", scut=4, normalize=FALSE)
nets_gt_remodel <- remodel.cna(nets_gt, member=membership_gt, method="sum", col.edge="feature", scut=4, normalize=FALSE)

w1_ras = (E(nets_ras_remodel[[1]]$community.network)$weight)*1
w2_ras = (E(nets_ras_remodel[[2]]$community.network)$weight)*1
w1_gt = (E(nets_gt_remodel[[1]]$community.network)$weight)*1
w2_gt = (E(nets_gt_remodel[[2]]$community.network)$weight)*1


par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_ras_remodel[[1]], layout=layout_2d[1:9,], weights = w1_ras, vertex.label=NA, main="ras_gtp")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
plot.cna(nets_ras_remodel[[2]], layout=layout_2d[1:9,], weights = w2_ras, vertex.label=NA, main="ras_gdp")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_ras.pdf")

par(mfrow=c(1,2), mar = c(2,2,4,2), pty="s")
plot.cna(nets_gt_remodel[[1]], layout=layout_2d, weights = w1_gt, vertex.label=NA, main="transducin_gtp")
for (i in 1:length(name_gt)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_gt[i], col="black", cex=1)
}
plot.cna(nets_gt_remodel[[2]], layout=layout_2d, weights = w2_gt, vertex.label=NA, main="transducin_gdp")
for (i in 1:length(name_gt)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_gt[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_transducin.pdf")

save(nets_ras, nets_gt,
     nets_ras_remodel, nets_gt_remodel,
     membership_ras, membership_gt,
     name_ras, name_gt,
     layout_2d,
     file="network_comparison.RData")

