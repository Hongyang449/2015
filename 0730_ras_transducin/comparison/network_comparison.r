## name: network_comparison.r
## date: 08/18/2015

# use sd instead of var to calculate edge colors between networks
# color the edge "separately"

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/gt_networks.RData")
load("/Users/hyangl/Desktop/0105_manuscript/dccm/cij_consensus_2grps.RData")

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout.RData")
load("/Users/hyangl/Desktop/0105_manuscript/pca/sse_for_plot.RData")

source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/draw.box.R")
source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/label_ras.R")

pdb_ras <- read.pdb("/Users/hyangl/Desktop/0105_manuscript/dccm/pdb/gtp_plot.pdb")

# cna; 1:166
cna_ras_gtp <- cna(lmi_consensus_gtp[1:166,1:166], cutoff.cij=0)
cna_ras_gdp <- cna(lmi_consensus_gdp[1:166,1:166], cutoff.cij=0)

nets_ras <- list(gtp=cna_ras_gtp,gdp=cna_ras_gdp)
nets_gt <- nets_gt_3grps[1:2]

#### ras 9 grps, gt 11 grps

## manually input membership in ras/gt
membership_ras <- rep(0,166)
membership_ras[c(1:9,18:24,38:56)] <- 1
membership_ras[10:17] <- 2
membership_ras[25:37] <- 3
membership_ras[57:75] <- 4
membership_ras[c(76:85,110:121,140:153)] <- 5
membership_ras[86:109] <- 6
membership_ras[122:130] <- 9
membership_ras[131:139] <- 7
membership_ras[154:166] <- 8

name_ras <- c("b1-b3,a1","p-loop","SI","SII","b4-b6","a3","a4","a5","loop8")
names(name_ras) <- 1:length(name_ras)

# membership_gt only cover residues 31-339; residues 108,109,308,313 are not in the 305*305 cij_gt matrix

membership_gt <- rep(0,339)
membership_gt[c(31:35,44:53,180:195)] <- 1
membership_gt[36:43] <- 2
membership_gt[c(54:72,148:172)] <- 10
membership_gt[c(73:107,110:147)] <- 11
membership_gt[173:179] <- 3
membership_gt[196:214] <- 4
membership_gt[c(215:226,259:274,316:329)] <- 5
membership_gt[227:258] <- 6
membership_gt[275:299] <- 9
membership_gt[c(300:307,309:312,314:315)] <- 7
membership_gt[330:339] <- 8

name_gt <- c("b1-b3,a1","p-loop","SI","SII","b4-b6", "SIII,a3","a4","a5","loop8","H1","H2")
names(name_gt) <- 1:length(name_gt)

layout_ras <- layout[1:9,]
layout_ras[4,] <- layout_ras[4,] - c(3,0)
layout_gt <- layout
              
source("/Users/hyangl/Desktop/0730_ras_transducin/transducin/remodel.cna.03.R")

nets_ras_remodel <- remodel.cna(nets_ras, member=membership_ras, method="mean", ne=10, scut=4)
nets_gt_remodel <- remodel.cna(nets_gt, member=membership_gt[membership_gt!=0], method="mean", ne=10, scut=4)

w1_ras = exp(-E(nets_ras_remodel[[1]]$community.network)$weight)*20
w2_ras = exp(-E(nets_ras_remodel[[2]]$community.network)$weight)*20
w1_gt = exp(-E(nets_gt_remodel[[1]]$community.network)$weight)*20
w2_gt = exp(-E(nets_gt_remodel[[2]]$community.network)$weight)*20


par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_ras_remodel[[1]], layout=layout_ras, weights = w1_ras, vertex.label=NA, main="ras_gtp")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
plot.cna(nets_ras_remodel[[2]], layout=layout_ras, weights = w2_ras, vertex.label=NA, main="ras_gdp")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_ras.pdf")

par(mfrow=c(1,2), mar = c(2,2,4,2), pty="s")
plot.cna(nets_gt_remodel[[1]], layout=layout_gt, weights = w1_gt, vertex.label=NA, main="transducin_gtp")
for (i in 1:length(name_gt)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_gt[i], col="black", cex=1)
}
plot.cna(nets_gt_remodel[[2]], layout=layout_gt, weights = w2_gt, vertex.label=NA, main="transducin_gdp")
for (i in 1:length(name_gt)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_gt[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_transducin.pdf")

save(nets_ras_remodel, nets_gt_remodel,
     membership_ras, membership_gt,
     name_ras, name_gt,
     layout_ras, layout_gt,
     file="network_comparison.RData")

