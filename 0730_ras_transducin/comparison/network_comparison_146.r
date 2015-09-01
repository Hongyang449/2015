## name: network_comparison_146.r
## date: 08/18/2015

# only use the aligned 146 residues to build networks.
# use sd instead of var to calculate edge colors between networks


load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/gt_networks.RData")
load("/Users/hyangl/Desktop/0105_manuscript/dccm/cij_consensus_2grps.RData")

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0105_manuscript/pca/sse_for_plot.RData")

source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/draw.box.R")
source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/label_ras.R")

pdb_ras <- read.pdb("/Users/hyangl/Desktop/0105_manuscript/dccm/pdb/gtp_plot.pdb")

# cij
inds1 <- as.numeric(ali["resno_ras",ali["aligned",]==1])

inds2 <- as.numeric(ali["resno_Gt",ali["aligned",]==1])
inds2 <- which(as.numeric(ali["resno_Gt",ali["membership_xinqiu",]!="0"]) %in% inds2)

cij_ras_gtp_146 <- matrix(-1,nrow=169,ncol=169)
cij_ras_gdp_146 <- matrix(-1,nrow=169,ncol=169)
cij_gt_gtp_146 <- matrix(-1,nrow=169,ncol=169)
cij_gt_gdp_146 <- matrix(-1,nrow=169,ncol=169)

cij_ras_gtp_146[inds1,inds1] <- lmi_consensus_gtp[inds1,inds1]
cij_ras_gdp_146[inds1,inds1] <- lmi_consensus_gdp[inds1,inds1]

cij <- nets_gt$gtp$cij[inds2,inds2]
cij[cij>0] = exp(-cij[cij>0])
cij_gt_gtp_146[inds1,inds1] <- cij
cij <- nets_gt$gdp$cij[inds2,inds2]
cij[cij>0] = exp(-cij[cij>0])
cij_gt_gdp_146[inds1,inds1] <- cij

# 2D cij plot: gtp vs gdp
cij_ras_gtp_vs_gdp_146 <- cij_ras_gtp_146
cij_ras_gtp_vs_gdp_146[lower.tri(cij_ras_gtp_vs_gdp_146)] <- cij_ras_gdp_146[lower.tri(cij_ras_gdp_146)]
cij_gt_gtp_vs_gdp_146 <- cij_gt_gtp_146
cij_gt_gtp_vs_gdp_146[lower.tri(cij_gt_gtp_vs_gdp_146)] <- cij_gt_gdp_146[lower.tri(cij_gt_gdp_146)]
which(abs(cij_ras_gtp_146-cij_ras_gdp_146) > 0.2)

x11(type="cairo")
plot.new()
plot.dccm(cij_ras_gtp_vs_gdp_146,helix.col="gray20",sheet.col="gray80",
          scales=list(x=list(at=label_ab, cex=0.8, rot=90), y=list(at=label_ab, cex=0.8)),
          main="Ras",sse=sse, at=c(-1,-0.7,-0.4,0.4,0.7,1))
draw.box(sse)
label_ras(c("GTP","GDP"))
savePlot(filename="figures/cij_ras_gtp_vs_gdp_146.png", type="png")

x11(type="cairo")
plot.new()
plot.dccm(cij_gt_gtp_vs_gdp_146,helix.col="gray20",sheet.col="gray80",
          scales=list(x=list(at=label_ab, cex=0.8, rot=90), y=list(at=label_ab, cex=0.8)),
          main="Transducin",sse=sse, at=c(-1,-0.7,-0.4,0.4,0.7,1))
draw.box(sse)
label_ras(c("GTP","GDP"))
savePlot(filename="figures/cij_gt_gtp_vs_gdp_146.png", type="png")

# cna
cna_ras_gtp_146 <- cna(cij_ras_gtp_146[inds1,inds1], cutoff.cij=0)
cna_ras_gdp_146 <- cna(cij_ras_gdp_146[inds1,inds1], cutoff.cij=0)

cna_gt_gtp_146 <- cna(cij_gt_gtp_146[inds1,inds1], cutoff.cij=0)
cna_gt_gdp_146 <- cna(cij_gt_gdp_146[inds1,inds1], cutoff.cij=0)

nets_ras_146 <- list(gtp=cna_ras_gtp_146,gdp=cna_ras_gdp_146)
nets_gt_146 <- list(gtp=cna_gt_gtp_146,gdp=cna_gt_gdp_146)

# membership
membership_ras_146 <- rep(0,166)
membership_ras_146[c(5:9,18:26,29:30,38:46,50:56)] <- 1
membership_ras_146[10:17] <- 2
membership_ras_146[c(31:37)] <- 3
membership_ras_146[c(57:75)] <- 4
membership_ras_146[c(76:85,110:121,140:153)] <- 5
membership_ras_146[c(86:103,108:109)] <- 6
membership_ras_146[c(126:139)] <- 7
membership_ras_146[c(154:163)] <- 8

membership_146 <- membership_ras_146[membership_ras_146!=0]

name_ras <- c("b1-b3,a1","p-loop","SI","SII","b4-b6","a3","a4","a5")
names(name_ras) <- 1:length(name_ras)

nets_ras_146_remodel <- remodel.cna(nets_ras_146, member = membership_146, method="mean", ne=10, scut=4)
nets_gt_146_remodel <- remodel.cna(nets_gt_146, member = membership_146, method="mean", ne=10, scut=4)

w1_ras = exp(-E(nets_ras_146_remodel[[1]]$community.network)$weight)*20
w2_ras = exp(-E(nets_ras_146_remodel[[2]]$community.network)$weight)*20
w1_gt = exp(-E(nets_gt_146_remodel[[1]]$community.network)$weight)*20
w2_gt = exp(-E(nets_gt_146_remodel[[2]]$community.network)$weight)*20

layout_146 <- layout[1:8,]
layout_146[4,] <- layout_146[4,] - c(3,0)

par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_ras_146_remodel[[1]], layout=layout_146, weights = w1_ras, vertex.label=NA, main="ras_gtp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
plot.cna(nets_ras_146_remodel[[2]], layout=layout_146, weights = w2_ras, vertex.label=NA, main="ras_gdp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_ras_146.pdf")

par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_gt_146_remodel[[1]], layout=layout_146, weights = w1_gt, vertex.label=NA, main="transducin_gtp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
plot.cna(nets_gt_146_remodel[[2]], layout=layout_146, weights = w2_gt, vertex.label=NA, main="transducin_gdp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_transducin_146.pdf")

save(inds1, inds2, 
     cij_ras_gtp_146, cij_ras_gdp_146,
     cij_gt_gtp_146, cij_gt_gdp_146,
     nets_ras_146_remodel, nets_gt_146_remodel,
     membership_ras_146, membership_146, layout_146,
     file="network_comparison_146.RData")

