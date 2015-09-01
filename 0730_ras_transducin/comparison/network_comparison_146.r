## name: network_comparison_146.r
## date: 08/18/2015

# only use the aligned 146 residues to build networks.
# use sd instead of var to calculate edge colors between networks


load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/gt_networks.RData")
load("/Users/hyangl/Desktop/0105_manuscript/dccm/cij_consensus_2grps.RData")

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/0105_manuscript/pca/sse_for_plot.RData")

source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/draw.box.R")
source("/Users/hyangl/Desktop/0105_manuscript/dccm/functions/label_ras.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.cna2.R")

pdb_ras <- read.pdb("/Users/hyangl/Desktop/0105_manuscript/dccm/pdb/gtp_plot.pdb")

# cij
inds1 <- as.numeric(ali["resno_ras",ali["membership_146",]!="0"])

inds2 <- as.numeric(ali["resno_gt",ali["membership_146",]!="0"])
inds2 <- which(as.numeric(ali["resno_gt",ali["membership_gt",]!="0"]) %in% inds2)

cij_ras_gtp_146 <- matrix(-1,nrow=169,ncol=169)
cij_ras_gdp_146 <- matrix(-1,nrow=169,ncol=169)
cij_gt_gtp_146 <- matrix(-1,nrow=169,ncol=169)
cij_gt_gdp_146 <- matrix(-1,nrow=169,ncol=169)

cij_ras_gtp_146[inds1,inds1] <- lmi_consensus_gtp[inds1,inds1]
cij_ras_gdp_146[inds1,inds1] <- lmi_consensus_gdp[inds1,inds1]

cij <- nets_gt_3grps$gtp$cij[inds2,inds2]
cij[cij>0] = exp(-cij[cij>0])
cij_gt_gtp_146[inds1,inds1] <- cij
cij <- nets_gt_3grps$gdp$cij[inds2,inds2]
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
          main="Ras",sse=sse)
draw.box(sse)
label_ras(c("GTP","GDP"))
savePlot(filename="figures/cij_ras_gtp_vs_gdp_146.png", type="png")

x11(type="cairo")
plot.new()
plot.dccm(cij_gt_gtp_vs_gdp_146,helix.col="gray20",sheet.col="gray80",
          scales=list(x=list(at=label_ab, cex=0.8, rot=90), y=list(at=label_ab, cex=0.8)),
          main="Transducin",sse=sse)
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
membership_146 <- as.numeric(ali["membership_146",ali["membership_146",]!="0"])

name_ras <- c("b1-b3,a1","p-loop","SI","SII","b4-b6","a3","a4","a5")
names(name_ras) <- 1:length(name_ras)

nets_ras_146_remodel <- remodel.cna(nets_ras_146, member = membership_146, method="sum", col.edge="feature", scut=4)
nets_gt_146_remodel <- remodel.cna(nets_gt_146, member = membership_146, method="sum", col.edge="feature", scut=4)

w1_ras = exp(-E(nets_ras_146_remodel[[1]]$community.network)$weight)*20
w2_ras = exp(-E(nets_ras_146_remodel[[2]]$community.network)$weight)*20
w1_gt = exp(-E(nets_gt_146_remodel[[1]]$community.network)$weight)*20
w2_gt = exp(-E(nets_gt_146_remodel[[2]]$community.network)$weight)*20

par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_ras_146_remodel[[1]], layout=layout_2d[1:8,], weights = w1_ras, vertex.label=NA, main="ras_gtp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
plot.cna(nets_ras_146_remodel[[2]], layout=layout_2d[1:8,], weights = w2_ras, vertex.label=NA, main="ras_gdp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_ras_146.pdf")

par(mfrow=c(1,2), mar = c(6,6,12,6), pty="s")
plot.cna(nets_gt_146_remodel[[1]], layout=layout_2d[1:8,], weights = w1_gt, vertex.label=NA, main="transducin_gtp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
plot.cna(nets_gt_146_remodel[[2]], layout=layout_2d[1:8,], weights = w2_gt, vertex.label=NA, main="transducin_gdp_146")
for (i in 1:length(name_ras)) {
  x=locator(n=1)
  text(x=x$x[1], y=x$y[1], name_ras[i], col="black", cex=1)
}
dev.copy2pdf(file="figures/2D_cna_transducin_146.pdf")

save(inds1, inds2, 
     cij_ras_gtp_146, cij_ras_gdp_146,
     cij_gt_gtp_146, cij_gt_gdp_146,
     nets_ras_146_remodel, nets_gt_146_remodel,
     ali, membership_146, layout_2d,
     file="network_comparison_146.RData")

