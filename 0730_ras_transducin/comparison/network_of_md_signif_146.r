## name: network_of_md_signif_146.r
## date: 09/05/2015

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna3.R")
source("/Users/hyangl/Desktop/2015/functions/get.signif.R")
source("/Users/hyangl/Desktop/2015/functions/get.community.cij.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.nets.R")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")

p.cutoff <- 0.05
cutoff <- seq(0.0,0.7,by=0.05)
layout <- layout_2d
layout_2d <- layout[1:9,]

membership_146 <- as.numeric(ali["membership_146",ali["membership_146",]!="0"])

inds1 <- as.numeric(ali["resno_ras",ali["membership_146",]!="0"])
inds2 <- as.numeric(ali["resno_gt",ali["membership_146",]!="0"])
inds2 <- which(as.numeric(ali["resno_gt",ali["membership_gt",]!="0"]) %in% inds2)

# prepare the dummy nets!
cij_gtp <- filter.dccm(array_cij_ca_pearson_gtp[inds1,inds1,1], cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp, cutoff.cij=0)
cna_gdp <- cna_gtp
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

cij_pearson <- list(gtp=array_cij_ca_pearson_gtp[inds1,inds1,],gdp=array_cij_ca_pearson_gdp[inds1,inds1,])
cij_lmi <- list(gtp=array_cij_ca_lmi_gtp[inds1,inds1,],gdp=array_cij_ca_lmi_gdp[inds1,inds1,])
cmap <- list(gtp=cmap_ca_gtp_dist10[inds1,inds1],gdp=cmap_ca_gdp_dist10[inds1,inds1])

########## pearson ###########

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_ras_146_pearson <- lapply(as.list(cutoff), function(x) {
  get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=x, membership=membership_146)
  })
names(community_cij_md_ras_146_pearson) <- cutoff

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_146_pearson <- lapply(community_cij_md_ras_146_pearson, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_md_ras_146_pearson) <- cutoff

# remodel nets
nets_md_ras_146_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_146, signif=p_community_cij_md_ras_146_pearson[[x]],
    layout_2d=layout_2d, p.cutoff=p.cutoff)
  })
names(nets_md_ras_146_pearson_remodel) <- cutoff

########## lmi ###########

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_ras_146_lmi <- lapply(as.list(cutoff), function(x) {
  get.community.cij(cij=cij_lmi, cmap=cmap, cutoff.cij=x, membership=membership_146)
  })
names(community_cij_md_ras_146_lmi) <- cutoff

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_146_lmi <- lapply(community_cij_md_ras_146_lmi, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_md_ras_146_lmi) <- cutoff

# remodel nets
nets_md_ras_146_lmi_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_lmi, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_146, signif=p_community_cij_md_ras_146_lmi[[x]],
    layout_2d=layout_2d, p.cutoff=p.cutoff)
  })
names(nets_md_ras_146_lmi_remodel) <- cutoff

############### Gt ################
layout_2d <- layout

cij_gtp <- filter.dccm(net.md$gtp$cijs0[inds2,inds2,], cmap=cm$GTP[inds2,inds2], cutoff.cij=0.6)
cij_gdp <- filter.dccm(net.md$gdp$cijs0[inds2,inds2,], cmap=cm$GDP[inds2,inds2], cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp, cutoff.cij=0)
cna_gdp <- cna(cij_gdp, cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

cij_pearson <- list(gtp=net.md$gtp$cijs0[inds2,inds2,], gdp=net.md$gdp$cijs0[inds2,inds2,])
cmap <- list(gtp=cm$GTP[inds2,inds2], gdp=cm$GDP[inds2,inds2])

########## pearson ###########

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_gt_146_pearson <- lapply(as.list(cutoff), function(x) {
  get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=x, membership=membership_146)
  })
names(community_cij_md_gt_146_pearson) <- cutoff

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_gt_146_pearson <- lapply(community_cij_md_gt_146_pearson, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_md_gt_146_pearson) <- cutoff

# remodel nets
nets_md_gt_146_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_146, signif=p_community_cij_md_gt_146_pearson[[x]],
    layout_2d=layout_2d, p.cutoff=p.cutoff)
  })
names(nets_md_gt_146_pearson_remodel) <- cutoff




save(community_cij_md_ras_146_pearson, community_cij_md_ras_146_lmi,
     p_community_cij_md_ras_146_pearson, p_community_cij_md_ras_146_lmi,
     nets_md_ras_146_pearson_remodel, nets_md_ras_146_lmi_remodel,
     community_cij_md_gt_146_pearson,
     p_community_cij_md_gt_146_pearson,
     nets_md_gt_146_pearson_remodel,
     file="network_of_md_signif_146.RData")


# plot!
layout_2d <- layout[1:8,]
plot.nets(nets_md_ras_146_pearson_remodel[as.character(seq(0.2,0.35,by=0.05))], layout_2d=layout_2d)
mtext("Ras_146_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_146_pearson_signif_0.2_0.35.pdf")

plot.nets(nets_md_ras_146_pearson_remodel[as.character(seq(0.4,0.55,by=0.05))], layout_2d=layout_2d)
mtext("Ras_146_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_146_pearson_signif_0.4_0.55.pdf")

plot.nets(nets_md_ras_146_lmi_remodel[as.character(seq(0.3,0.45,by=0.05))], layout_2d=layout_2d)
mtext("Ras_146_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_146_lmi_signif_0.3_0.45.pdf")

plot.nets(nets_md_ras_146_lmi_remodel[as.character(seq(0.5,0.65,by=0.05))], layout_2d=layout_2d)
mtext("Ras_146_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_146_lmi_signif_0.5_0.65.pdf")

plot.nets(nets_md_gt_146_pearson_remodel[as.character(seq(0.4,0.6,by=0.05))], layout_2d=layout_2d)
mtext("Gt_146_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_gt_146_pearson_signif_0.4_0.6.pdf")





