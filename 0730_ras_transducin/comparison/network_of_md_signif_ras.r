## name: network_of_md_signif_ras.r
## date: 09/05/2015

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna3.R")
source("/Users/hyangl/Desktop/2015/functions/get.signif.R")
source("/Users/hyangl/Desktop/2015/functions/get.community.cij.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.nets.R")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")

p.cutoff <- 0.05
cutoff <- seq(0.0,0.7,by=0.05)
layout_2d <- layout_2d[1:9,]

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

# prepare the dummy nets!
cij_gtp <- filter.dccm(array_cij_ca_pearson_gtp[,,1], cutoff.cij=0.6)
cij_gdp <- filter.dccm(array_cij_ca_pearson_gdp[,,1], cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp[1:166,1:166], cutoff.cij=0)
cna_gdp <- cna(cij_gdp[1:166,1:166], cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

cij_pearson <- list(gtp=array_cij_ca_pearson_gtp[1:166,1:166,],gdp=array_cij_ca_pearson_gdp[1:166,1:166,])
cij_lmi <- list(gtp=array_cij_ca_lmi_gtp[1:166,1:166,],gdp=array_cij_ca_lmi_gdp[1:166,1:166,])
cmap <- list(gtp=cmap_ca_gtp_dist10[1:166,1:166],gdp=cmap_ca_gdp_dist10[1:166,1:166])

########## pearson ###########

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_ras_pearson <- lapply(as.list(cutoff), function(x) {
  get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=x, membership=membership_ras)
  })
names(community_cij_md_ras_pearson) <- cutoff

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_pearson <- lapply(community_cij_md_ras_pearson, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_md_ras_pearson) <- cutoff

# remodel nets
nets_md_ras_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_ras, signif=p_community_cij_md_ras_pearson[[x]],
    layout_2d=layout_2d, p.cutoff=p.cutoff)
  })
names(nets_md_ras_pearson_remodel) <- cutoff

########## lmi ###########

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_ras_lmi <- lapply(as.list(cutoff), function(x) {
  get.community.cij(cij=cij_lmi, cmap=cmap, cutoff.cij=x, membership=membership_ras)
  })
names(community_cij_md_ras_lmi) <- cutoff

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_lmi <- lapply(community_cij_md_ras_lmi, function(x) {
  get.signif(x[[1]],x[[2]])
  })
names(p_community_cij_md_ras_lmi) <- cutoff

# remodel nets
nets_md_ras_lmi_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_lmi, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_ras, signif=p_community_cij_md_ras_lmi[[x]],
    layout_2d=layout_2d, p.cutoff=p.cutoff)
  })
names(nets_md_ras_lmi_remodel) <- cutoff

save(community_cij_md_ras_pearson, community_cij_md_ras_lmi,
     p_community_cij_md_ras_pearson, p_community_cij_md_ras_lmi,
     nets_md_ras_pearson_remodel, nets_md_ras_lmi_remodel,
     file="network_of_md_signif_ras.RData")


# plot!
library(bio3d)
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_signif_ras.RData")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")
layout_2d <- layout_2d[1:9,]

plot.nets(nets_md_ras_pearson_remodel[as.character(seq(0.2,0.35,by=0.05))], layout_2d=layout_2d)
mtext("Ras_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson_signif_0.2_0.35.pdf")

plot.nets(nets_md_ras_pearson_remodel[as.character(seq(0.4,0.55,by=0.05))], layout_2d=layout_2d)
mtext("Ras_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson_signif_0.4_0.55.pdf")

plot.nets(nets_md_ras_lmi_remodel[as.character(seq(0.3,0.45,by=0.05))], layout_2d=layout_2d)
mtext("Ras_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_lmi_signif_0.3_0.45.pdf")

plot.nets(nets_md_ras_lmi_remodel[as.character(seq(0.5,0.65,by=0.05))], layout_2d=layout_2d)
mtext("Ras_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_lmi_signif_0.5_0.65.pdf")

plot.nets(nets_md_ras_pearson_remodel[as.character(seq(0.3,0.6,by=0.1))], layout_2d=layout_2d)
mtext("Ras_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson_signif_0.3_0.6.pdf")

plot.nets(nets_md_ras_lmi_remodel[as.character(seq(0.3,0.6,by=0.1))], layout_2d=layout_2d)
mtext("Ras_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_lmi_signif_0.3_0.6.pdf")

# raw plot - no cutoff
require(igraph)
layout(matrix(1:2, nrow=1))
width = 0.2
w1 <- (E(nets_md_ras_pearson_remodel[[1]][[1]]$community.network)$weight) * width
w2 <- (E(nets_md_ras_pearson_remodel[[1]][[2]]$community.network)$weight) * width
plot.cna(nets_md_ras_pearson_remodel[[1]][[1]], layout=layout_2d, weights = w1, 
  vertex.label=NA, main="gtp_cutoff.cij=0")
plot.cna(nets_md_ras_pearson_remodel[[1]][[2]], layout=layout_2d, weights = w1,
  vertex.label=NA, main="gdp_cutoff.cij=0")
mtext("Ras_MD_networks(pearson;signif)", line=-3, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson_signif_0.pdf")









