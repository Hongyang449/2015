## name: plot_signif_mu.r
## date: 09/14/2015

library(bio3d)
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_feature_ras.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna3.R")
source("/Users/hyangl/Desktop/2015/functions/get.signif.R")
source("/Users/hyangl/Desktop/2015/functions/remodel.nets.R")
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")

mu <- c(0,0.5,1)
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

for (i in 1:length(mu)) {
# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_pearson <- lapply(community_cij_md_ras_pearson, function(x) {
  get.signif(x[[1]],x[[2]],mu=mu[i])
  })
names(p_community_cij_md_ras_pearson) <- cutoff

# remodel nets
nets_md_ras_pearson_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_ras, signif=p_community_cij_md_ras_pearson[[x]],
    layout_2d=layout_2d[1:9,], p.cutoff=p.cutoff)
  })
names(nets_md_ras_pearson_remodel) <- cutoff

plot.nets(nets_md_ras_pearson_remodel[as.character(seq(0.3,0.6,by=0.1))], layout_2d=layout_2d)
mtext(paste0("Ras_MD_networks(pearson;signif;mu=", mu[i], ")"), line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file=paste0("figures/mu/ras_pearson_mu_", mu[i], ".pdf"))
}

for (i in 1:length(mu)) {
# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_lmi <- lapply(community_cij_md_ras_lmi, function(x) {
  get.signif(x[[1]],x[[2]],mu=mu[i])
  })
names(p_community_cij_md_ras_lmi) <- cutoff

# remodel nets
nets_md_ras_lmi_remodel <- lapply(as.list(1:length(cutoff)), function(x) {
  remodel.nets(cij=cij_lmi, cmap=cmap, cutoff.cij=cutoff[x], nets_dummy=nets_dummy,
    membership=membership_ras, signif=p_community_cij_md_ras_lmi[[x]],
    layout_2d=layout_2d[1:9,], p.cutoff=p.cutoff)
  })
names(nets_md_ras_lmi_remodel) <- cutoff

plot.nets(nets_md_ras_lmi_remodel[as.character(seq(0.3,0.6,by=0.1))], layout_2d=layout_2d)
mtext(paste0("Ras_MD_networks(lmi;signif;mu=", mu[i], ")"), line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file=paste0("figures/mu/ras_lmi_mu_", mu[i], ".pdf"))
}

