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
source("/Users/hyangl/Desktop/2015/functions/plot.nets.R")

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

# prepare the dummy nets!
cij_gtp <- filter.dccm(array_cij_ca_pearson_gtp[,,1], cutoff.cij=0.6)
cij_gdp <- filter.dccm(array_cij_ca_pearson_gdp[,,1], cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp[1:166,1:166], cutoff.cij=0)
cna_gdp <- cna(cij_gdp[1:166,1:166], cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

######### pearson #########

cij_pearson <- list(array_cij_ca_pearson_gtp[1:166,1:166,],array_cij_ca_pearson_gdp[1:166,1:166,])
cmap <- list(cmap_ca_gtp_dist10[1:166,1:166],cmap_ca_gdp_dist10[1:166,1:166])

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_ras_pearson <- 
  lapply(as.list(seq(0.0,0.7,by=0.1)), function(x) {
    get.community.cij(cij=cij_pearson, cmap=cmap, cutoff.cij=x, membership=membership_ras)
    })

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_pearson <- 
  lapply(community_cij_md_ras_pearson, function(x) {
    get.signif(x[[1]],x[[2]])
    })

# plot and save remodel nets
nets_md_ras_pearson_remodel <- NULL
i = 1
layout(matrix(1:8,nrow=2))
for (j in seq(0.0,0.3,by=0.1)) {
  # calculate and plot networks
  nets_remodel <- plot.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=j, 
    nets_dummy=nets_dummy, membership=membership_ras, 
    signif=p_community_cij_md_ras_pearson[[i]], layout_2d=layout_2d[1:9,])

  # list and save networks
  nets_md_ras_pearson_remodel <- c(nets_md_ras_pearson_remodel,
    list(nets_remodel))

  i <- i+1
}
mtext("Ras_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson_signif_0.0_0.3.pdf")

layout(matrix(1:8,nrow=2))
for (j in seq(0.4,0.7,by=0.1)) {
  # calculate and plot networks
  nets_remodel <- plot.nets(cij=cij_pearson, cmap=cmap, cutoff.cij=j, 
    nets_dummy=nets_dummy, membership=membership_ras, 
    signif=p_community_cij_md_ras_pearson[[i]], layout_2d=layout_2d[1:9,])

  # list and save networks
  nets_md_ras_pearson_remodel <- c(nets_md_ras_pearson_remodel,
    list(nets_remodel))

  i <- i+1
}
mtext("Ras_MD_networks(pearson;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_pearson_signif_0.4_0.7.pdf")

######### lmi #########

cij_lmi <- list(array_cij_ca_lmi_gtp[1:166,1:166,],array_cij_ca_lmi_gdp[1:166,1:166,])
cmap <- list(cmap_ca_gtp_dist10[1:166,1:166],cmap_ca_gdp_dist10[1:166,1:166])

# community_cij: list 0.0-0.7 -> list gtp/gdp -> array 9*9*10
community_cij_md_ras_lmi <-
  lapply(as.list(seq(0.0,0.7,by=0.1)), function(x) {
    get.community.cij(cij=cij_lmi, cmap=cmap, cutoff.cij=x, membership=membership_ras)
    })

# p_community_cij: list 0.0-0.7 -> matrix p.value
p_community_cij_md_ras_lmi <-
  lapply(community_cij_md_ras_lmi, function(x) {
    get.signif(x[[1]],x[[2]])
    })

# plot and save remodel nets
nets_md_ras_lmi_remodel <- NULL
i = 1
layout(matrix(1:8,nrow=2))
for (j in seq(0.0,0.3,by=0.1)) {
  # calculate and plot networks
  nets_remodel <- plot.nets(cij=cij_lmi, cmap=cmap, cutoff.cij=j,
    nets_dummy=nets_dummy, membership=membership_ras,
    signif=p_community_cij_md_ras_lmi[[i]], layout_2d=layout_2d[1:9,])

  # list and save networks
  nets_md_ras_lmi_remodel <- c(nets_md_ras_lmi_remodel,
    list(nets_remodel))

  i <- i+1
}
mtext("Ras_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_lmi_signif_0.0_0.3.pdf")

layout(matrix(1:8,nrow=2))
for (j in seq(0.4,0.7,by=0.1)) {
  # calculate and plot networks
  nets_remodel <- plot.nets(cij=cij_lmi, cmap=cmap, cutoff.cij=j,
    nets_dummy=nets_dummy, membership=membership_ras,
    signif=p_community_cij_md_ras_lmi[[i]], layout_2d=layout_2d[1:9,])

  # list and save networks
  nets_md_ras_lmi_remodel <- c(nets_md_ras_lmi_remodel,
    list(nets_remodel))

  i <- i+1
}
mtext("Ras_MD_networks(lmi;signif)", line=-50, outer=TRUE, cex=1.5)
dev.copy2pdf(file="figures/cna_md_ras_lmi_signif_0.4_0.7.pdf")

save(community_cij_md_ras_pearson, community_cij_md_ras_lmi,
     p_community_cij_md_ras_pearson, p_community_cij_md_ras_lmi,
     nets_md_ras_pearson_remodel, nets_md_ras_lmi_remodel,
     file="network_of_md_signif_ras.RData")
  
