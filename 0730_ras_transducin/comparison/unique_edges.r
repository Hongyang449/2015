## name: unique_edges.r
## date: 09/21/2015

library(bio3d)

load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")

source("/Users/hyangl/Desktop/2015/functions/get.unique.key.R")

load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_signif_ras.RData")

# add the "sse" colored by community colors
add.sse <- function(ylim,membership) {
  for (i in na.omit(unique(membership))) {
    pos <- bounds(which(membership == i))

    bo <- max(ylim) + (diff(ylim) * 0.001)
    to <- max(ylim) + (diff(ylim) * 0.04)
    rect(xleft=pos[,"start"], xright=pos[,"end"], ybottom=bo, ytop=to, col=vmd.colors()[i])
    bo <- min(ylim) - (diff(ylim) * 0.001)
    to <- min(ylim) - (diff(ylim) * 0.04)
    rect(xleft=pos[,"start"], xright=pos[,"end"], ybottom=bo, ytop=to, col=vmd.colors()[i])
  }
}

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])

attach(transducin)
aaseq_ras <- paste0(ali["1QRA_A",ali["membership_ras",]!="0"], c(1:166))
aaseq_gt <- paste0(ali["1TND_A",ali["membership_gt",]!="0"],pdbs$resno[1,!is.gap(pdbs)][1:305])

########### ras #########
cij_pearson <- list(gtp=array_cij_ca_pearson_gtp[1:166,1:166,],gdp=array_cij_ca_pearson_gdp[1:166,1:166,])
cij_lmi <- list(gtp=array_cij_ca_lmi_gtp[1:166,1:166,],gdp=array_cij_ca_lmi_gdp[1:166,1:166,])
cmap <- list(gtp=cmap_ca_gtp_dist10[1:166,1:166],gdp=cmap_ca_gdp_dist10[1:166,1:166])


cutoff <- seq(0.3,0.6,by=0.1)
# all couplings
num_key_ras_all_pearson <- lapply(as.list(cutoff), function(x) {
  get.unique.key(cij=cij_pearson, cmap=cmap, cutoff.cij=x, aaseq=aaseq_ras)
  })
names(num_key_ras_all_pearson) <- cutoff
# only inter couplings - remove intra couplings
num_key_ras_inter_pearson <- lapply(as.list(cutoff), function(x) {
  get.unique.key(cij=cij_pearson, cmap=cmap, cutoff.cij=x, 
    membership=membership_ras, intra=FALSE, aaseq=aaseq_ras)
  })
names(num_key_ras_inter_pearson) <- cutoff
# only intra couplings - remove inter couplings
num_key_ras_intra_pearson <- lapply(as.list(cutoff), function(x) {
  get.unique.key(cij=cij_pearson, cmap=cmap, cutoff.cij=x, 
    membership=membership_ras, inter=FALSE, aaseq=aaseq_ras)
  })
names(num_key_ras_intra_pearson) <- cutoff

########### gt #########
cij_pearson <- list(gtp=net.md$gtp$cijs0, gdp=net.md$gdp$cijs0)
cmap <- list(gtp=cm$GTP, gdp=cm$GDP)

cutoff <- seq(0.4,0.7,by=0.1)
# all couplings
num_key_gt_all_pearson <- lapply(as.list(cutoff), function(x) {
  get.unique.key(cij=cij_pearson, cmap=cmap, cutoff.cij=x, aaseq=aaseq_gt)
  })
names(num_key_gt_all_pearson) <- cutoff
# only inter couplings - remove intra couplings
num_key_gt_inter_pearson <- lapply(as.list(cutoff), function(x) {
  get.unique.key(cij=cij_pearson, cmap=cmap, cutoff.cij=x,
    membership=membership_gt, intra=FALSE, aaseq=aaseq_gt)
  })
names(num_key_gt_inter_pearson) <- cutoff
# only intra couplings - remove inter couplings
num_key_gt_intra_pearson <- lapply(as.list(cutoff), function(x) {
  get.unique.key(cij=cij_pearson, cmap=cmap, cutoff.cij=x,
    membership=membership_gt, inter=FALSE, aaseq=aaseq_gt)
  })
names(num_key_gt_intra_pearson) <- cutoff

save(num_key_ras_all_pearson, num_key_ras_inter_pearson, num_key_ras_intra_pearson,
     num_key_gt_all_pearson, num_key_gt_inter_pearson, num_key_gt_intra_pearson,
     file="unique_edges.RData")

# plot!
cutoff <- seq(0.3,0.6,by=0.1)
for (i in 2) {
  x <- num_key_ras_all_pearson[[i]]$num
  ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
  plot(x["gtp",], type="h",col="red", ylim=ylim,
    main=paste0("unique_edges(ras;all;pearson;cutoff.cij=",cutoff[i],")"),
    xlab="residues", ylab="number of unique edges")
  points(-x["gdp",], type="h",col="green")
  text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
    pos=3, col="red", cex=0.7)
  text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
    pos=1, col="green", cex=0.7)
  add.sse(ylim, membership_ras)
  legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

  dev.copy2pdf(file=paste0("figures/num_unique_edges_ras_all_pearson_",cutoff[i],".pdf"))
}

for (i in 2) {
  x <- num_key_ras_inter_pearson[[i]]$num
#  ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
  plot(x["gtp",], type="h",col="red", ylim=ylim,
    main=paste0("unique_edges(ras;inter;pearson;cutoff.cij=",cutoff[i],")"),
    xlab="residues", ylab="number of unique edges")
  points(-x["gdp",], type="h",col="green")
  text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
    pos=3, col="red", cex=0.7)
  text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
    pos=1, col="green", cex=0.7)
  add.sse(ylim, membership_ras)
  legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

  dev.copy2pdf(file=paste0("figures/num_unique_edges_ras_inter_pearson_",cutoff[i],".pdf"))
}

for (i in 2) {
  x <- num_key_ras_intra_pearson[[i]]$num
#  ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
  plot(x["gtp",], type="h",col="red", ylim=ylim,
    main=paste0("unique_edges(ras;intra;pearson;cutoff.cij=",cutoff[i],")"),
    xlab="residues", ylab="number of unique edges")
  points(-x["gdp",], type="h",col="green")
  text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
    pos=3, col="red", cex=0.7)
  text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
    pos=1, col="green", cex=0.7)
  add.sse(ylim, membership_ras)
  legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

  dev.copy2pdf(file=paste0("figures/num_unique_edges_ras_intra_pearson_",cutoff[i],".pdf"))
}

## Gt
membership_gt_339 <- rep(NA,339)
membership_gt_339[pdbs$resno[1,!is.gap(pdbs)]] <- membership_gt
cutoff <- seq(0.4,0.7,by=0.1)

for (i in 3) {
  x <- num_key_gt_all_pearson[[i]]$num
  y <- pdbs$resno[1,!is.gap(pdbs)][1:305]
  ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)

  plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
    main=paste0("unique_edges(gt;all;pearson;cutoff.cij=",cutoff[i],")"),
    xlab="residues", ylab="number of unique edges")
  points(x=y, y=-x["gdp",], type="h",col="green")
  text(x=y[which(x["gtp",]>=3)], y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
    pos=3, col="red", cex=0.7)
  text(x=y[which(x["gdp",]>=3)], y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
    pos=1, col="green", cex=0.7)
  add.sse(ylim, membership_gt_339)
  legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

  dev.copy2pdf(file=paste0("figures/num_unique_edges_gt_all_pearson_",cutoff[i],".pdf"))
}

for (i in 3) {
  x <- num_key_gt_inter_pearson[[i]]$num
  y <- pdbs$resno[1,!is.gap(pdbs)][1:305]
#  ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)

  plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
    main=paste0("unique_edges(gt;inter;pearson;cutoff.cij=",cutoff[i],")"),
    xlab="residues", ylab="number of unique edges")
  points(x=y, y=-x["gdp",], type="h",col="green")
  text(x=y[which(x["gtp",]>=3)], y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
    pos=3, col="red", cex=0.7)
  text(x=y[which(x["gdp",]>=3)], y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
    pos=1, col="green", cex=0.7)
  add.sse(ylim, membership_gt_339)
  legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

  dev.copy2pdf(file=paste0("figures/num_unique_edges_gt_inter_pearson_",cutoff[i],".pdf"))
}

for (i in 3) {
  x <- num_key_gt_intra_pearson[[i]]$num
  y <- pdbs$resno[1,!is.gap(pdbs)][1:305]
#  ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)

  plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
    main=paste0("unique_edges(gt;intra;pearson;cutoff.cij=",cutoff[i],")"),
    xlab="residues", ylab="number of unique edges")
  points(x=y, y=-x["gdp",], type="h",col="green")
  text(x=y[which(x["gtp",]>=3)], y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
    pos=3, col="red", cex=0.7)
  text(x=y[which(x["gdp",]>=3)], y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
    pos=1, col="green", cex=0.7)
  add.sse(ylim, membership_gt_339)
  legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

  dev.copy2pdf(file=paste0("figures/num_unique_edges_gt_intra_pearson_",cutoff[i],".pdf"))
}





