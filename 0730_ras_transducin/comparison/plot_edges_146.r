## name: plot_edge_146.r
## date: 10/01/2015

# generate the plot of 146 aligned residues

library(bio3d)

load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/unique_edges.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/weighted_edges.RData")

source("/Users/hyangl/Desktop/2015/functions/add.community.rect.R")

membership_146 <- as.numeric(ali["membership_146",ali["membership_146",]!="0"])

inds1 <- as.numeric(ali["resno_ras",ali["membership_146",]!="0"])
inds2 <- as.numeric(ali["resno_gt",ali["membership_146",]!="0"])
inds2 <- which(as.numeric(ali["resno_gt",ali["membership_gt",]!="0"]) %in% inds2)

############ ras ##########
cutoff <- seq(0.3,0.6,by=0.1)
i <- 2

## unique edges
x <- num_key_ras_all_pearson[[i]]$num[,inds1]
ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("unique_edges(ras_146;all;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="number of unique edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/num_unique_edges_ras_146_all_pearson_",cutoff[i],".pdf"))

x <- num_key_ras_inter_pearson[[i]]$num[,inds1]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("unique_edges(ras_146;inter;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="number of unique edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/num_unique_edges_ras_146_inter_pearson_",cutoff[i],".pdf"))

x <- num_key_ras_intra_pearson[[i]]$num[,inds1]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("unique_edges(ras_146;intra;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="number of unique edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/num_unique_edges_ras_146_intra_pearson_",cutoff[i],".pdf"))


## sum weighted edges
x <- weighted_key_ras_all_pearson[[i]]$sum[,inds1]
ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("weighted_edges(ras_146;all;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="sum weighted edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/weighted_edges_ras_146_all_pearson_",cutoff[i],".pdf"))

x <- weighted_key_ras_inter_pearson[[i]]$sum[,inds1]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("weighted_edges(ras_146;inter;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="sum weighted edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/weighted_edges_ras_146_inter_pearson_",cutoff[i],".pdf"))

x <- weighted_key_ras_intra_pearson[[i]]$sum[,inds1]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("weighted_edges(ras_146;intra;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="sum weighted edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/weighted_edges_ras_146_intra_pearson_",cutoff[i],".pdf"))



############ gt ##########
cutoff <- seq(0.4,0.7,by=0.1)
i <- 3

## unique edges
x <- num_key_gt_all_pearson[[i]]$num[,inds2]
ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("unique_edges(gt_146;all;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="number of unique edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/num_unique_edges_gt_146_all_pearson_",cutoff[i],".pdf"))

x <- num_key_gt_inter_pearson[[i]]$num[,inds2]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("unique_edges(gt_146;inter;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="number of unique edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/num_unique_edges_gt_146_inter_pearson_",cutoff[i],".pdf"))

x <- num_key_gt_intra_pearson[[i]]$num[,inds2]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("unique_edges(gt_146;intra;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="number of unique edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/num_unique_edges_gt_146_intra_pearson_",cutoff[i],".pdf"))

## sum weighted edges
x <- weighted_key_gt_all_pearson[[i]]$sum[,inds2]
ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("weighted_edges(gt_146;all;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="sum weighted edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/weighted_edges_gt_146_all_pearson_",cutoff[i],".pdf"))

x <- weighted_key_gt_inter_pearson[[i]]$sum[,inds2]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("weighted_edges(gt_146;inter;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="sum weighted edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/weighted_edges_gt_146_inter_pearson_",cutoff[i],".pdf"))

x <- weighted_key_gt_intra_pearson[[i]]$sum[,inds2]
#ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main=paste0("weighted_edges(gt_146;intra;pearson;cutoff.cij=",cutoff[i],")"),
  xlab="residues", ylab="sum weighted edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=4), y=x["gtp",which(x["gtp",]>=4)], labels=names(which(x["gtp",]>=4)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=4), y=-x["gdp",which(x["gdp",]>=4)], labels=names(which(x["gdp",]>=4)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_146)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)
dev.copy2pdf(file=paste0("figures/weighted_edges_gt_146_intra_pearson_",cutoff[i],".pdf"))

