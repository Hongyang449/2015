## name: sum_edges_of_ttest.r
## date: 09/22/2015

library(bio3d)

load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/ttest_of_cij.RData")

source("/Users/hyangl/Desktop/2015/functions/get.unique.key.R")
source("/Users/hyangl/Desktop/2015/functions/add.community.rect.R")

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])
attach(transducin)
aaseq_ras <- paste0(ali["1QRA_A",ali["membership_ras",]!="0"], c(1:166))
aaseq_gt <- paste0(ali["1TND_A",ali["membership_gt",]!="0"],pdbs$resno[1,!is.gap(pdbs)][1:305])

cutoff = 0.01

##### ras #####
cij_pearson <- list(gtp=array_cij_ca_pearson_gtp[1:166,1:166,],gdp=array_cij_ca_pearson_gdp[1:166,1:166,])
cij_lmi <- list(gtp=array_cij_ca_lmi_gtp[1:166,1:166,],gdp=array_cij_ca_lmi_gdp[1:166,1:166,])

cij_filter <- matrix(0, nrow=dim(p_ras_pearson), ncol=dim(p_ras_pearson))
cij_filter[p_ras_pearson < cutoff] <- 1

sum_key_ras_all_pearson <- get.unique.key(cij_pearson, extra.filter=cij_filter,
  aaseq=aaseq_ras, count=FALSE)
sum_key_ras_inter_pearson <- get.unique.key(cij_pearson, extra.filter=cij_filter,
  aaseq=aaseq_ras, count=FALSE, membership=membership_ras, intra=FALSE)
sum_key_ras_intra_pearson <- get.unique.key(cij_pearson, extra.filter=cij_filter,
  aaseq=aaseq_ras, count=FALSE, membership=membership_ras, inter=FALSE)

##### gt #####
cij_pearson <- list(gtp=net.md$gtp$cijs0, gdp=net.md$gdp$cijs0)
cmap <- list(gtp=cm$GTP, gdp=cm$GDP)

cij_filter <- matrix(0, nrow=dim(p_gt_pearson), ncol=dim(p_gt_pearson))
cij_filter[p_gt_pearson < cutoff] <- 1

sum_key_gt_all_pearson <- get.unique.key(cij_pearson, extra.filter=cij_filter,
  aaseq=aaseq_gt, count=FALSE)
sum_key_gt_inter_pearson <- get.unique.key(cij_pearson, extra.filter=cij_filter,
  aaseq=aaseq_gt, count=FALSE, membership=membership_gt, intra=FALSE)
sum_key_gt_intra_pearson <- get.unique.key(cij_pearson, extra.filter=cij_filter,
  aaseq=aaseq_gt, count=FALSE, membership=membership_gt, inter=FALSE)

save(sum_key_ras_all_pearson, sum_key_ras_inter_pearson, sum_key_ras_intra_pearson,
     sum_key_gt_all_pearson, sum_key_gt_inter_pearson, sum_key_gt_intra_pearson,
     file="sum_edges_of_ttest.RData")

# plot!
#### ras ####

## all couplings
x <- sum_key_ras_all_pearson$sum
ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(ras;all;pearson)",
  xlab="residues", ylab="sum of edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=3), y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=3), y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_ras)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_edges_ras_all_pearson.pdf")

## inter couplings only
x <- sum_key_ras_inter_pearson$sum
# ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(ras;inter;pearson)",
  xlab="residues", ylab="sum of edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=3), y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=3), y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_ras)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_edges_ras_inter_pearson.pdf")

## intra couplings only
x <- sum_key_ras_intra_pearson$sum
# ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(ras;intra;pearson)",
  xlab="residues", ylab="sum of edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=3), y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=3), y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_ras)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_edges_ras_intra_pearson.pdf")

## diff couplings
x <- sum_key_ras_all_pearson$sum
x["gtp",] <- x["gtp",] - x["gdp",]
x["gdp",] <- -x["gtp",]
x[x<0] <- 0
# ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(ras;all;pearson;diff!)",
  xlab="residues", ylab="sum of edges")
points(-x["gdp",], type="h",col="green")
text(x=which(x["gtp",]>=1), y=x["gtp",which(x["gtp",]>=1)], labels=names(which(x["gtp",]>=1)),
  pos=3, col="red", cex=0.7)
text(x=which(x["gdp",]>=1), y=-x["gdp",which(x["gdp",]>=1)], labels=names(which(x["gdp",]>=1)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_ras)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_diff_edges_ras_all_pearson.pdf")

#### gt ####
membership_gt_339 <- rep(NA,339)
membership_gt_339[pdbs$resno[1,!is.gap(pdbs)]] <- membership_gt
y <- pdbs$resno[1,!is.gap(pdbs)][1:305]

## all couplings
x <- sum_key_gt_all_pearson$sum
ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(gt;all;pearson)",
  xlab="residues", ylab="sum of edges")
points(x=y, y=-x["gdp",], type="h",col="green")
text(x=y[which(x["gtp",]>=3)], y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
  pos=3, col="red", cex=0.7)
text(x=y[which(x["gdp",]>=3)], y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_gt_339)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_edges_gt_all_pearson.pdf")

## inter couplings
x <- sum_key_gt_inter_pearson$sum
# ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(gt;inter;pearson)",
  xlab="residues", ylab="sum of edges")
points(x=y, y=-x["gdp",], type="h",col="green")
text(x=y[which(x["gtp",]>=3)], y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
  pos=3, col="red", cex=0.7)
text(x=y[which(x["gdp",]>=3)], y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_gt_339)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_edges_gt_inter_pearson.pdf")

## intra couplings
x <- sum_key_gt_intra_pearson$sum
# ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(gt;intra;pearson)",
  xlab="residues", ylab="sum of edges")
points(x=y, y=-x["gdp",], type="h",col="green")
text(x=y[which(x["gtp",]>=3)], y=x["gtp",which(x["gtp",]>=3)], labels=names(which(x["gtp",]>=3)),
  pos=3, col="red", cex=0.7)
text(x=y[which(x["gdp",]>=3)], y=-x["gdp",which(x["gdp",]>=3)], labels=names(which(x["gdp",]>=3)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_gt_339)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_edges_gt_intra_pearson.pdf")

## all couplings
x <- sum_key_gt_all_pearson$sum
x["gtp",] <- x["gtp",] - x["gdp",]
x["gdp",] <- -x["gtp",]
x[x<0] <- 0
# ylim <- range(c(range(x["gtp",]), range(-x["gdp",]))) + c(-1,1)
plot(x=y, y=x["gtp",], type="h",col="red", ylim=ylim,
  main="sum_edges_of_ttest(gt;all;pearson;diff!)",
  xlab="residues", ylab="sum of edges")
points(x=y, y=-x["gdp",], type="h",col="green")
text(x=y[which(x["gtp",]>=1)], y=x["gtp",which(x["gtp",]>=1)], labels=names(which(x["gtp",]>=1)),
  pos=3, col="red", cex=0.7)
text(x=y[which(x["gdp",]>=1)], y=-x["gdp",which(x["gdp",]>=1)], labels=names(which(x["gdp",]>=1)),
  pos=1, col="green", cex=0.7)
add.community.rect(ylim, membership_gt_339)
legend("topright", legend=c("GTP","GDP"), col=c("red","green"), pch="|", cex=1)

dev.copy2pdf(file="figures/sum_diff_edges_gt_all_pearson.pdf")


