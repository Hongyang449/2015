## name: key_edges_md.r
## date: 09/08/0215

# key_pearson: list 0.0-0.7 -> list gtp/gdp -> key edges matrix n*4

library(bio3d)

load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_signif_ras.RData")
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_signif_gt.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cij.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")

source("/Users/hyangl/Desktop/2015/functions/key.extract.R")
source("/Users/hyangl/Desktop/2015/functions/key.reorder.R")

dim_ras <- 166
dim_gt <- 305

# attach pdbs of gt to renumber the keys (the original number is for the 305 non-gap positions)

key.renumber <- function(key) {
  attach(transducin)
  lapply(key, function(x) {
    lapply(x, function(y) {
      if(is.matrix(y)) { y[,1:2] <- pdbs$resno[1,!is.gap(pdbs)][as.numeric(y[,1:2])] }
      return(y)
      })
    })
}

# add mean, sd, aa names to the key matrix
key.addvalue <- function(key, mat_mean, mat_sd, pos2aa) {
  if(is.matrix(key)) {
    key <- cbind(key, matrix(NA, nrow=nrow(key), ncol=4))
    key[,5] <- apply(key, 1, function(x) { mat_mean[x[1],x[2]]} )
    key[,6] <- apply(key, 1, function(x) { mat_sd[x[1],x[2]]} )
    key[,7:8] <- names(pos2aa[key[,1:2]])
  }
  return(key)
}

# cij/cmap to calculate mean/sd matrices
cij_gt_pearson <- list(gtp=net.md$gtp$cijs0, gdp=net.md$gdp$cijs0)
cmap_gt <- list(gtp=cm$GTP, gdp=cm$GDP)
cij_ras_pearson <- list(gtp=array_cij_ca_pearson_gtp[1:166,1:166,],
  gdp=array_cij_ca_pearson_gdp[1:166,1:166,])
cmap_ras <- list(gtp=cmap_ca_gtp_dist10[1:166,1:166],gdp=cmap_ca_gdp_dist10[1:166,1:166])

# The mean/sd matrix list of gtp/gdp
mean_cij_ras_pearson <- list(gtp=apply(array_cij_ca_pearson_gtp[1:166,1:166,],c(1:2),mean),
  gdp=apply(array_cij_ca_pearson_gdp[1:166,1:166,],c(1:2),mean))
sd_cij_ras_pearson <- list(gtp=apply(array_cij_ca_pearson_gtp[1:166,1:166,],c(1:2),sd),
  gdp=apply(array_cij_ca_pearson_gdp[1:166,1:166,],c(1:2),sd))
mean_cij_gt_pearson <- list(gtp=apply(net.md$gtp$cijs0,c(1:2),mean),
  gdp=apply(net.md$gdp$cijs0,c(1:2),mean))
sd_cij_gt_pearson <- list(gtp=apply(net.md$gtp$cijs0,c(1:2),sd),
  gdp=apply(net.md$gdp$cijs0,c(1:2),sd))

# generate the position vecktor with name of aa
pos2aa_ras <- 1:166
names(pos2aa_ras) <- ali["1QRA_A",ali["membership_ras",]!="0"]
pos2aa_ras <- list(pos2aa_ras,pos2aa_ras)
# gt I need the map of 305 position (mean/sd matrices are 305*305), then renumber them
pos2aa_gt <- 1:305
names(pos2aa_gt) <- ali["1TND_A",ali["membership_gt",]!="0"]
pos2aa_gt <- list(pos2aa_gt,pos2aa_gt)

# postion mapping between ras and gt
pos_ras2gt <- matrix(NA, nrow=2, ncol=146)
pos_ras2gt[1,] <- paste0(ali["resno_ras",ali["membership_146",]!="0"],
  ali["1QRA_A",ali["membership_146",]!="0"])
pos_ras2gt[2,] <- paste0(ali["resno_gt",ali["membership_146",]!="0"],
  ali["1TND_A",ali["membership_146",]!="0"])
rownames(pos_ras2gt) <- c("ras","gt")

## switch2 - alpha3 ##
edge <- c(4,6)
ras_04_06 <- key.extract(nets_md_ras_pearson_remodel, edge)
ras_04_06<- key.reorder(ras_04_06, dim=dim_ras)
ras_04_06 <- lapply(ras_04_06, function(x) {
    Map(key.addvalue, x, mean_cij_ras_pearson, sd_cij_ras_pearson, pos2aa_ras)
    })
gt_04_06 <- key.extract(nets_md_gt_pearson_remodel, edge)
gt_04_06<- key.reorder(gt_04_06, dim=dim_gt)
gt_04_06 <- lapply(gt_04_06, function(x) {
    Map(key.addvalue, x, mean_cij_gt_pearson, sd_cij_gt_pearson, pos2aa_gt)
    })
gt_04_06 <- key.renumber(gt_04_06)

## lobe2 - alpha4 ##
edge <- c(5,7)
ras_05_07 <- key.extract(nets_md_ras_pearson_remodel, edge)
ras_05_07<- key.reorder(ras_05_07, dim=dim_ras)
ras_05_07 <- lapply(ras_05_07, function(x) {
    Map(key.addvalue, x, mean_cij_ras_pearson, sd_cij_ras_pearson, pos2aa_ras)
    })
gt_05_07 <- key.extract(nets_md_gt_pearson_remodel, edge)
gt_05_07<- key.reorder(gt_05_07, dim=dim_gt)
gt_05_07 <- lapply(gt_05_07, function(x) {
    Map(key.addvalue, x, mean_cij_gt_pearson, sd_cij_gt_pearson, pos2aa_gt)
    })
gt_05_07 <- key.renumber(gt_05_07)

## alpha3 - loop8 ##
edge <- c(6,9)
ras_06_09 <- key.extract(nets_md_ras_pearson_remodel, edge)
ras_06_09<- key.reorder(ras_06_09, dim=dim_ras)
ras_06_09 <- lapply(ras_06_09, function(x) {
    Map(key.addvalue, x, mean_cij_ras_pearson, sd_cij_ras_pearson, pos2aa_ras)
    })
gt_06_09 <- key.extract(nets_md_gt_pearson_remodel, edge)
gt_06_09<- key.reorder(gt_06_09, dim=dim_gt)
gt_06_09 <- lapply(gt_06_09, function(x) {
    Map(key.addvalue, x, mean_cij_gt_pearson, sd_cij_gt_pearson, pos2aa_gt)
    })
gt_06_09 <- key.renumber(gt_06_09)

# beta1-3/alpha1 - alpha5
edge <- c(1,8)
ras_01_08 <- key.extract(nets_md_ras_pearson_remodel, edge)
ras_01_08<- key.reorder(ras_01_08, dim=dim_ras)
ras_01_08 <- lapply(ras_01_08, function(x) {
    Map(key.addvalue, x, mean_cij_ras_pearson, sd_cij_ras_pearson, pos2aa_ras)
    })
gt_01_08 <- key.extract(nets_md_gt_pearson_remodel, edge)
gt_01_08<- key.reorder(gt_01_08, dim=dim_gt)
gt_01_08 <- lapply(gt_01_08, function(x) {
    Map(key.addvalue, x, mean_cij_gt_pearson, sd_cij_gt_pearson, pos2aa_gt)
    })
gt_01_08 <- key.renumber(gt_01_08)


save(ras_04_06, gt_04_06,
     ras_05_07, gt_05_07,
     ras_06_09, gt_06_09,
     ras_01_08, gt_01_08,
     file="key_edges_md.RData")

# example: how to generate .csv then convert it to .xlsx
mat1 <- ras_04_06["0.3"][[1]]$gtp
mat1[,1] <- paste0(mat1[,1],mat1[,7])
mat1[,2] <- paste0(mat1[,2],mat1[,8])
mat1[,7:8] <- pos_ras2gt[2,match(mat1[,1:2], pos_ras2gt[1,])]
write.csv(mat1, file="ras_04_06.csv")

mat2 <- ras_04_06["0.3"][[1]]$gdp
mat2[,1] <- paste0(mat2[,1],mat2[,7])
mat2[,2] <- paste0(mat2[,2],mat2[,8])
mat2[,7:8] <- pos_ras2gt[2,match(mat2[,1:2], pos_ras2gt[1,])]
write.csv(mat2, file="tmp.csv")

ras_04_06[as.character(seq(0.3,0.6,by=0.1))]

mat1 <- gt_04_06["0.4"][[1]]$gtp
mat1[,1] <- paste0(mat1[,1],mat1[,7])
mat1[,2] <- paste0(mat1[,2],mat1[,8])
mat1[,7:8] <- pos_ras2gt[1,match(mat1[,1:2], pos_ras2gt[2,])]
write.csv(mat1, file="gt_04_06.csv")

mat2 <- gt_04_06["0.4"][[1]]$gdp
mat2[,1] <- paste0(mat2[,1],mat2[,7])
mat2[,2] <- paste0(mat2[,2],mat2[,8])
mat2[,7:8] <- pos_ras2gt[1,match(mat2[,1:2], pos_ras2gt[2,])]
write.csv(mat2, file="tmp.csv")

gt_04_06[as.character(seq(0.4,0.7,by=0.1))]




