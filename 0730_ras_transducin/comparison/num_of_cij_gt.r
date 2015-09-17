## name: num_of_cij_gt.r
## date: 09/09/2015

load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/net.md.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/transducin/md.cmap.RData")

cij_pearson <- list(gtp=net.md$gtp$cijs0, gdp=net.md$gdp$cijs0)
cmap <- list(gtp=cm$GTP,gdp=cm$GDP)

cutoff <- seq(0.05,0.95,by=0.05)

num_cij <- matrix(0, nrow=8, ncol=19)
rownames(num_cij) <- c("pearson_gtp_cmap", "pearson_gtp_nocmap",
                       "pearson_gdp_cmap", "pearson_gdp_nocmap",
                       "lmi_gtp_cmap", "lmi_gtp_nocmap",
                       "lmi_gdp_cmap", "lmi_gdp_nocmap")
colnames(num_cij) <- cutoff

# different cutoff with/without cmap filter
for (j in 1:length(cutoff)) {
  cij <- filter.dccm(cij_pearson$gtp, cmap=cmap$gtp, cutoff.cij=cutoff[j])
  num_cij["pearson_gtp_cmap",j] <- length(which(cij[upper.tri(cij)]!=0))
  cij <- filter.dccm(cij_pearson$gtp, cutoff.cij=cutoff[j])
  num_cij["pearson_gtp_nocmap",j] <- length(which(cij[upper.tri(cij)]!=0))

  cij <- filter.dccm(cij_pearson$gdp, cmap=cmap$gdp, cutoff.cij=cutoff[j])
  num_cij["pearson_gdp_cmap",j] <- length(which(cij[upper.tri(cij)]!=0))
  cij <- filter.dccm(cij_pearson$gdp, cutoff.cij=cutoff[j])
  num_cij["pearson_gdp_nocmap",j] <- length(which(cij[upper.tri(cij)]!=0))

#  cij <- filter.dccm(cij_lmi$gtp, cmap=cmap$gtp, cutoff.cij=cutoff[j])
#  num_cij["lmi_gtp_cmap",j] <- length(which(cij[upper.tri(cij)]!=0))
#  cij <- filter.dccm(cij_lmi$gtp, cutoff.cij=cutoff[j])
#  num_cij["lmi_gtp_nocmap",j] <- length(which(cij[upper.tri(cij)]!=0))

#  cij <- filter.dccm(cij_lmi$gdp, cmap=cmap$gdp, cutoff.cij=cutoff[j])
#  num_cij["lmi_gdp_cmap",j] <- length(which(cij[upper.tri(cij)]!=0))
#  cij <- filter.dccm(cij_lmi$gdp, cutoff.cij=cutoff[j])
#  num_cij["lmi_gdp_nocmap",j] <- length(which(cij[upper.tri(cij)]!=0))
}

# num of edge in cmap
num_gtp_cmap <- length(which(cmap$gtp[upper.tri(cmap$gtp)]!=0))
num_gdp_cmap <- length(which(cmap$gdp[upper.tri(cmap$gdp)]!=0))

# plot
plot(x=cutoff,y=num_cij["pearson_gtp_cmap",], xlab="cutoff.cij",
  ylab="number of edges", main="gt_edges_vs_cutoff_pearson", col="red",pch=1)
points(x=cutoff,y=num_cij["pearson_gtp_nocmap",],col="red",pch=4)
points(x=cutoff,y=num_cij["pearson_gdp_cmap",],col="green",pch=1)
points(x=cutoff,y=num_cij["pearson_gdp_nocmap",],col="green",pch=4)

# lines of cmap
segments(0,num_gtp_cmap,1,num_gtp_cmap,col="red",lwd=1,lty=2)
segments(0,num_gdp_cmap,1,num_gdp_cmap,col="green",lwd=1,lty=2)

legend("topright", c("gtp_cmap","gtp_nocmap","gdp_cmap","gdp_nocmap"),
  col=c("red","red","green","green"), pch=c(1,4,1,4),cex=0.7)

dev.copy2pdf(file="figures/gt_edges_vs_cutoff_pearson.pdf")

plot(x=cutoff[6:19],y=num_cij["pearson_gtp_cmap",6:19], xlab="cutoff.cij",
  ylab="number of edges", main="gt_edges_vs_cutoff_pearson", col="red",pch=1,
  ylim=c(0,4000))
points(x=cutoff[6:19],y=num_cij["pearson_gtp_nocmap",6:19],col="red",pch=4)
points(x=cutoff[6:19],y=num_cij["pearson_gdp_cmap",6:19],col="green",pch=1)
points(x=cutoff[6:19],y=num_cij["pearson_gdp_nocmap",6:19],col="green",pch=4)

# lines of cmap
segments(0,num_gtp_cmap,1,num_gtp_cmap,col="red",lwd=1,lty=2)
segments(0,num_gdp_cmap,1,num_gdp_cmap,col="green",lwd=1,lty=2)

legend("topright", c("gtp_cmap","gtp_nocmap","gdp_cmap","gdp_nocmap"),
  col=c("red","red","green","green"), pch=c(1,4,1,4),cex=0.7)

dev.copy2pdf(file="figures/gt_edges_vs_cutoff_pearson_0.30_0.95.pdf")

#plot(x=cutoff[5:19],y=num_cij["lmi_gtp_cmap",5:19], xlab="cutoff.cij",
#  ylab="number of edges", main="gt_edges_vs_cutoff_lmi", col="red",pch=1,
#  ylim=c(0,3000))
#points(x=cutoff[5:19],y=num_cij["lmi_gtp_nocmap",5:19],col="red",pch=4)
#points(x=cutoff[5:19],y=num_cij["lmi_gdp_cmap",5:19],col="green",pch=1)
#points(x=cutoff[5:19],y=num_cij["lmi_gdp_nocmap",5:19],col="green",pch=4)
#
## lines of cmap
#segments(0,num_gtp_cmap,1,num_gtp_cmap,col="red",lwd=1,lty=2)
#segments(0,num_gdp_cmap,1,num_gdp_cmap,col="green",lwd=1,lty=2)
#
#legend("topright", c("gtp_cmap","gtp_nocmap","gdp_cmap","gdp_nocmap"),
#  col=c("red","red","green","green"), pch=c(1,4,1,4),cex=0.7)
#
#dev.copy2pdf(file="figures/gt_edges_vs_cutoff_lmi_0.25_0.95.pdf")

save(num_cij, num_gtp_cmap, num_gdp_cmap,
     file="num_of_cij_gt.RData")

