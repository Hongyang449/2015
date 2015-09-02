## name: network_of_nma_ras.r 
## date: 09/02/2015

load("/Users/hyangl/Desktop/2015/0105_manuscript/pca/pdbs_ras.RData")
load("/Users/hyangl/Desktop/0105_manuscript/pca/ras_annotate.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/Desktop/0730_ras_transducin/comparison/layout_2d.RData")
load("/Users/hyangl/Desktop/2015/0105_manuscript/dccm/cmap.RData")

source("/Users/hyangl/Desktop/2015/functions/remodel.cna2.R")


# NMA
modes <- nma.pdbs(pdbs, rm.gaps=TRUE, ncore=4)

# load ras_annotate - here I cluster gtp/gdp based on auto pdb annotate
# need to further remove some (e.g. state1/GEF/inhibitor) structures
ras_annotate <- as.matrix(ras_annotate)
ids <- substr(basename(pdbs$id),1,6)
colors_gxp <- ras_annotate[, "ligandcolors"][ids]

# gtp
files <- pdbs$id[colors_gxp=="red"]
pdbs_gtp <- pdbaln(files=files, outfile=NULL, ncore=10)
modes_gtp <- nma.pdbs(pdbs_gtp, rm.gaps=TRUE, ncore=4)

# gdp
files <- pdbs$id[colors_gxp=="green"]
pdbs_gdp <- pdbaln(files=files, outfile=NULL, ncore=10)
modes_gdp <- nma.pdbs(pdbs_gdp, rm.gaps=TRUE, ncore=4)

# calculate pearson correlation
cij_nma_pearson_gtp <- dccm(modes_gtp)
cij_nma_pearson_gdp <- dccm(modes_gdp)

# membership
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

# remodel.cna() only use the cij to remodel; I don't need to calculate cna but pass in the cij
# What I did is provide dummy nets_pearson but change the cij (minus.log!)
# This is much faster when j is very small

# prepare the dummy nets!
cij_gtp <- filter.dccm(cij_nma_pearson_gtp$all.dccm, 
  cmap=cmap_ca_gtp_dist10, cutoff.cij=0.6)
cij_gdp <- filter.dccm(cij_nma_pearson_gdp$all.dccm, 
  cmap=cmap_ca_gdp_dist10, cutoff.cij=0.6)
cna_gtp <- cna(cij_gtp[1:166,1:166], cutoff.cij=0)
cna_gdp <- cna(cij_gdp[1:166,1:166], cutoff.cij=0)
nets_dummy <- list(gtp=cna_gtp,gdp=cna_gdp)

save(pdbs, pdbs_gtp, pdbs_gdp,
     modes, modes_gtp, modes_gdp,
     cij_nma_pearson_gtp, cij_nma_pearson_gdp,
     nets_nma_ras_pearson_remodel,
     file="network_of_nma_ras.RData")


