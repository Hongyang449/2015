
## name: vmd_of_networks.r
## date: 09/11/2015

## manually copy file: /Users/hyangl/Desktop/0105_manuscript/dccm/pdb/gtp_plot.pdb to view/
## vim .pdb rm ligands and rename it to ras_166.pdb

library(bio3d)
load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_signif_ras.RData")

pdb <- read.pdb("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/view/ras_166.pdb")
x <- nets_md_ras_pearson_remodel[[1]]$gtp
vmd.cna(x, pdb, vmdfile="view/cna_ras_166.vmd", pdbfile="view/cna_ras_166.pdb",launch=T)



