###### find signif pairs between communities

## name: key_md_ras.r
## date: 09/07/2015

# which(p<0.05) -> vec2mat -> list pairs
# list 0.0-0.7 -> list signif_edges -> mat key_residue_pairs

library(bio3d)

load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/network_of_md_signif_ras.RData")

source("/Users/hyangl/Desktop/2015/functions/position_vec2mat.R")

############## pearson ###############

# !! bug
# nets_md_ras_pearson_remodel[[8]]$gdp$community.key.cij
#      [,1] [,2] [,3]
# [1,]   18    1    5
# [2,]  146    1    5

# list of key residue pairs
key_pearson_gtp <- NULL
key_pearson_gdp <- NULL
for (i in 1:(length(p_community_cij_md_ras_pearson)-1)) {
  # p.value matrix of edges
  mat_p <- p_community_cij_md_ras_pearson[[i]]

  # use position_vec2mat() to return the pairs of communities that distinguish two states
  # use upper.tri then small community goes first
  mat_pairs <- position_vec2mat(which(mat_p[upper.tri(mat_p)]<0.05), 
    dim=dim(mat_p)[1], method="upper")

  # extract all key residue pairs
  mat_key_gtp <- nets_md_ras_pearson_remodel[[i]]$gtp$community.key.cij
  mat_key_gdp <- nets_md_ras_pearson_remodel[[i]]$gdp$community.key.cij

  # list key_gtp for a certain cutoff.cij
  key_gtp <- NULL
  key_gdp <- NULL 

  for(j in 1:dim(mat_pairs)[1]) {
    # inds for a certain pair of communities
    inds <- intersect(which(mat_key_gtp[,3]==mat_pairs[j,1]),
      which(mat_key_gtp[,4]==mat_pairs[j,2]))
    # list the key residue pairs contributing to the edge
    if (length(inds)==1) x=t(as.matrix(mat_key_gtp[inds,]))
    else x=mat_key_gtp[inds, ]
    key_gtp <- c(key_gtp, list(x))

    inds <- intersect(which(mat_key_gdp[,3]==mat_pairs[j,1]),
      which(mat_key_gdp[,4]==mat_pairs[j,2]))
    if (length(inds)==1) x=t(as.matrix(mat_key_gdp[inds,]))
    else x=mat_key_gdp[inds, ]
    key_gdp <- c(key_gdp, list(x))
  }

  key_pearson_gtp <- c(key_pearson_gtp, list(key_gtp))
  key_pearson_gdp <- c(key_pearson_gdp, list(key_gdp))
}


############## lmi ###############

# list of key residue pairs
key_lmi_gtp <- NULL
key_lmi_gdp <- NULL
for (i in 1:length(p_community_cij_md_ras_lmi)) {
  # p.value matrix of edges
  mat_p <- p_community_cij_md_ras_lmi[[i]]

  # use position_vec2mat() to return the pairs of communities that distinguish two states
  # use upper.tri then small community goes first
  mat_pairs <- position_vec2mat(which(mat_p[upper.tri(mat_p)]<0.05),
    dim=dim(mat_p)[1], method="upper")

  # extract all key residue pairs
  mat_key_gtp <- nets_md_ras_lmi_remodel[[i]]$gtp$community.key.cij
  mat_key_gdp <- nets_md_ras_lmi_remodel[[i]]$gdp$community.key.cij

  # list key_gtp for a certain cutoff.cij
  key_gtp <- NULL
  key_gdp <- NULL

  for(j in 1:dim(mat_pairs)[1]) {
    # inds for a certain pair of communities
    inds <- intersect(which(mat_key_gtp[,3]==mat_pairs[j,1]),
      which(mat_key_gtp[,4]==mat_pairs[j,2]))
    # list the key residue pairs contributing to the edge
    if (length(inds)==1) x=t(as.matrix(mat_key_gtp[inds,]))
    else x=mat_key_gtp[inds, ]
    key_gtp <- c(key_gtp, list(x))

    inds <- intersect(which(mat_key_gdp[,3]==mat_pairs[j,1]),
      which(mat_key_gdp[,4]==mat_pairs[j,2]))
    if (length(inds)==1) x=t(as.matrix(mat_key_gdp[inds,]))
    else x=mat_key_gdp[inds, ]
    key_gdp <- c(key_gdp, list(x))
  }

  key_lmi_gtp <- c(key_lmi_gtp, list(key_gtp))
  key_lmi_gdp <- c(key_lmi_gdp, list(key_gdp))
}


save(key_pearson_gtp, key_pearson_gdp,
     key_lmi_gtp, key_lmi_gdp,
     p_community_cij_md_ras_pearson, p_community_cij_md_ras_lmi,
     nets_md_ras_pearson_remodel, nets_md_ras_lmi_remodel,
     file="key_md_ras.RData")



