## date: 09/10/2015
## Given a certain community edge, this function extract all the key residue edges from a list of networks
## e.g. nets/key: list 0.0-0.7 -> list gtp/gdp -> key edges matrix n*4

key.extract <- function(nets, edge) {
  key <- lapply(nets, function(x) {
    lapply(x, function(y) {
      inds <- intersect(which(y$community.key.cij[,3]==edge[1]),
        which(y$community.key.cij[,4]==edge[2]))
      if (length(inds)==1) return(t(as.matrix(y$community.key.cij[inds,])))
      else return(y$community.key.cij[inds,])
      })
    })
  names(key) <- names(nets)
  return(key)
}

# key.extract() prototype
#key_pearson_gtp_04_06 <- NULL
#key_pearson_gdp_04_06 <- NULL
#for (i in 1:length(nets_md_ras_pearson_remodel)) {
#
#  # extract all key residue pairs
#  mat_key_gtp <- nets_md_ras_pearson_remodel[[i]]$gtp$community.key.cij
#  mat_key_gdp <- nets_md_ras_pearson_remodel[[i]]$gdp$community.key.cij
#
#  inds <- intersect(which(mat_key_gtp[,3]==edge[1]), 
#    which(mat_key_gtp[,4]==edge[2]))
#  if (length(inds)==1) key_gtp=t(as.matrix(mat_key_gtp[inds,]))
#  else key_gtp=mat_key_gtp[inds, ]  
#  inds <- intersect(which(mat_key_gdp[,3]==edge[1]), 
#    which(mat_key_gdp[,4]==edge[2]))
#  if (length(inds)==1) key_gdp=t(as.matrix(mat_key_gdp[inds,]))
#  else key_gdp=mat_key_gdp[inds, ]  
#
#  key_pearson_gtp_04_06 <- c(key_pearson_gtp_04_06, list(key_gtp))
#  key_pearson_gdp_04_06 <- c(key_pearson_gdp_04_06, list(key_gdp))
#}

