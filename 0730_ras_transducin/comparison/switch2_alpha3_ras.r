## name: switch2_alpha3_pairs.r
## date: 09/08/0215

library(bio3d)

load("/Users/hyangl/Desktop/2015/0730_ras_transducin/comparison/key_md_ras.RData")

edge_in <- c(4,6)

########### pearson ############

# list 0.0-0.7 -> mat n*4 unique pairs for edge_in
key_unique_pearson_gtp <- NULL
key_unique_pearson_gdp <- NULL
for (i in 1:length(key_pearson_gtp)) {
  for (j in 1:length(key_pearson_gtp[[i]])) {

    # in case of dim=0*4 error
    if (dim(key_pearson_gtp[[i]][[j]])[1]!=0) edge <- key_pearson_gtp[[i]][[j]][1,3:4]
    else edge <- key_pearson_gdp[[i]][[j]][1,3:4]

    # find residue pairs of this edge
    if (all(edge==edge_in)) {
      key_gtp <- key_pearson_gtp[[i]][[j]]
      key_gdp <- key_pearson_gdp[[i]][[j]]
    }
  }

  # find unique residue pairs that contribute to this edge
  inds_gtp <- which(!(apply(key_gtp,1,list) %in% apply(key_gdp,1,list)))
  inds_gdp <- which(!(apply(key_gdp,1,list) %in% apply(key_gtp,1,list)))

  # length==1, convert the vector to matrix
  if (length(inds_gtp)==1) x=t(as.matrix(key_gtp[inds_gtp,]))
  else x=key_gtp[inds_gtp, ]
  if (length(inds_gdp)==1) y=t(as.matrix(key_gdp[inds_gdp,]))
  else y=key_gdp[inds_gdp, ]

  key_unique_pearson_gtp <- c(key_unique_pearson_gtp, list(x))
  key_unique_pearson_gdp <- c(key_unique_pearson_gdp, list(y))
}

########### lmi ############

# list 0.0-0.7 -> mat n*4 unique pairs for edge_in
key_unique_lmi_gtp <- NULL
key_unique_lmi_gdp <- NULL
for (i in 1:length(key_lmi_gtp)) {
  for (j in 1:length(key_lmi_gtp[[i]])) {

    # in case of dim=0*4 error
    if (dim(key_lmi_gtp[[i]][[j]])[1]!=0) edge <- key_lmi_gtp[[i]][[j]][1,3:4]
    else edge <- key_lmi_gdp[[i]][[j]][1,3:4]

    # find residue pairs of this edge
    if (all(edge==edge_in)) {
      key_gtp <- key_lmi_gtp[[i]][[j]]
      key_gdp <- key_lmi_gdp[[i]][[j]]
    }
  }

  # find unique residue pairs that contribute to this edge
  inds_gtp <- which(!(apply(key_gtp,1,list) %in% apply(key_gdp,1,list)))
  inds_gdp <- which(!(apply(key_gdp,1,list) %in% apply(key_gtp,1,list)))

  # length==1, convert the vector to matrix
  if (length(inds_gtp)==1) x=t(as.matrix(key_gtp[inds_gtp,]))
  else x=key_gtp[inds_gtp, ]
  if (length(inds_gdp)==1) y=t(as.matrix(key_gdp[inds_gdp,]))
  else y=key_gdp[inds_gdp, ]

  key_unique_lmi_gtp <- c(key_unique_lmi_gtp, list(x))
  key_unique_lmi_gdp <- c(key_unique_lmi_gdp, list(y))
}

########### reorder residue pairs ##############
## Here I map the coordinates to matrix then order() as a inds.

######### pearson ########

key_unique_pearson_gtp_reorder <- key_unique_pearson_gtp
for (i in 1:length(key_unique_pearson_gtp)) {
  inds <- rep(NA, length=dim(key_unique_pearson_gtp[[i]])[1])
  if(length(inds)==0) next
  for (j in 1:dim(key_unique_pearson_gtp[[i]])[1]) {
    mat <- matrix(NA, nrow=166, ncol=166)
    mat[key_unique_pearson_gtp[[i]][j,1], key_unique_pearson_gtp[[i]][j,2]] <- 1
    inds[j] <- which(!is.na(t(mat)))
  }
  if(length(inds)==1) key_unique_pearson_gtp_reorder[[i]] <- key_unique_pearson_gtp[[i]]
  else key_unique_pearson_gtp_reorder[[i]] <- key_unique_pearson_gtp[[i]][order(inds),]
}
key_unique_pearson_gdp_reorder <- key_unique_pearson_gdp
for (i in 1:length(key_unique_pearson_gdp)) {
  inds <- rep(NA, length=dim(key_unique_pearson_gdp[[i]])[1])
  if(length(inds)==0) next
  for (j in 1:dim(key_unique_pearson_gdp[[i]])[1]) {
    mat <- matrix(NA, nrow=166, ncol=166)
    mat[key_unique_pearson_gdp[[i]][j,1], key_unique_pearson_gdp[[i]][j,2]] <- 1
    inds[j] <- which(!is.na(t(mat)))
  }
  if(length(inds)==1) key_unique_pearson_gdp_reorder[[i]] <- key_unique_pearson_gdp[[i]]
  else key_unique_pearson_gdp_reorder[[i]] <- key_unique_pearson_gdp[[i]][order(inds),]
}

######### lmi ########
key_unique_lmi_gtp_reorder <- key_unique_lmi_gtp
for (i in 1:length(key_unique_lmi_gtp)) {
  inds <- rep(NA, length=dim(key_unique_lmi_gtp[[i]])[1])
  if(length(inds)==0) next
  for (j in 1:dim(key_unique_lmi_gtp[[i]])[1]) {
    mat <- matrix(NA, nrow=166, ncol=166)
    mat[key_unique_lmi_gtp[[i]][j,1], key_unique_lmi_gtp[[i]][j,2]] <- 1
    inds[j] <- which(!is.na(t(mat)))
  }
  if(length(inds)==1) key_unique_lmi_gtp_reorder[[i]] <- key_unique_lmi_gtp[[i]]
  else key_unique_lmi_gtp_reorder[[i]] <- key_unique_lmi_gtp[[i]][order(inds),]
}
key_unique_lmi_gdp_reorder <- key_unique_lmi_gdp
for (i in 1:length(key_unique_lmi_gdp)) {
  inds <- rep(NA, length=dim(key_unique_lmi_gdp[[i]])[1])
  if(length(inds)==0) next
  for (j in 1:dim(key_unique_lmi_gdp[[i]])[1]) {
    mat <- matrix(NA, nrow=166, ncol=166)
    mat[key_unique_lmi_gdp[[i]][j,1], key_unique_lmi_gdp[[i]][j,2]] <- 1
    inds[j] <- which(!is.na(t(mat)))
  }
  if(length(inds)==1) key_unique_lmi_gdp_reorder[[i]] <- key_unique_lmi_gdp[[i]]
  else key_unique_lmi_gdp_reorder[[i]] <- key_unique_lmi_gdp[[i]][order(inds),]
}

save(key_unique_pearson_gtp, key_unique_pearson_gdp,
     key_unique_lmi_gtp, key_unique_lmi_gdp,
     key_unique_pearson_gtp_reorder, key_unique_pearson_gdp_reorder,
     key_unique_lmi_gtp_reorder, key_unique_lmi_gdp_reorder,
     file="switch2_alpha3_pairs.RData")


