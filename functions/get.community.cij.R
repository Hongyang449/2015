# get.community.cij() calculate community.cij of multiple single cij matrix
# input cij, cmap are lists!
# date: 09/04/2015

get.community.cij <- function(cij, cmap, cutoff.cij, membership) {

  # find the dim of community.cij
  l <- length(unique(membership))

  # generate the dummy networks
  # only pass cij to nets_dummy then use remodel.cna to calculate community.cij
  nets_dummy <- NULL
  for (j in 1:length(cij)) {
    nets_dummy <- c(nets_dummy, list(cna(cij[[j]][,,1], cutoff.cij=0.5)))
  }

  # calculate community cij
  array_community_cij <- NULL
  for (j in 1:length(cij)) {

    # use the consensus cij (set values to 1) as the cmap/extra.filter
    community_cij <- array(0, dim=c(l,l,dim(cij[[j]])[3]))
    cij_filter <- filter.dccm(cij[[j]], cmap=cmap[[j]], cutoff.cij=cutoff.cij)
    cij_filter[cij_filter!=0] <- 1

    for (i in 1:dim(cij[[j]])[3]) {

      # filter cij; minus.log
      cij1 <- filter.dccm(cij[[j]][,,i], cutoff.cij=0, extra.filter=cij_filter)
      cij1 <- abs(cij1)
      cij1[cij1>=1] <- 0.9999
      cij1[cij1>0] <- -log(cij1[cij1>0])

      # pass cij into nets_dummy
      nets <- nets_dummy
      nets[[j]]$cij <- cij1

      # use remoel.cna() to calculate community_cij
      nets_remodel <- remodel.cna(nets, member=membership, method="sum",
        col.edge="feature", scut=4, normalize=FALSE)

      # collect community.cij
      community_cij[,,i] <- nets_remodel[[j]]$community.cij
    }

    array_community_cij <- c(array_community_cij, list(community_cij))
  }

  return (array_community_cij)
}


