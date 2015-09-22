## date: 09/20/2015
## This function calculate the number of unique edges of a residue in a certain state
## It also return a -1/0/1 matrix in which 1 means unique edge in the first state

## If membership is provided, intra-community couplings are removed/only saved.
get.unique.key <- function(cij, cmap, scut=4, cutoff.cij=0, extra.filter=NULL, membership=NULL, 
  intra=TRUE, inter=TRUE, aaseq=NULL) {
  cij1 <- filter.dccm(cij[[1]], cmap=cmap[[1]], cutoff.cij=cutoff.cij, extra.filter=extra.filter)
  cij2 <- filter.dccm(cij[[2]], cmap=cmap[[2]], cutoff.cij=cutoff.cij, extra.filter=extra.filter)

  cij1[cij1!=0] <- 1
  cij2[cij2!=0] <- 1

  cij1[diag.ind(cij1,n=scut) | t(diag.ind(cij1,n=scut))] <- 0
  cij2[diag.ind(cij2,n=scut) | t(diag.ind(cij2,n=scut))] <- 0
  
  # remove intra couplings
  if(!is.null(membership) & !intra) {
    for(i in unique(membership)) {
      inds <- which(membership==i)
      cij1[inds,inds] <- 0
      cij2[inds,inds] <- 0
    }
  }

  # remove inter couplings
  if(!is.null(membership) & !inter) {
    cij_filter <- matrix(0, nrow=nrow(cij1), ncol=ncol(cij1))
    for(i in unique(membership)) {
      inds <- which(membership==i)
      cij_filter[inds,inds] <- 1
    }
    cij1 <- cij1 * cij_filter
    cij2 <- cij2 * cij_filter
  }

  # cij diff matrix
  cij_diff <- cij1 - cij2
  num_unique <- apply(cij_diff, 1, function(x) {
    c(length(which(x==1)),length(which(x==-1)))
    })
  colnames(num_unique) <- aaseq
  rownames(num_unique) <- names(cij)
  return(list(cij=cij_diff, num=num_unique))
}


