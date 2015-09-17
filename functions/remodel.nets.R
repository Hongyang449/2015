## date: 09/06/2015

remodel.nets <- function(cij, cmap, cutoff.cij, nets_dummy, membership, signif, layout_2d, p.cutoff) {

  # calculate consensus cij
  cij1 <- filter.dccm(cij[[1]], cmap=cmap[[1]], cutoff.cij=cutoff.cij)
  cij2 <- filter.dccm(cij[[2]], cmap=cmap[[2]], cutoff.cij=cutoff.cij)

  # abs!
  cij1 <- abs(cij1)
  cij2 <- abs(cij2)

  # minus.log cij
  cij1[cij1>=1] <- 0.9999
  cij1[cij1>0] <- -log(cij1[cij1>0])
  cij2[cij2>=1] <- 0.9999
  cij2[cij2>0] <- -log(cij2[cij2>0])

  # change the cij of dummy networks
  nets <- nets_dummy
  nets[[1]]$cij <- cij1
  nets[[2]]$cij <- cij2

  # remodel networks
  nets_remodel <- remodel.cna(nets, member=membership, method="sum",
    col.edge="significance", scut=4, normalize=FALSE, signif=signif, p.cutoff=p.cutoff)

  return(nets_remodel)
}

