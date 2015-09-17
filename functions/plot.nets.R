## date: 09/10/2015
## This function plot a list of remodeled networks returned from remodel.cna()

plot.nets <- function(nets, layout_2d) {
  layout(matrix(1:(2*length(nets)),nrow=2))
  for(i in 1:length(nets)) {
    w1 <- E(nets[[i]][[1]]$community.network)$weight
    w2 <- E(nets[[i]][[2]]$community.network)$weight

    plot.cna(nets[[i]][[1]], layout=layout_2d, weights = w1, vertex.label=NA,
      main=paste0(names(nets[[i]])[1],"_cutoff.cij=",names(nets)[i]))
    plot.cna(nets[[i]][[2]], layout=layout_2d, weights = w2, vertex.label=NA,
      main=paste0(names(nets[[i]])[2],"_cutoff.cij=",names(nets)[i]))
  }
}

