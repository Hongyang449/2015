## date: 09/21/2015
# add the "sse" colored by community colors
add.community.rect <- function(ylim,membership) {
  for (i in na.omit(unique(membership))) {
    pos <- bounds(which(membership == i))

    bo <- max(ylim) + (diff(ylim) * 0.001)
    to <- max(ylim) + (diff(ylim) * 0.04)
    rect(xleft=pos[,"start"], xright=pos[,"end"], ybottom=bo, ytop=to, col=vmd.colors()[i])
    bo <- min(ylim) - (diff(ylim) * 0.001)
    to <- min(ylim) - (diff(ylim) * 0.04)
    rect(xleft=pos[,"start"], xright=pos[,"end"], ybottom=bo, ytop=to, col=vmd.colors()[i])
  }
}


