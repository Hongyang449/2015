## date: 09/10/2015
## Given a list of key residue edges returned from key.extract(), this function reorder it
## e.g. key: list 0.0-0.7 -> list gtp/gdp -> key edges matrix n*4

key.reorder <- function(key, dim) {
  key_reorder <- lapply(key, function(x) {
    lapply(x, function(y) {
      inds <- rep(NA, length=nrow(y))
      if(length(inds)==0 || length(inds)==1) return(y)
      for (j in 1:nrow(y)) {
        mat <- matrix(NA, nrow=dim, ncol=dim)
        mat[y[j,1],y[j,2]] <- 1
        inds[j] <- which(!is.na(t(mat)))
      }
      return(y[order(inds),])
      })
    })
  names(key_reorder) <- names(key)
  return(key_reorder)
}
