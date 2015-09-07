# date: 09/06/2015

# position_vec2mat() map the position vector of a matrix to a X*2 matrix, 
# in which the 1st column is the row and the 2nd column is the column position 
# It can map the position vector to whole/upper/tri matrix coordinates

position_vec2mat <- function(x, dim, method=c('whole', 'upper', 'lower')) {

  method <- match.arg(method)

  position_vec2mat_whole <- function(x, dim) {
    mat <- matrix(0, nrow=length(x),ncol=2)
    for(i in 1:length(x)) {
      if((x[i] %% dim) == 0) {
        mat[i,1] <- x[i] - (x[i]/dim - 1) * dim
        mat[i,2] <- x[i]/dim
      }
      else {
        mat[i,1] <- x[i] - as.integer(x[i]/dim) * dim
        mat[i,2] <- as.integer(x[i]/dim) + 1
      }    
    }
    return(mat)
  }

  mat <- switch(method,
    whole = {
      position_vec2mat_whole(x, dim)
    },
    upper = {
      mat_tri <- matrix(0,nrow=dim,ncol=dim)
      mat_tri[upper.tri(mat_tri)] <- 1:(dim*(dim-1)/2)
      y <- x
      for(i in 1:length(x)) { y[i] <- which(mat_tri==x[i]) }
      mat <- position_vec2mat_whole(y,dim)
    },
    lower = {
      mat_tri <- matrix(0,nrow=dim,ncol=dim)
      mat_tri[lower.tri(mat_tri)] <- 1:(dim*(dim-1)/2)
      y <- x
      for(i in 1:length(x)) { y[i] <- which(mat_tri==x[i]) }
      mat <- position_vec2mat_whole(y,dim)
    })

  return(mat)
}


