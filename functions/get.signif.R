# get.signif() to calculate t.test values
# input two sets of cij arrays, this function returns the t.test p-value of each edge
# date: 09/04/2015

get.signif <- function(cij1, cij2, mu=0, upper=FALSE) {
  # mat_p restore the p.value of t.test of a certain community edge
  # diag value are NA (t.test error for identical values) 
  mat_p <- matrix(NA, nrow=dim(cij1)[1], ncol=dim(cij1)[2])

  # use '[' to extract the 1:2 dimensions of a array!
  # then we can iterate the matrix useing for loop!
  # e.g. if the array cij1 is 9*9*10; the corresponding mat_1 is 81*10
  mat_1 <- apply(cij1, 3, '[')
  mat_2 <- apply(cij2, 3, '[')

  # inds: whether upper or undiaganol
  if(upper) inds <- which(upper.tri(cij1[,,1]))
  else inds <- which(!diag.ind(cij1[,,1]))
  for (i in inds) {
    if (abs(mean(mat_1[i,]) - mean(mat_2[i,])) > mu) {
      if (mean(mat_1[i,]) >= mean(mat_2[i,]))
      mat_p[i] <- t.test(mat_1[i,], mat_2[i,], mu=mu)$p.value
      else mat_p[i] <- t.test(mat_2[i,], mat_1[i,], mu=mu)$p.value
    }
    else mat_p[i] <- NaN
  }
  return( mat_p )
}

