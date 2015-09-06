vec2color <- function(vec, pal=c("blue", "green", "red"), n) {
  ##-- Define a color scale from a numeric vector
  require(classInt)
  return( findColours(classIntervals(vec, n=n, style="equal"), pal) )
}

get.edgecolor <- function(cij, method, colmap, cutoff, n = 30, ...) {

   if(length(cij) == 1 || method == 'none') return(NULL)

   require(classInt)

   pcij <- do.call(rbind, lapply(cij, function(x) x[lower.tri(x)]))

   color.code <- switch(method,  
      variance = {
         vars = apply(pcij, 2, sd)
         rep(vars, each=nrow(pcij))
      }, 
      feature = {
         as.vector( apply(pcij, 2, function(x) {
            if(all(x==0) || ((max(x) - min(x))/max(x)) < cutoff)
               return(rep(0, length(x)))
            class <- suppressWarnings(classIntervals(x, n=2, style='equal'))
            flag <- as.integer(x >= class$brks[2])
            flag <- flag * (1:length(flag))
         }) )
      })
  
   colors <- rep(NA, length(color.code)) 
   color.code <- color.code[as.vector(pcij > 0)] # exclude pairs that have no edge 
   if(length(unique(color.code)) == 1)
      colors[as.vector(pcij > 0)] <- rep(colmap[1], length(color.code))
   else 
      colors[as.vector(pcij>0)] <- switch(method, 
         variance = suppressWarnings( vec2color(color.code, pal = colmap, n = n) ),
         feature = suppressWarnings( vec2color(color.code, pal = colmap[sort(unique(color.code)+1)], 
                       n = length(unique(color.code))) )
      )
   colors <- split(colors, f = rep(1:nrow(pcij), ncol(pcij)))
   colors <- lapply(colors, function(x) x[!is.na(x)])
   return(colors)
}

normalize.cij <- function(cij, ...) {
   pcij <- do.call(rbind, lapply(cij, function(x) x[lower.tri(x)]))
   pcij[pcij <= 0] <- NA
   for(i in 1:nrow(pcij)) {
      vmin <- min(pcij[i, ], na.rm = TRUE)
      vmax <- max(pcij[i, ], na.rm = TRUE)
      if(vmax > vmin)
         pcij[i, ] <- (pcij[i, ] - vmin) / (vmax - vmin) * 0.9 + 0.05
      else
         pcij[i, !is.na(pcij[i, ])] <- 0.95
   }
   pcij[is.na(pcij)] <- 0

   ncij <- apply(pcij, 1, function(x) {
      xx = cij[[1]]
      xx[lower.tri(xx)] <- x
      xx[upper.tri(xx)] <- t(xx)[upper.tri(xx)]
      return(list(xx))
   } )
   do.call("c", ncij)
}

# Methods to calculate community cijs:
#    - Mean, use the mean of top 'ne' inter-community cijs as the cij of the CG community network
#    - Max,  use the maximum inter-community cij as the cij of the community network
#    - Sum,  use the sum of all inter-community cijs as the cij of the community network
remodel.cna <- function(x, member = NULL, col = NULL, minus.log = TRUE,
       method = c('none', 'sum', 'mean', 'max'), ne=3, scut=4, normalize = TRUE, 
       vmd.color = TRUE, col.edge=c('none', 'variance', 'feature'),
       colmap.edge = NULL,  coledge.cutoff = 0.5, ...) {

   require(igraph)
   method <- match.arg(method)
   col.edge <- match.arg(col.edge)

   # assume input 'x' is a network ensemble 
   if(inherits(x, "cna")) x = list(x)
   
   # will return with network names
   net.names <- names(x)

   # if not provided, the membership is extracted from the network
   if(is.null(member)) {
      member = lapply(x, function(y) y$communities$membership)
   } else {
      member = rep(list(member), length(x))
      if(method == 'none') {
         warning('method is "none" but member is not NULL. Set method to "sum"')
         method = 'sum'
      }
   }
   
   if(method != 'none') {
      # calculate cijs for CG networks with defined membership
      cg.cij <- lapply(1:length(x), function(i) {
         cij <- x[[i]]$cij
         member <- member[[i]]
         n.mem <- length(unique(member))
         n2 <- pairwise(n.mem)
 
         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         cij[diag.ind(cij, n=scut)] <- 0
         cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]
   
         w <- apply(n2, 1, function(i) {
           ind1 <- which(member %in% i[1])
           ind2 <- which(member %in% i[2])
           switch(method, 
              max = max(cij[ind1, ind2]),
              mean= mean(sort(cij[ind1, ind2], decreasing=TRUE)[1:ne]),  
              sum = sum(cij[ind1, ind2]) 
           )
         } )
         cg.cij <- matrix(1, nrow=n.mem, ncol=n.mem)
         cg.cij[lower.tri(cg.cij)] <- w
         cg.cij[upper.tri(cg.cij)] <- t(cg.cij)[upper.tri(cg.cij)]
         cg.cij
      } )
   } else {
      cg.cij <- lapply(x, function(y) {
         cij <- y$community.cij
         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         cij
      })
   }

   if(method != 'none' || 
      any(sapply(x, function(y) is.null(y$community.key.cij))) ) {
      ###############################################################
      # Get the indices of cijs that contribute to the community cijs.
      # This will be stored in the returned networks and will be used
      # with another application related to plot 3D networks with VMD.
      # For 2D plot, this variable is not used.
      key <- lapply(1:length(x), function(i) {
         cij <- x[[i]]$cij
         member <- member[[i]]
         n.mem <- length(unique(member))
         n2 <- pairwise(n.mem)

         if(minus.log) cij[cij>0] = exp(-cij[cij>0])
         cij[diag.ind(cij, n=scut)] <- 0
         cij[lower.tri(cij)] <- t(cij)[lower.tri(cij)]
   
         val <- apply(n2, 1, function(i) {
           ind1 <- which(member %in% i[1])
           ind2 <- which(member %in% i[2])
           k <- switch(method, 
              max = which.max(cij[ind1, ind2]),
              mean= order(cij[ind1, ind2], decreasing=TRUE)[1:ne],
              sum = 1:length(cij[ind1, ind2])
           )
           k2 <- floor((k-1)/length(ind1)) + 1
           k1 <- (k-1) %% length(ind1) + 1
           list(i=ind1[k1], j=ind2[k2])
         })
         cbind(unlist(lapply(val, '[[', 'i')), unlist(lapply(val, '[[', 'j')))
      })
      #################################################################
   } else {
      key <- lapply(x, '[[', 'community.key.cij')
   }

   if(length(x) == 1) {
      edge.color = NULL
   } else {
      if(col.edge == 'none') {
         edge.color = NULL
      } else {
         check <- TRUE
         if(length(unique(sapply(member, length))) != 1) check = FALSE
         else 
            check <- all(apply(do.call(rbind, member), 2, function(x) length(unique(x))==1))
         if(!check) edge.color=NULL
         else {
            # set default colormap according to the color method for edges
            if(is.null(colmap.edge)) {
               colmap.edge <- switch(col.edge, 
                  variance = c('blue', 'red'),
                  feature = {
                     tcol = c('gray', 'red', 'darkgreen', 'blue') 
                     tcol = union(tcol, colors()[-1])
                     tcol[1:(length(x)+1)]
               })
            }
            if(col.edge=='feature' && length(colmap.edge) != (length(x)+1))
               stop('Number of colors does not match input number of networks')

            # Calculate the variance of CG cijs across networks.
            # The values will be used to color the edges of CG networks
            edge.color <- get.edgecolor(cij = cg.cij, method = col.edge,
                   colmap = colmap.edge, cutoff = coledge.cutoff, ...)
         }
      }
   }
   if(is.null(edge.color)) 
      edge.color <- lapply(x, function(y) 
          get.edge.attribute(y$community.network, name='color') )

   # Normalize CG cijs for each network
   if(normalize) cg.cij <- normalize.cij(cg.cij, ...)

   # calculate node sizes and colors
   cg.node.size <- lapply(member, table)
   n.mem <- sapply(member, function(y) length(unique(y)))
   if(is.null(col)) {
      check <- all(sapply(1:length(n.mem), 
                  function(i) n.mem[i] == length(V(x[[i]]$community.network))))
      if(check) 
         col <- lapply(x, function(y) get.vertex.attribute(y$community.network, 'color'))
      else 
         col = lapply(member, function(y) {
                    col <- 1:length(unique(y))
                    if(vmd.color) vmd.colors()[col]
                    else col
               } )
   } else {
      if(!all(sapply(n.mem, '==', length(col))))
         stop("Length of color vector doesn't match number of communities")
      if(is.numeric(col) && vmd.color) col = vmd.colors()[col]
      col = rep(list(col), length(member))
   }
   # update network components
   x <- lapply(1:length(x), function(i) {
      y = x[[i]]
      cij <- cg.cij[[i]]
      if(minus.log && (method != 'sum' || normalize)) {
         cij[cij>=1] <- 0.9999
         cij[cij>0] <- -log(cij[cij>0])
      }
      y$community.network <- graph.adjacency(cij, 
                                 mode = "undirected",
                                 weighted = TRUE,
                                 diag = FALSE)
      y$community.network <- set.vertex.attribute(y$community.network, "size", value=cg.node.size[[i]])

      y$community.network <- set.edge.attribute(y$community.network, "color", value=edge.color[[i]])
      y$communities$membership <- member[[i]]
      y$community.cij <- cij
      inds = which(abs(y$cij[key[[i]]]) > 0)
      y$community.key.cij <- key[[i]][inds, ]

      y$network <- set.vertex.attribute(y$network, "color", value= col[[i]][member[[i]]])
      y$community.network <- set.vertex.attribute(y$community.network, "color", value = col[[i]])
      if(!is.null(y$community.reindex)) {
         if(vmd.color) y$community.reindex = match(col[[i]], vmd.colors())
         else y$community.reindex = 1:length(unique(member[[i]]))
      }
      y
   } )

   names(x) <- net.names
   return(x)
}
