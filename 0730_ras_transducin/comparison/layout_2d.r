## name: layout_2d.r
## date: 09/01/2015

# This is the layout_2d from Xinqiu

layout_2d <- matrix(0,nrow=11,ncol=2)
layout_2d[1,] <- c(15.964872, 76.19087)
layout_2d[2,] <- c(-1.961059, 85.52253)
layout_2d[3,] <- c(-1.706022, 76.64920)
layout_2d[4,] <- c(10.038941, 82.52253)
layout_2d[5,] <- c(20.331522, 87.59111)
layout_2d[6,] <- c(11.590258, 96.27745)
layout_2d[7,] <- c(22.588111, 96.51467)
layout_2d[8,] <- c(27.513000, 76.73330)
layout_2d[9,] <- c(29.588111, 92.51467)
layout_2d[10,] <- c(-21.961059, 81.52253)
layout_2d[11,] <- c(-21.961059, 91.52253)

# modify it for ras

layout_2d[3,] <- layout_2d[3,] + c(3,0)
layout_2d[4,] <- layout_2d[4,] + c(0,2)
layout_2d[7,] <- layout_2d[7,] + c(2,-1)
layout_2d[9,] <- layout_2d[9,] + c(0,-3)

save(layout_2d, file="layout_2d.RData")

