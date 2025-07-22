#' Data structure with the  projected and boundary by node and class.
#' 
#' @param ppf is a PPforest object
#' @param id is a vector with the selected projection directions
#' @param sp is the simplex dimensions, if k is the number of classes sp = k - 1
#' @param dx first direction included in id
#' @param dy second direction included in id
#' @return Data frame needed to visualize a ternary plot
#' @export
#' @importFrom magrittr %>%
#' @references da da Silva, N., Cook, D. & Lee, EK. Interactive graphics for visually diagnosing forest classifiers in R. Comput Stat 40, 3105–3125 (2025). https://doi.org/10.1007/s00180-023-01323-x
#' @examples
#' #crab data set with all the observations used as training
#' pprf.crab <- PPforest(data = crab, std ='min-max', y = "Type",
#'  size.tr = 1, m = 100, size.p = .5, PPmethod = 'LDA')
#'  require(dplyr)
#' pl_ter <- function(dat, dx, dy ){
#'   p1  <- dat[[1]] %>% dplyr::filter(pair %in% paste(dx, dy, sep = "-") ) %>%
#'     dplyr::select(Class, x, y) %>%
#'     ggplot2::ggplot(ggplot2::aes(x, y, color = Class)) +
#'     ggplot2::geom_segment(data = dat[[2]], ggplot2::aes(x = x1, xend = x2,
#'                                                y = y1, yend = y2), color = "black" ) +
#'     ggplot2::geom_point(size = I(3), alpha = .5) +
#'     ggplot2::labs(y = " ",  x = " ") +
#'     ggplot2::theme(legend.position = "none", aspect.ratio = 1) +
#'     ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2") +
#'     ggplot2::labs(x = paste0("T", dx, " "), y = paste0("T", dy, " ")) +
#'     ggplot2::theme(aspect.ratio = 1)
#'   p1
#' }
#' #ternary plot in tree different selected dierections
#'  pl_ter(ternary_str(pprf.crab, id = c(1, 2, 3), sp = 3, dx = 1, dy = 2), 1, 2 )
#' 
ternary_str <-  function(ppf, id, sp, dx, dy){
  x <- NULL
  y <- NULL
  
f.helmert <- function(d)
{
  helmert <- rep(1 / sqrt(d), d)
  for (i in 1:(d - 1))
  {
    x <- rep(1 / sqrt(i * (i + 1)), i)
    x <- c(x, -i / sqrt(i * (i + 1)))
    x <- c(x, rep(0, d - i - 1))
    helmert <- rbind(helmert, x)
  }
  
  return(helmert)
}

makePairs <- function(dat, id) {
  aux <- dat[ ,-c(1, 2) ]
  
  d <- aux[, id]
  grid <- expand.grid(x = id, y = id)
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(
      Class = dat[, 1],
      ids = dat[, 2],
      x = dat[, xcol+2],
      y = dat[, ycol+2],
      pair = paste(grid[i, ], collapse = '-')
    )
  }))
  
  all
}

#ppf PPforest object
#id select proj directions
ternarydata <- function(ppf, id){
  n.class <- ppf$train %>% dplyr::select(tidyselect::all_of(ppf$class.var) )%>% unique() %>% nrow()
  projct <- t(f.helmert(nrow(unique(data.frame(ppf$train[, ppf$class.var]))))[-1,])
  
  dat3 <-
    data.frame(
      Class = ppf$train[, ppf$class.var],
      ids = 1:nrow(ppf$train),
      proj.vote = as.matrix(ppf$votes) %*% projct
    )
  
  ##with 3 or less classes
  empt <- rep(1:nrow(dat3), 3)
 
  if (n.class > 3) {
    gg1 <-  makePairs(dat3,id)
  }
  
  gg1 <-  makePairs(dat3, id )
  
  return(gg1)
}


f_composition <- function(data) {
  d <- dim(data)[2]
  hm <- f.helmert(d)
  x <- data - matrix(1 / d, dim(data)[1], d)
  return((x %*% t(hm))[,-1])
}

simplex <- function(sp) {
  vert <- f_composition(diag(sp + 1))
  colnames(vert) <- paste0("d", 1:ncol(vert))
  
  wires <-
    do.call(expand.grid, list(c(1:nrow(vert)), c(1:nrow(vert))))
  
  structure(list(points = vert,
                 edges = wires[!(wires[, 1] == wires[, 2]),]))
}

##ternary plot str

  s <- simplex(sp)
  pts <- data.frame(s$points)

  gg1 <- ternarydata(ppf,id)
  
  edg <- data.frame(x1=pts[,dx][s$edges[,1]], x2=pts[,dx][s$edg[,2]],
 
                                       y1=pts[,dy][s$edg[,1]], y2=pts[,dy][s$edg[,2]])


return( list(gg1, edg) )  
}


