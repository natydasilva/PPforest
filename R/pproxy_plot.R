#' Proximity matrix visualization
#'
#' @param ppf a PPforest object
#' @param type is MDS or heat
#' @param k number of dimensions of the MDS layout 
#' @return proximity matrix plot, if type is MDS then a MDS plot using proximity matrix information is shown and if type is heat a heat map of the proximity matrix is shown
#' @export
#' @examples
#' #leukemia data set with all the observations used as training
#' pprf.leukemia <- PPforest(data = leukemia, class = 'Type',
#' size.tr = 1, m = 70, size.p = .4, PPmethod = 'PDA', strata = TRUE)
#' pproxy_plot(pprf.leukemia, type= "MDS", k = 3)
pproxy_plot <- function(ppf, type = "heat", k) {
  if (type == "heat") {
    value <- NULL
    Var1 <- NULL
    Var2 <- NULL
    
    id <- diag(dim(ppf$train)[1])
    id[lower.tri(id, diag = TRUE)] <- ppf[[9]]$proxi
    id[upper.tri(id)] <- t(id)[upper.tri(id)]
    m.prox <- reshape2::melt(id)
    m.prox$Var2 <- as.factor(m.prox$Var2)
    m.prox$Var2 <-
      factor(m.prox$Var2, levels = rev(levels(m.prox$Var2)))
    m.prox$Var1 <- as.factor(m.prox$Var1)
    
    a <-
      ggplot2::ggplot(m.prox, ggplot2::aes(Var1, Var2)) + ggplot2::xlab("") +
      ggplot2::ylab("") + ggplot2::geom_tile(ggplot2::aes(fill = value)) + ggplot2::scale_fill_gradient(high = "#132B43",
                                                                                                        low = "#56B1F7", name = "Proximity") + ggplot2::theme(aspect.ratio = 1)
    
    plotly::ggplotly(a)
    
  } else {
    value <- NULL
    Var1 <- NULL
    Var2 <- NULL
    MDS1 <- NULL
    MDS2 <- NULL
    Class <- NULL
    x <- NULL
    y <- NULL
    Type <- NULL
    
    
    id <- diag(dim(ppf$train)[1])
    id[lower.tri(id, diag = TRUE)] <- 1 - ppf[[9]]$proxi
    rf.mds <-
      stats::cmdscale(d = stats::as.dist(id), eig = TRUE, k = k)
    colnames(rf.mds$points) <- paste("MDS", 1:k, sep = "")
    nlevs <- nlevels(ppf$train[, 1])
    
    if (k == 2) {
      df <- data.frame(Class = ppf$train[, 1], rf.mds$points)
      a <-
        ggplot2::ggplot(data = df) + ggplot2::geom_point(
          ggplot2::aes(x = MDS1, y = MDS2, color = Class),size = I(3),alpha = .5
        ) + ggplot2::theme(aspect.ratio = 1) +
        ggplot2::scale_colour_brewer(type = "qual",palette = "Dark2",name = "Class")
      
      plotly::ggplotly(a)
      
    } else {
      df <- data.frame(Class = ppf$train[, 1], rf.mds$points)
      makePairs <- function(data)
      {
        grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
        grid <- subset(grid, x != y)
        all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
          xcol <- grid[i, "x"]
          ycol <- grid[i, "y"]
          data.frame(
            xvar = names(data)[ycol], yvar = names(data)[xcol],
            x = data[, xcol], y = data[, ycol], data
          )
        }))
        all$xvar <- factor(all$xvar, levels = names(data))
        all$yvar <- factor(all$yvar, levels = names(data))
        densities <-
          do.call("rbind", lapply(1:ncol(data), function(i) {
            data.frame(
              xvar = names(data)[i], yvar = names(data)[i], x = data[, i]
            )
          }))
        list(all = all, densities = densities)
      }
      
      # expand data frame for pairs plot
      gg1 = makePairs(df[,-1])
      
      # new data frame mega iris
      mega_data = data.frame(gg1$all, Class = rep(df$Class, length = nrow(gg1$all)))
      
      #pairs plot
      
      
      a <-
        ggplot2::ggplot(mega_data, ggplot2::aes_string(x = "x", y = "y")) +  
        ggplot2::theme(legend.position =  'bottom', aspect.ratio = 1) +
        ggplot2::facet_grid(xvar ~ yvar, scales = "free") +
        ggplot2::geom_point( ggplot2::aes(colour = Class), na.rm = TRUE, alpha = 0.5, size = I(3)) +
        ggplot2::stat_density( ggplot2::aes_string(x = "x", y = "..scaled.."),
          data = gg1$densities, position = "identity",
          colour = "grey20", geom = "line",  size =I(.5)) +
        ggplot2::scale_colour_brewer(type = "qual",palette = "Dark2", name =
                                       "Class")
      
      plotly::ggplotly(a)
      
    }
    
  }
}
