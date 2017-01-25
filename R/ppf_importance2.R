#' Importance variable visualization for each tree in the forest
#' 
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param ppf is a PPforest object
#' @param tree.id id number of the tree we want to see the specific importance measure will be in red
#' @return A side-by-side jittered dotplot with the absolute importance measure for each tree in the forest
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' pprf.crab <- PPforest2(data = crab, class = "Type",
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA', strata = TRUE)
#' ppf_importance2(crab, "Type", pprf.crab, tree.id = 1)
ppf_importance2 <- function(data, class, ppf, tree.id = 1) {
  var <- NULL
  value <- NULL
  ids <- NULL
  node <- NULL
  Variables <- NULL
  Abs.importance <- NULL
  
  
  
  bnf <- function(x){
    bn <- abs(x$projbest.node)
    bn[bn == 0] <- NA
    data.frame(node = 1:dim(x$projbest.node)[1], bn)
  }
  bestnode <- ppf[[8]] %>%  lapply(bnf) %>%dplyr::bind_rows()
  
  
  colnames(bestnode)[-1] <-
    colnames(data[,-as.numeric(1:ncol(data) %*% as.numeric(colnames(data) ==
                                                             class))])
  bestnode$node <- as.factor(bestnode$node)
  
  
  
  
  aux <-
    bestnode %>% dplyr::mutate(ids = rep(1:ppf$n.tree,each = dim(ppf[[8]][[2]]$projbest.node)[1])) %>% tidyr::gather(var,value,-ids,-node)
  aux$Variables <- as.numeric(as.factor(aux$var))
  aux$Abs.importance <- round(aux$value,2)
  
  p <-
    ggplot2::ggplot(
      dplyr::filter(aux,!ids %in% 1), ggplot2::aes(
        x = Variables , y = Abs.importance ,group = ids,key = ids
      )
    ) +
    ggplot2::geom_jitter(height = 0,size = I(2),alpha = 0.3) + ggplot2::facet_grid(node ~
                                                                                     .) + ggplot2::scale_x_discrete(limits = levels(as.factor(aux$var))) + ggplot2::ggtitle("Importance variable for each tree") +
    ggplot2::theme(legend.position = "none",axis.text.x = ggplot2::element_text(angle = 90))
  
  p <-
    p + ggplot2::geom_jitter(
      data = dplyr::filter(aux,ids %in% tree.id), ggplot2::aes(
        x = Variables , y = Abs.importance ,group = ids,colour = "red"
      ),height = 0,size = I(2)
    ) 
  
  
  plotly::ggplotly(p,tooltip = c("y","key"))
  
  
}