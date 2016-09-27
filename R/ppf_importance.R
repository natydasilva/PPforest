#' Importance variable visualization 
#' 
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param ppf is a PPforest object
#' @param global is TRUE if we want to see the global importance of the forest
#' @param weight is TRUE if we want to see a weighted mesure of the forest importance based on out of bag trees errors
#' @return A dotplot with a global measure of importance  variables in the PPforest.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #leukemia data set with all the observations used as training
#' pprf.leukemia <- PPforest(data = leukemia, class = "Type",
#'  size.tr = 1, m = 70, size.p = .4, PPmethod = 'PDA', strata = TRUE)
#' ppf_importance(data = leukemia, class = "Type", pprf.leukemia, global = TRUE, weight = TRUE) 
ppf_importance <- function(data , class, ppf, global = TRUE, weight = TRUE) {
  x <- data %>% dplyr::select(-get(class)) %>%
    apply(2,FUN=scale)
  y <- data %>% dplyr::select(get(class))
  value <- NULL
  variable <- NULL
  node <- NULL
  
  mat.proj <- abs(plyr::ldply(ppf[[8]][[2]], function(x) data.frame(node = 1:dim(x$projbest.node)[1], x$projbest.node)))
  colnames(mat.proj)[-1] <- colnames(dplyr::select(data,-get(class)))
  
  index.part <- plyr::ldply(ppf[[8]][[2]], function(x) data.frame(index = x$Tree.Struct[, 5][x$Tree.Struct[, 
                                                                                                           5] != 0]))
  n.vars <- ncol( mat.proj[, -1] )
  index.mat <- matrix(rep(index.part[, 1], ncol(x) ), ncol = ncol(x), byrow = F)
  
  oob.error.tree <- rep(ppf[[6]], each = length(unique(mat.proj$node)))
  imp.weight <- mat.proj[, -1] * index.mat * (1 - oob.error.tree)
  
  
  mmat.vi <- reshape2::melt(mat.proj, id.vars = "node")
  mat.vi.w <- data.frame(node = mat.proj$node, imp.weight)
  colnames(mat.vi.w)[-1] <- colnames(x)
  mmat.vi.w <- reshape2::melt(mat.vi.w, id.vars = "node")
  
  if(global ){
    if(weight){
      import.vi.wg <- mmat.vi.w %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value)) %>% dplyr::arrange(dplyr::desc(mean))
      import.vi.wg$variable <- with(import.vi.wg, reorder(variable, mean))
      
      a <- ggplot2::ggplot(import.vi.wg, ggplot2::aes(x = mean, y = variable)) + ggplot2::geom_point()+ggplot2::theme(aspect.ratio=1)
      print(import.vi.wg)
      
    }else{
      import.vi <- mmat.vi %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value)) %>% dplyr::arrange(dplyr::desc(mean))
      import.vi$variable <- with(import.vi, reorder(variable, mean))
      
      a <- ggplot2::ggplot(import.vi, ggplot2::aes(x = mean, y = variable)) + ggplot2::geom_point() + 
        ggplot2::facet_wrap(~node)+ ggplot2::theme(aspect.ratio=1)
      print(import.vi) 
    }
  }else{
    if(weight){
      import.vi.w <- mmat.vi.w%>% dplyr::filter( value != 0) %>% dplyr::group_by(node, variable) %>% dplyr::summarise(mean = mean(value)) %>% 
        dplyr::arrange(dplyr::desc(mean))
      
      import.vi.w$variable <- with(import.vi.w, reorder(variable, mean))
      a <- ggplot2::ggplot(import.vi.w, ggplot2::aes(x = mean, y = variable)) + ggplot2::geom_point() + 
        ggplot2::facet_grid(~node) +ggplot2::theme(aspect.ratio=1)
      print(import.vi.w)
      
    }else{
      
      import.vi <- mmat.vi%>% dplyr::filter(value!= 0)%>%dplyr::group_by(node, variable) %>% dplyr::summarise(mean = mean(value)) %>% dplyr::arrange(dplyr::desc(mean))
      
      import.vi$variable <- with(import.vi, reorder(variable, mean))
      
      a <- ggplot2::ggplot(import.vi, ggplot2::aes(x = mean, y = variable)) + ggplot2::geom_point() + 
        ggplot2::facet_grid(~node) +ggplot2::theme(aspect.ratio=1)
      print(import.vi) 
    }
  
  }
  
  plotly::ggplotly(a)
  
}



