#' Importance variable visualization 
#' 
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param ppf is a PPforest object
#' @param global is TRUE if we want to see the global importance of the forest
#' @param weight is TRUE if we want to see a weighted mesure of the forest importance based on out of bag trees errors
#' @return A dotplot with a global measure of importance  variables in the PPforest.
#' @param interactive if is TRUE will use plotly to translate an static ggplot object
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' pprf.crab <- PPforest2(data = crab, class = "Type",
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA', strata = TRUE)
#' ppf_importance(data = crab, class = "Type", pprf.crab, global = TRUE,
#'  weight = FALSE, interactive = TRUE) 
ppf_importance <- function(data , class, ppf, global = TRUE, weight = TRUE, interactive) {
  x <- data %>% dplyr::select(-get(class)) %>%
    apply(2, FUN = scale)
  y <- data %>% dplyr::select(get(class))
  value <- NULL
  variable <- NULL
  node <- NULL
  
  mat.proj <- lapply(ppf[[8]], function(x) data.frame(node = 1:dim(x[[2]])[1], abs(x[[2]]))) %>% dplyr::bind_rows()
  colnames(mat.proj)[-1] <- colnames(dplyr::select(data,-get(class)))
  
  index.part <- lapply(ppf[[8]], function(x) data.frame(index = x$Tree.Struct[, 5][x$Tree.Struct[, 5] != 0]))%>%dplyr::bind_rows()
  n.vars <- ncol( mat.proj[, -1] )
  index.mat <- matrix(rep(index.part[, 1], ncol(x)), ncol =ncol(x), byrow = F)
  
  oob.error.tree <- rep(ppf[[6]], each = length(unique(mat.proj$node)))
  imp.weight <- mat.proj[, -1] * index.mat * (1 - oob.error.tree)
  
  
  #mmat.vi <- reshape2::melt(mat.proj, id.vars = "node")
  mmat.vi <- mat.proj %>% tidyr::gather(variable, value, -node )
  mat.vi.w <- data.frame(node = mat.proj$node, imp.weight)
  colnames(mat.vi.w)[-1] <- colnames(x)
  mat.proj %>% tidyr::gather(variable, value, -node )
  mmat.vi.w <- mat.vi.w %>% tidyr::gather(variable, value, -node )
  
  if(global ){
    if(weight){
      import.vi.wg <- mmat.vi.w %>% dplyr::group_by(variable) %>% dplyr::summarise(mean = mean(value)) %>% dplyr::arrange(dplyr::desc(mean))
      import.vi.wg$variable <- with(import.vi.wg, reorder(variable, mean))
      
      a <- ggplot2::ggplot(import.vi.wg, ggplot2::aes(x = mean, y = variable)) + ggplot2::geom_point() + ggplot2::theme(aspect.ratio=1)
      print(import.vi.wg)
      
    }else{
      import.vi <- mmat.vi %>% dplyr::group_by(variable) %>% dplyr::mutate(mean = mean(value)) %>% dplyr::arrange(dplyr::desc(mean))
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
        ggplot2::facet_grid(~node) + ggplot2::theme(aspect.ratio=1)
      print(import.vi) 
    }
  
  }
  
  if(interactive){
    plotly::ggplotly(a)
  }else{
   
    print(a)
    
    }
 
}



