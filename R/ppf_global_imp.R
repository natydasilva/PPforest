#' Global importance measure for a PPforest object
#' 
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param ppf is a PPforest object
#' @return Data frame with the global importance measure
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' 
#' pprf.crab <- PPforest(data = crab, std = TRUE, class = 'Type',
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA', parallel = TRUE, cores = 2)
#'  
#' ppf_global_imp(data = crab, class = 'Type', pprf.crab) 
#' 
ppf_global_imp <- function(data, class, ppf) {
    
    value <- NULL
    variable <- NULL
    node <- NULL
    
    y <- data[, class]
    
    mat.proj <- lapply(ppf[["output.trees"]], function(x) {
        data.frame(node = 1:nrow(x[[2]]), abs(x[[2]]))
    }) %>% dplyr::bind_rows()
    
    colnames(mat.proj)[-1] <- colnames(dplyr::select(data, -get(class)))
    
    index <- lapply(ppf[["output.trees"]], function(x) {
        data.frame(index = x$Tree.Struct[, "Index"][x$Tree.Struct[, "Index"] != 0])
    }) %>% dplyr::bind_rows()
    
    n.vars <- ncol(mat.proj) - 1
    index.mat <- matrix(rep(index[, 1], n.vars), ncol = n.vars, byrow = F)
    
    oob.error.tree <- rep(ppf[["oob.error.tree"]], each = length(unique(mat.proj$node)))
    imp.weight <- mat.proj[, -1] * index.mat * (1 - oob.error.tree)
    
    mat.vi.w <- data.frame(node = mat.proj$node, imp.weight)
    colnames(mat.vi.w)[-1] <- colnames(mat.proj[, -1])
    mmat.vi.w <- mat.vi.w %>% tidyr::gather(variable, value, -node)
    
    import.vi.wg <- mmat.vi.w %>% dplyr::group_by(variable) %>% dplyr::filter(value != 0) %>% 
        dplyr::summarise(mean = mean(value)) %>% dplyr::arrange(dplyr::desc(mean)) %>% dplyr::mutate(variable = stats::reorder(variable, 
        mean))
    
    
    import.vi.wg
}



