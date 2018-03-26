#' Global importance measure for a PPforest object as the average IMP PPtree measure over all the trees 
#' in the forest
#' 
#' @param ppf is a PPforest object
#' @param class A character with the name of the class variable. 
#' @return Data frame with the global importance measure
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' 
#' pprf.crab <- PPforest(data = crab, std =TRUE, class = 'Type',
#'  size.tr = 1, m = 100, size.p = .5, PPmethod = 'LDA')
#'  ppf_avg_imp(pprf.crab, 'Type') 
#'  
ppf_avg_imp <- function(ppf, class) {
    node.id <- NULL
    nodecl <- NULL
    node <- NULL
    clnd <- NULL
    impaux <- NULL
    Class <- NULL
    variable <- NULL
    value <- NULL
    tr <- NULL
    
    nn <- data.frame(nn = 1:sum(ppf[["output.trees"]][[1]]$Tree.Struct[, 4] != 0))
    nodecl <- function(x) {
        aux <- node_data(ppf = ppf, x)
        aux$node.id <- as.factor(aux$node.id)
        aux %>% dplyr::group_by(node.id) %>% dplyr::summarise(clt = length(unique(Class)))
    }
    
    
    mat.proj <- lapply(ppf[["output.trees"]], function(x) {
        data.frame(node = 1:sum(x$Tree.Struct[, 4] != 0), abs(x[[2]]))
    }) %>% dplyr::bind_rows()
    
    
    infond <- apply(data.frame(1:ppf$n.tree), 1, function(x) nodecl(x)$clt)  #info to weight importance
    info <- data.frame(clnd = matrix(infond, ncol = 1, nrow = ppf$n.tree * nrow(infond), byrow = T))
    colnames(mat.proj)[-1] <- colnames(dplyr::select(ppf$train, -class))
    
    mat.proj %>% dplyr::bind_cols(clnd = info) %>% dplyr::mutate(tr = rep(1:ppf$n.tree, dim(nn)[1])) %>% 
        tidyr::gather(variable, value, -node, -tr, -clnd) %>% dplyr::mutate(impaux = value/clnd) %>% 
        dplyr::group_by(variable, tr) %>% dplyr::summarise(mean = sum(impaux)) %>% dplyr::group_by(variable) %>% 
        dplyr::summarise(mean = mean(mean)) %>% dplyr::arrange(dplyr::desc(mean)) %>% dplyr::mutate(variable = stats::reorder(variable, 
        mean))
    
}
