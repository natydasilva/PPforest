#' Grow a PPtree_split for each bootstrap sample
#' 
#' For each bootstrap sample grow a projection persuit tree (PPtree object).
#' @importFrom magrittr %>%
#' @param data.b are the bootstrap samples from training set.
#' @param size.p proportion of random sample variables in each split.
#' @param PPmethod is the projection pursuit index to be optimized, options LDA or PDA, by default it is LDA.
#' @param lambda a parameter for PDA index
#' @param ... arguments to be passed to methods
#' @return data frame with trees_pp output for all the bootstraps samples.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #leukemia data set
#' leukemia.b <- ppf_bootstrap(data = leukemia, class = "Type", m = 70) 
#' leukemia.trees <- trees_pp(data.b = leukemia.b, size.p = .4, PPmethod = 'PDA', lambda = .1) 
#' str(leukemia.trees,max.level=1)
trees_pp <- function(data.b, size.p = 0.9, PPmethod = "LDA", lambda = 0.1, ...) {
    . <- NULL

    names(data.b)[1] <- "class"
    if (PPmethod == "LDA") {
        trees <- data.b %>% dplyr::do(tr = PPtree_split(class~., data = .,  PPmethod = "LDA",  size.p = size.p, 
            ...))
        
    } else {
        trees <- data.b %>% dplyr::do(tr = PPtree_split(class~. ,  data = ., PPmethod = "PDA", size.p = size.p, 
            lambda, ...))
    }
    trees
} 
