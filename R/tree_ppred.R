#' Obtain predicted class for new data using PPforest 
#' 
#' Vector with predicted values from a PPforest.
#' @param xnew data frame with explicative variables used to get new predicted values.
#' @param output.tree trees classifiers from trees_pp function or PPforest object
#' @param ... arguments to be passed to methods
#' @return predicted values from PPforest  
#' @export
#' @importFrom magrittr %>%
#' @examples 
#' training.id <- train_fn(data = leukemia, class = "Type", size.p = 0.9)
#' leukemia.b <- ppf_bootstrap(data = leukemia[training.id$id,], class = "Type", m = 70) 
#' leukemia.trees <- trees_pp(data.b = leukemia.b, size.p = .4, PPmethod = 'LDA')
#' pr <- tree_ppred( xnew = leukemia[-training.id$id, -1] , leukemia.trees)
tree_ppred <- function(xnew, output.tree, ...) {
    . <- NULL
    
    votes <- output.tree %>% dplyr::do(tr = PPtreeViz::PP.classify(test.data = xnew, Tree.result = .$tr, Rule = 1))
    
    out <- votes %>% dplyr::do(pred = .$tr[[2]])
    
    vote.mat <- matrix(unlist(out$pred), ncol = dim(xnew)[[1]], byrow = T)
    
    
    max.vote <- apply(vote.mat, 2, function(x) {
        t1 <- table(x)
        names(t1)[which.max(t1)]
    })
    
    
    return(list(out, vote.mat, max.vote))
} 
