#' Obtain predicted class for new data using PPforest 
#' 
#' Vector with predicted values from a PPforest.
#' @param object trees classifiers from trees_pp function or PPforest object
#' @param xnew data frame with explicative variables used to get new predicted values.
#' @param parallel if TRUE, apply function in parallel
#' @param cores The number of cores to use for parallel execution. By default is 2 cores.
#' @param ... arguments to be passed to methods
#' @return predicted values form PPforest
#' @export
#' @importFrom magrittr %>%
#' @examples 
#' crab.trees <- baggtree(data = crab, class = "Type", 
#' m =  200, PPmethod = 'PDA', lambda = .1, size.p = 0.4 ) 
#' pr <- trees_pred(  crab.trees,xnew = crab[, -1] )
trees_pred <- function( object, xnew, parallel = FALSE, cores = 2, ...) {
  
  doMC::registerDoMC(cores)
                                                                     
  votes <- plyr::ldply(object, function(x) as.numeric(PPclassify2(Tree.result = x[[1]], test.data = xnew, Rule = 1)[[2]]) ,  .parallel = parallel)
  
 
  max.vote <- mvote(as.matrix((votes[ , -1])))
  
  colnames(votes) <- NULL
  vote.mat <- as.matrix(votes[,-1], ncol = dim(xnew)[[1]], byrow = T)
  

  return(list(vote.mat, max.vote))
} 


