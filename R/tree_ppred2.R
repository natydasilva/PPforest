#' Obtain predicted class for new data using PPforest 
#' 
#' Vector with predicted values from a PPforest.
#' @param xnew data frame with explicative variables used to get new predicted values.
#' @param output.tree trees classifiers from trees_pp function or PPforest object
#' @param parallel if TRUE, apply function in parallel
#' @param cores The number of cores to use for parallel execution. By default is 2 cores.
#' @return predicted values from PPforest  
#' @export
#' @importFrom magrittr %>%
#' @examples 
#' training.id <- train_fn(data = leukemia, class = "Type", size.p = 0.9)
#' leukemia.trees <- baggtree(data = leukemia[training.id$id,], class = "Type", 
#' m =  70, PPmethod = 'PDA', lambda = .1, size.p = 0.4 ) 
#' pr <- tree_ppred2( xnew = leukemia[-training.id$id, -1] , leukemia.trees)
tree_ppred2 <- function(xnew, output.tree, parallel = FALSE, cores = 2) {
  
  doMC::registerDoMC(cores)
                                                                     
  votes <- plyr::ldply(output.tree, function(x) as.numeric(PPclassify2(Tree.result = x[[1]], test.data = xnew, Rule = 1)[[2]]) ,  .parallel = parallel)
  
 
  max.vote <- mvote(as.matrix((votes[ , -1])))
  
  colnames(votes) <- NULL
  vote.mat <- as.matrix(votes[,-1], ncol = dim(xnew)[[1]], byrow = T)
  

  return(list(vote.mat, max.vote))
} 


