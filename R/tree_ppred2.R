#' Obtain predicted class for new data using PPforest 
#' 
#' Vector with predicted values from a PPforest.
#' @param xnew data frame with explicative variables used to get new predicted values.
#' @param output.tree trees classifiers from trees_pp function or PPforest object
#' @param cores The number of cores to use for parallel execution. By default is 2 cores.
#' @return predicted values from PPforest  
#' @export
#' @importFrom magrittr %>%
#' @examples 
#' training.id <- train_fn(data = leukemia, class = "Type", size.p = 0.9)
#' leukemia.trees <- baggtree(data = leukemia[training.id$id,], class = "Type", 
#' m =  70, PPmethod = 'PDA', lambda = .1, size.p = 0.4 ) 
#' pr <- tree_ppred2( xnew = leukemia[-training.id$id, -1] , leukemia.trees)
tree_ppred2 <- function(xnew, output.tree, cores = 2) {
  
  doMC::registerDoMC(cores)
  
  votes <- plyr::ldply(output.tree, function(x) as.numeric(PPtreeViz::PP.classify(test.data = xnew, Tree.result = x[[1]], Rule = 1)[[2]]) ,  .parallel = TRUE)
  
  
 
  max.vote <-  as.numeric((votes %>% tidyr::gather(vars, votes,-bootsam) %>% 
                            plyr::ddply( plyr::.(vars), function(x){
                              t1 <- table(x$votes)
                              names(t1)[which.max(t1)]
                              }, .parallel = TRUE))[,2])

  colnames(votes) <- NULL
  vote.mat <- as.matrix(votes[,-1], ncol = dim(xnew)[[1]], byrow = T)
  
  # max.vote <- as.numeric(apply(vote.mat, 2, function(x){
  #   
  #   t1 <- table(x)
  #   
  #   names(t1)[which.max(t1)]
  # }))
  # 
  
  return(list(vote.mat, max.vote))
} 


