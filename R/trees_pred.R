#' Obtain predicted class for new data from baggtree function or PPforest 
#' 
#' @param object  Projection pursuit classification forest structure from PPforest or baggtree
#' @param xnew data frame with explicative variables used to get new predicted values.
#' @param parallel logical condition, if it is TRUE then  parallelize the function
#' @param cores number of cores used in the parallelization
#' @param rule split rule 1: mean of two group means 2: weighted mean of two group means - weight with group size 3: weighted mean of two group means - weight with group sd 4: weighted mean of two group means - weight with group se 5: mean of two group medians 6: weighted mean of two group medians - weight with group size 7: weighted mean of two group median - weight with group IQR 8: weighted mean of two group median - weight with group IQR and size
#' @return predicted values from PPforest or baggtree
#' @export
#' @importFrom magrittr %>%
#' @examples 
#' \dontrun{
#' crab.trees <- baggtree(data = crab, class = 'Type', 
#' m =  200, PPmethod = 'LDA', lambda = .1, size.p = 0.4 )
#'  
#' pr <- trees_pred(  crab.trees,xnew = crab[, -1], parallel= FALSE, cores = 2)
#' 
#' pprf.crab <- PPforest(data = crab, class = 'Type',
#'  std = FALSE, size.tr = 2/3, m = 100, size.p = .4, PPmethod = 'LDA', parallel = TRUE )
#'  
#' trees_pred(pprf.crab, xnew = pprf.crab$test ,parallel = TRUE)
#' }
trees_pred <- function(object, xnew, parallel = FALSE, cores = 2, rule = 1) {

        if (parallel) {
         
         doParallel::registerDoParallel(cores)
          
        }
        if (inherits(object,"PPforest")) {
            votes <-  plyr::ldply(object[[8]], function(x) as.numeric(PPforest::PPclassify2(Tree.result = x, 
                test.data = xnew, Rule = rule)[[2]]), .parallel = parallel)[, -1]
         
          
        } else {
            votes <- plyr::ldply(object, function(x) as.numeric(PPclassify2(Tree.result = x[[1]], 
                test.data = xnew, Rule = rule)[[2]]), .parallel = parallel)[, -1]
            
            
        }
  
 if(parallel){
   doParallel::stopImplicitCluster()
 }
    max.vote <- mvote(as.matrix((votes)))
    
    colnames(votes) <- NULL
    vote.mat <- as.matrix(votes, ncol = dim(xnew)[[1]], byrow = T)
    result <- list(vote.mat, max.vote)
    names(result) <- c("predtree", "predforest")
    return(result)
}


