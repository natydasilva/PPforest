#' Obtain predicted class for new data from baggtree function 
#' 
#' Vector with predicted values from a baggtree function.
#' @param object trees classifiers from trees_pp function or PPforest object
#' @param xnew data frame with explicative variables used to get new predicted values.
#' @param ... arguments to be passed to methods
#' @return predicted values form PPforest
#' @export
#' @importFrom magrittr %>%
#' @examples 
#' crab.trees <- baggtree(data = crab, class = "Type", 
#' m =  200, PPmethod = 'LDA', lambda = .1, size.p = 0.4 )
#'  
#' pr <- trees_pred(  crab.trees,xnew = crab[, -1] )
#' 
#' 
#' pprf.crab <- PPforest(data = crab, class = "Type",
#'  std = TRUE, size.tr = 1, m = 100, size.p = .4, PPmethod = 'LDA' )
#'  
#' trees_pred(pprf.crab, xnew = pprf.crab$train[, -1] )
#' 
trees_pred <- function( object, xnew, parallel = FALSE, cores = 2, ...) {

       if(class(object) == "PPforest"){
          # votes <- plyr::ldply(object[[8]], function(x)
          #   as.numeric(PPclassify2(Tree.result = x, test.data = xnew, Rule = 1)[[2]]) )

         votes <- object[[8]] %>% purrr::map_df(.f = function(x)
           as.numeric(PPclassify2(Tree.result = x, test.data = xnew, Rule = 1)[[2]])) %>%t() 

        }else{
 #votes <- plyr::ldply(object, function(x) as.numeric(PPclassify2(Tree.result = x[[1]], test.data = xnew, Rule = 1)[[2]]) )
          votes <- object %>% purrr::map_df(.f = function(x)
            as.numeric(PPclassify2(Tree.result = x[[1]], test.data = xnew, Rule = 1)[[2]])) %>%t()

          
       }
  max.vote <- mvote(as.matrix((votes)))
  
  colnames(votes) <- NULL
  vote.mat <- as.matrix(votes, ncol = dim(xnew)[[1]], byrow = T)
  
  return(list(vote.mat, max.vote))
} 


