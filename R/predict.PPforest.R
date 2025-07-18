#' Predict method for PPforest objects
#'
#' @param object A fitted PPforest object
#' @param newdata A data frame with predictors (same structure as training data)
#' @param rule Split rule used in classification (integer from 1 to 8)
#'  1: mean of two group means 
#'  2: weighted mean of two group means - weight with group size
#'  3: weighted mean of two group means - weight with group sd 
#'  4: weighted mean of two group means - weight with group se 
#'  5: mean of two group medians 
#'  6: weighted mean of two group medians - weight with group size 
#'  7: weighted mean of two group median - weight with group IQR 
#'  8: weighted mean of two group median - weight with group IQR and size
#' @param parallel Logical, whether to use parallel processing
#' @param cores Number of cores to use if parallel = TRUE
#' @param ... Additional arguments (ignored)
#'
#' @return A list with:
#' \describe{
#'   \item{predtree}{Matrix with individual tree predictions}
#'   \item{predforest}{Final predicted classes based on majority vote}
#' }
#' @export
#' @examples 
#' \dontrun{
#' pprf.crab <- PPforest(data = crab, class = 'Type',
#'  xstd = 'min-max', size.tr = 2/3, m = 100, size.p = .4, PPmethod = 'LDA', parallel = TRUE )
#'  
#' predict(pprf.crab, newdata = pprf.crab$test, parallel = TRUE) 
#' 
#' }
predict.PPforest <- function(object, newdata, rule = 1, parallel = TRUE, cores = 2, ...) {
  
  if (!inherits(object, "PPforest")) {
    stop("Object must be of class 'PPforest'")
  }
  
  if (missing(newdata)) {
    stop("'newdata' must be provided for prediction")
  }
  
  if (parallel) {
    doParallel::registerDoParallel(cores)
  }
  
  votes <- plyr::ldply(
    object[[8]],  #  contains the trees in a list
    function(x) {
     pred <- as.numeric(PPforest::PPclassify(Tree.result = x, test.data = newdata, Rule = rule)[[2]])
    },
    .parallel = parallel
  )[, -1]
  
  if (parallel) {
    doParallel::stopImplicitCluster()
  }
  colnames(votes) <- NULL
  
  vote.mat <- as.matrix(votes, ncol = nrow(newdata), byrow = TRUE)
  max.vote <- mvote(as.matrix(vote.mat) ) 
  max.vote_cl <- as.factor(max.vote)
  levels(max.vote_cl) <- object$vote.mat_cl
  
  return(list(predtree = vote.mat, predforest = max.vote, predforest_cl <- max.vote_cl ))
}