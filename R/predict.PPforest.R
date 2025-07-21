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
#' set.seed(1883)
#' train <- sample(1:nrow(crab), nrow(crab)*.7)
#' crab_train <- data.frame(crab[train, ])
#' crab_test <- data.frame(crab[-train, ])
#' pprf.crab <- PPforest(data = crab_train, class = 'Type',
#'  xstd = 'min-max', size.tr = 0.7, m = 200, size.p = .4, PPmethod = 'LDA', parallel = TRUE )
#'  
#' pred <- predict(pprf.crab, newdata = crab_test[,-1], parallel = TRUE) 
#' 1 - sum(as.numeric(as.factor(crab_test[,1])) == as.numeric(as.factor(pred[[2]])))/length(pred[[2]])
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
  
   if(object$xstd != "no"){
   if (object$xstd == "min-max") {
      dataux <- newdata %>%  tibble::as_tibble()
      mincol <- dataux |> apply( 2, min)
      maxmincol <- dataux |> apply(2, function(x) {
        max(x) - min(x)
      })
    }
  #   if (object$xstd  == "quant") {
  #     dataux <- xnew |>  tibble::as_tibble()
  #     mincol <- dataux |> apply( 2, stats::quantile, 0.05) 
  #     maxmincol <- dataux |> apply(2, function(x) {
  #       stats::quantile(x, 0.95) - stats::quantile(x, 0.05)
  #     }) 
  #   }
  #   
  #   if (object$xstd  == 'scale') {
  #     dataux <- xnew %>% 
  #       apply(2, FUN = scale) %>%
  #       tibble::as_tibble()
  #     
  #   }
  #   if(object$xstd  %in% c('min-max', 'quant')) {
     datauxscale <- data.frame((dataux  - matrix(mincol, nrow(dataux), ncol(dataux), byrow = T)) / matrix(maxmincol, nrow(dataux), ncol(dataux), byrow = T))
    newdata <- datauxscale
  #     }else{
  #     xnew <- dataux
  #     #colnames(data)[1] <- object$class.var
  #     } 
}
  votes <- plyr::ldply(
    object$output.trees,  #  contains the trees in a list
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