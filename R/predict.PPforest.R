#' Obtain predicted class for new data using PPforest 
#' 
#' Vector with predicted values from a PPforest.
#' @param object  an object of class PPforest
#' @param class A character with the name of the class variable. 
#' @param test.data the test dataset
#' @param Rule split rule \itemize{ \item{1: mean of two group means} \item{ 2: weighted mean of two group means - weight with group size} \item{3: weighted mean of two group means - weight with group sd} \item{ 4: weighted mean of two group means - weight with group se} \item{5: mean of two group medians} \item{ 6: weighted mean of two group medians - weight with group size} \item{ 7: weighted mean of two group median - weight with group IQR} \item{8: weighted mean of two group median - weight with group IQR and size} }
#' @param true.class true class of test dataset if available
#' @param ... arguments to be passed to methods
#' @return predicted values form PPforest  
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #use 2/3 of the observations as training data and 1/3 as test
#' pprf.leukemia <- PPforest(data = leukemia, class = "Type",
#'  size.tr = 2/3, m = 70, size.p = .4, PPmethod = 'PDA', strata = TRUE)
#' pr.leukemia <- predict( object = pprf.leukemia, class= "Type", test.data = pprf.leukemia$test)
#' pr.leukemia
predict.PPforest <- function(object, class, test.data, Rule, true.class = NULL, ...) {
    . <- NULL
    Var2 <- NULL
    votes <- object[[8]] %>% dplyr::do(tr = PPtreeViz::PP.classify(test.data = test.data, Tree.result = .$tr, Rule = 1, 
        ...))
    
    out <- votes %>% dplyr::do(pred = .$tr[[2]])
    
    
    vote.mat <- matrix(unlist(out$pred), ncol = length(votes[[1]][[1]][[2]]), byrow = T)
    max.vote <- apply(vote.mat, 2, function(x) {
        t1 <- table(x)
        names(t1)[which.max(t1)]
    })
    
    votes2 <- matrix(0, ncol = length(unique(object$train[, class])), nrow = nrow(test.data))
    colnames(votes2) <- levels(object$train[, class])
    dt <- reshape2::melt(vote.mat)
    aux <- plyr::ddply(dt, plyr::.(Var2), function(x) table(x$value))
    
    
    vote.matrix.prop <- aux[,-1]/rowSums(aux[,-1])
    
    list(vote.matrix.prop, max.vote)
} 



