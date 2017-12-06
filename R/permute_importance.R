#' Obtain the permuted importance variable measure
#' 
#' @param ppf is a PPforest object
#' @return A data frame with permuted importance measures, imp is the permuted importance measure defined in Brieman paper,
#' imp2 is the permuted importance measure defined in randomForest package, the standard deviation (sd.im and sd.imp2) for each measure is computed and the 
#' also the standardized mesure. 
#' @useDynLib PPforest
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' pprf.crab <- PPforest(data = crab, class = 'Type',
#' std = TRUE, size.tr = 1, m = 100, size.p = .4, PPmethod = 'LDA', parallel = TRUE, core = 2)
#' permute_importance(ppf = pprf.crab) 
#' 
permute_importance <- function(ppf) {
    sd <- NULL
    imp <- NULL
    imp2 <- NULL
    sd.imp <- NULL
    sd.imp2 <- NULL
    
    train <- as.matrix(ppf$train[, -which(colnames(ppf$train) == ppf$class.var)])
    classes <- as.integer(unlist(ppf$train[, which(colnames(ppf$train) == ppf$class.var)]))
    oobid <- apply(ppf$oob.obs, 1, function(x) which(x == 1) - 1)
    
    permute <- oobid %>% lapply(function(x) sample(x, length(x)))
    trees <- ppf[[8]]
    noob <- as.integer(lapply(oobid, length))
    TRstrL <- trees %>% lapply(function(x) as.matrix(x[[1]]))
    
    TRsplL <- trees %>% lapply(function(x) as.matrix(x[[3]]))
    TRprnodeL <- trees %>% lapply(function(x) as.matrix(x[[2]]))
    
    
    corr.oob.per <- imposoon(train, classes, oobid, permute, trees, noob, TRstrL, TRsplL, TRprnodeL)
    
    rank.var <- t(apply(corr.oob.per, 1, rank, ties.method = "random"))
    corr.oob <- (1 - ppf$oob.error.tree) * unlist(lapply(oobid, length))
    n.oob <- unlist(lapply(permute, length))
    
    # imp is the permuted importance using accuracy in the randomforest paper imp2 is the
    # permuted importance using the error like in randomForest package
    imp.pl <- data.frame(nm = colnames(ppf$train)[colnames(ppf$train) != ppf$class.var], imp = apply(corr.oob.per, 
        2, function(x) mean((corr.oob - x))), sd.imp = apply(corr.oob.per, 2, function(x) sd(corr.oob - 
        x)), imp2 = apply(corr.oob.per, 2, function(x) mean(((1 - x/n.oob) - ppf$oob.error.tree))), 
        sd.imp2 = apply(corr.oob.per, 2, function(x) sd(((1 - x/n.oob) - ppf$oob.error.tree)))) %>% 
        
    dplyr::mutate(imp2.std = imp2/sd.imp2, imp.std = imp/sd.imp) %>% dplyr::arrange(imp)
    imp.pl$nm <- factor(imp.pl$nm, levels = imp.pl[order(imp.pl$imp2), "nm"])
    
    imp.pl
}


