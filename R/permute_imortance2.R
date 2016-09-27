#' Obtain the importance variable measure based in permutation
#' 
#' @param ppf is a PPforest object
#' @return permuted importance measure
#' @useDynLib PPforest2
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' pprf.crab<- PPforest(data = crab, class = "Type",
#'  size.tr = 1, m = 100, size.p = .4, PPmethod = 'LDA', strata = TRUE)
#' permute_importance2(ppf = pprf.crab) 
permute_importance2 <- function(ppf){
  train <- as.matrix(ppf$train[, -which(colnames(ppf$train) == ppf$class.var)])
  classes <- as.integer(ppf$train[,which(colnames(ppf$train) == ppf$class.var)])
  oobid <- apply(ppf$oob.obs, 1, function(x) which(x) - 1) 
  permute <- plyr:: llply(oobid , function(x) sample(x,length(x))) 
  trees <- ppf[[8]][[2]] #list of trees
  noob <- as.integer( lapply(oobid,length) ) 
  TRstrL <- plyr::llply(trees, function(x) as.matrix(x$Tree.Struct))
  TRsplL <- plyr::llply(trees, function(x) as.matrix(x$splitCutoff.node))
  TRprnodeL <- plyr::llply(trees, function(x) as.matrix(x$projbest.node))
  
  
  corr.oob.per <- imposoon(train, classes, oobid, permute, trees, noob , TRstrL,TRsplL,TRprnodeL)
  
  rank.var <- t(apply(corr.oob.per, 1, rank, ties.method = "random"))
  corr.oob <- (1-ppf$oob.error.tree) * unlist(lapply(oobid, length))
  n.oob <-  unlist(lapply(permute, length))
  imp.pl <- data.frame( nm = colnames(ppf$train)[colnames(ppf$train)!=ppf$class.var], 
                        imp=apply(corr.oob.per, 2, function(x) mean((corr.oob-x))), 
                        imp2=apply(corr.oob.per, 2, function(x) mean(((1-x/n.oob)-ppf$oob.error.tree))), 
                        sd =apply(corr.oob.per, 2, function(x) sd( ((1-x/n.oob)-ppf$oob.error.tree)))) %>%
    dplyr::mutate(imp.sd =imp2/sd)%>% dplyr::arrange(imp)
  
  imp.pl$nm <- factor(imp.pl$nm, levels= imp.pl[order( imp.pl$imp2), "nm"])
  
  imp.pl
}


