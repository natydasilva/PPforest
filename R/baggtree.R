#' Grow a PPtree_split for each bootstrap sample
#' 
#' For each bootstrap sample grow a projection persuit tree (PPtree object).
#' @importFrom magrittr %>%
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable.
#' @param m is the number of bootstrap replicates, this corresponds with the number of trees to grow. To ensure that each observation is predicted a few times we have to select this number no too small. \code{m = 500} is by default.
#' @param size.p proportion of random sample variables in each split.
#' @param PPmethod is the projection pursuit index to be optimized, options LDA or PDA, by default it is LDA.
#' @param lambda a parameter for PDA index
#' @param ... arguments to be passed to methods
#' @return data frame with trees_pp output for all the bootstraps samples.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #leukemia data set
#' leukemia.trees <- baggtree(data = leukemia, class = "Type", 
#'  m =  70, PPmethod = 'PDA', lambda = .1, size.p = 0.4 ) 
#' str(leukemia.trees, max.level = 1)

baggtree <- function(data , class , m = 500, PPmethod = "LDA", lambda = 0.1, size.p = 1){
   bootsam <- NULL
   . <- NULL
boottree <- function(data , class , PPmethod, lambda, size.p ){
  
  origclass <- data[,class]
  origdata <- data[ , setdiff( colnames( data ), class )]
   bt1 <- boot(as.matrix(as.numeric(as.factor( origclass ) ) ), as.matrix( origdata ) )
 
   tree <- PPtree_split2( class, data = data[ ( bt1 + 1 ), ] , PPmethod, lambda ,size.p )
   tree
   
}

doMC::registerDoMC(2)

# getDoParWorkers()
## [1] 4

data.frame(bootsam = 1:m) %>% 
    plyr::dlply( plyr::.(bootsam), function(x) boottree(data , class, PPmethod , lambda , size.p ) ,  .parallel = TRUE)


}







