#' For each bootstrap sample grow a projection pursuit tree (PPtree object).
#' 
#' @importFrom magrittr %>%
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable.
#' @param m is the number of bootstrap replicates, this corresponds with the number of trees to grow. To ensure that each observation is predicted a few times we have to select this number no too small. \code{m = 500} is by default.
#' @param size.p proportion of random sample variables in each split.
#' @param PPmethod is the projection pursuit index to be optimized, options LDA or PDA, by default it is LDA.
#' @param lambda a parameter for PDA index
#' @param parallel logical condition, if it is TRUE then  parallelize the function
#' @param cores number of cores used in the parallelization
#' @return data frame with trees_pp output for all the bootstraps samples.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set
#' crab.trees <- baggtree(data = crab, class = 'Type',
#' m =  200, PPmethod = 'LDA', lambda = .1, size.p = 0.5 , parallel = TRUE, cores = 2)
#' str(crab.trees, max.level = 1)
#' 
baggtree <- function(data, class, m = 500, PPmethod = "LDA", lambda = 0.1, size.p = 1, parallel = FALSE, 
    cores = 2) {
    
    bootsam <- NULL
    . <- NULL
    
    boottree <- function(data, class, PPmethod, lambda, size.p) {
        
        origclass <- data[, class]
        origdata <- data[, setdiff(colnames(data), class)]
        origdata <- as.matrix(origdata)
        origclass <- as.numeric(as.factor(unlist(origclass)))
        
        bt1 <- boot(as.matrix(origclass), origdata)
        
        f <- stats::as.formula(paste(class, "~.", sep = ""))
        tree <- PPtree_split(f, data = data[(bt1 + 1), ], PPmethod, lambda, size.p = size.p)
        
        list(tree, bt1)
        
    }
    
  
        if(parallel) {
    data <- data
    
    doParallel::registerDoParallel(cores)
            
           pp <- plyr::dlply(dplyr::data_frame(bootsam = 1:m), plyr::.(bootsam), function(x) boottree(data, 
                class, PPmethod, lambda, size.p), .parallel = parallel)
          
           doParallel::stopImplicitCluster()
        }else{
            pp <- plyr::dlply(dplyr::data_frame(bootsam = 1:m), plyr::.(bootsam), function(x) boottree(data, 
                class, PPmethod, lambda, size.p), .parallel = parallel)
            
        }
    return(pp)
}

