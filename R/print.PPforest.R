#' Print PPforest object
#'
#' @param x is a PPforest class object
#' @param ... additional parameter
#' @return printed results for PPforest object
#' @export
#' 
print.PPforest <- function(x, ...) {
    cat("\nCall:\n", deparse(x$call), "\n")
    cat("               Type of random forest: ", x$type, "\n", sep = "")
    cat("                     Number of trees: ", x$n.tree, "\n", sep = "")
    cat("No. of variables tried at each split: ", x$n.var, "\n\n", sep = "")
    
    cat("        OOB estimate of  error rate: ", round(x$oob.error.forest * 100, digits = 2), 
        "%\n", sep = "")
    cat("Confusion matrix:\n")
    print(x$confusion)
    
}
