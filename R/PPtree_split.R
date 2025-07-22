#' Projection pursuit classification tree with random variable selection in each split
#' 
#' Find tree structure using  projection pursuit indices of classification in each split.
#' @usage PPtree_split(form, data, PPmethod='LDA', 
#' size.p=1,  lambda = 0.1,...) 
#' @param form A character with the name of the class variable.
#' @param data Data frame with the complete data set.
#' @param PPmethod index to use for projection pursuit: 'LDA' and 'PDA'
#' @param size.p proportion of variables randomly sampled in each split, default is 1, returns a PPtree.
#' @param lambda penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.
#' @param ... arguments to be passed to methods
#' @return An object of class \code{PPtreeclass} with components
#' \item{Tree.Struct}{Tree structure of projection pursuit classification tree}
#' \item{projbest.node}{1-dim optimal projections of each split node}
#' \item{splitCutoff.node}{cutoff values of each split node}
#' \item{origclass_num}{original class numeric} 
#' \item{origdata}{original data}
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @useDynLib PPforest
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords tree
#' @examples
#' #crab data set
#' 
#' Tree.crab <- PPtree_split('Type~.', data = crab, PPmethod = 'LDA', size.p = 0.5)
#' Tree.crab
#'
PPtree_split <- function(form, data, PPmethod = "LDA", size.p = 1, lambda = 0.1, ...) {
    
    formula <- stats::as.formula(form)
    mf <- stats::model.frame(formula, data = data)
    origclass <- stats::model.response(mf)
    
    cls <- all.vars(formula)[[1]]
    
    origdata <- data[, -which(colnames(data) %in% cls)]
    origdata <- as.matrix(origdata)
    pp <- ncol(origdata)
    origclass_num <- as.numeric(as.factor(origclass))
    
    
    g <- table(origclass_num)
    G <- length(g)
    
    Tree.final <- treeconstruct(origclass_num, origdata, Treestruct = cbind(1:(2 * G - 1), matrix(0, 
        ncol = 4, nrow = 2 * G - 1)), id = 0, rep = 1, rep1 = 2, rep2 = 1, projbestnode = matrix(0, 
        ncol = pp, nrow = 1), splitCutoffnode = matrix(0, ncol = 8, nrow = 1), PPmethod, lambda, 
        size.p)
    
    
    
    Tree.Struct <- Tree.final$Treestruct
    colnames(Tree.Struct) <- c("id", "L.node.ID", "R.F.node.ID", "Coef.ID", "Index")
    # projbest.node <- Tree.final$projbestnode[-1, ]
    
    
    if (nrow(Tree.final$splitCutoffnode) == 2) {
        splitCutoff.node <- data.frame(splitCutoffnode = t(Tree.final$splitCutoffnode[-1, ]))
        colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
        projbest.node <- t(as.matrix(Tree.final$projbestnode[-1, ]))
        
    } else {
        splitCutoff.node <- data.frame(splitCutoffnode = Tree.final$splitCutoffnode[-1, ])
        colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
        projbest.node <- Tree.final$projbestnode[-1, ]
    }
    treeobj <- list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, splitCutoff.node = splitCutoff.node, 
        origclass = origclass, origclass_num = origclass_num, origdata = origdata)
    
    class(treeobj) <- append(class(treeobj), "PPtreeclass")
    
    return(treeobj)
}
