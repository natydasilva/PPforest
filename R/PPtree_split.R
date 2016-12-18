#' Projection pursuit classification tree with random variable selection in each split
#' 
#' Find tree structure using various projection pursuit indices of classification in each split.
#' @usage PPtree_split(form, data, PPmethod='LDA', 
#' size.p=1,  lambda=0.1,...) 
#' @param form A character with the name of the class variable.
#' @param data Data frame with the complete data set.
#' @param PPmethod index to use for projection pursuit: 'LDA', 'PDA'
#' @param size.p proportion of variables randomly sampled in each split, default is 1, returns a PPtree.
#' @param lambda penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.
#' @param ... arguments to be passed to methods
#' @return An object of class \code{PPtreeclass} with components
#' \item{Tree.Struct}{Tree structure of projection pursuit classification tree}
#' \item{projbest.node}{1-dim optimal projections of each split node}
#' \item{splitCutoff.node}{cutoff values of each split node}
#' \item{origclass}{original class} 
#' \item{origdata}{original data}
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK (2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @useDynLib PPforest2
#' @importFrom Rcpp evalCpp
#' @export
#' @keywords tree
#' @examples
#' #leukemia data set
#' Tree.leukemia <- PPtree_split("Type~.", data = leukemia,PPmethod = "PDA", size.p = 0.9)
#' Tree.leukemia
#' #crab data set
#' Tree.crab <- PPtree_split("Type~.", data = crab, PPmethod = "LDA", size.p = 1)
#' Tree.crab
PPtree_split <- function(form, data,  PPmethod = "LDA", size.p = 1,  lambda = 0.1, ...) {
     TOL <- NULL
     formula <- stats::as.formula(form)
     mf <- stats::model.frame(formula, data = data)
     origclass <- stats::model.response(mf)
    
     
 
    cls <- all.vars(formula)[[1]]
    
    origdata <- data[,-which(colnames(data)%in%cls)]
    origdata <- as.matrix(origdata)
    pp <- ncol(origdata)
    origclass <- as.numeric(as.factor(origclass))

    Find.proj <- function(origclass, origdata, PPmethod,size.p,lambda,...) {
    findprojwrap(origclass,origdata, PPmethod,
                sizep= size.p, lambda  )
    }
    

    Tree.construct <- function(origclass, origdata, Tree.Struct, 
                               id, rep, rep1, rep2, projbest.node, splitCutoff.node, 
                               PPmethod, lambda, size.p,...) {
      origclass <- as.integer(origclass)
      n <- nrow(origdata)
      g <- table(origclass)
      G <- length(g)
      if (length(Tree.Struct) == 0) {
        Tree.Struct <- matrix(1:(2 * G - 1), ncol = 1)
        Tree.Struct <- cbind(Tree.Struct, 0, 0, 0, 0)
      }
      if (G == 1) {
        Tree.Struct[id, 3] <- as.numeric(names(g))
        list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
             splitCutoff.node = splitCutoff.node, rep = rep, 
             rep1 = rep1, rep2 = rep2)
      }else{
        Tree.Struct[id, 2] <- rep1
        rep1 <- rep1 + 1
        Tree.Struct[id, 3] <- rep1
        rep1 <- rep1 + 1
        Tree.Struct[id, 4] <- rep2
        rep2 <- rep2 + 1
        
        a <- findprojwrap(origclass,origdata, PPmethod,
                          sizep= size.p, lambda  )
        splitCutoff.node <- rbind(splitCutoff.node, t(a$C))
        Tree.Struct[id, 5] <- a$Index
        projbest.node <- rbind(projbest.node, t(a$Alpha))
        t.class <- origclass
        t.data <- origdata
        t.class <- t.class * a$IOindexL
        t.n <- length(t.class[t.class == 0])
        t.index <- sort.list(t.class)
        t.index <- sort(t.index[-(1:t.n)])
        t.class <- t.class[t.index]
        t.data <- origdata[t.index, ]
        
        print(Tree.Struct[id, 2])
        
        b <- Tree.construct(t.class, t.data, Tree.Struct, 
                            id=Tree.Struct[id, 2], rep, rep1, rep2, projbest.node, 
                            splitCutoff.node, PPmethod, lambda,size.p, 
                            ...)
        Tree.Struct <- b$Tree.Struct
        projbest.node <- b$projbest.node
        splitCutoff.node <- b$splitCutoff.node
        rep <- b$rep
        rep1 <- b$rep1
        rep2 <- b$rep2
        t.class <- origclass
        t.data <- origdata
        t.class <- (t.class * a$IOindexR)
        t.n <- length(t.class[t.class == 0])
        t.index <- sort.list(t.class)
        t.index <- sort(t.index[-(1:t.n)])
        t.class <- t.class[t.index]
        t.data <- origdata[t.index, ]
        n <- nrow(t.data)
        G <- length(table(t.class))
        print(Tree.Struct[id, 3])
        b <- Tree.construct(t.class, t.data, Tree.Struct, 
                            Tree.Struct[id, 3], rep, rep1, rep2, projbest.node, 
                            splitCutoff.node, PPmethod,  lambda,size.p, 
                            ...)
        
        Tree.Struct <- b$Tree.Struct
        projbest.node <- b$projbest.node
        splitCutoff.node <- b$splitCutoff.node
        rep <- b$rep
        rep1 <- b$rep1
        rep2 <- b$rep2
      }
      list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
           splitCutoff.node = splitCutoff.node, rep = rep, rep1 = rep1, 
           rep2 = rep2)
    }
    
    splitCutoff.node <- NULL
    projbest.node <- NULL
    Tree.Struct <- NULL
    id <- 1
    rep1 <- 2
    rep2 <- 1
    rep <- 1
    Tree.final <- Tree.construct(origclass, origdata, Tree.Struct, 
                                 id, rep, rep1, rep2, projbest.node, splitCutoff.node, 
                                 PPmethod,lambda=lambda,size.p=size.p, ...)
    Tree.Struct <- Tree.final$Tree.Struct
    colnames(Tree.Struct) <- c("id", "L.node.ID", "R.F.node.ID", 
                               "Coef.ID", "Index")
    projbest.node <- Tree.final$projbest.node
    splitCutoff.node <- Tree.final$splitCutoff.node
    colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
    treeobj <- list(Tree.Struct = Tree.Struct, projbest.node = projbest.node, 
                    splitCutoff.node = splitCutoff.node, origclass = origclass, 
                    origdata = origdata)
    class(treeobj) <- append(class(treeobj), "PPtreeclass")
    return(treeobj)
}
