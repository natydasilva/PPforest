#' Projection Pursuit Random Forest
#'
#'\code{PPforest} implements a random forest using projection pursuit trees algorithm (based on PPtreeViz package).
#' @usage PPforest(data, class, std = TRUE, size.tr, m, PPmethod, size.p,
#'  lambda = .1, parallel = FALSE, cores = 2)
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable.
#' @param std if TRUE standardize the data set, needed to compute global importance measure. 
#' @param size.tr is the size proportion of the training if we want to split the data in training and test.
#' @param m is the number of bootstrap replicates, this corresponds with the number of trees to grow. To ensure that each observation is predicted a few times we have to select this number no too small. \code{m = 500} is by default.
#' @param PPmethod is the projection pursuit index to optimize in each classification tree. The options are \code{LDA} and \code{PDA}, linear discriminant and penalized linear discriminant. By default it is \code{LDA}.
#' @param size.p proportion of variables randomly sampled in each split.
#' @param lambda penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.
#' @param parallel logical condition, if it is TRUE then  parallelize the function
#' @param cores number of cores used in the parallelization
#' @return An object of class \code{PPforest} with components.
#' \item{prediction.training}{predicted values for training data set.}
#' \item{training.error}{error of the training data set.}
#' \item{prediction.test}{predicted values for the test data set if \code{testap = TRUE}(default).}
#' \item{error.test}{error of the test data set if \code{testap = TRUE}(default).}
#' \item{oob.error.forest}{out of bag error in the forest.}
#' \item{oob.error.tree}{out of bag error for each tree in the forest.}
#' \item{boot.samp}{information of bootrap samples.}
#' \item{output.trees}{output from a \code{trees_pp} for each bootrap sample.}
#' \item{proximity}{Proximity matrix, if two cases are classified in the same terminal node then the proximity matrix is increased by one in \code{PPforest} there are one terminal node per class.}
#' \item{votes}{ a matrix with one row for each input data point and one column for each class, giving the fraction of (OOB) votes from the \code{PPforest}.}
#' \item{n.tree}{number of trees grown in \code{PPforest}.}
#' \item{n.var}{number of predictor variables selected to use for spliting at each node.}
#' \item{type}{classification.}
#' \item{confusion}{confusion matrix of the prediction (based on OOB data).}
#' \item{call}{the original call to \code{PPforest}.}
#' \item{train}{is the training data based on \code{size.tr} sample proportion}
#' \item{test}{is the test data based on \code{1-size.tr} sample proportion}
#' @export
#' @examples
#' #crab example with all the observations used as training
#' 
#'pprf.crab <- PPforest(data = crab, class = 'Type',
#'  std = FALSE, size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA' , parallel = TRUE, cores = 2)
#' pprf.crab
#' 
PPforest <- function(data, class, std = TRUE, size.tr = 2/3, m = 500, PPmethod, size.p, lambda = 0.1, 
    parallel = FALSE, cores = 2) {
    
    Var1 <- NULL
    tree <- NULL
    pred <- NULL
    id <- NULL
    if (std) {
        dataux <- data %>% dplyr::select(-get(class)) %>% apply(2, FUN = scale) %>% dplyr::as_data_frame()
        data <- data.frame(data[, class], dataux)
        colnames(data)[1] <- class
    }
    
    clnum <- as.numeric(as.factor(data[, class]))
    tr.index <- trainfn(as.matrix(clnum), as.matrix(data[, setdiff(colnames(data), class)]), 
        sizetr = size.tr) + 1
    train <- data %>% dplyr::slice(tr.index)
    
    type = "Classification"
    var.sel <- round((ncol(train) - 1) * size.p)
    
    outputaux <- baggtree(data = train, class = class, m = m, PPmethod = PPmethod, lambda = lambda, 
        size.p = size.p, parallel = parallel, cores = cores)
    
    output <- lapply(outputaux, function(x) x[[1]])
    
    data.b <- lapply(outputaux, function(x) x[[2]])
    pred.tr <- trees_pred(outputaux, xnew = dplyr::select(train, -get(class)), parallel, cores = cores)
    
    expand.grid.ef <- function(seq1, seq2) {
        data.frame(a = rep.int(seq1, length(seq2)), b = rep.int(seq2, rep.int(length(seq1), 
            length(seq2))))
    }
    pos <- expand.grid.ef(1:dim(train)[1], 1:dim(train)[1])
    
    proximity <- proximi((pred.tr[[1]]), m)
    
    index <- oobindex(data.b, m)
    
    oob.obs <- oobobs(index)
    
    
    mvote.oob <- mvoteoob(pred.tr[[1]], oob.obs)
    oob.pred <- mvote.oob[, length(unique(clnum)) + 1]
    
    votes <- mvote.oob[, -(length(unique(clnum)) + 1)]
    colnames(votes) <- levels(train[, class])
    
    vote.matrix.prop <- votes/rowSums(votes)
    
    oob.error <- 1 - sum(diag(table(oob.pred, unlist(train[, class]))))/length(unlist(train[, 
        class]))
    
    oob.err.tree <- ooberrortree(pred.tr[[1]], oob.obs, as.numeric(as.factor(unlist(train[, 
        class]))), m)
    
    error.tr <- 1 - sum(as.numeric(as.factor(unlist(train[, class]))) == pred.tr[[2]])/length(pred.tr[[2]])
    test <- data[-tr.index, ] %>% dplyr::select(-get(class)) %>% dplyr::filter_()
    
    if (dim(test)[1] != 0) {
        pred.test <- trees_pred(outputaux, xnew = test, parallel, cores = cores)
        error.test <- 1 - sum(as.numeric(as.factor(data[-tr.index, class])) == pred.test[[2]])/length(pred.test[[2]])
        pred.test = as.factor(pred.test[[2]])
        levels(pred.test) <- levels(unlist(train[, class]))
    } else {
        pred.test <- NULL
        error.test <- NULL
        test <- NULL
    }
    
    oob.pred <- as.factor(oob.pred)
    if (is.factor(unlist(train[, class]))) {
        levels(oob.pred) <- levels(unlist(train[, class]))
    } else {
        levels(oob.pred) <- levels(as.factor(unlist(train[, class])))
    }
    
    prediction.training <- as.factor(pred.tr[[2]])
    if (is.factor(unlist(train[, class]))) {
        levels(prediction.training) <- levels(unlist(train[, class]))
    } else {
        levels(prediction.training) <- levels(as.factor(unlist(train[, class])))
    }
    
    
    tab.tr <- table(Observed = unlist(train[, class]), Predicted = oob.pred)
    
    class.error <- 1 - diag(tab.tr)/((stats::addmargins(tab.tr, 2))[, "Sum"])
    confusion <- cbind(tab.tr, class.error = round(class.error, 2))
    
    results <- list(predicting.training = prediction.training, training.error = error.tr, prediction.test = pred.test, 
        error.test = error.test, oob.error.forest = oob.error, oob.error.tree = oob.err.tree, 
        boot.samp = data.b, output.trees = output, proximity = proximity, votes = vote.matrix.prop, 
        prediction.oob = oob.pred, n.tree = m, n.var = var.sel, type = "Classification", confusion = confusion, 
        call = match.call(), train = train, test = test, vote.mat = pred.tr[[1]], class.var = class, 
        oob.obs = oob.obs)
    
    class(results) <- "PPforest"
    
    return(results)
    
}
