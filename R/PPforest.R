#' Projection Pursuit Random Forest
#'
#'\code{PPforest} implements a random forest using projection pursuit trees algorithm (based on PPtreeViz package).
#' @usage PPforest(data, class, size.tr, m, PPmethod, size.p, strata = TRUE, lambda=.1)
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param size.tr is the size proportion of the training if we want to split the data in training and test.
#' @param m is the number of bootstrap replicates, this corresponds with the number of trees to grow. To ensure that each observation is predicted a few times we have to select this nunber no too small. \code{m = 500} is by default.
#' @param PPmethod is the projection pursuit index to optimize in each classification tree. The options are \code{LDA} and \code{PDA}, linear discriminant and penalized linear discriminant. By default it is \code{LDA}.
#' @param size.p proportion of variables randomly sampled in each split.
#' @param strata if set \code{TRUE} then the bootrap samples are stratifyed by class variable.
#' @param lambda penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.
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
#' \item{confusion}{confusion matrix of the prediction (based on OOb data).}
#' \item{call}{the original call to \code{PPforest}.}
#' \item{train}{is the training data based on \code{size.tr} sample proportion}
#' \item{test}{is the test data based on \code{1-size.tr} sample proportion}
#' @export
#' @examples
#' #leukemia data set with all the observations used as training
#' pprf.leukemia <- PPforest(data = leukemia, class = "Type",
#'  size.tr = 1, m = 70, size.p = .4, PPmethod = 'PDA', strata = TRUE)
#' pprf.leukemia
PPforest <- function(data, class,  size.tr = 2/3, m = 500, PPmethod, size.p, strata = TRUE, lambda = 0.1) {
   
    Var1 <- NULL
    tr.index <- train_fn(data, class,  size.tr)
    train <-  data %>% 
              dplyr::slice(as.numeric(tr.index$id))
    
    
    type = "Classification"
    var.sel <- floor((ncol(train)-1) * size.p)
    
    if (strata == TRUE) {
        data.b <- ppf_bootstrap(data = train, class, m, strata)  
        output <- data.b %>% trees_pp(size.p, PPmethod, lambda = 0.1)
    } else {
        data.b <- ppf_bootstrap(data = train, class, m, strata = FALSE)
        output <- data.b %>% trees_pp( size.p, PPmethod, lambda = 0.1)
    }
    
 
    pred.tr <- tree_ppred(xnew = dplyr::select(train,-get(class)), output)
 
    pos <- expand.grid(a = 1:dim(train)[1], b = 1:dim(train)[1])
    tri.low <- pos %>% dplyr::filter(pos[, 1] >= pos[, 2])
    
    same.node <- data.frame(tri.low, dif = apply(t(pred.tr[[2]]), 2, function(x) x[tri.low[, 1]] == x[tri.low[, 
        2]]))
    proximity <- data.frame(same.node[, c(1:2)], proxi = apply(same.node[, -c(1:2)], 1, function(x) sum(x == 1))/dim((pred.tr[[2]]))[1])
    
    l.train <- 1:nrow(train)
    index <- lapply(attributes(data.b)$indices, function(x) x + 1)
    
    oob.obs <- plyr::ldply(index, function(x) (!l.train %in% x))
    
    oob.pred <- sapply(X = 1:nrow(train), FUN = function(i) {
      if(sum(oob.obs[, i]>0)){
        t1 <- table(pred.tr[[2]][oob.obs[, i] == TRUE, i])
        names(t1)[which.max(t1)]
      }else{
        print("More trees are needed to get the oob predictions")
      }
    })
    

    votes <- plyr::mdply( dplyr::data_frame( id = 1:nrow(train)) , function(id) {
      x <- pred.tr[[2]][oob.obs[, id] == TRUE, id] 
      table(factor(x, levels=levels(train[, class])) )  
      })[,-1]

#     votes <- matrix(0, ncol = length(unique(train[, class])), nrow = nrow(train))
#     colnames(votes) <- levels(train[, class])
#     
    
#     for (i in 1:nrow(train)) {
#         cond <- colnames(votes) %in% names(oob.mat[[i]])
#         votes[i, cond] <- oob.mat[[i]]
#     }
#     

    
    vote.matrix.prop <-votes/rowSums(votes)

    oob.error <- 1 - sum(diag(table(oob.pred, train[, class])))/length(train[, class])
    
    
    m.pred.tr <- reshape2::melt(pred.tr[[2]])
    m.oob.obs <- reshape2::melt(as.matrix(oob.obs))
    m.pred.tr$oob <- m.oob.obs$value
    m.pred.tr$class <- rep(train[, class], each = m)
    
    
    oob.err.tree <- plyr::ddply(m.pred.tr[m.pred.tr$oob, ], plyr::.(Var1), function(x) {
        dd <- diag(table(x$value, x$class))
        1 - sum(dd)/length(x$value)
    })$V1
    
    
    
    error.tr <- 1 - sum(train[, class] == pred.tr[[3]])/length(pred.tr[[3]])
    test <- data[-tr.index$id,]%>%  
            dplyr::select(-get(class))%>%
            dplyr::filter_()
    
    if (dim(test)[1] != 0) {
        pred.test <- tree_ppred(xnew = test, output)
        error.test <- 1 - sum(data[-tr.index$id, class] == pred.test[[3]])/length(pred.test[[3]])
    } else {
        pred.test <- NULL
        error.test <- NULL
        test <- NULL
    }
    
    oob.pred <- as.factor(oob.pred)
    levels(oob.pred) <- levels(train[, class])
    
    tab.tr <- table(Observed = train[, class], Predicted = oob.pred)
    
    class.error <- 1 - diag(tab.tr)/((stats::addmargins(tab.tr, 2))[, "Sum"])
    confusion <- cbind(tab.tr, class.error = round(class.error, 2))
    
    results <- list(prediction.training = pred.tr[[3]], training.error = error.tr, prediction.test = pred.test[[3]], 
        error.test = error.test, oob.error.forest = oob.error, oob.error.tree = oob.err.tree, boot.samp = data.b, 
        output.trees = output, proximity = proximity, votes = vote.matrix.prop, prediction.oob = oob.pred, n.tree = m, 
        n.var = var.sel, type = "Classification", confusion = confusion, call = match.call(), train = train, test = test, 
        vote.mat = pred.tr[[2]],class.var=class, oob.obs = oob.obs)
    class(results) <- "PPforest"
    
    return(results)
    
} 
