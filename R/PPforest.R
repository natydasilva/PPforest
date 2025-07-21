#' Projection Pursuit Random Forest
#'
#'\code{PPforest} implements a random forest using projection pursuit trees algorithm (based on PPtreeViz package).
#' @usage PPforest(data, class, xstd = 'min-max', size.tr, m, PPmethod, size.p,
#'  lambda = .1, parallel = FALSE, cores = 2, rule = 1)
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable.
#' @param xstd if TRUE standardize the data set, needed to compute global importance measure. 
#' @param size.tr is the size proportion of the training if we want to split the data in training and test.
#' @param m is the number of bootstrap replicates, this corresponds with the number of trees to grow. To ensure that each observation is predicted a few times we have to select this number no too small. \code{m = 500} is by default.
#' @param PPmethod is the projection pursuit index to optimize in each classification tree. The options are \code{LDA} and \code{PDA}, linear discriminant and penalized linear discriminant. By default it is \code{LDA}.
#' @param size.p proportion of variables randomly sampled in each split.
#' @param lambda penalty parameter in PDA index and is between 0 to 1 . If \code{lambda = 0}, no penalty parameter is added and the PDA index is the same as LDA index. If \code{lambda = 1} all variables are treated as uncorrelated. The default value is \code{lambda = 0.1}.
#' @param parallel logical condition, if it is TRUE then  parallelize the function
#' @param cores number of cores used in the parallelization
#' @param rule split rule 1: mean of two group means 2: weighted mean of two group means - weight with group size 3: weighted mean of two group means - weight with group sd 4: weighted mean of two group means - weight with group se 5: mean of two group medians 6: weighted mean of two group medians - weight with group size 7: weighted mean of two group median - weight with group IQR 8: weighted mean of two group median - weight with group IQR and size
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
#'@references da Silva, N., Cook, D., & Lee, E. K. (2021). A projection pursuit forest 
#'algorithm for supervised classification. Journal of Computational and Graphical Statistics,
#' 30(4), 1168-1180.
#' @export
#' @examples
#' #crab example with all the observations used as training
#' set.seed(123)
#'pprf.crab <- PPforest(data = crab, class = 'Type',
#'  xstd = 'no', size.tr = 0.7, m = 200, size.p = .4, 
#'  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
#' pprf.crab
#' 
PPforest <- function(data, class, xstd = 'min-max', size.tr = 2/3, m = 500, PPmethod, size.p, lambda = 0.1, 
    parallel = FALSE, cores = 2, rule = 1) {
    
    Var1 <- NULL
    tree <- NULL
    pred <- NULL
    id <- NULL

    cllev <- levels(as.factor(data[, class]))
    clnum <- as.numeric(as.factor(data[, class]))
    tr.index <- trainfn(as.matrix(clnum), as.matrix(data[, setdiff(colnames(data), class)]), 
                        sizetr = size.tr) + 1
    tr.index <- as.vector(tr.index)
    train <- data %>% dplyr::slice(tr.index)
    test <- data[-tr.index, ] 
    #%>% dplyr::select(-(!!class)) 
    
    type = "Classification"
    
    # Variable scaling.
    if (xstd != "no") {
      
      if (xstd == "min-max") {
        dataux <- train %>% 
          dplyr::select(-(!!class)) |>  tibble::as_tibble()
        mincol <- dataux |> apply( 2, min)
        maxmincol <- dataux |> apply(2, function(x) {
          max(x) - min(x)
        })
      }
      if (xstd == "quant") {
        dataux <- train %>% 
          dplyr::select(-(!!class)) |>  tibble::as_tibble()
        mincol <- dataux |> apply( 2, stats::quantile, 0.05) 
        maxmincol <- dataux |> apply(2, function(x) {
          stats::quantile(x, 0.95) - stats::quantile(x, 0.05)
        }) 
      }
    
      if (xstd == 'scale') {
        # dataux1 <- train %>% 
        #   dplyr::select(-(!!class)) %>% 
        #   apply(2, FUN = scale) %>%
        #   tibble::as_tibble()
        meanx <- train |>  dplyr::select(-(!!class)) |>  dplyr::summarise(dplyr::across(dplyr::everything(), mean)) |> unlist()
        sdx   <- train |>  dplyr::select(-(!!class)) |>  dplyr::summarise(dplyr::across(dplyr::everything(), stats::sd))   |> unlist()
        dataux <- train |> 
          dplyr::select(-(!!class)) |>
          dplyr::mutate( dplyr::across( dplyr::everything(), 
                        ~ (.x - meanx[cur_column()]) / sdx[cur_column()]))
    }
    if(xstd %in% c('min-max', 'quant')) {
    trainscale <- (dataux  - matrix(mincol, nrow(dataux), ncol(dataux), byrow = T)) / matrix(maxmincol, nrow(dataux), ncol(dataux), byrow = T)
    train <- data.frame(train[, class], trainscale)
    colnames(train)[1] <- class
    if (dim(test)[1] != 0){
      testscale <- (test |> dplyr::select(-(!!class))  - matrix(mincol, nrow(test), ncol(test), byrow = T)) / matrix(maxmincol, nrow(test), ncol(test), byrow = T)
      test <-  data.frame(test[, class], testscale)
      colnames(test)[1] <- class
      }
    
    }else{
      train <- data.frame(data[, class], dataux)
      colnames(train)[1] <- class
      if (dim(test)[1] != 0){
        testscale <- test |> 
          dplyr::select(-(!!class)) |>
          dplyr::mutate(dplyr::across(dplyr::everything(), 
                        ~ (.x - meanx[cur_column()]) / sdx[cur_column()]))
        test <- data.frame(test[, class], testscale)
        colnames(test)[1] <- class
        
      }
      }
    }  
    
    outputaux <- baggtree(data = train, class = class, m = m, PPmethod = PPmethod, lambda = lambda, 
        size.p = size.p, parallel = parallel, cores = cores)
    
    var.sel <- sum(outputaux[[1]][[1]]$projbest.node[1,]!=0)
    output <- lapply(outputaux, function(x) x[[1]])
    
    data.b <- lapply(outputaux, function(x) x[[2]])
    pred.tr <- trees_pred(outputaux, xnew = dplyr::select(train, -(!!class)), parallel, cores = cores, rule = rule)
    
    expand.grid.ef <- function(seq1, seq2) {
        data.frame(a = rep.int(seq1, length(seq2)), b = rep.int(seq2, rep.int(length(seq1), 
            length(seq2))))
    }
    pos <- expand.grid.ef(1:dim(train)[1], 1:dim(train)[1])
    
    proximity <- proximi((pred.tr$predtree), m)
    
    index <- oobindex(data.b, m)
    
    oob.obs <- oobobs(index)
    
    mvote.oob <- mvoteoob(pred.tr$predtree, oob.obs)
    oob.pred <- mvote.oob[, length(unique(clnum)) + 1]
    
    votes <- mvote.oob[, -(length(unique(clnum)) + 1)]
    colnames(votes) <- levels(train[, class])
    
    vote.matrix.prop <- votes/rowSums(votes)
    
    oob.error <- 1 - sum(diag(table(oob.pred, unlist(train[, class]))))/length(unlist(train[, 
        class]))
    
    oob.err.tree <- ooberrortree(pred.tr$predtree, oob.obs, as.numeric(as.factor(unlist(train[, 
        class]))), m)
    
  
    error.tr <- 1 - sum(as.numeric(as.factor(unlist(train[, class]))) == pred.tr$predforest)/length(pred.tr$predforest)
   
    
    if (dim(test)[1] != 0){
      
        pred.test <- trees_pred(outputaux, xnew = dplyr::select(test, -(!!class)), parallel, cores = cores, rule = rule)
        error.test <- 1 - sum(as.numeric(as.factor(data[-tr.index, class])) == pred.test[[2]])/length(pred.test[[2]])
        pred.test = as.factor(pred.test[[2]])
        levels(pred.test) <- levels(as.factor(train[, class]))
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
    
    prediction.training <- as.factor(pred.tr$predforest)
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
        call = match.call(), train = train, test = test, vote.mat = pred.tr$predtree, vote.mat_cl= cllev, class.var = class,
        oob.obs = oob.obs, xstd = xstd)
  
    
    class(results) <- "PPforest"
    
    return(results)
    
}
