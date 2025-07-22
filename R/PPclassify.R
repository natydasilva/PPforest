
#' Predict class for the test set and calculate prediction error after finding the PPtree structure, .
#' 
#' @usage PPclassify( Tree.result, test.data = NULL, Rule = 1, true.class = NULL)  
#' @param Tree.result the result from PPtree_split
#' @param test.data  the test dataset
#' @param Rule split rule 1:mean of two group means, 2:weighted mean, 3: mean of max(left group) and min(right group), 4: weighted mean of max(left group) and min(right group)
#' @param true.class true class of test dataset if available
#' @return predict.class predicted class
#' @return predict.error prediction error
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection pursuit classification tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @keywords tree
PPclassify<- function(Tree.result, test.data = NULL, Rule = 1, true.class = NULL) {
   cllev <- levels(as.factor(Tree.result[[4]]) )
   
   if (is.null(test.data)){
    test.data <- Tree.result$origdata
   }
   test.data <- as.matrix(test.data)
   
    if (!is.null(true.class)) {
        true.class <- as.matrix(true.class)
        if (nrow(true.class) == 1) 
            true.class <- t(true.class)
        if (!is.numeric(true.class)) {
            class.name <- names(table(true.class))
            temp <- rep(0, nrow(true.class))
            for (i in 1:length(class.name)) temp <- temp + (true.class == class.name[i]) * i
            true.class <- temp
        }
    }
    
    
    n <- nrow(test.data)
    class.temp <- rep(1, n)
    test.class.index <- matrix(rep(0, n), ncol = n)
    temp <- PPclassindex(class.temp, test.class.index, as.matrix(test.data), as.matrix(Tree.result$Tree.Struct), 
        as.matrix(Tree.result$projbest.node), as.matrix(Tree.result$splitCutoff.node), 0, Rule)
    
    
    test.class <- rep(0, n)
    IOindex <- rep(1, n)
    if (dim(as.matrix(temp$testclassindex[-1, ]))[2] == 1) {
        temp <- PPclassification(as.matrix(Tree.result$Tree.Struct), t(as.matrix(temp$testclassindex[-1, 
            ])), as.vector(IOindex), as.vector(test.class), 0, 0)
    } else {
        temp <- PPclassification(as.matrix(Tree.result$Tree.Struct), as.matrix(temp$testclassindex[-1, 
            ]), as.vector(IOindex), as.vector(test.class), 0, 0)
    }
   
    
    if (!is.null(true.class)) {
        predict.error <- sum(true.class != temp$testclass)/dim(true.class)[1]
    } else {
        predict.error <- NA
    }
    
    temp$testclass <-as.factor(temp$testclass)
    levels(temp$testclass) <- cllev
    
    list(predict.error = predict.error, predict.class = temp$testclass)
    
}
