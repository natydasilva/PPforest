#' OOB error visualization
#' 
#'\code{ppf_oob_error} Plot the cummulative oob error as a function of number of trees 
#' @usage ppf_oob_error(ppf, nsplit1, nsplit2)
#' @param ppf a PPforest object
#' @param nsplit1 number, increment of the sequence where cumulative oob error rate is computed in the  1/3 trees
#' @param nsplit2 number, increment of the sequence where cumulative oob error rate is computed in the  2/3 trees
#' @return a plot with the cumulative oob error rate
#' @export
#' @examples
#' pprf.leukemia <- PPforest(data = leukemia, class = "Type", 
#' size.tr = 1, m = 70, size.p = .4, PPmethod = 'PDA', strata = TRUE)
#' ppf_oob_error(ppf = pprf.leukemia, nsplit1 = 5, nsplit2 = 15)
ppf_oob_error <- function(ppf, nsplit1, nsplit2) {
    ntree <- NULL
    value <- NULL
    variable <- NULL
    error.cum <- function(ppf, m) {
        l.train <- 1:nrow(ppf$train)
        index2 <- lapply(attributes(ppf$boot.samp)$indices[1:m], function(x) x + 1)
        
        oob.obs <- plyr::ldply(index2, function(x) (!l.train %in% x))
        pred.mtree <- ppf$vote.mat[1:m, ]
        
        oob.pred <- sapply(X = 1:nrow(ppf$train), FUN = function(i) {
            t1 <- table(pred.mtree[oob.obs[, i] == TRUE, i])
            names(t1)[which.max(t1)]
        })
        
        
        oob.mat <- sapply(X = 1:nrow(ppf$train), FUN = function(i) {
            table(pred.mtree[oob.obs[, i] == TRUE, i])
        })
        
        aux <- unlist(lapply(oob.pred, is.null))
        oob.all <- 1 - sum(diag(table(unlist(oob.pred[!aux]), ppf$train[!aux, 1])))/length(ppf$train[!aux, 1])
        tab.err <- table(unlist(oob.pred[!aux]), ppf$train[!aux, 1])
        oob.class <- 1 - diag(tab.err)/apply(tab.err, 2, sum)
        c(oob.all, oob.class)
    }
    
    oob.err.sp <- plyr::mdply(data.frame(m = round(c(seq(2, round(ppf$n.tree*1/3), nsplit1),seq( round(ppf$n.tree*1/3)+1, ppf$n.tree,nsplit2)))), error.cum, ppf = ppf)
    
    names(oob.err.sp)[1:2] <- c("ntree", "all")
    
    oob.pl <- reshape2::melt(oob.err.sp, id.vars = "ntree")
    
    myColors <- c("#000000", RColorBrewer::brewer.pal(length(unique((oob.pl$variable)))-1,"Dark2"))
    names(myColors) <- levels(oob.pl$Class)
    
    if (max(oob.pl$value) < 0.5) {
        p1 <- ggplot2::ggplot(data = oob.pl, ggplot2::aes(x = ntree, y = value, color = variable)) + ggplot2::geom_point() + 
            ggplot2::geom_line(size=I(0.5)) + ggplot2::scale_y_continuous(name = "OOB error rate", limits = c(-0.1, max(oob.pl$value)+.1)) + ggplot2::scale_x_continuous(name = "Number of trees")
       a <- p1 + ggplot2::theme(legend.position = "bottom",aspect.ratio=1) + ggplot2::scale_color_manual(values = myColors, name =
                                                                          "Class") + ggplot2::ggtitle("Cumulative OOB error")
      
         plotly::ggplotly(a)
     
        
    } else {
        p1 <- ggplot2::ggplot(data = oob.pl, ggplot2::aes(x = ntree, y = value, color = variable)) + ggplot2::geom_point() + 
            ggplot2::geom_line(size=I(.5)) + ggplot2::scale_y_continuous(name = "OOB error rate", limits = c(-0.1, 1.1)) + ggplot2::scale_x_continuous(name = "Number of trees")
        a <- p1 +ggplot2::theme(legend.position = "bottom",aspect.ratio=1) + ggplot2::scale_color_manual(values = myColors, name =
                                                                                                           "Class") + ggplot2::ggtitle("Cumulative OOB error")
        
         
         plotly::ggplotly(a)
      
    }
} 
