#' OOB error visualization
#' 
#'\code{ppf_oob_error} Plot the cummulative oob error as a function of number of trees 
#' @usage ppf_oob_error(ppf, nsplit1, interactive)
#' @param ppf a PPforest object
#' @param nsplit1 number, increment of the sequence where cumulative oob error rate is computed in the  1/3 trees
#' @param interactive if is TRUE will use plotly to translate an static ggplot object
#' @return a plot with the cumulative oob error rate
#' @export
#' @examples
#' pprf.crab<- PPforest(data = crab, class = "Type", 
#' size.tr = 1, m = 200, size.p = .4, PPmethod = 'LDA')
#' ppf_oob_error(ppf = pprf.crab, nsplit1 = 15, interactive = TRUE)
ppf_oob_error <- function(ppf, nsplit1,interactive){
  ntree <- NULL
  value <- NULL
  variable <- NULL
  tree.id <- NULL
  OOB.error <- NULL
  Class <- NULL
   ppf_oob_aux<-function(ppf, nsplit1){

  error.cum <- function(ppf, m) {
    l.train <- 1:nrow(ppf$train)
    index2 <- lapply(as.numeric(attributes(ppf$boot.samp)$names[1:m]), function(x)
        x + 1)
    
    
    oob.obs <- index2 %>%  lapply(function(x)
      data.frame(obs=!l.train %in% x)) %>% dplyr::bind_cols() %>%t()
    pred.mtree <- ppf$vote.mat[1:m,]
    
    
    
    oob.pred <-
      sapply(
        X = 1:nrow(ppf$train), FUN = function(i) {
          t1 <- table(pred.mtree[oob.obs[, i] == TRUE, i])
          names(t1)[which.max(t1)]
        }
      )
    
    
    oob.mat <- sapply(
      X = 1:nrow(ppf$train), FUN = function(i) {
        table(pred.mtree[oob.obs[, i] == TRUE, i])
      }
    )
    
    aux <- unlist(lapply(oob.pred, is.null))
    oob.all <-
      1 - sum(diag(table(unlist(oob.pred[!aux]), ppf$train[!aux, 1]))) / length(ppf$train[!aux, 1])
    tab.err <- table(unlist(oob.pred[!aux]), ppf$train[!aux, 1])
    oob.class <- 1 - diag(tab.err) / apply(tab.err, 2, sum)
    c(oob.all, oob.class)
  }
  
  
  mm <- data.frame(m = round(seq(
    2, round(ppf$n.tree),nsplit1)))
  
  errcfun <- function(x){
    error.cum(ppf,x)
  }
  
  oob.err.sp <- data.frame(mm, apply(mm, 1,errcfun)  %>% t() )
  
  
  names(oob.err.sp)[1:2] <- c("tree.id", "All")
  
  #oob.pl <- reshape2::melt(oob.err.sp, id.vars = "tree.id")
  oob.pl <- oob.err.sp %>% tidyr::gather(variable, value, -tree.id)
  colnames(oob.pl)[2:3] <- c("Class", "OOB.error")
  oob.pl
}

oob.pl <- ppf_oob_aux(ppf, nsplit1 = round(nrow(ppf$train)/13) )



myColors <- c("#000000", RColorBrewer::brewer.pal(length(unique(ppf$train[, ppf$class.var])),"Dark2"))
names(myColors) <- levels(oob.pl$Class)


p1 <- oob.pl %>% ggplot2::ggplot(ggplot2::aes( x = tree.id, y = OOB.error , colour = Class)) + 
  ggplot2::geom_point(alpha = .5) + ggplot2::geom_line(size = I(0.5), alpha = .5) + ggplot2::labs(y = "OOB error rate", 
    x = "Number of trees", title = "Cumulative OOB error") + ggplot2::ylim(c(0,1)) +
  ggplot2::theme(legend.position = "none", aspect.ratio = 1) + ggplot2::scale_color_manual(values = myColors)

if(interactive){
  plotly::ggplotly(p1,tooltip = c("colour","y","x") )
  
}else{
  
  print(p1)
  
}



}
