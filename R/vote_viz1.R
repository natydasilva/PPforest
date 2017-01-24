#' Vote matrix visualization
#'  
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param ppf is a PPforest object
#' @return A side-by-side jittered dotplot with the probability for each class to be classify in the corresponding class column.
#' Color by true observed class for each observation.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' pprf.crab <- PPforest2(data = crab, class = "Type",
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA', strata = TRUE)
#' vote_viz1(crab,"Type",pprf.crab)
vote_viz1 <- function(data, class, ppf){
  Class <- NULL
  Probability <- NULL
  Type <- NULL
  reddata <- data.frame(
      ids = 1:nrow(data), Type = ppf$train[,1], ppf$votes, pred = ppf$prediction.oob
    ) %>%  tidyr::gather(Class,Probability,-dplyr::one_of("pred","ids",class))
  
  
  
  p <-
    ggplot2::ggplot(data = reddata, ggplot2::aes(Class, Probability, color = Type) )+ ggplot2::geom_jitter(height =
                                                                                          0,size = I(3), alpha = .5) +
    ggplot2::theme(legend.position = "none",axis.text.x  = ggplot2::element_text(angle = 45, vjust = 0.5), aspect.ratio = 1) +
    ggplot2::labs(x = "Class", title = "Side by side plot") + ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2")
  
  
  plotly::ggplotly(p,tooltip = c("x","y")) 
  
  
}
