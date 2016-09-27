#' Parallel plot for the data set colored by class
#'
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @return A parallel coordinate plot with the data selt collored by class
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #leukemia data set with all the observations used as training
#' pprf.leukemia <- PPforest(data = leukemia, class = "Type",
#'  size.tr = 1, m = 70, size.p = .4, PPmethod = 'PDA', strata = TRUE)
#' parallel_dat(leukemia,"Type")
parallel_dat <- function(data, class){

var <-NULL
Value <- NULL
Variables <- NULL
ids <- NULL
Class <- NULL
  aux <- data %>% dplyr::arrange_(class) %>% dplyr::mutate(ids = 1:nrow(data))  %>% 
    tidyr::gather(var, Value, -dplyr::one_of(class, 'ids')  )
  aux$Variables <- as.factor(aux$var)
  colnames(aux)[which(colnames(data)==class)] <- "Class"
  p <-
    ggplot2::ggplot(aux, ggplot2::aes(
      x = Variables , y = Value ,group = ids, colour = Class
    )) +
    ggplot2::geom_line(alpha = 0.3) + ggplot2::scale_x_discrete(limits = levels(as.factor(aux$var))) + ggplot2::ggtitle("Parallel plot CRAB data") +
    ggplot2::theme(legend.position = "none",axis.text.x = ggplot2::element_text(angle = 90)) + ggplot2::scale_colour_brewer(type = "qual",palette =
                                                            "Dark2")
  

  
  plotly::ggplotly(p,tooltip = c("colour","y","key"))
}