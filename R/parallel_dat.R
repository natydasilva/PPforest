#' Parallel plot for the data set colored by class, the data are standarized in the function
#'
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @return A parallel coordinate plot with the data selt collored by class
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' pprf.crab <- PPforest2(data =crab, class = "Type",
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA', strata = TRUE)
#' parallel_dat(crab,"Type")
parallel_dat <- function(data, class){

var <-NULL
Value <- NULL
Variables <- NULL
ids <- NULL
Class <- NULL
sd <- NULL
Type <-NULL


myscale <- function(x) (x - mean(x)) / sd(x)
scale.dat <- data %>% dplyr::mutate_each(dplyr::funs(myscale),-dplyr::matches(class)) 
scale.dat.melt <- scale.dat %>%  dplyr::mutate(ids = 1:nrow(data)) %>% tidyr::gather(var,Value,-Type,-ids)
scale.dat.melt$Variables <- as.numeric(as.factor(scale.dat.melt$var))
colnames(scale.dat.melt)[1] <- "Class"

  p <- ggplot2::ggplot(scale.dat.melt, ggplot2::aes(x = Variables, y = Value, 
                                  group = ids, key = ids, colour = Class, var = var)) +
    ggplot2::geom_line(alpha = 0.3) + ggplot2::scale_x_discrete(limits = levels(as.factor(scale.dat.melt$var)), expand = c(0.01,0.01)) +
    ggplot2::ggtitle("Data parallel plot ") + ggplot2::theme(legend.position = "none", axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5), aspect.ratio = 1) + 
    ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2")
  
  

  plotly::ggplotly(p,tooltip = c("var","colour","y","key"))
}

