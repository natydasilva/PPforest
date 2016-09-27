#' Draws bootstrap samples with strata option.
#'
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param m is the number of bootstrap replicates we want to sample.
#' @param strata is TRUE if the bootstrap samples are stratified by class variable.
#' @return data frame object with m bootstrap samples
#' @export
#' @importFrom magrittr %>%
#' @examples
#'leukemia.b <- ppf_bootstrap(data = leukemia, class = "Type",  m = 200) 
#'index <- lapply(attributes(leukemia.b)$indices, function(x) x + 1)
ppf_bootstrap <- function(data, class, m = 500, strata = TRUE) {
    . <- NULL
    samp <- NULL
    n <- nrow(data)
    class.id <- data %>%
      dplyr::select_(class) %>%
      dplyr::mutate(id = 1:n) 
    
   # class.id <- dplyr::data_frame_(id = 1:n, group = class)
    if (strata == TRUE) {
        samp.g <- replicate(m, class.id %>%
                              dplyr::group_by_(class) %>% 
                              dplyr::do(samp = sort(sample(.$id, replace = TRUE))) %>% 
                              dplyr::ungroup() %>% 
                              dplyr::select(samp), simplify = FALSE)
        
        attr(data, "indices") <- lapply(samp.g, function(x) as.numeric(sort(unlist(x))) - 1)
    } else {
        attr(data, "indices") <- replicate(m, sort(sample(n, replace = TRUE) - 1), simplify = FALSE)
    }
    
    attr(data, "drop") <- TRUE
    attr(data, "group_sizes") <- rep(n, m)
    attr(data, "biggest_group_size") <- n
    attr(data, "labels") <- data.frame(replicate = 1:m)
    attr(data, "vars") <- list(quote(replicate))
    class(data) <- c("grouped_df", "tbl_df", "tbl", "data.frame")
    
    data
} 
