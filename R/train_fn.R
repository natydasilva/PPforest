#' Obtain stratified sample of the data 
#' 
#' Index id for training set, sample in each class with the same proportion.
#' @param data Data frame with the complete data set.
#' @param class A character with the name of the class variable. 
#' @param size.p Sample size proportion in each class group.
#' @return numeric vector  with the sample indexes for training data
#' @export
#' @importFrom magrittr %>%
#' @examples
#' training.id <- train_fn(data = leukemia, class = "Type", size.p = 0.7)
#' training.id
train_fn <- function(data, class, size.p = 0.7) {
    id <- NULL
    n <- nrow(data)
    class.id <- data %>%
      dplyr::select_(class) %>%
      dplyr::mutate(id = 1:n) 
     
    class.id %>% 
      dplyr::group_by_(class) %>% 
      dplyr::sample_frac(size.p) %>% 
      dplyr::arrange(id) %>%
      dplyr::ungroup() %>% 
      dplyr::select(id)
} 
