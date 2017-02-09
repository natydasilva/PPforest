#' Gobal importance measure for a PPforest object as the average IMP PPtree measure over all the trees 
#' in the forest
#' 
#' @param ppf is a PPforest object
#' @return Data frame with the global importance measure
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#' pprf.crab <- PPforest(data = crab, std =TRUE, class = "Type",
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA')
#'  #ppf_avg_pptree_imp(pprf.crab) 
ppf_avg_pptree_imp <- function(ppf) {
node.id <- NULL
nodecl <- NULL
data<-NULL
node<- NULL
Class<- NULL
variable <-NULL
value<-NULL
tr<-NULL

  nn <- data.frame(nn = 1:sum(ppf[["output.trees"]][[1]]$Tree.Struct[, 4]!=0))
  nodecl <- function(x) {
   aux <- node_data(ppf = ppf, x ) 
   aux$node.id <- as.factor( aux$node.id)
     aux %>% dplyr::group_by(node.id) %>% 
      dplyr::summarise(clt=length(unique(Class)))  
  }
  
  #dat_pl <- apply(nn, 1,  nodecl) %>% lapply(data.frame) %>% dplyr::bind_rows()

  mat.proj <- lapply(ppf[["output.trees"]], function(x){
    data.frame( node = 1:sum(x$Tree.Struct[, 4]!=0), abs(x[[2]]))
  }) %>% dplyr::bind_rows() 
  
  info <- apply(data.frame(1:ppf$n.tree),1, function(x) nodecl(x))
  # lapply(info, function(x) {
  #   if(node.id%in%nn){
  #   node.id)!=nn
  # }
  #mat.proj%>% group_by()lapply(info, function(x), merge())
  colnames(mat.proj)[-1] <- colnames(dplyr::select(data,-get(class)))
  
  mat.proj %>% dplyr::mutate(tr = rep(1:nrow(ppf$train), each = length(nn[,1]))) %>% 
    tidyr::gather(variable, value, -node, -tr ) %>% 
    dplyr::group_by(variable, tr) %>% 
    dplyr::summarise(mean = sum(value)) %>%
    dplyr::group_by(variable) %>% 
    dplyr::summarise(mean = mean(mean)) %>%
    dplyr::arrange(dplyr::desc(mean) ) 
  #weight for the number of classes'
  }
