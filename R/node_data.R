#' Data structure with the  projected and boundary by node and class.
#' @param ppf is a PPforest object
#' @param tr numerical value to  identify a tree
#' @param Rule split rule 1:mean of two group means, 2:weighted mean, 3: mean of max(left group) and min(right group), 4: weighted mean of max(left group) and min(right group)
#' @return Data frame with projected data for each class and node id and the boundaries
#' @export
#' @importFrom magrittr %>%
#' @examples
#' #crab data set with all the observations used as training
#'
#' pprf.crab <- PPforest(data = crab, std = 'min-max', y = 'Type',
#'  size.tr = 1, m = 200, size.p = .5, PPmethod = 'LDA')
#' node_data(ppf = pprf.crab, tr = 1) 
#' 
node_data <- function(ppf, tr, Rule = 1) {
    ind_node <- function(PPclassOBJ, node.id, Rule) {
        searchGroup <- function(node.id, TS, gName) {
            flag <- TRUE
            sel.id <- TS[node.id, 2:3]
            LR.id <- c(TRUE, FALSE)
            sel.group <- NULL
            i <- 1
            while ((sel.id[i] != 0) && (i < length(sel.id))) {
                if (TS[sel.id[i], 2] != 0) {
                  sel.id <- c(sel.id, TS[sel.id[i], 2:3])
                  if (LR.id[i]) 
                    LR.id <- c(LR.id, c(TRUE, TRUE)) else LR.id <- c(LR.id, c(FALSE, FALSE))
                }
                if (TS[sel.id[i + 1], 2] != 0) {
                  sel.id <- c(sel.id, TS[sel.id[i + 1], 2:3])
                  if (LR.id[i + 1]) 
                    LR.id <- c(LR.id, c(TRUE, TRUE)) else LR.id <- c(LR.id, c(FALSE, FALSE))
                }
                i <- i + 2
            }
            sel.Name <- TS[sel.id[which(TS[sel.id, 2] == 0)], 3]
            selName <- sort(gName[sel.Name])
            L.list <- sort(gName[sel.Name[LR.id[which(TS[sel.id, 2] == 0)]]])
            R.list <- sort(gName[sel.Name[!LR.id[which(TS[sel.id, 2] == 0)]]])
            return(list(selName = selName, Llist = L.list, Rlist = R.list))
        }
        
        TS <- PPclassOBJ$Tree.Struct
        Alpha <- PPclassOBJ$projbest.node
        cut.off <- PPclassOBJ$splitCutoff.node
        origdata <- PPclassOBJ$origdata
        origclass <- PPclassOBJ$origclass
        p <- ncol(origdata)
        gName <- names(table(origclass))
        
        if (TS[node.id, 2] != 0) {
            SG.result <- searchGroup(node.id, TS, gName)
            selG <- SG.result$selName
            selL <- SG.result$Llist
            selR <- SG.result$Rlist
            sel.id <- NULL
            LR.class <- NULL
            for (i in 1:length(selG)) {
                sel.id <- c(sel.id, which(origclass == selG[i]))
                LR.class <- c(LR.class, rep(ifelse(sum(selL == selG[i]) != 0, "L", "R"), length(which(origclass == 
                  selG[i]))))
            }
            proj.data <- c(as.matrix(origdata) %*% as.matrix(Alpha[TS[node.id, 4], ]))[sel.id]
            
            
            proj.class <- origclass[sel.id]
            plot.data <- data.frame(proj.data = proj.data, origclass = proj.class, cut = cut.off[TS[node.id, 
                4], Rule], node.id = node.id, LR.class, Dir = as.factor(proj.data > cut.off[TS[node.id, 
                4], Rule]))
            colnames(plot.data)[2] <- "Class"
            
            plot.data
            
        }
    }
    
   
    nn <- data.frame(nn = ppf[["output.trees"]][[tr]]$Tree.Struct[ppf[["output.trees"]][[tr]]$Tree.Struct[, 
        4] != 0, 1])
    densf <- function(x) {
        ind_node(PPclassOBJ = ppf[["output.trees"]][[tr]], node.id = x, Rule = 1)
    }
    
    dat_pl <- lapply(nn[, 1], densf) %>% lapply(data.frame) %>% dplyr::bind_rows()
    
    dat_pl
    
}
