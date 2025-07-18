## ----libraries, cache = FALSE, echo = FALSE, message = FALSE, warning = FALSE----
require(PPforest)
require(dplyr)
require(RColorBrewer)
require(GGally)
require(gridExtra)
require(PPtreeViz)
require(ggplot2)
require(knitr)
set.seed(310756) #reproducibility

## ----hooks, echo = FALSE------------------------------------------------------

knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE, autodep=TRUE, cache.lazy=FALSE )
opts_knit$set(eval.after = 'fig.cap')
theme_set(theme_bw(base_family="serif"))

## ----descri, fig.align="center", fig.cap=capmatrix,  fig.show='hold', fig.height = 5, fig.width = 5, echo=FALSE, eval=FALSE----
# 
# a <- GGally::ggpairs(PPforest::crab,
#     columns = 2:6,
#     ggplot2::aes(colour = Type, alpha=.1),
#     lower = list(continuous = 'points'),
#     axisLabels='none',
#     upper=list(continuous='blank')
#      , legend = NULL)
# 
# capmatrix<-"Scatter plot matrix of crab data "
# a

