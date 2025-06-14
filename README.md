
`PPforest` package
======================
Natalia da Silva, Dianne Cook & Eun-Kyung Lee 



<img src="man/figures/PPforest.png" align="right" alt="" width="160" />
\vspace{1cm}

[![CRAN Status](https://www.r-pkg.org/badges/version/PPforest)]( https://CRAN.R-project.org/package=PPforest)[![CRAN\_Download\_Badge](https://cranlogs.r-pkg.org/badges/grand-total/PPforest)](https://cran.r-project.org/package=PPforest) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/PPforest)](https://www.r-pkg.org/pkg/PPforest)[![Travis-CI Build Status](https://travis-ci.org/natydasilva/PPforest.svg?branch=master)](https://travis-ci.org/natydasilva/PPforest)


Introduction
============


The `PPforest` package (projection pursuit random forest) contains functions to run a projection pursuit random forest for classification problems. This method utilize combinations of variables in each tree construction.  In a random forest each split is based on a single variable, chosen from a subset of predictors. In the `PPforest`, each split is based on a linear combination of randomly chosen variables. The linear combination is computed by optimizing a projection pursuit index, to get a projection of the variables that best separates the classes. The `PPforest` uses the `PPtree` algorithm, which fits a single tree to the data. Utilizing linear combinations of variables to separate classes takes the correlation between variables into account, and can outperform the basic forest when separations between groups occurs on combinations of variables. Two projection pursuit indexes, LDA and PDA, are used for `PPforest`.

To improve the speed performance `PPforest` package, `PPtree` algorithm was translated to Rcpp. 
`PPforest` package utilizes a number of R packages some of them included in "suggests" not to load them all at package start-up.

The development version of`PPforest` can be installed from github using:

```r
library(devtools)
install_github("natydasilva/PPforest")
library(PPforest)
```


Overview PPforest package
-------------------------

`PPforest` package implements a classification random forest using projection pursuit classification trees. The following table present all the functions in `PPforest` package.

| Function |Description |
| ----------------- | --------------------------------------------------------------  | 
|baggtree|For each bootstrap sample grow a projection pursuit tree (PPtree object).|
|node_data|Data structure with the  projected and boundary by node and class|
|permute_importance|Obtain the permuted importance variable measure|
|ppf_avg_imp| Computes a global importance measure for a PPforest object, average importance measure for a pptree over all the trees.| 
|PPclassify| Predict class for the test set and calculate prediction error after finding the PPtree structure|
|ppf_global_imp| Computes a global importance measure for a PPforest object|
|PPforest|Runs a Projection pursuit random forest|
|PPtree_split|Projection pursuit classification tree with random variable selection in each split|
|print.PPforest| Print PPforest object|
|predict.PPforest|Predict a PPforest object for newdata|
|ternary_str|Data structure with the  projected and boundary by node and class|
|tree_pred|Obtain predicted class for new data using PPforest t.|

Also `PPforest` package includes some data set that were used to test the predictive performance of our method. The data sets included are: crab, fishcatch, glass, image, leukemia, lymphoma NCI60, parkinson and wine.

 Example
------------
Australian crab data set will be used as example. This data contains measurements on rock crabs of the genus Leptograpsus. There are 200 observations from two species (blue and orange) and for each specie (50 in each one) there are 50 males and 50 females. Class variable has 4 classes with the combinations of specie and sex (BlueMale, BlueFemale, OrangeMale and OrangeFemale). The data were collected on site at Fremantle, Western Australia. For each specimen, five measurements were made, using vernier calipers.

1. FL the size of the frontal lobe length, in mm
2. RW rear width, in mm
3. CL length of mid line of the carapace, in mm
4. CW maximum width of carapace, in mm
5. BD depth of the body; for females, measured after displacement of the abdomen, in mm


```PPforest``` function runs a projection pursuit random forest.  The arguments are a data frame with the data information, class with the name of the class variable argument.  size.tr to specify the proportion of observations using in the training. Using this function we have the option to split the data in training and test using size.tr directly. `size.tr` is the proportion of data used in the training and the test proportion will be 1- `size.tr`.
The number of trees in the forest is specified using the argument `m`. The argument size.p is the sample proportion of the variables used in each node split, `PPmethod` is the projection pursuit index to be optimized,  two options LDA and PDA are available.

```r 

pprf.crab <- PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7, m = 200,
                                size.p =  .5,  PPmethod = 'LDA',  parallel =TRUE, cores = 2)

pprf.crab

Call:
 PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7,      m = 200, PPmethod = "LDA", size.p = 0.5, parallel = TRUE,      cores = 2) 
               Type of random forest: Classification
                     Number of trees: 200
No. of variables tried at each split: 3

        OOB estimate of  error rate: 4.29%
Confusion matrix:
             BlueFemale BlueMale OrangeFemale OrangeMale class.error
BlueFemale           35        0            0          0        0.00
BlueMale              3       32            0          0        0.09
OrangeFemale          0        0           32          3        0.09
OrangeMale            0        0            0         35        0.00

 

```

`PPforest` print a summary result from the model with the confusion matrix information and the oob-error rate in a similar way randomForest packages does.

This function returns the predicted values of the training data, training error, test error and predicted test values. Also there is the information about out of bag error for the forest and also for each tree in the forest. Bootstrap samples, output of all the trees in the forest from , proximity matrix and vote matrix, number of trees grown in the forest, number of predictor variables selected to use for splitting at each node. Confusion matrix of the prediction (based on OOb data), the training data and test data and vote matrix are also returned.

The printed version of a `PPforest` object follows the `randomForest` printed version to make them comparable. Based on confusion matrix, we can observe that the biggest error is for BlueMale class. Most of the wrong classified values are between BlueFemale and BlueMale.

```r
str(pprf.crab, max.level = 1)
List of 22
 $ predicting.training: Factor w/ 4 levels "BlueFemale","BlueMale",..: 2 1 2 2 2 2 1 2 1 2 ...
 $ training.error     : num 0.0429
 $ prediction.test    : Factor w/ 4 levels "BlueFemale","BlueMale",..: 1 2 2 2 1 2 2 2 1 2 ...
 $ error.test         : num 0.05
 $ oob.error.forest   : num 0.0429
 $ oob.error.tree     : num [1:200, 1] 0.3158 0.1404 0.16 0.102 0.0769 ...
 $ boot.samp          :List of 200
 $ output.trees       :List of 200
 $ proximity          : num [1:140, 1:140] 0 0.68 0.855 0.845 0.785 0.87 0.35 0.71 0.395 0.48 ...
 $ votes              : num [1:140, 1:4] 0.415 0.802 0.481 0.191 0.2 ...
  ..- attr(*, "dimnames")=List of 2
 $ prediction.oob     : Factor w/ 4 levels "BlueFemale","BlueMale",..: 2 1 2 2 2 2 1 2 1 2 ...
 $ n.tree             : num 200
 $ n.var              : int 3
 $ type               : chr "Classification"
 $ confusion          : num [1:4, 1:5] 35 3 0 0 0 32 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 2
 $ call               : language PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7, m = 200, PPmethod = "LDA",      size.p = 0.5, para| __truncated__
 $ train              :'data.frame':	140 obs. of  6 variables:
 $ test               :'data.frame':	60 obs. of  5 variables:
 $ vote.mat           : num [1:200, 1:140] 1 1 2 2 2 2 1 2 1 1 ...
  ..- attr(*, "dimnames")=List of 2
 $ vote.mat_cl        : chr [1:4] "BlueFemale" "BlueMale" "OrangeFemale" "OrangeMale"
 $ class.var          : chr "Type"
 $ oob.obs            : num [1:200, 1:140] 0 1 1 0 1 0 0 0 0 1 ...
 - attr(*, "class")= chr "PPforest"
```

The `PPforest` object can be used to predict new data using the `predict` function. The predicted values are returned as a factor with the class levels.
```{r}

pred.crab <- predict(pprf.crab, newdata = crab[1:10,-1 ])

pred.crab[[3]] 
 [1] BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale
 [9] BlueFemale BlueFemale
Levels: BlueFemale BlueMale OrangeFemale OrangeMale
```

