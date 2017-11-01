
The `PPforest` package
======================
N. da Silva, D. Cook & E. Lee 
2017-11-01


<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/natydasilva/PPforest.svg?branch=master)](https://travis-ci.org/natydasilva/PPforest)


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
|predict.PPforest|Predict class for the test set and calculate prediction error|
|ternary_str|Data structure with the  projected and boundary by node and class|
|tree_pred|Obtain predicted class for new data using PPforest t.|

Also `PPforest` package includes some data set that were used to test the predictive performance of our method. The data sets included are: crab, fishcatch, glass, image, leukemia, lymphoma NCI60, parkinson and wine.

```PPforest``` function runs a projection pursuit random forest.  The arguments are a data frame with the data information, class with the name of the class variable argument.  size.tr to specify the proportion of observations using in the training. Using this function we have the option to split the data in training and test using size.tr directly. `size.tr` is the proportion of data used in the training and the test proportion will be 1- `size.tr`.
The number of trees in the forest is specified using the argument `m`. The argument size.p is the sample proportion of the variables used in each node split, `PPmethod` is the projection pursuit index to be optimized,  two options LDA and PDA are available.

```{r }

pprf.crab <- PPforest::PPforest(data = crab, class = "Type", size.tr = 1, m = 200,
                                size.p =  .5,  PPmethod = 'LDA',  parallel =TRUE, cores = 2)

pprf.crab
```

`PPforest` print a summary result from the model with the confusion matrix information and the oob-error rate in a similar way randomForest packages does.

This function returns the predicted values of the training data, training error, test error and predicted test values. Also there is the information about out of bag error for the forest and also for each tree in the forest. Bootstrap samples, output of all the trees in the forest from , proximity matrix and vote matrix, number of trees grown in the forest, number of predictor variables selected to use for splitting at each node. Confusion matrix of the prediction (based on OOb data), the training data and test data and vote matrix are also returned.

The printed version of a `PPforest` object follows the `randomForest` printed version to make them comparable. Based on confusion matrix, we can observe that the biggest error is for BlueMale class. Most of the wrong classified values are between BlueFemale and BlueMale.

```{r}
str(pprf.crab,max.level=1 )
```
