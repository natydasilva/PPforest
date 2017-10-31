
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/natydasilva/PPforest.svg?branch=master)](https://travis-ci.org/natydasilva/PPforest)


##Introduction

The `PPforest` package (projection pursuit random forest) contains functions to run a projection pursuit random forest for classification problems. This method utilize combinations of variables in each tree construction.  In a random forest each split is based on a single variable, chosen from a subset of predictors. In the `PPforest`, each split is based on a linear combination of randomly chosen variables. The linear combination is computed by optimizing a projection pursuit index, to get a projection of the variables that best separates the classes. The `PPforest` uses the `PPtree` algorithm, which fits a single tree to the data. Utilizing linear combinations of variables to separate classes takes the correlation between variables into account, and can outperform the basic forest when separations between groups occurs on combinations of variables. Two projection pursuit indexes, LDA and PDA, are used for `PPforest`.

To improve the speed performance `PPforest` package, `PPtree` algorithm was translated to Rcpp. 
`PPforest` package utilizes a number of R packages some of them included in "suggests" not to load them all at package start-up.
The development version of`PPforest` can be installed from github using:

```{r,echo=FALSE, message=FALSE, warning=FALSE,eval=FALSE}
library(devtools)
#install_github("natydasilva/PPforest")
library(PPforest)
```

##Projection pursuit classification forest

In `PPforest`, projection pursuit classification trees  are used as the individual model to be combined in the forest. The original algorithm is in `PPtreeViz` package,  we translate the original tree algorithm into `Rcpp` to improve the speed performance to run the forest.  

One important characteristic of PPtree is that treats the data always as a two-class system,  when the classes are more than two the algorithm uses a two step  projection pursuits optimization in every node split.
Let  $(X_i,y_i)$ the data set, $X_i$ is a  p-dimensional vector of explanatory variables and  $y_i\in {1,2,\ldots G}$ represents class information with $i=1,\ldots n$.

In the first step optimize a projection pursuit index to find an optimal one-dimension projection $\alpha^*$ for separating all classes in the current data. With the projected data redefine the problem in a two class problem by comparing means, and assign a new label $G1$ or $G2$ to each observation, a new variable $y_i^*$ is created.  The new groups $G1$ and $G2$ can contain more than one original classes. Next step is to find an optimal one-dimensional projection $\alpha$, using $(X_i,y_i^*)$ to separate the two class problem $G1$ and $G2$. The best separation of $G1$ and $G2$ is determine in this step and the decision rule is defined for the current node, if $\sum_{i=1}^p \alpha_i M1< c$ then assign $G1$ to the left node else assign $G2$ to the right node, where $M1$ is the mean of $G1$.
For each groups we can repeat all the previous steps until $G1$ and $G2$ have only one class from the original classes. Base on this process to grow the tree, the depth of PPtree is at most the number of classes because one class is assigned only to one final node.

Trees from `PPtree` algorithm are simple, they use the association between variables to find separation. If a linear boundary exists, `PPtree` produces a tree without misclassification.

Projection pursuit random forest algorithm description


1. Let N the number of cases in the training set $\Theta=(X,Y)$, $B$ bootstrap samples from the training set are taking (samples of size N with replacement).

2. For each bootstrap sample a \verb PPtree  is grown to the largest extent possible $h(x, {\Theta_k})$. No pruning. This tree is grown using step 3 modification.

3. Let M the number of input variables, a number of $m<<M$ variables are selected at random at each node and the best split based on a linear combination of these randomly chosen variables. The linear combination is computed by optimizing a projection pursuit index, to get a projection of the variables that best separates the classes.

4.  Predict the classes of each case not included in the bootstrap sample and compute oob error.

5.  Based on majority vote predict the class for new data.

###Overview PPforest package
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

