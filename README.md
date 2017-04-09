
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/natydasilva/PPforest.svg?branch=master)](https://travis-ci.org/natydasilva/PPforest)

# PPforest
A random forest is an ensemble learning method, build on bagged trees with random feature selection. 
The bagging provides power for classification because it yields information about variable importance,
predictive error and proximity of observations. This research adapts the random forest to utilize combinations 
of variables in the tree construction, which we call the projection pursuit classification random forest (PPF). 
In a random forest each split is based on a single variable, chosen from a subset of predictors. In the PPF, each 
split is based on a linear combination of randomly chosen variables from a subset of predictors. The linear
combination is computed by optimizing a projection pursuit index, to get a projection of the variables that best 
separates the classes. The PPF uses the PPtree algorithm (Lee, Cook, Park, Lee et al. 2013) , which fits a single tree to the data. 
Utilizing linear combinations of variables to separate classes takes the correlation between variables into account, 
and can outperform the basic forest when separations between groups occur on combinations of variables. Two projection 
pursuit indexes, LDA and PDA, are used for PPF. The methods are implemented into PPforest R package available in this repository.
