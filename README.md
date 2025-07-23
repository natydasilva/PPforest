
`PPforest` package
======================
Natalia da Silva, Dianne Cook & Eun-Kyung Lee 



<img src="man/figures/PPforest.png" align="right" alt="" width="160" />


[![CRAN Status](https://www.r-pkg.org/badges/version/PPforest)]( https://CRAN.R-project.org/package=PPforest)[![CRAN\_Download\_Badge](https://cranlogs.r-pkg.org/badges/grand-total/PPforest)](https://cran.r-project.org/package=PPforest) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/PPforest)](https://www.r-pkg.org/pkg/PPforest)


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
|ternary_str|Data structure with the projected and boundary by node and class|
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
set.seed(123)
pprf.crab <- PPforest::PPforest(data = crab, y = "Type", size.tr = 0.7, m = 200,
                                size.p =  .5,  PPmethod = 'LDA',  parallel =TRUE, cores = 2)

pprf.crab

The downloaded binary packages are in
	/var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//RtmpirXASo/downloaded_packages
── R CMD build ────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/Users/nataliadasilva/Documents/Research/R_packages/PPforest/DESCRIPTION’ ...
─  preparing ‘PPforest’: (623ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘PPforest_0.1.4.tar.gz’
   
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
  /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//RtmpirXASo/PPforest_0.1.4.tar.gz \
  --install-tests 
* installing to library ‘/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library’
* installing *source* package ‘PPforest’ ...
** this is package ‘PPforest’ version ‘0.1.4’
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 17.0.0 (clang-1700.0.13.5)’
using SDK: ‘’
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c optim_index.cpp -o optim_index.o
g++ -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -L/opt/homebrew/opt/libomp/lib -lomp -o PPforest.so RcppExports.o optim_index.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5 -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/usr/local/Cellar/gcc/5.3.0/lib/gcc/5' not found
installing to /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/00LOCK-PPforest/00new/PPforest/libs
** R
** data
*** moving datasets to lazyload DB
** byte-compile and prepare package for lazy loading
Note: possible error in 'baggtree(data = train, ': unused argument (class = y) 
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (PPforest)
> a <- GGally::ggpairs(PPforest::crab,
+     columns = 2:6,
+     ggplot2::aes(colour = Type, alpha=.1),
+     lower = list(continuous = 'points'),
+     axisLabels='none',
+     upper=list(continuous='blank')
+      , legend = NULL)
> capmatrix<-"Scatter plot matrix of crab data "
> a
> Tree.crab <- PPforest::PPtree_split("Type~.", data = crab, PPmethod = "LDA", size.p = 0.6)
Warning: internal error 1 in R_decompress1 with libdeflate
Error: lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/PPforest/R/PPforest.rdb' is corrupt
Error during wrapup: not that many frames on the stack
Error: no more error handlers available (recursive errors?); invoking 'abort' restart

> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Error in ns_s3_methods(package) : 
  lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/PPforest/R/PPforest.rdb' is corrupt
In addition: Warning message:
In ns_s3_methods(package) :
  internal error 1 in R_decompress1 with libdeflate

Restarting R session...
> library(devtools)
Loading required package: usethis
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Restarting R session...
> library(PPforest)
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

Restarting R session...
> library(PPforest)
> 'pprf.crab <- PPforest(data = crab, y = 'Type',
Error: unexpected symbol in "'pprf.crab <- PPforest(data = crab, y = 'Type"
In addition: Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no',  size.tr = 1, m = 200, size.p = .5, 
+  PPmethod = 'LDA', parallel = TRUE, cores = 2)
> devtools:document()
Error: object 'devtools' not found

> library(devtools)
Loading required package: usethis
> devtools:document()
Error: object 'devtools' not found

> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Writing ppf_global_imp.Rd
Restarting R session...
> library(PPforest)
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Writing ternary_str.Rd
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

Restarting R session...
> library(PPforest)
> devtools::install()
── R CMD build ────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/Users/nataliadasilva/Documents/Research/R_packages/PPforest/DESCRIPTION’ ...
─  preparing ‘PPforest’: (626ms)
✔  checking DESCRIPTION meta-information
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘PPforest_0.1.4.tar.gz’
   
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
  /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//RtmpwLV23j/PPforest_0.1.4.tar.gz \
  --install-tests 
* installing to library ‘/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library’
* installing *source* package ‘PPforest’ ...
** this is package ‘PPforest’ version ‘0.1.4’
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 17.0.0 (clang-1700.0.13.5)’
using SDK: ‘’
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c optim_index.cpp -o optim_index.o
g++ -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -L/opt/homebrew/opt/libomp/lib -lomp -o PPforest.so RcppExports.o optim_index.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5 -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/usr/local/Cellar/gcc/5.3.0/lib/gcc/5' not found
installing to /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/00LOCK-PPforest/00new/PPforest/libs
** R
** data
*** moving datasets to lazyload DB
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (PPforest)
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

> vignette("PPforest")
Warning messages:
1: In get(name, envir = env) :
  internal error 1 in R_decompress1 with libdeflate
2: In get(name, envir = env) :
  internal error 1 in R_decompress1 with libdeflate
3: In get(name, envir = env) :
  internal error 1 in R_decompress1 with libdeflate
4: In get(name, envir = env) : restarting interrupted promise evaluation
5: In get(name, envir = env) :
  internal error 1 in R_decompress1 with libdeflate
6: vignette ‘PPforest’ not found 

Restarting R session...
> library(PPforest)
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> devtools::install()
── R CMD build ────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/Users/nataliadasilva/Documents/Research/R_packages/PPforest/DESCRIPTION’ ...
─  preparing ‘PPforest’: (622ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘PPforest_0.1.4.tar.gz’
   
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
  /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmprtnczk/PPforest_0.1.4.tar.gz \
  --install-tests 
* installing to library ‘/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library’
* installing *source* package ‘PPforest’ ...
** this is package ‘PPforest’ version ‘0.1.4’
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 17.0.0 (clang-1700.0.13.5)’
using SDK: ‘’
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c optim_index.cpp -o optim_index.o
g++ -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -L/opt/homebrew/opt/libomp/lib -lomp -o PPforest.so RcppExports.o optim_index.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5 -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/usr/local/Cellar/gcc/5.3.0/lib/gcc/5' not found
installing to /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/00LOCK-PPforest/00new/PPforest/libs
** R
** data
*** moving datasets to lazyload DB
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (PPforest)
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 1, m = 3, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
Error: lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/PPforest/R/PPforest.rdb' is corrupt
In addition: Warning message:
internal error 1 in R_decompress1 with libdeflate 

> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 1, m = 3, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
Error: lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/PPforest/R/PPforest.rdb' is corrupt
In addition: Warning messages:
1: restarting interrupted promise evaluation 
2: internal error 1 in R_decompress1 with libdeflate 

Restarting R session...
> library(devtools)
Loading required package: usethis
> devtools:document()
Error: object 'devtools' not found

> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Restarting R session...
> library(PPforest)
> pprf.crab <- PPforest(data = crab, y = 'Type',
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

+ pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 1, m = 3, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
+ 
> set.seed(123)
> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 1, m = 3, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
> 'pprf.crab 
+ 
> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 1, m = 3, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
> pprf.crab

Call:
 PPforest(data = crab, y = "Type", xstd = "no", size.tr = 1, m = 3,      PPmethod = "LDA", size.p = 1, parallel = TRUE, cores = 2,      rule = 1) 
               Type of random forest: Classification
                     Number of trees: 3
No. of variables tried at each split: 5

        OOB estimate of  error rate: 26.5%
Confusion matrix:
             BlueFemale BlueMale OrangeFemale OrangeMale class.error
BlueFemale           48        2            0          0        0.04
BlueMale             21       29            0          0        0.42
OrangeFemale         14        0           33          3        0.34
OrangeMale           13        0            0         37        0.26
> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 1, m = 100, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
> pprf.crab

Call:
 PPforest(data = crab, y = "Type", xstd = "no", size.tr = 1, m = 100,      PPmethod = "LDA", size.p = 1, parallel = TRUE, cores = 2,      rule = 1) 
               Type of random forest: Classification
                     Number of trees: 100
No. of variables tried at each split: 5

        OOB estimate of  error rate: 6%
Confusion matrix:
             BlueFemale BlueMale OrangeFemale OrangeMale class.error
BlueFemale           47        3            0          0        0.06
BlueMale              6       44            0          0        0.12
OrangeFemale          0        0           47          3        0.06
OrangeMale            0        0            0         50        0.00
> pprf.crab <- PPforest(data = crab, y = 'Type',
+  xstd = 'no', size.tr = 0.8, m = 100, size.p = 1, 
+  PPmethod = 'LDA' , parallel = TRUE, cores = 2, rule = 1)
> pprf.crab

Call:
 PPforest(data = crab, y = "Type", xstd = "no", size.tr = 0.8,      m = 100, PPmethod = "LDA", size.p = 1, parallel = TRUE, cores = 2,      rule = 1) 
               Type of random forest: Classification
                     Number of trees: 100
No. of variables tried at each split: 5

        OOB estimate of  error rate: 6.25%
Confusion matrix:
             BlueFemale BlueMale OrangeFemale OrangeMale class.error
BlueFemale           37        3            0          0        0.07
BlueMale              4       36            0          0        0.10
OrangeFemale          0        0           37          3        0.07
OrangeMale            0        0            0         40        0.00
> pprf.crab$error.test
[1] 0.05
> devtools:document()
Error: object 'devtools' not found

> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Writing PPforest.Rd
Writing node_data.Rd
Writing permute_importance.Rd
Writing ppf_avg_imp.Rd
Writing ppf_global_imp.Rd
Writing predict.PPforest.Rd
Writing ternary_str.Rd
Restarting R session...
> library(PPforest)
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

Restarting R session...
> library(PPforest)
> vignette("PPforest")
Warning messages:
1: In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/
2: vignette ‘PPforest’ not found 

> devtools::build_vignettes()
ℹ Installing PPforest in temporary library
ℹ Building vignettes for PPforest
--- re-building ‘PPforest-vignette.Rmd’ using rmarkdown


processing file: PPforest-vignette.Rmd
                                                                                                      
output file: PPforest-vignette.knit.md

/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc +RTS -K512m -RTS PPforest-vignette.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /Users/nataliadasilva/Documents/Research/R_packages/PPforest/vignettes/PPforest-vignette.html --lua-filter /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/lua/latex-div.lua --lua-filter /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/lua/table-classes.lua --embed-resources --standalone --section-divs --template /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmd/h/default.html --highlight-style pygments --css /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/templates/html_vignette/resources/vignette.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmp31wW3F/rmarkdown-strd8ab67ca318b.html --citeproc 

Output created: PPforest-vignette.html
--- finished re-building ‘PPforest-vignette.Rmd’

ℹ Copying vignettes
ℹ Moving PPforest-vignette.html and PPforest-vignette.R to doc/
ℹ Copying PPforest-vignette.Rmd to doc/
ℹ Building vignette index
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

Restarting R session...
> library(PPforest)
> remotes::install_local(build_vignettes = TRUE)
Skipping install of 'PPforest' from a local remote, the SHA1 (0.1.4) has not changed since last install.
  Use `force = TRUE` to force installation
> remotes::install_local(PPforest, build_vignettes = TRUE)
Error in path.expand(path) : invalid 'path' argument

> remotes::install_local(build_vignettes = TRUE)
Skipping install of 'PPforest' from a local remote, the SHA1 (0.1.4) has not changed since last install.
  Use `force = TRUE` to force installation
> devtools::install(PPforest,build_vignettes = TRUE)
Error in `package_file()`:
! `path` must be a string.
Run `rlang::last_trace()` to see where the error occurred.

> vignette(package="ggplot2") 
> vignette(package="ggplot2-specs") 
Error in find.package(package, lib.loc) : 
  there is no package called ‘ggplot2-specs’

> vignette(package="ggplot2") 
> vignette(package="ggplot2/ggplot2-specs") 
Error in find.package(package, lib.loc) : 
  there is no package called ‘ggplot2/ggplot2-specs’

> vignette( "ggplot2-specs", package="ggplot2")
> vignette( "PPforest", package="PPforest")
Warning message:
vignette ‘PPforest’ not found 

> install.packages('PPforest')
Error in install.packages : Updating loaded packages

Restarting R session...
> install.packages("PPforest")
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/PPforest_0.1.3.tgz'
Content type 'application/x-gzip' length 2033156 bytes (1.9 MB)
==================================================
downloaded 1.9 MB


The downloaded binary packages are in
	/var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmpi4FUi6/downloaded_packages
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

Restarting R session...
> library(PPforest)
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> cran_stats("ODRF")
Error in cran_stats("ODRF") : could not find function "cran_stats"

> install.packages('dlstats')
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/dlstats_0.1.7.tgz'
Content type 'application/x-gzip' length 366165 bytes (357 KB)
==================================================
downloaded 357 KB


The downloaded binary packages are in
	/var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmp52HY9a/downloaded_packages
> library(dlstats)
> cran_stats("ODRF")
        start        end downloads package
1  2023-02-01 2023-02-28         3    ODRF
2  2023-03-01 2023-03-31       321    ODRF
3  2023-04-01 2023-04-30       206    ODRF
4  2023-05-01 2023-05-31       247    ODRF
5  2023-06-01 2023-06-30       306    ODRF
6  2023-07-01 2023-07-31       199    ODRF
7  2023-08-01 2023-08-31       234    ODRF
8  2023-09-01 2023-09-30       190    ODRF
9  2023-10-01 2023-10-31       181    ODRF
10 2023-11-01 2023-11-30       219    ODRF
11 2023-12-01 2023-12-31       190    ODRF
12 2024-01-01 2024-01-31       407    ODRF
13 2024-02-01 2024-02-29       165    ODRF
14 2024-03-01 2024-03-31       198    ODRF
15 2024-04-01 2024-04-30       230    ODRF
16 2024-05-01 2024-05-31       222    ODRF
17 2024-06-01 2024-06-30       175    ODRF
18 2024-07-01 2024-07-31       217    ODRF
19 2024-08-01 2024-08-31       233    ODRF
20 2024-09-01 2024-09-30       193    ODRF
21 2024-10-01 2024-10-31       184    ODRF
22 2024-11-01 2024-11-30       208    ODRF
23 2024-12-01 2024-12-31       201    ODRF
24 2025-01-01 2025-01-31       269    ODRF
25 2025-02-01 2025-02-28       225    ODRF
26 2025-03-01 2025-03-31       207    ODRF
27 2025-04-01 2025-04-30       337    ODRF
28 2025-05-01 2025-05-31       231    ODRF
29 2025-06-01 2025-06-30       247    ODRF
30 2025-07-01 2025-07-22       178    ODRF
> cran_stats("PPforest")
        start        end downloads  package
1  2020-01-01 2020-01-31       518 PPforest
2  2020-02-01 2020-02-29       446 PPforest
3  2020-03-01 2020-03-31       485 PPforest
4  2020-04-01 2020-04-30       530 PPforest
5  2020-05-01 2020-05-31       676 PPforest
6  2020-06-01 2020-06-30       576 PPforest
7  2020-07-01 2020-07-31       484 PPforest
8  2020-08-01 2020-08-31       481 PPforest
9  2020-09-01 2020-09-30       606 PPforest
10 2020-10-01 2020-10-31       480 PPforest
11 2020-11-01 2020-11-30       574 PPforest
12 2020-12-01 2020-12-31       530 PPforest
13 2021-01-01 2021-01-31       426 PPforest
14 2021-02-01 2021-02-28       491 PPforest
15 2021-03-01 2021-03-31       478 PPforest
16 2021-04-01 2021-04-30       420 PPforest
17 2021-05-01 2021-05-31       438 PPforest
18 2021-06-01 2021-06-30       299 PPforest
19 2021-07-01 2021-07-31       367 PPforest
20 2021-08-01 2021-08-31       465 PPforest
21 2021-09-01 2021-09-30       329 PPforest
22 2021-10-01 2021-10-31       696 PPforest
23 2021-11-01 2021-11-30       805 PPforest
24 2021-12-01 2021-12-31       490 PPforest
25 2022-01-01 2022-01-31       301 PPforest
26 2022-02-01 2022-02-28       318 PPforest
27 2022-03-01 2022-03-31       451 PPforest
28 2022-04-01 2022-04-30       363 PPforest
29 2022-05-01 2022-05-31       403 PPforest
30 2022-06-01 2022-06-30       303 PPforest
31 2022-07-01 2022-07-31       386 PPforest
32 2022-08-01 2022-08-31       298 PPforest
33 2022-09-01 2022-09-30       439 PPforest
34 2022-10-01 2022-10-31       400 PPforest
35 2022-11-01 2022-11-30       358 PPforest
36 2022-12-01 2022-12-31       227 PPforest
37 2023-01-01 2023-01-31       217 PPforest
38 2023-02-01 2023-02-28       227 PPforest
39 2023-03-01 2023-03-31       327 PPforest
40 2023-04-01 2023-04-30       281 PPforest
41 2023-05-01 2023-05-31       307 PPforest
42 2023-06-01 2023-06-30       217 PPforest
43 2023-07-01 2023-07-31       223 PPforest
44 2023-08-01 2023-08-31       249 PPforest
45 2023-09-01 2023-09-30       206 PPforest
46 2023-10-01 2023-10-31       211 PPforest
47 2023-11-01 2023-11-30       229 PPforest
48 2023-12-01 2023-12-31       213 PPforest
49 2024-01-01 2024-01-31       466 PPforest
50 2024-02-01 2024-02-29       180 PPforest
51 2024-03-01 2024-03-31       225 PPforest
52 2024-04-01 2024-04-30       241 PPforest
53 2024-05-01 2024-05-31       245 PPforest
54 2024-06-01 2024-06-30       178 PPforest
55 2024-07-01 2024-07-31       280 PPforest
56 2024-08-01 2024-08-31       257 PPforest
57 2024-09-01 2024-09-30       249 PPforest
58 2024-10-01 2024-10-31       211 PPforest
59 2024-11-01 2024-11-30       192 PPforest
60 2024-12-01 2024-12-31       199 PPforest
61 2025-01-01 2025-01-31       267 PPforest
62 2025-02-01 2025-02-28       232 PPforest
63 2025-03-01 2025-03-31       171 PPforest
64 2025-04-01 2025-04-30       218 PPforest
65 2025-05-01 2025-05-31       175 PPforest
66 2025-06-01 2025-06-30       200 PPforest
67 2025-07-01 2025-07-22       159 PPforest
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Restarting R session...
> library(PPforest)
> vignette("PPforest")
Warning messages:
1: In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/
2: vignette ‘PPforest’ not found 

> devtools::install(build_vignettes = TRUE)
These packages have more recent versions available.
It is recommended to update all of them.
Which would you like to update?

1: All                                
2: CRAN packages only                 
3: None                               
4: rprojroot  (2.0.4 -> 2.1.0 ) [CRAN]
5: commonmark (1.9.5 -> 2.0.0 ) [CRAN]
6: ggstats    (0.9.0 -> 0.10.0) [CRAN]
7: GGally     (2.2.1 -> 2.3.0 ) [CRAN]

Enter one or more numbers, or an empty line to skip updates: 1
rprojroot  (2.0.4 -> 2.1.0 ) [CRAN]
commonmark (1.9.5 -> 2.0.0 ) [CRAN]
S7         (NA    -> 0.2.0 ) [CRAN]
ggstats    (0.9.0 -> 0.10.0) [CRAN]
GGally     (2.2.1 -> 2.3.0 ) [CRAN]
Installing 5 packages: rprojroot, commonmark, S7, ggstats, GGally
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/rprojroot_2.1.0.tgz'
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/commonmark_2.0.0.tgz'
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/S7_0.2.0.tgz'
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/ggstats_0.10.0.tgz'
trying URL 'https://cran.rstudio.com/bin/macosx/big-sur-arm64/contrib/4.5/GGally_2.2.1.tgz'

The downloaded binary packages are in
	/var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmpah1M5n/downloaded_packages
── R CMD build ────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/Users/nataliadasilva/Documents/Research/R_packages/PPforest/DESCRIPTION’ ...
─  preparing ‘PPforest’: (667ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  installing the package to build vignettes
✔  creating vignettes (13.4s)
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts (424ms)
─  checking for empty or unneeded directories
─  building ‘PPforest_0.1.4.tar.gz’
   
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
  /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmpah1M5n/PPforest_0.1.4.tar.gz \
  --install-tests 
* installing to library ‘/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library’
* installing *source* package ‘PPforest’ ...
** this is package ‘PPforest’ version ‘0.1.4’
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 17.0.0 (clang-1700.0.13.5)’
using SDK: ‘’
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c optim_index.cpp -o optim_index.o
g++ -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -L/opt/homebrew/opt/libomp/lib -lomp -o PPforest.so RcppExports.o optim_index.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5 -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/usr/local/Cellar/gcc/5.3.0/lib/gcc/5' not found
installing to /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/00LOCK-PPforest/00new/PPforest/libs
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (PPforest)
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> library(PPforest)
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> browseVignettes()
> devtools::install(build_vignettes = TRUE)
Error in eapply(pkgload::ns_env(pkg$package), force, all.names = TRUE) : 
  lazy-load database '/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/PPforest/R/PPforest.rdb' is corrupt
In addition: Warning message:
In eapply(pkgload::ns_env(pkg$package), force, all.names = TRUE) :
  internal error 1 in R_decompress1 with libdeflate

Restarting R session...
> library(PPforest)
> browseVignettes().
Error: unexpected symbol in "browseVignettes()."

> browseVignettes()
> devtools::build_vignettes()
ℹ Installing PPforest in temporary library
ℹ Building vignettes for PPforest
--- re-building ‘PPforest-vignette.Rmd’ using rmarkdown


processing file: PPforest-vignette.Rmd
                                                                                                      
output file: PPforest-vignette.knit.md

/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/pandoc +RTS -K512m -RTS PPforest-vignette.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /Users/nataliadasilva/Documents/Research/R_packages/PPforest/vignettes/PPforest-vignette.html --lua-filter /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/lua/latex-div.lua --lua-filter /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/lua/table-classes.lua --embed-resources --standalone --section-divs --template /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmd/h/default.html --highlight-style pygments --css /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/rmarkdown/rmarkdown/templates/html_vignette/resources/vignette.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//RtmpKikogx/rmarkdown-strf1a86c149f95.html --citeproc 

Output created: PPforest-vignette.html
--- finished re-building ‘PPforest-vignette.Rmd’

ℹ Copying vignettes
ℹ Moving PPforest-vignette.html and PPforest-vignette.R to doc/
ℹ Copying PPforest-vignette.Rmd to doc/
ℹ Building vignette index
> browseVignettes()
> vignettes("PPforest")
Error in vignettes("PPforest") : could not find function "vignettes"

> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> devtools::install()
── R CMD build ────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/Users/nataliadasilva/Documents/Research/R_packages/PPforest/DESCRIPTION’ ...
─  preparing ‘PPforest’: (666ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts
─  checking for empty or unneeded directories
─  building ‘PPforest_0.1.4.tar.gz’
   
Running /Library/Frameworks/R.framework/Resources/bin/R CMD INSTALL \
  /var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T//Rtmpq7uxh5/PPforest_0.1.4.tar.gz \
  --install-tests 
* installing to library ‘/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library’
* installing *source* package ‘PPforest’ ...
** this is package ‘PPforest’ version ‘0.1.4’
** using staged installation
** libs
using C++ compiler: ‘Apple clang version 17.0.0 (clang-1700.0.13.5)’
using SDK: ‘’
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c RcppExports.cpp -o RcppExports.o
g++ -std=gnu++14 -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG  -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/Rcpp/include' -I'/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/RcppArmadillo/include' -I/opt/R/arm64/include -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp    -fPIC  -mtune=native -g -O2 -Wall -pedantic   -c optim_index.cpp -o optim_index.o
g++ -std=gnu++14 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -L/opt/homebrew/opt/libomp/lib -lomp -o PPforest.so RcppExports.o optim_index.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5 -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/usr/local/Cellar/gcc/5.3.0/lib/gcc/5' not found
installing to /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/00LOCK-PPforest/00new/PPforest/libs
** R
** data
*** moving datasets to lazyload DB
** byte-compile and prepare package for lazy loading
** help
*** installing help indices
*** copying figures
** building package indices
** installing vignettes
** testing if installed package can be loaded from temporary location
** checking absolute paths in shared objects and dynamic libraries
** testing if installed package can be loaded from final location
** testing if installed package keeps a record of temporary installation path
* DONE (PPforest)
> browseVignettes()
> vignette("PPforest")
Warning message:
vignette ‘PPforest’ not found 

> library(help="ODRF")
> library(ggplot2)
> library(dplyr)
Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union
> library(patchwork)  # para mostrar dos gráficos lado a lado
> 
> # Simular datos: relación no lineal con crecimiento exponencial
> set.seed(123)
> df <- data.frame(
+     x = 10^(runif(100, 0, 4)),  # valores entre 1 y 10000
+     ruido = rnorm(100, 0, 0.1)
+ ) |> mutate(y = log10(x) + ruido)
> 
> # Gráfico sin transformación (escala lineal)
> p1 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     labs(title = "Escala lineal", x = "x", y = "y") +
+     theme_minimal()
> 
> # Gráfico con eje X en escala logarítmica
> p2 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     scale_x_log10() +
+     labs(title = "Escala logarítmica (x)", x = "log10(x)", y = "y") +
+     theme_minimal()
> 
> # Mostrar lado a lado
> p1 + p2
> library(ggplot2)
> library(dplyr)
> library(tidyr)
> 
> # Simular datos
> set.seed(123)
> df <- tibble(
+     x = 10^(runif(100, 0, 4)),
+     ruido = rnorm(100, 0, 0.1)
+ ) |>
+     mutate(y = log10(x) + ruido)
> 
> # Reestructurar para comparar
> df_long <- df |>
+     mutate(x_log = log10(x)) |>
+     pivot_longer(cols = c(x, x_log), names_to = "escala", values_to = "x_val") |>
+     mutate(escala = ifelse(escala == "x", "Escala lineal", "Escala log10"))
> 
> # Graficar con facet
> ggplot(df_long, aes(x = x_val, y = y)) +
+     geom_point() +
+     facet_wrap(~ escala, scales = "free_x") +
+     labs(x = "x", y = "y") +
+     theme_minimal()
> ?scale_y_log10
> p1 + p2
> library(ggplot2)
> library(dplyr)
> library(patchwork)  # para mostrar dos gráficos lado a lado
> 
> # Simular datos: relación no lineal con crecimiento exponencial
> set.seed(123)
> df <- data.frame(
+     x = 10^(runif(100, 0, 4)),  # valores entre 1 y 10000
+     ruido = rnorm(100, 0, 0.1)
+ ) |> mutate(y = log10(x) + ruido)
> 
> # Gráfico sin transformación (escala lineal)
> p1 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     labs(title = "Escala lineal", x = "x", y = "y") +
+     theme_minimal()
> 
> # Gráfico con eje X en escala logarítmica
> p2 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     scale_x_log10() +
+     labs(title = "Escala logarítmica (x)", y = "log10(y)", x = "x") +
+     theme_minimal()
> 
> # Mostrar lado a lado
> p1 + p2
> library(ggplot2)
> library(dplyr)
> library(patchwork)  # para mostrar dos gráficos lado a lado
> 
> # Simular datos: relación no lineal con crecimiento exponencial
> set.seed(123)
> df <- data.frame(
+     x = 10^(runif(100, 0, 4)),  # valores entre 1 y 10000
+     ruido = rnorm(100, 0, 0.1)
+ ) |> mutate(y = log10(x) + ruido)
> 
> # Gráfico sin transformación (escala lineal)
> p1 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     labs(title = "Escala lineal", x = "x", y = "y") +
+     theme_minimal()
> 
> # Gráfico con eje X en escala logarítmica
> p2 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     scale_y_log10() +
+     labs(title = "Escala logarítmica (x)", y = "log10(y)", x = "x") +
+     theme_minimal()
> 
> # Mostrar lado a lado
> p1 + p2
Warning messages:
1: In transformation$transform(x) : NaNs produced
2: In scale_y_log10() :
  log-10 transformation introduced infinite values.
3: Removed 1 row containing missing values or values outside the scale range (`geom_point()`). 

> library(ggplot2)
> library(dplyr)
> library(patchwork)  # para mostrar dos gráficos lado a lado
> 
> # Simular datos: relación no lineal con crecimiento exponencial
> set.seed(123)
> df <- data.frame(
+     x = 10^(runif(100, 0, 4)),  # valores entre 1 y 10000
+     ruido = rnorm(100, 0, 0.1)
+ ) |> mutate(y = log10(x) + ruido)
> 
> # Gráfico sin transformación (escala lineal)
> p1 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     labs(title = "Escala lineal", x = "x", y = "y") +
+     theme_minimal()
> 
> # Gráfico con eje X en escala logarítmica
> p2 <- ggplot(df, aes(x = x, y = y)) +
+     geom_point() +
+     scale_x_log10() +
+     labs(title = "Escala logarítmica (x)", x = "log10(x)", y = "y") +
+     theme_minimal()
> 
> # Mostrar lado a lado
> p1 + p2
> usethis::use_release_issue()
✔ Setting active project to "/Users/nataliadasilva/Documents/Research/R_packages/PPforest".
Current version is 0.1.4.
What should the release version be? (0 to exit) 

1: major --> 1.0.0
2: minor --> 0.2.0
3: patch --> 0.1.5

Selection: 2
✔ Opening URL <https://github.com/natydasilva/PPforest/issues/3>.
> usethis::use_github_links()
✔ Adding "https://github.com/natydasilva/PPforest/issues" to BugReports.
ℹ There is 1 uncommitted file:
• DESCRIPTION
! Is it ok to commit it?

1: Not now
2: Definitely
3: No way

Selection: 3
Restarting R session...
> library(PPforest)
> usethis::use_release_issue()
✔ Setting active project to "/Users/nataliadasilva/Documents/Research/R_packages/PPforest".
Current version is 0.1.4.
What should the release version be? (0 to exit) 

1: major --> 1.0.0
2: minor --> 0.2.0
3: patch --> 0.1.5

Selection: 2
✔ Opening URL <https://github.com/natydasilva/PPforest/issues/4>.
> 
> 
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Restarting R session...
> library(PPforest)
> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Writing NCI60.Rd
Writing crab.Rd
Writing fishcatch.Rd
Writing glass.Rd
Writing image.Rd
Writing leukemia.Rd
Writing lymphoma.Rd
Writing olive.Rd
Writing parkinson.Rd
Writing wine.Rd
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

Restarting R session...
> library(PPforest)
> usethis::use_github_links()
✔ Setting active project to "/Users/nataliadasilva/Documents/Research/R_packages/PPforest".
ℹ There is 1 uncommitted file:
• DESCRIPTION
! Is it ok to commit it?

1: Nope
2: Negative
3: Definitely

Selection: 2
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Restarting R session...
> library(PPforest)
> usethis::use_github_links()
✔ Setting active project to "/Users/nataliadasilva/Documents/Research/R_packages/PPforest".
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

> usethis::use_release_issue()
Current version is 0.1.4.
What should the release version be? (0 to exit) 

1: major --> 1.0.0
2: minor --> 0.2.0
3: patch --> 0.1.5

Selection: 2
✔ Opening URL <https://github.com/natydasilva/PPforest/issues/5>.
> urlchecker::url_check()
! Warning: README.md:11:379 
Moved
[![CRAN Status](https://www.r-pkg.org/badges/version/PPforest)]( https://CRAN.R-project.org/package=PPforest)[![CRAN\_Download\_Badge](https://cranlogs.r-pkg.org/badges/grand-total/PPforest)](https://cran.r-project.org/package=PPforest) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/PPforest)](https://www.r-pkg.org/pkg/PPforest)[![Travis-CI Build Status](https://travis-ci.org/natydasilva/PPforest.svg?branch=master)](https://travis-ci.org/natydasilva/PPforest)
                                                                                                                                                                                                                                                                                                                                                                                          ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                                                                                                                                                                                                                                                                                                                                          https://app.travis-ci.com/natydasilva/PPforest
> devtools::check(remote = TRUE, manual = TRUE)
══ Documenting ════════════════════════════════════════════════════════════════════════════════════
ℹ Updating PPforest documentation
ℹ Loading PPforest

══ Building ═══════════════════════════════════════════════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ────────────────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/Users/nataliadasilva/Documents/Research/R_packages/PPforest/DESCRIPTION’
─  preparing ‘PPforest’: (706ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  installing the package to build vignettes
✔  creating vignettes (13.8s)
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts (456ms)
─  checking for empty or unneeded directories
─  building ‘PPforest_0.1.4.tar.gz’
   
══ Checking ═══════════════════════════════════════════════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_REMOTE_               : TRUE
• _R_CHECK_CRAN_INCOMING_                      : TRUE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ────────────────────────────────────────────────────────────────────────────────────
─  using log directory ‘/private/var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T/RtmpErGO1y/file10ae51f300630/PPforest.Rcheck’
─  using R version 4.5.0 (2025-04-11)
─  using platform: aarch64-apple-darwin20
─  R was compiled by
       Apple clang version 14.0.0 (clang-1400.0.29.202)
       GNU Fortran (GCC) 14.2.0
─  running under: macOS Sequoia 15.5
─  using session charset: UTF-8
─  using option ‘--as-cran’
✔  checking for file ‘PPforest/DESCRIPTION’
─  checking extension type ... Package
─  this is package ‘PPforest’ version ‘0.1.4’
─  package encoding: UTF-8
─  checking CRAN incoming feasibility ... [4s/19s] NOTE (19s)
   Maintainer: ‘Natalia da Silva <natalia.dasilva@fcea.edu.uy>’
   
   Found the following (possibly) invalid URLs:
     URL: https://travis-ci.org/natydasilva/PPforest (moved to https://app.travis-ci.com/natydasilva/PPforest)
       From: README.md
       Status: 200
       Message: OK
✔  checking package namespace information
✔  checking package dependencies (606ms)
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files ...
✔  checking for hidden files and directories ...
✔  checking for portable file names
✔  checking for sufficient/correct file permissions
✔  checking whether package ‘PPforest’ can be installed (8.5s)
─  used C++ compiler: ‘Apple clang version 17.0.0 (clang-1700.0.13.5)’
─  used SDK: ‘’
✔  checking installed package size ...
✔  checking package directory ...
N  checking for future file timestamps (1m 0.6s)
   unable to verify current time
✔  checking ‘build’ directory
✔  checking DESCRIPTION meta-information ...
✔  checking top-level files ...
✔  checking for left-over files ...
✔  checking index information ...
✔  checking package subdirectories ...
✔  checking code files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded ...
✔  checking whether the package can be loaded with stated dependencies ...
✔  checking whether the package can be unloaded cleanly ...
✔  checking whether the namespace can be loaded with stated dependencies ...
✔  checking whether the namespace can be unloaded cleanly ...
✔  checking use of S3 registration (2.6s)
✔  checking dependencies in R code ...
✔  checking S3 generic/method consistency ...
✔  checking replacement functions ...
✔  checking foreign function calls ...
✔  checking R code for possible problems (1.8s)
✔  checking Rd files ...
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references ...
✔  checking for missing documentation entries ...
✔  checking for code/documentation mismatches (491ms)
✔  checking Rd \usage sections ...
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples ...
✔  checking contents of ‘data’ directory ...
✔  checking data for non-ASCII characters ...
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves ...
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking line endings in Makefiles
✔  checking compilation flags in Makevars ...
✔  checking for GNU extensions in Makefiles ...
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles ...
✔  checking pragmas in C/C++ headers and code
✔  checking compilation flags used
✔  checking compiled code ...
✔  checking installed files from ‘inst/doc’ ...
✔  checking files in ‘vignettes’ ...
✔  checking examples (3s)
✔  checking for unstated dependencies in vignettes ...
✔  checking package vignettes ...
✔  checking re-building of vignette outputs (5.1s)
✔  checking PDF version of manual (2.3s)
N  checking HTML version of manual ...
   Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
   Please obtain a recent version of HTML Tidy by downloading a binary
   release or compiling the source code from <https://www.html-tidy.org/>.
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   See
     ‘/private/var/folders/sp/wknmqx111dgdcwbczth6r4540000gn/T/RtmpErGO1y/file10ae51f300630/PPforest.Rcheck/00check.log’
   for details.
   
── R CMD check results ──────────────────────────────────────────────────────── PPforest 0.1.4 ────
Duration: 1m 49.7s

❯ checking CRAN incoming feasibility ... [4s/19s] NOTE
  Maintainer: ‘Natalia da Silva <natalia.dasilva@fcea.edu.uy>’
  
  Found the following (possibly) invalid URLs:
    URL: https://travis-ci.org/natydasilva/PPforest (moved to https://app.travis-ci.com/natydasilva/PPforest)
      From: README.md
      Status: 200
      Message: OK

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
  Please obtain a recent version of HTML Tidy by downloading a binary
  release or compiling the source code from <https://www.html-tidy.org/>.

0 errors ✔ | 0 warnings ✔ | 3 notes ✖
> install.packages('tidy')
Warning in install.packages :
  package ‘tidy’ is not available for this version of R

A version of this package for your version of R might be available elsewhere,
see the ideas at
https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

> devtools::document()
ℹ Updating PPforest documentation
ℹ Loading PPforest
Restarting R session...
> library(PPforest)
> usethis::use_github_links()
✔ Setting active project to "/Users/nataliadasilva/Documents/Research/R_packages/PPforest".
Warning message:
In loadNamespace(package, lib.loc = libLoc) :
  package ‘PPforest’ has no 'package.rds' in Meta/

> usethis::use_release_issue()
Current version is 0.1.4.
What should the release version be? (0 to exit) 

1: major --> 1.0.0
2: minor --> 0.2.0
3: patch --> 0.1.5

Selection: 2
✔ Opening URL <https://github.com/natydasilva/PPforest/issues/6>.
> urlchecker::url_check()
✔ All URLs are correct!
Error in `challenge_uncommitted_changes()`:
✖ Uncommitted changes. Please commit before continuing.
Run `rlang::last_trace()` to see where the error occurred.

> usethis::use_version('minor')
✔ Adding "0.2.0" to Version.
✔ Adding new heading to NEWS.md.
ℹ There are 2 uncommitted files:
• DESCRIPTION
• NEWS.md
ℹ There are 2 uncommitted files:
Selection: 
> pprf.crab <- PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7, m = 200,
+                                 size.p =  .5,  PPmethod = 'LDA',  parallel =TRUE, cores = 2)
Error in PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7,  : 
  unused argument (class = "Type")

> 
> pprf.crab <- PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7, m = 200,
+                                 size.p =  .5,  PPmethod = 'LDA',  parallel =TRUE, cores = 2)
Error in PPforest::PPforest(data = crab, class = "Type", size.tr = 0.7,  : 
  unused argument (class = "Type")

> 
> pprf.crab <- PPforest::PPforest(data = crab, y = "Type", size.tr = 0.7, m = 200,
+                                 size.p =  .5,  PPmethod = 'LDA',  parallel =TRUE, cores = 2)
> pprf.crab

Call:
 PPforest::PPforest(data = crab, y = "Type", size.tr = 0.7, m = 200,      PPmethod = "LDA", size.p = 0.5, parallel = TRUE, cores = 2) 
               Type of random forest: Classification
                     Number of trees: 200
No. of variables tried at each split: 3

        OOB estimate of  error rate: 6.43%
Confusion matrix:
             BlueFemale BlueMale OrangeFemale OrangeMale class.error
BlueFemale           33        2            0          0        0.06
BlueMale              4       31            0          0        0.11
OrangeFemale          0        0           33          2        0.06
OrangeMale            0        1            0         34        0.03
> ```r 
Error: attempt to use zero-length variable name

> pprf.crab

Call:
 PPforest::PPforest(data = crab, y = "Type", size.tr = 0.7, m = 200,      PPmethod = "LDA", size.p = 0.5, parallel = TRUE, cores = 2) 
               Type of random forest: Classification
                     Number of trees: 200
No. of variables tried at each split: 3

        OOB estimate of  error rate: 6.43%
Confusion matrix:
             BlueFemale BlueMale OrangeFemale OrangeMale class.error
BlueFemale           33        2            0          0        0.06
BlueMale              4       31            0          0        0.11
OrangeFemale          0        0           33          2        0.06
OrangeMale            0        1            0         34        0.03
 

```

`PPforest` print a summary result from the model with the confusion matrix information and the oob-error rate in a similar way randomForest packages does.

This function returns the predicted values of the training data, training error, test error and predicted test values. Also there is the information about out of bag error for the forest and also for each tree in the forest. Bootstrap samples, output of all the trees in the forest from , proximity matrix and vote matrix, number of trees grown in the forest, number of predictor variables selected to use for splitting at each node. Confusion matrix of the prediction (based on OOb data), the training data and test data and vote matrix are also returned.

The printed version of a `PPforest` object follows the `randomForest` printed version to make them comparable. Based on confusion matrix, we can observe that the biggest error is for BlueMale class. Most of the wrong classified values are between BlueFemale and BlueMale.

```r
str(pprf.crab, max.level = 1)
List of 28
 $ predicting.training: Factor w/ 4 levels "BlueFemale","BlueMale",..: 2 2 2 2 1 2 2 1 2 1 ...
 $ training.error     : num 0.0571
 $ prediction.test    : Factor w/ 4 levels "BlueFemale","BlueMale",..: 2 2 2 2 2 2 2 2 2 2 ...
 $ error.test         : num 0.0667
 $ oob.error.forest   : num 0.0643
 $ oob.error.tree     : num [1:200, 1] 0.1961 0.0625 0.1429 0.0408 0.1923 ...
 $ boot.samp          :List of 200
 $ output.trees       :List of 200
 $ proximity          : num [1:140, 1:140] 0 0.79 0.735 0.8 0.36 0.825 0.71 0.265 0.46 0.36 ...
 $ votes              : num [1:140, 1:4] 0.243 0.288 0.183 0.235 0.835 ...
  ..- attr(*, "dimnames")=List of 2
 $ prediction.oob     : Factor w/ 4 levels "BlueFemale","BlueMale",..: 2 2 2 2 1 2 2 1 2 1 ...
 $ n.tree             : num 200
 $ n.var              : int 3
 $ type               : chr "Classification"
 $ confusion          : num [1:4, 1:5] 33 4 0 0 2 31 0 1 0 0 ...
  ..- attr(*, "dimnames")=List of 2
 $ call               : language PPforest::PPforest(data = crab, y = "Type", size.tr = 0.7, m = 200, PPmethod = "LDA", size.p = 0.5,      parallel| __truncated__
 $ train              :'data.frame':	140 obs. of  6 variables:
 $ test               :'data.frame':	60 obs. of  6 variables:
 $ vote.mat           : num [1:200, 1:140] 1 2 4 2 2 2 2 4 1 1 ...
  ..- attr(*, "dimnames")=List of 2
 $ vote.mat_cl        : chr [1:4] "BlueFemale" "BlueMale" "OrangeFemale" "OrangeMale"
 $ class.var          : chr "Type"
 $ oob.obs            : num [1:200, 1:140] 0 0 0 0 0 0 1 1 0 1 ...
 $ std                : chr "scale"
 $ dataux             : num [1:140, 1:5] -2.14 -1.71 -1.65 -1.36 -1.28 ...
  ..- attr(*, "dimnames")=List of 2
  ..- attr(*, "scaled:center")= Named num [1:5] 15.5 12.8 32 36.3 14
  .. ..- attr(*, "names")= chr [1:5] "FL" "RW" "CL" "CW" ...
  ..- attr(*, "scaled:scale")= Named num [1:5] 3.47 2.57 7.06 7.84 3.4
  .. ..- attr(*, "names")= chr [1:5] "FL" "RW" "CL" "CW" ...
 $ mincol             : NULL
 $ maxmincol          : NULL
 $ train_mean         : Named num [1:5] 15.5 12.8 32 36.3 14
  ..- attr(*, "names")= chr [1:5] "FL" "RW" "CL" "CW" ...
 $ train_sd           : Named num [1:5] 3.47 2.57 7.06 7.84 3.4
  ..- attr(*, "names")= chr [1:5] "FL" "RW" "CL" "CW" ...
 - attr(*, "class")= chr "PPforest"
```

The `PPforest` object can be used to predict new data using the `predict` function. The predicted values are returned as a factor with the class levels.

```r

pred.crab <- predict(pprf.crab, newdata = crab[1:10,-1 ])

pred.crab[[3]] 
 [1] BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale BlueFemale
 [9] BlueFemale BlueFemale
 
Levels: BlueFemale BlueMale OrangeFemale OrangeMale
```

