#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// 
// // [[Rcpp::export]]
// CharacterVector ppcl(  List tree, NumericMatrix x, Function f) {
//   CharacterVector res = f( tree,x);
//   return res;
// }


// RObject callFunction(  IntegerVector cl,NumericMatrix x,  Function f) {
//   RObject res = f(cl,x);
//   return res;
// }


double innerProduct(NumericVector x, NumericVector y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
}

// NumericVector PPpred( List tree, NumericMatrix testdata ) {
//   NumericMatrix TRstr = tree[1];
//   NumericMatrix TRprnode = tree[2];
//   NumericMatrix TRspl = tree[3];

// [[Rcpp::export]]
NumericVector PPpred( NumericMatrix TRstr, NumericMatrix TRprnode, NumericMatrix TRspl, NumericMatrix testdata) {
  int testn = testdata.nrow();
  NumericVector cl = testdata.nrow();
  
  //NumericVector terminal = TRstr( TRstr(,4) == 0, 1)
  for (int i = 0; i<testn; i++) {
    
    int nd = 0;
    double pp = 0.0;
    
    while ( TRstr(nd, 3) >0 ) {
      int p = (TRstr(nd , 3) - 1);
      
      pp = innerProduct(TRprnode(p, _), testdata(i,_));
      
      if( pp < TRspl(p, 0) ) {
        nd = (TRstr(nd, 1) - 1);
      } else {
        nd = (TRstr(nd, 2) -1);
      }
    }
    cl[i] = TRstr(nd, 2);
  }
  return cl;
}

 // [[Rcpp::export]]
NumericMatrix imposoon(NumericMatrix train,
                   NumericVector classes, List oobid, List permute, List trees , IntegerVector noob,
                     List TRstrL, List TRsplL, List TRprnodeL) {
  int ntree = trees.size();
  int p = train.ncol();
  NumericMatrix pp( ntree, p);
  
    //corr.oob.per <- matrix(0,ncol=ncol(ppf$train)-1,nrow=length(permute))
     for(int i = 0; i<p; i++) {
      for (int j = 0; j<ntree; j++) {
        //for (int k=0; k<p; k++){
        //indi = oobid[j];
        NumericMatrix d(noob[j],p);
        NumericVector id = oobid[j];
        NumericVector idp = permute[j];
        NumericVector pred( noob[j] );
        NumericVector cloob( noob[j]);
        
        for( int k = 0; k < noob[j]; k++) {
          cloob(k) = classes(id[k]);
          for (int l = 0; l < p; l++) {
            if (l == i)  {
              d(k, l) = train( idp[k], l);
            } else {
              d(k, l) = train( id[k] , l);
              }
            }
        }
        NumericMatrix TRstr = TRstrL[j];
        NumericMatrix TRspl = TRsplL[j];
        NumericMatrix TRprnode = TRprnodeL[j];
  
        pred =  PPpred(TRstr,TRprnode, TRspl, d);
        double M = 0.0;
        for(int m = 0; m<pred.size(); m++){
          if(pred[m] == cloob[m]) {
            M ++;
            }
        }
        pp(j,i) =  M;
      }
    }

  return pp;
  }

