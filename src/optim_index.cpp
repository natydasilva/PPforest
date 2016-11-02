#include <RcppArmadillo.h>  
// [[Rcpp::depends(RcppArmadillo)]] 
using namespace arma; 
using namespace Rcpp;


std::map<int, int> tableC(IntegerVector x) {
  std::map<int, int> counts;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    counts[x[i]]++;
  }
  
  return counts;
}  

// [[Rcpp::export]]
double LDAindex(IntegerVector origclass, NumericMatrix origdata, 
                NumericMatrix proj=NumericMatrix(0),bool weight=true){
  double index;
  int n=origdata.nrow(),p=origdata.ncol(),q=proj.ncol(),p1=proj.nrow();
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(origclass);
  int g=gn.size();  
  if(p1!=p)
    q=p;
  NumericVector allmean(q);
  NumericMatrix W(q,q),WB(q,q),gsum(q,g);        
  NumericMatrix projdata(n,q);
  if(p1!=p||p1==1){
    projdata=origdata; 
  } else{
    for(int i=0;i<n;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          projdata(i,j)+=origdata(i,k)*proj(k,j);
        }
      }
    }
  }
  for(int i=0;i<n;i++){
    for(int k=0;k<q;k++){
      allmean(k)+=projdata(i,k)/n;
      gsum(k,(origclass(i)-1))+=projdata(i,k);
    }
  }
  for(int i=0;i<n;i++){
    int l=origclass[i]-1;
    double gn1;
    if(weight){
      gn1=gn(l);
    } else{
      gn1=n/g; 
    }
    for(int j1=0;j1<q;j1++){
      for(int j2=0;j2<=j1;j2++){
        W(j1,j2)+=((projdata(i,j1)-gsum(j1,l)/gn(l))*
          (projdata(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
        W(j2,j1)=W(j1,j2);
        double temp=((projdata(i,j1)-gsum(j1,l)/gn(l))*
                     (projdata(i,j2)-gsum(j2,l)/gn(l))+
                     (gsum(j1,l)/gn(l)-allmean(j1))*
                     (gsum(j2,l)/gn(l)-allmean(j2)))/gn(l)*gn1;
        WB(j1,j2)+=temp;
        WB(j2,j1)=WB(j1,j2);
      }
    }
  }
  Function det=base["det"];
  index=1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
  return index;
}


   
// [[Rcpp::export]] 
  List LDAopt(IntegerVector origclass,arma::mat origdata,int q=1, 
              std::string PPmethod="LDA",bool weight=true){
    
  int n=origdata.n_rows, p=origdata.n_cols;
  
  // group totals, group means and overall mean
  std::map<int, int> ng = tableC(origclass); 
  int g = ng.size();
  
   arma::colvec mean_all(p);
   arma::mat  mean_g(g, p);
   
    for (int i=0; i < p; i++){
      mean_all[i] = mean( origdata.col(i) ) ;
      for (int k=0; k < g; k++) {
        NumericVector tot(g);
        for (int j=0; j<n; j++) {
          if (origclass[j] == (k+1) ) tot(k) += origdata(j,i) ;  
        }
        mean_g(k, i) = tot(k)/ng[k+1] ; 
      }
    }  
    
    // between SS matrix
    arma::mat  B(p, p, arma::fill::zeros);
    arma::mat  W(p, p, arma::fill::zeros);
   for (int k=0; k < g; k++) {
     arma::vec  temp1(p);
     arma::mat Temp1(p,p, arma::fill::zeros);
      double gn1 = 0.0;
      
        if(weight) {
          gn1 = ng[k+1];
        } else {
          gn1 = n/g;
        }
        temp1 = arma::trans(mean_g.row(k)) - mean_all;
        Temp1 = gn1 * temp1 * arma::trans(temp1); 
        for (int i=0; i < p; i++) {
          for (int j=0; j < p; j++) {
            B(i,j) += Temp1(i,j);
          }
        }
   }
   
   // Whithin SS matrix
   arma::mat rr(n, p, arma::fill::zeros); 
        for (int i=0; i < n; i++) {
          for (int k=0; k < g; k++) {    
            if ( origclass[i] == (k+1) ) {
              for (int j=0; j<p; j++)
              rr(i,j) = origdata(i,j) - mean_g(k,j); 
            }
          }
        }
   W = arma::trans(rr) * rr;
        
// eigen decomposition and compute projection
  W = arma::inv(W+B);
  B = W * B;
  cx_vec eigval;
  cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, B);
  int mm = eigval.index_max();
  
//   arma::vec temp_proj = real(eigvec.col(mm));
//   NumericMatrix proj(p,1);  
//   for (int i=0; i < p; i++)  proj(i,1) =  temp_proj[i];
   // PAra poder usar LDAindex, hay que transformar de rcpp a rcpp armadillo
   //proj = origdata * proj;
//   double index = LDAindex(origclass, origdata, proj, weight);
   // Rcpp::Named("indexbest")=index,

   return Rcpp::List::create( Rcpp::Named("projbest")= real(eigvec.col(mm)),
                              Rcpp::Named("projdata") = origdata*real(eigvec.col(mm)));                            
}

// [[Rcpp::export]]
double PDAindex(IntegerVector origclass, NumericMatrix origdata,
                NumericMatrix proj=NumericMatrix(0),bool weight=true,
                double lambda=0.1){
  double index;
  int n=origdata.nrow(),p=origdata.ncol(),q=proj.ncol(),p1=proj.nrow();
  Environment base("package:base");
  Function table=base["table"];
  NumericVector gn=table(origclass);
  int g=gn.size();  
  NumericMatrix W(p,p),WB(p,p),gsum(p,g);
  NumericVector allmean(p);
  if(p1!=p)
    q=p;
  for(int i=0;i<n;i++){
    for(int k=0;k<p;k++){
      allmean(k)+=origdata(i,k)/n;
      gsum(k,(origclass(i)-1))+=origdata(i,k);
    }
  }
  for(int i=0;i<n;i++){
    int l=origclass[i]-1;
    double gn1;
    if(weight){
      gn1=gn(l);
    } else{
      gn1=n/g; 
    }
    for(int j1=0;j1<p;j1++){
      for(int j2=0; j2<=j1; j2++){
        double temp1,temp2;
        if(j1!=j2){
          temp1=(1-lambda)*((origdata(i,j1)-gsum(j1,l)/gn(l))*
            (origdata(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
          
          temp2=(1-lambda)*((origdata(i,j1)-gsum(j1,l)/gn(l))*
            (origdata(i,j2)-gsum(j2,l)/gn(l)))+
            (gsum(j1,l)/gn(l)-allmean(j1))*
            (gsum(j2,l)/gn(l)-allmean(j2))/gn(l)*gn1;
        } else{
          temp1=((origdata(i,j1)-gsum(j1,l)/gn(l))*
            (origdata(i,j2)-gsum(j2,l)/gn(l)))/gn(l)*gn1;
          
          temp2=((origdata(i,j1)-gsum(j1,l)/gn(l))*
            (origdata(i,j2)-gsum(j2,l)/gn(l))+
            (gsum(j1,l)/gn(l)-allmean(j1))*
            (gsum(j2,l)/gn(l)-allmean(j2)))/gn(l)*gn1;              
        }  
        W(j1,j2)+=temp1;  
        WB(j1,j2)+=temp2;
        W(j2,j1)=W(j1,j2);            
        WB(j2,j1)=WB(j1,j2);
      }
    }      
  }
  NumericMatrix Wt(q,p),WBt(q,p);    
  NumericMatrix Wtt(q,q),WBtt(q,q); 
  if(p1!=p||p1==1){
    Wtt=W;
    WBtt=WB;
  } else{
    for(int i=0;i<p;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          Wt(j,i)+=W(k,i)*proj(k,j);
          WBt(j,i)+=WB(k,i)*proj(k,j);               
        }
      }
    }    
    for(int i=0;i<q;i++){
      for(int j=0;j<q;j++){
        for(int k=0;k<p;k++){
          Wtt(i,j)+=Wt(i,k)*proj(k,j);
          WBtt(i,j)+=WBt(i,k)*proj(k,j);               
        }
      }
    }      
  }   
  Function det=base["det"];
  index=1.0-as<double>(det(wrap(Wtt)))/as<double>(det(wrap(WBtt)));   
  return index;
}



// [[Rcpp::export]]
List PDAopt(IntegerVector origclass,arma::mat origdata,int q=1, 
            std::string PPmethod="PDA",bool weight=true, double lambda=0.1){
  
  int n=origdata.n_rows, p=origdata.n_cols;
  
  // group totals, group means and overall mean
  std::map<int, int> ng = tableC(origclass); 
  int g = ng.size();
  
  arma::colvec mean_all(p);
  arma::mat  mean_g(g, p);
  
  //data.std <- as.matrix(origdata)
  //class.table<-table(origclass)
  // g<-length(class.table)
  // class.name<-names(class.table)
  //p<-ncol(data.std)
  // n<-nrow(data.std)
  //mean.g<-matrix(apply(data.std,2,function(x) 
  //    tapply(x,origclass,mean,na.rm=TRUE)),ncol=p)
  //  mean.all<-matrix(apply(data.std,2,mean),ncol=p)           
  for (int i=0; i < p; i++){
    mean_all[i] = mean( origdata.col(i) ) ;
    for (int k=0; k < g; k++) {
      NumericVector tot(g);
      for (int j=0; j<n; j++) {
        if (origclass[j] == (k+1) ) tot(k) += origdata(j,i) ;  
      }
      mean_g(k, i) = tot(k)/ng[k+1] ; 
    }
  }  
  
  // between SS matrix
  arma::mat  B(p, p, arma::fill::zeros);
  arma::mat  W(p, p, arma::fill::zeros);
  
  //B<-matrix(0,ncol=p,nrow=p)
  //W<-matrix(0,ncol=p,nrow=p) 
  
  for (int k=0; k < g; k++) {
    arma::vec  temp1(p);
    arma::mat Temp1(p,p, arma::fill::zeros);
    double gn1 = 0.0;
    
    if(weight) {
      gn1 = ng[k+1];
    } else {
      gn1 = n/g;
    }
    temp1 = arma::trans(mean_g.row(k)) - mean_all;
    Temp1 = gn1 * temp1 * arma::trans(temp1); 
    for (int i=0; i < p; i++) {
      for (int j=0; j < p; j++) {
        B(i,j) += Temp1(i,j);
      }
    }
  }
  
  // Whithin SS matrix
  arma::mat rr(n, p, arma::fill::zeros); 
  for (int i=0; i < n; i++) {
    for (int k=0; k < g; k++) {    
      if ( origclass[i] == (k+1) ) {
        for (int j=0; j<p; j++)
          rr(i,j) = origdata(i,j) - mean_g(k,j); 
      }
    }
  }
  W = arma::trans(rr) * rr;
  
  //arma::mat  Wt(p, p, arma::fill::zeros);
  arma::mat Wt = (1-lambda)*W;
  Wt.diag() = W.diag() ;
  
  
  
  // eigen decomposition and compute projection
  Wt = arma::inv(Wt+B);
  B = Wt * B;
  cx_vec eigval;
  cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, B);
  int mm = eigval.index_max();
  
  return Rcpp::List::create( Rcpp::Named("projbest") = real(eigvec.col(mm)),
                             Rcpp::Named("projdata") = origdata*real(eigvec.col(mm)));                            
}

// ===================================== 
// 
// 
// // [[Rcpp::export]]  
// extern "C" SEXP fastLm(SEXP ys, SEXP Xs) {
//   
//   Rcpp::NumericVector yr(ys);                 // creates Rcpp vector from SEXP
//   Rcpp::NumericMatrix Xr(Xs);                 // creates Rcpp matrix from SEXP
//   int n = Xr.nrow(), k = Xr.ncol();
//   
//   arma::mat X(Xr.begin(), n, k, false);       // reuses memory and avoids extra copy
//   arma::colvec y(yr.begin(), yr.size(), false);
//   
//   arma::colvec coef = arma::solve(X, y);      // fit model y ~ X
//   arma::colvec resid = y - X*coef;            // residuals
//   
//   double sig2 = arma::as_scalar( arma::trans(resid)*resid/(n-k) );
//   // std.error of estimate
//   arma::colvec stderrest = arma::sqrt( sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
//   
//   return Rcpp::List::create(
//     Rcpp::Named("coefficients") = coef,
//     Rcpp::Named("stderr")       = stderrest
//   ) ;
//   
// }
// 
// // [[Rcpp::export]]
// arma::mat outerC2(arma::colvec a, arma::rowvec b) {
//   return a*b;
// }
// 
// // [[Rcpp::export]]
// arma::mat cube_sum(Rcpp::NumericVector vx) {
//   
//   Rcpp::IntegerVector x_dims = vx.attr("dim");
//   arma::cube x(vx.begin(), x_dims[0], x_dims[1], x_dims[2], false);
//   
//   arma::mat result(x.n_cols, x.n_slices);
//   for (unsigned int i = 0; i < x.n_slices; i++) {
//     result.col(i) = arma::conv_to<arma::colvec>::from(arma::sum(x.slice(i)));  
//   }
//   
//   return result;
// }
// 
// // [[Rcpp::export]]
// arma::mat cube_sum1(arma::cube c) {
//   arma::mat ss(c.n_rows, c.n_cols, arma::fill::zeros);
//   for(int i = 0; i < c.n_slices; i++) {
//     //ss += c.tube(arma::span::all, arma::span(i));
//     ss += c.slice(i);
//   }
//   return ss;
// }
// 
// 
// // [[Rcpp::export]]
// arma::mat prueba(arma::mat a, arma::vec b) {
//   arma::mat c(a.n_rows, a.n_cols) ;
//   int n = b.size();
//   
//   arma:cube d(a.n_rows, a.n_cols, n);
//   for(int i =0; i <n; i++) {
//     d.slice(i) =   b[i] * a;
//   }
//   c = cube_sum1(d);
//   return c;
// }

