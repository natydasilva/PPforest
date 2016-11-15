#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]] 
using namespace arma; 
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec tableC(arma::vec x) {
  arma::vec values = arma::unique(x);
    arma::vec counts(values.size());
    for(int i = 0; i < values.size(); i++) {
      for(int j=0; j<x.size();j++){
        if(x(j)==values(i)){
          counts(i)++;
            }
      }

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
double signC(double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return -1;
  }
}
   
// [[Rcpp::export]] 
  arma::vec LDAopt(arma::vec origclass,arma::mat origdata,int q=1, 
              std::string PPmethod="LDA",bool weight=true){
    
  int n=origdata.n_rows, p=origdata.n_cols;
  
  // group totals, group means and overall mean
  arma::vec ng = tableC(origclass); 
  int g = ng.size();
  arma::vec clval = arma::unique(origclass);
  
  
   arma::colvec mean_all(p);
   arma::mat  mean_g(g, p);
   
    for (int i=0; i < p; i++){
      mean_all[i] = mean( origdata.col(i) ) ;
      for (int k=0; k < g; k++) {
        double tot = 0.0;
        for (int j=0; j<n; j++) {
          if (origclass(j) == clval(k) ) tot += origdata(j,i) ;  
        }
        mean_g(k, i) = tot/ng(k) ; 
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
          gn1 = ng(k);
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
            if ( origclass(i) == clval(k) ) {
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
  
   arma::vec bestproj = real(eigvec.col(mm));
//   NumericMatrix proj(p,1);  
//   for (int i=0; i < p; i++)  proj(i,1) =  temp_proj[i];
   // PAra poder usar LDAindex, hay que transformar de rcpp a rcpp armadillo
   //proj = origdata * proj;
//   double index = LDAindex(origclass, origdata, proj, weight);
   // Rcpp::Named("indexbest")=index,

   // return Rcpp::List::create( Rcpp::Named("projbest")= temp_proj,
   //                            Rcpp::Named("projdata") = origdata*temp_proj
   //                           );    
   return bestproj;
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
arma::vec PDAopt(arma::vec origclass,arma::mat origdata,int q=1, 
            std::string PPmethod="PDA",bool weight=true, double lambda=0.1){
  
  int n=origdata.n_rows, p=origdata.n_cols;
  
  // group totals, group means and overall mean
  arma::vec ng = tableC(origclass); 
  int g = ng.size();
  arma::vec clval = arma::unique(origclass);
  arma::colvec mean_all(p);
  arma::mat  mean_g(g, p);
  
          
  for (int i=0; i < p; i++){
    mean_all[i] = mean( origdata.col(i) ) ;
    for (int k=0; k < g; k++) {
      double tot = 0.0;
      for (int j=0; j<n; j++) {
        if (origclass(j) == clval(k) ) tot += origdata(j,i) ;  
      }
      mean_g(k, i) = tot/ng(k) ; 
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
      gn1 = ng(k);
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
      if ( origclass(i) == clval(k) ) {
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
  
  arma::vec bestproj = real(eigvec.col(mm));
  
  // return Rcpp::List::create( Rcpp::Named("projbest") = temp_proj,
  //                            Rcpp::Named("projdata") = origdata*temp_proj);    
  
  return bestproj;
}



// [[Rcpp::export]]
arma::uvec varselect(int p, int s){//p number of variables, s number of selected variables
  arma::uvec id = arma::linspace<uvec>(0, p-1,p); //integer sequence from 0 to p-1
  arma:: uvec idx = arma::linspace<uvec>(0,s-1,s);
  arma::uvec idsuff = arma::shuffle<uvec>(id);// suffle the indexes
  arma::uvec ids = idsuff.elem(idx);//select only s indexes
  return arma::sort(ids);//sort the selected variables
}

//Subset of random selected variables from origdata used in each node partition,
//remove variables with 0 variance
// [[Rcpp::export]]
List datanode(arma::mat origdata, double sizep){
  //remove the variable with zero variance
  arma::mat sdcol = arma::stddev(origdata); //sd by column
  arma::uvec idx = find(sdcol.row(0) > 0);//find indices of non-zero elements, or elements satisfying a relational condition
  arma::mat redudata = origdata.cols(idx);//select the columns with 0 variance
  int p = redudata.n_cols;
  int sp = std::round(sizep*p);
  
  arma::uvec vrnd = varselect(p, sp);//variable selection for each node partition in Rcpp
  arma::mat datanode = redudata.cols(vrnd);//data with the selected variables for each node partition
  
  return Rcpp::List::create( Rcpp::Named("data") = datanode,
                             Rcpp::Named("varselected") = vrnd);
}

//--------First split, redefine the problem in a two class problem,this is used inside findproj
// [[Rcpp::export]]
arma::vec split_rel(arma::vec origclass,arma::mat origdata, arma::colvec  projdata){
//origdata here are after variable selection (varselect and datenode functions)  
  
int n = origdata.n_rows;
  
// group totals, group means and overall mean
  arma::vec ng = tableC(origclass); 
  int g = ng.size();
  arma::vec clval = arma::unique(origclass);
  arma::colvec mean_g(g);
  arma::vec newclass(n);
  
  if (g==2) {
    //IntegerVector class_rel = origclass;
    newclass = origclass;
    
  } else { 
    for (int k=0; k < g; k++) {
      double tot=0.0;
      for (int j=0; j<n; j++) {
        if (origclass(j) == clval(k) ) tot += projdata(j);  
      }
      mean_g(k) = tot/ng(k) ; 
    }

  arma::uvec mlist = sort_index(mean_g);
  arma::vec sm = sort(mean_g);
  arma::vec msort = diff(sm);
  
  int mm = msort.index_max();
  double pm = (sm(mm)+sm((mm+1)))/2.0;
  //arma::vec newclass(n);
  
  for (int k=0; k < g; k++) {
    for(int i=0; i<n;i++){
      if (origclass(i) == clval(k) ) {
        if((mean_g(k)-pm)>0){
          newclass(i)=2;
        }else{
          newclass(i)=1;
        }
      }
    }
  }
  }
  
  return newclass;
  
}


//Finds the 1D projection for separating all classes
// [[Rcpp::export]]
List findproj(arma::vec origclass,
              arma::mat origdata, std::string PPmethod="LDA", 
              double lambda=0.1){
  
  arma::vec a1(origdata.n_cols) , a2(origdata.n_cols), a(origdata.n_cols);
  arma::vec ng = tableC(origclass); 
  int g = ng.size();
  arma::vec projdata(origclass.size());
 // arma::vec a( a1.size() );
 // double indexbest=0.0;
 
  if (PPmethod =="LDA"){
    a1 = LDAopt(origclass,origdata,1, "LDA",true);
   //double indexbest = LDAindex(origclass,origdata);
   } else {
    a1 = PDAopt(origclass, origdata, 1,"PDA",true, lambda);
  //double indexbest = PDAindex(origclass,origdata)
  }

  int  index = arma::index_max(arma::abs(a1));
  double sign = signC(a1(index));
  arma::vec classe = split_rel(origclass, origdata,  origdata*a1); 

  if (g > 2) {
    
    if (PPmethod =="LDA"){
      a2 = LDAopt(classe,origdata,1, "LDA",true);
      } else {
      a2 = PDAopt(classe, origdata, 1,"PDA",true, lambda);
    }
      
    double sign2 = signC(a2(index));
      if (sign != sign2) {
        a  = -1*a2;
      } else {
        a = a2;
      }
  } else {
    a = a1;
  }
 
 
  projdata = origdata*a;
  return Rcpp::List::create(Rcpp::Named("projdata") = projdata,
                           Rcpp::Named("projbest") = a,
                           Rcpp::Named("class")= classe);
}



// [[Rcpp::export]]
List findprojPDA(arma::vec origclass,
              arma::mat origdata, 
              double lambda=0.1){
  
  arma::vec a = PDAopt(origclass, origdata, lambda);
  arma::vec classe = split_rel(origclass, origdata,  origdata*a );
  
  arma::vec ng = tableC(origclass);
  int g = ng.size();
   
   int  index = arma::index_max(arma::abs(a));
   double sign = signC(a(index));

  arma::vec  out;
  
  if (g > 2) {
    arma::vec a2 = PDAopt(classe, origdata, lambda);
    double sign2 = signC(a2(index));
    if (sign != sign2) {
       out  = (-1)*a2;
    } else {
      out = a2;
    }
  } else {
    out = a;
  }

    arma::vec projdata = origdata*a;
    return Rcpp::List::create(Rcpp::Named("projdata") = projdata,
                              Rcpp::Named("projbest") = out,
                              Rcpp::Named("class")= classe);
}
  


  // [[Rcpp::export]]
  List findprojLDA(arma::vec origclass,
                   arma::mat origdata){
    
    arma::vec a = LDAopt(origclass, origdata);
    arma::vec classe = split_rel(origclass, origdata,  origdata*a );
    
    arma::vec ng = tableC(origclass);
    int g = ng.size();
    
    int  index = arma::index_max(arma::abs(a));
    double sign = signC(a(index));
    
    arma::vec  out;
    
    if (g > 2) {
      arma::vec a2 = LDAopt(classe, origdata);
      double sign2 = signC(a2(index));
      if (sign != sign2) {
        out  = (-1)*a2;
      } else {
        out = a2;
      }
    } else {
      out = a;
    }
    
    arma::vec projdata = origdata*a;
    return Rcpp::List::create(Rcpp::Named("projdata") = projdata,
                              Rcpp::Named("projbest") = out,
                              Rcpp::Named("class")= classe);
  }




//   n <- nrow(origdata)
//     p <- ncol(origdata)
//     g <- table(origclass)
//     g.name <- names(g)
//     G <- length(g) 
//     
//     proj.data <-  a[[2]] #projected data
// # sign <- sign(a$projbest[abs(a$projbest) == max(abs(a$projbest))])
// # index <- (1:p) * (abs(a$projbest) == max(abs(a$projbest)))
// # index <- index[index > 0] #index of the biggest coefficient
//     sign <-signC(max(abs(a$projbest)))
//     index <- which.max(abs(a$projbest))
//     
//     
//     if (G == 2) { # only two classes not to relabel the classes in to 
//       class <- split_rel(origclass ,origdata ,proj.data )  
//     } else { # need to transform the problem in a two class problem
// # m <- tapply(c(proj.data), origclass, mean) # by class mean
//       sd <- tapply(c(proj.data), origclass, sd) # by class sd
//       sd.sort <- sort.list(sd)
//       class <- split_rel(origclass ,origdata ,proj.data )  
// #redefine the problem in a two class problem if there are more than 2 classes
//       
// # redo the 1D projection
//       
//       if (PPmethod == "LDA") {
//         a <- LDAopt(class,origdata, q=1, 
//                     PPmethod="LDA",weight=TRUE)
//         indexbest <- LDAindex(class,origdata, 
//                               proj=as.matrix(a[[2]]), weight=TRUE)
//       } else if (PPmethod == "PDA") {
//         a <- PDAopt(class, origdata, weight=TRUE, q = 1, lambda = lambda)
//         indexbest <- PDAindex(class,origdata, 
//                               proj=as.matrix(a[[2]]), weight=TRUE,lambda = lambda)
//       } 
//       if (sign != signC(a$projbest[index])) 
//         a$projbest <- -a$projbest
//         proj.data <- as.matrix(origdata) %*% a$projbest
//     }
//   

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


