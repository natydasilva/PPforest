#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]] 
using namespace arma; 
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec tableC(arma::vec x) {
  arma::vec values = arma::unique(x);
  arma::vec counts(values.size(), fill::zeros);
  
  for(int i = 0; i < values.size(); i++) {
    for(int j = 0; j < x.size(); j++){
      if(x(j) == values(i)){
        counts(i)++;
      }
    }
    
  }
  
  return counts;
}  


// [[Rcpp::export]]
double roundme(double x)
{
  if (x < 0.0)
    return (int)(x - 0.5);
  else
    return (int)(x + 0.5);
}

// [[Rcpp::export]]
arma::vec LDAindex2(arma::vec origclass, arma::mat origdata, 
                 arma::mat proj, bool weight=true){
  double index=0.0;
  int n = origdata.n_rows, p = origdata.n_cols; 
  int q = proj.n_cols, p1 = proj.n_rows;

  
   
   arma::vec clval = arma::unique(origclass);
   //----
   arma::vec newclass(n, fill::zeros);
   for(int i = 0; i < n; i++){
   for (int k = 0; k < clval.size(); k++) {
       if (origclass(i) == clval(k) ){
         newclass(i) = k + 1;
       }
     }

   }
   origclass = newclass;
   arma::vec gn = tableC(origclass);
   int g = gn.size(); 

   //----
   
   
   if(p1!= p){
     q = p;
     }
    arma::vec allmean(q, fill::zeros);
   arma::mat W(q, q, fill::zeros), WB(q, q, fill::zeros ),gsum(q, g, fill::zeros );        
   arma::mat projdata(n, q, fill::zeros);
   
  if(p1!= p||p1 == 1){
    projdata = origdata;
  }else{
    for(int i = 0; i < n; i++ ){
      for(int j = 0; j < q; j++ ){
        for(int k = 0; k < p ; k++ ){
          projdata(i, j )+= origdata( i, k )*proj(k, j );
        }
      }
    }
  }
  for(int i = 0; i < n; i++){
    for(int k = 0; k < q; k++){
      allmean(k) += projdata(i, k) / n;
      gsum(k, (origclass(i) - 1) ) += projdata(i, k);
    }
  }
  double temp = 0.0;
  for(int i = 0; i < n; i++){
    int l = origclass(i) - 1;
    double gn1;
    if(weight){
      gn1 = gn(l);
    } else{
      gn1 = n / g;
    }
    
    for(int j1 = 0; j1 < q; j1++){
      for(int j2 = 0; j2 <= j1; j2++){
        W(j1, j2) += ( (projdata(i, j1) - gsum(j1, l) / gn(l) )*
          (projdata(i, j2) - gsum(j2, l) / gn(l) ) ) / gn(l)*gn1;
        W(j2, j1) = W(j1, j2);
        temp = ( (projdata(i, j1) - gsum(j1, l) / gn(l) )*
                     (projdata(i,j2)-gsum(j2,l)/gn(l) ) +
                     (gsum(j1, l) / gn(l) - allmean(j1) )*
                     (gsum(j2, l) / gn(l) - allmean(j2))) / gn(l)*gn1;
        WB(j1, j2)+= temp;
        WB(j2, j1)= WB(j1, j2);
      }
    }
  }

  index = 1.0 - det(W) / det(WB);
  
   arma::vec out(q, fill::zeros);
   out(0)  = index; 
    return out;
}


// [[Rcpp::export]] 
double signC(double x) {
  if (x > 0) {
    return 1;
  } else if (x == 0) {
    return 0;
  } else {
    return - 1;
  }
}
   
// [[Rcpp::export]] 
  arma::vec LDAopt(arma::vec origclass, arma::mat origdata, int q = 1, 
              std::string PPmethod = "LDA", bool weight = true){
    
  int n = origdata.n_rows, p = origdata.n_cols;
  
  // group totals, group means and overall mean
  arma::vec clval = arma::unique(origclass);
  arma::vec ng(clval.size(), fill::zeros); 
  ng = tableC(origclass); 

  int g = ng.size();
   arma::colvec mean_all(p);
   arma::mat  mean_g(g, p);
   
    for (int i = 0; i < p; i++){
      mean_all[i] = mean( origdata.col(i) ) ;// variable mean
      for (int k = 0; k < g; k++) { // in groups
        double tot = 0.0;
        for (int j = 0; j < n; j++) {//in observations
          if (origclass(j) == clval(k) ) tot += origdata(j,i) ; //total by variable and group 
        }
        mean_g(k, i) = tot / ng(k) ; //mean by variable and group
      }
    }  
    
    // between SS matrix
    arma::mat  B(p, p, arma::fill::zeros);
    arma::mat  W(p, p, arma::fill::zeros);
    double gdbl = ng.size();
    
   for (int k = 0; k < g; k++) {
     arma::vec  temp1(p);
     arma::mat Temp1(p, p, arma::fill::zeros);
      double gn1 = 0.0;
      
        if(weight) {
          gn1 = ng(k);
        } else {
          gn1 = n/gdbl;
        }
        temp1 = arma::trans(mean_g.row(k)) - mean_all;
        Temp1 = gn1 * temp1 * arma::trans(temp1); 
        //Rcout << gn1;
        for (int i = 0; i < p; i++) {
          for (int j = 0; j < p; j++) {
            B(i,j) += Temp1(i, j);
          }
        }
   }
   
   // Whithin SS matrix
   arma::mat rr(n, p, arma::fill::zeros); 
        for (int i = 0; i < n; i++) {
          for (int k = 0; k < g; k++) {    
            if ( origclass(i) == clval(k) ) {
              for (int j = 0; j < p; j++)
              rr(i, j) = origdata(i, j) - mean_g(k, j); 
            }
          }
        }
   W = arma::trans(rr) * rr;
        
// eigen decomposition and compute projection
  W = arma::inv(W + B);
  B = W * B;
  cx_vec eigval;
  cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, B);
  int mm = eigval.index_max();
  
   arma::vec bestproj = real(eigvec.col(mm));
 
   return bestproj;
}



// [[Rcpp::export]]
double PDAindex2(arma::vec origclass, arma::mat origdata,
                arma::mat proj,bool weight = true,
                double lambda = 0.1){

  double index = 0.0;
  int n = origdata.n_rows, p = origdata.n_cols; 
  
  int q = proj.n_cols, p1 = proj.n_rows;
  
  arma::vec gn = tableC(origclass);
  int g = gn.size();  
  
 arma::vec clval = arma::unique(origclass);
 //----
 arma::vec newclass(n, fill::zeros);
 for(int i = 0; i < n; i++){
   for (int k = 0; k < clval.size(); k++) {
     if (origclass(i) == clval(k) ){
       newclass(i) = k + 1;
     }
   }
   
 }
 origclass = newclass;
 
  arma::mat W(p, p, fill::zeros), WB(p, p, fill::zeros), gsum(p, g, fill::zeros);
  arma::vec allmean(p, fill::zeros);
  
  if(p1!= p)
    q = p;
  for(int i = 0; i < n; i++){
    for(int k = 0; k < p; k++){
      allmean(k)+= origdata(i, k) / n;
      gsum(k, (origclass(i) - 1 ) )+= origdata(i, k);
    }
  }
  for(int i=0;i<n;i++){
    int l=origclass(i)-1;
    double gn1;
    if(weight){
      gn1=gn(l);
    } else{
      gn1=n/g; 
    }
    for(int j1=0;j1<p;j1++){
      for(int j2=0; j2<=j1; j2++){
        double temp1 = 0.0; 
        double temp2 = 0.0;
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
  arma::mat Wt(q,p,fill::zeros),WBt(q,p,fill::zeros);    
  arma::mat Wtt(q,q,fill::zeros),WBtt(q,q,fill::zeros); 
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
  
  index = 1.0-det(Wtt)/det(WBtt);
  return index;
}




// [[Rcpp::export]]
arma::vec PDAopt(arma::vec origclass,arma::mat origdata,int q=1, 
            std::string PPmethod="PDA",bool weight=true, double lambda=0.1){
  
  int n=origdata.n_rows, p=origdata.n_cols;
  
  // group totals, group means and overall mean
  arma::vec clval = arma::unique(origclass);
  arma::vec ng(clval.size(), fill::zeros);
  ng = tableC(origclass); 
  int g = ng.size();
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
  
  //B=matrix(0,ncol=p,nrow=p)
  //W=matrix(0,ncol=p,nrow=p) 
  
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
  //int sp = std::round<int>(sizep*p);
  int sp = roundme(sizep*p);
  
  arma::uvec vrnd = varselect(p, sp);//variable selection for each node partition in Rcpp
  arma::mat datanode = redudata.cols(vrnd);//data with the selected variables for each node partition
  
  return Rcpp::List::create( Rcpp::Named("data") = datanode,
                             Rcpp::Named("varselected") = vrnd);
}

//--------First split, redefine the problem in a two class problem,this is used inside findproj
// [[Rcpp::export]]
arma::vec split_rel(arma::vec origclass, arma::mat origdata, arma::colvec  projdata){
//origdata here are after variable selection (varselect and datenode functions)  
  
int n = origdata.n_rows;
  
// group totals, group means and overall mean
  arma::vec ng = tableC(origclass); 
  int g = ng.size();
  arma::vec clval = arma::unique(origclass);
  arma::colvec mean_g(g);
  arma::vec newclass(n, fill::zeros);
  
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
              arma::mat origdata, std::string PPmethod, 
              double lambda=0.1){
  
  arma::vec a1(origdata.n_cols, fill::zeros) , a2(origdata.n_cols, fill::zeros), a(origdata.n_cols, fill::zeros);
  arma::vec ng = tableC(origclass); 
  int g = ng.size();
  arma::vec projdata(origclass.size());
 
 
  if(PPmethod.compare("LDA")==0){
    //Rcout << "LDA\n";
    a1 = LDAopt(origclass, origdata,  1, "LDA", true);
    //Rcout << a1;
   }else{
    //Rcout << "PDA\n";
     a1 = PDAopt(origclass, origdata, 1,"PDA", true, lambda);
     }

  
  int  index = arma::index_max(arma::abs(a1));
  double sign = signC(a1(index));
  arma::vec classe = split_rel(origclass, origdata,  origdata*a1); 

  if (g > 2) {
        if (PPmethod.compare("LDA")==0){
          //Rcout << "still LDA\n";
          
      a2 = LDAopt(classe,origdata,1, "LDA",true);
          //Rcout<< a2;
      } else {
        //Rcout << "still PDA\n";
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
                           Rcpp::Named("projbest") = a);
}


//subseting a vector based on a condition == val
// [[Rcpp::export]]
arma::uvec arma_sub_cond(arma::vec x, int val) {
  arma::uvec ids = find(x == val); // Find indices
  return ids;
}

//compute quantiles using R quantile function
// [[Rcpp::export]]
double quantileCpp(arma::vec x, double probs) {
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  double ans = 0.0;
  ans = as<double>(quantile(x, probs));

  return ans;
}

//compute quantiles, check translate to armadillo
// [[Rcpp::export]] 
NumericVector quant(NumericVector x, NumericVector q) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000000001)];
}
//Compute boundaries to use for future predictive observations. Use findproj output
// [[Rcpp::export]] 
arma::vec nodestr(arma::vec classe, arma::vec projdata){
int n = projdata.n_rows; 
 
 // group totals, group means and overall mean
 arma::vec clval = arma::unique(classe);
 arma::vec ng(clval.size(), fill::zeros);
 ng = tableC(classe); 
 int g = ng.size();
 //Rcout << ng; 
 //Rcout << '\n';

 arma::colvec mean_g(g);
 arma::colvec ss_g(g);
 arma::colvec sd(g);
 arma::colvec median_g(g);
 arma::colvec IQR_g(g);
   for (int k=0; k < g; k++) {
     double tot = 0.0, tot2 = 0.0;
     arma::uvec ids = arma_sub_cond(classe,clval(k));
      median_g(k) = median(projdata(ids));
      if(n>1){
      IQR_g(k) = quantileCpp(projdata(ids),0.75) - quantileCpp(projdata(ids),0.25);
      }else{
        IQR_g(k) = 0;
      }
     for (int j = 0; j<n; j++) {
       if (classe(j) == clval(k) ) {
         tot += projdata(j) ;
         tot2 += projdata(j)*projdata(j);
         
       }
     }
     
     //Rcout << ng(k);
      mean_g(k) = tot/ng(k) ; 
    
     if(n > 1){
     ss_g(k) =  tot2 -ng(k)*mean_g(k)*mean_g(k);
     sd(k) = sqrt(ss_g(k)/(ng(k)-1));// no prob but wr
     }else{
       sd(k) = 0;
     }
   }
   
    
    arma::uvec mlist = sort_index(mean_g);
    arma::colvec mLR = sort(mean_g);
    arma::vec sdLR = sd(mlist);
    arma::vec medianLR = median_g(mlist);
    arma::vec IQRLR = IQR_g(mlist);
   arma::vec nLR = tableC(classe)(mlist);
   arma::vec C(8);
   // // int mm = msort.index_max();
   C(0) = ( mLR(0) + mLR(1) )/2;
   C(1) = (mLR(0)*nLR(1) + mLR(1)*nLR(0))/sum(nLR);
 
   
   if(sum(sdLR == 0) != 0){
     C(2) =  C(0);
   }else{
     C(2) = (mLR(0)*sdLR(1) + mLR(1)*sdLR(0))/sum(sdLR);
   }
   
   if(sum(sdLR == 0) != 0){
     C(3) = C(1);
   }else{
     C(3) = (mLR(0)* sdLR(1)/sqrt(nLR(1)) + mLR(1)*sdLR(0)/sqrt(nLR(0))) /(sdLR(0)/sqrt(nLR(0)) +
        sdLR(1)/sqrt(nLR(1)));
   }
    
    C(4) = (medianLR(0) + medianLR(1))/2;
   C(5) = (medianLR(0)*nLR(1) + medianLR(1)*nLR(0))/sum(nLR);
   if(sum(IQRLR == 0) != 0){
     C(6) = C(4);
   }else{
     C(6) = (medianLR(0)*IQRLR(1) + medianLR(1)*IQRLR(0))/sum(IQRLR); 
   }
    
    if(sum(IQRLR == 0) != 0){
      C(7) = C(5);
    }else{
      C(7) = (medianLR(0)*(IQRLR(1)/sqrt(nLR(1))) + medianLR(1)*(IQRLR(0)/sqrt(nLR(0))))/((IQRLR(0)/sqrt(nLR(0))) + 
          (IQRLR(1)/sqrt(nLR(1))));
    
    }
   

   return C;
    
 }



//[[Rcpp::export]]
List findprojwrap(arma::vec origclass,arma::mat origdata, std::string PPmethod,
                  double sizep=1, double lambda =.1){

  int pp = origdata.n_cols;
  List dataspl = datanode(origdata, sizep );
  origdata = as<mat>(dataspl["data"]);
  arma::uvec vrnd = as<uvec>(dataspl["varselected"]);

  List oneDproj = findproj(origclass, origdata, PPmethod, lambda);
  arma::vec projdata = as<vec>(oneDproj["projdata"]);
  
  arma::vec classe = split_rel(origclass, origdata, projdata);
 
  arma::mat projbest = as<mat>(oneDproj["projbest"]);

       arma::vec C = nodestr(classe, projdata);

       arma::mat Alpha = zeros<mat>(pp);
      
       
       for(int i=0; i<vrnd.size(); i++){
            int v = vrnd(i);
           Alpha(v) = projbest(i,0);

       }

  arma::vec indexbest(1, fill::zeros);
if(PPmethod.compare("LDA")==0){
 //if(PPmethod == "LDA"){
 indexbest = LDAindex2(classe,origdata,projbest);
  }else{

  //if(PPmethod == "PDA"){
    indexbest = PDAindex2(classe,origdata,projbest);
  }

       int n = projdata.n_rows;

         // group totals, group means and overall mean
         arma::vec clval = arma::unique(classe);
         arma::vec ng(clval.size(), fill::zeros);
         ng = tableC(classe);
         int g = ng.size();
         arma::colvec mean_g(g,fill::zeros);
         
         for (int k=0; k < g; k++) {
           double tot = 0.0;
           for (int j = 0; j<n; j++) {
             if (classe(j) == clval(k) ) {
               tot += projdata(j) ;

               }
           }
          
          mean_g(k) = tot/ng(k);
        
          
          
         }
         
           arma::colvec mLR = sort(mean_g);
           arma::uvec sortLR  = sort_index(mean_g);
           
         arma::vec IOindexL(classe.size(), fill::zeros);
         arma::vec IOindexR(classe.size(),fill::zeros);

        
         for(int i=0; i<classe.size(); i++){
           if(classe(i)==(clval(sortLR(0)))){
           IOindexL(i) = true;
           }else{
             IOindexL(i) = false;
           }
           if(classe(i)==(clval(sortLR(1)))){
           IOindexR(i) = true;
           }else{
             IOindexR(i) = false;
            }
           }





return Rcpp::List::create(Rcpp::Named("Index") = indexbest,Rcpp::Named("Alpha") = Alpha, Rcpp::Named("C") = C,
                          Rcpp::Named("IOindexL")=IOindexL, Rcpp::Named("IOindexR") = IOindexR,Rcpp::Named("classe") = classe,
                                      Rcpp::Named("projdata") = projdata );

}

//tree structure
// [[Rcpp::export]]
List treeconstruct(arma::vec origclass, arma::mat origdata, arma::mat Treestruct, int id, int rep, int rep1, int rep2, arma::mat projbestnode, arma::mat  splitCutoffnode,
                   std::string PPmethod, double lambda = 0.1, double sizep = 1) {
  
  int n = origdata.n_rows;
  arma::vec cl2 = unique(origclass);
  arma::vec g(cl2.size(), fill::zeros); 
  g = tableC(origclass);

  int G = g.size();
  
  List a;
  List b;
  
  arma::vec C(8, fill::zeros);
  arma::vec classe(n, fill::zeros);


  if(G == 1){
    //check as.numeric group names
  
   Treestruct(id, 2) = cl2(0);

  return Rcpp::List::create( Rcpp::Named("Treestruct") = Treestruct,  Rcpp::Named("projbestnode") = projbestnode,
                               Rcpp::Named("splitCutoffnode")=splitCutoffnode, Rcpp::Named("rep")=rep,
                               Rcpp::Named("rep1") = rep1, Rcpp::Named("rep2") = rep2);
  }else{
    Treestruct(id, 1) = rep1;
    rep1 = rep1 + 1;
    Treestruct(id, 2) = rep1;
    rep1 = rep1 + 1;
    Treestruct(id, 3) = rep2;
    rep2 = rep2 + 1;
    
    // Rcout<< Treestruct;
    a = findprojwrap(origclass, origdata, PPmethod, sizep, lambda);
    classe = as<vec>(a["classe"]);
    C = nodestr(classe,as<vec>(a["projdata"]));
   
    splitCutoffnode.insert_rows( splitCutoffnode.n_rows,C.t());

  
  
  }
      Treestruct(id, 4) = as<double>(a["Index"]);

      arma::mat Alpha = as<mat>(a["Alpha"]);
      
      projbestnode.insert_rows(projbestnode.n_rows, Alpha.t());
       
        arma::vec tclass = origclass;
        arma::mat tdata = origdata;
        arma::vec IOindexL = as<vec>(a["IOindexL"]);
        tclass = tclass%IOindexL;

      arma::uvec tnaux= find(tclass>0);
      arma::uvec tindex = sort(tnaux);
         
      tclass = tclass(tindex);
      tdata = tdata.rows(tindex);
        
       b = treeconstruct(tclass, tdata,  Treestruct, Treestruct(id, 1)-1, rep,
                        rep1, rep2,  projbestnode,
                         splitCutoffnode, PPmethod,  lambda, sizep);


       Treestruct = as<mat>(b["Treestruct"]);

        projbestnode = as<mat>(b["projbestnode"]);
        splitCutoffnode = as<mat>(b["splitCutoffnode"]);
        rep = as<int>(b["rep"]);
        rep1 =  as<int>(b["rep1"]);
        rep2 =  as<int>(b["rep2"]);

        tclass = origclass;
        tdata =  origdata;
        arma::vec IOindexR = as<vec>(a["IOindexR"]);
        tclass = tclass%IOindexR;

        tnaux= find(tclass>0);
        tindex = sort(tnaux);
        tclass = tclass(tindex);
        tdata = tdata.rows(tindex);
          
        n =  tdata.n_rows;
        g.zeros();
        g = tableC(tclass);
        G = g.size();
        
        b = treeconstruct(tclass, tdata, Treestruct,
                            Treestruct(id, 2)-1, rep, rep1, rep2, projbestnode,
                            splitCutoffnode, PPmethod, lambda,sizep);

        Treestruct = as<mat>(b["Treestruct"]);
        projbestnode = as<mat>(b["projbestnode"]);
        splitCutoffnode = as<mat>(b["splitCutoffnode"]);
        rep = as<int>(b["rep"]);
        rep1 =  as<int>(b["rep1"]);
        rep2 =  as<int>(b["rep2"]);
        
       return b; 
        
  }



// [[Rcpp::export]]
arma::vec csample_num( arma::vec x,
                           int size,
                           bool replace,
                           arma::vec prob) {
arma::vec ret = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

// [[Rcpp::export]]
arma::vec boot( arma::mat origclass, arma::mat origdata) {
 //

   int n=origdata.n_rows;
   arma::vec id = arma::linspace<vec>(0, n-1,n); //integer sequence from 0 to n-1
    arma::vec clval = arma::unique(origclass.col(0));
    origclass.insert_cols(origclass.n_cols,id);
    
    arma::vec bootsel;
       for (int k = 0; k < clval.size(); k++) {
         arma::uvec idx = find(origclass.col(0) == clval(k));//find indices of non-zero elements, or elements satisfying a relational condition
         arma::mat reduclass = origclass.rows(idx);
         int nc = reduclass.n_rows;
           arma::vec samp = csample_num(reduclass.col(1), reduclass.n_rows, true,  ones<vec>(nc));
       bootsel = join_cols(bootsel,samp);
       }
  return bootsel;
}


//stratified training sample
// [[Rcpp::export]]
arma::vec trainfn( arma::mat origclass, arma::mat origdata, double sizetr){

  int n = origdata.n_rows;
  arma::vec id = arma::linspace<vec>( 0, n-1, n ); //integer sequence from 0 to n-1
  arma::vec clval = arma::unique(origclass.col(0));
  origclass.insert_cols(origclass.n_cols, id);
  
  arma::vec trainsel;
  for (int k = 0; k < clval.size(); k++) {
    arma::uvec idx = find(origclass.col(0) == clval(k));//find indices of non-zero elements, or elements satisfying a relational condition
    arma::mat reduclass = origclass.rows(idx);
    int nc = roundme(reduclass.n_rows*sizetr);
    arma::vec samp = csample_num(reduclass.col(1), nc, false,  ones<vec>(reduclass.n_rows));
    trainsel = join_cols(trainsel,samp);
  }
  return sort(trainsel);
}

//proximity matrix
// [[Rcpp::export]]
arma::mat proximi(arma::mat predtrnt, int m){
 arma::mat predtr = predtrnt.t();
arma::mat prox(predtr.n_rows, predtr.n_rows, fill::zeros);

for(int k = 0; k < predtr.n_cols; k++) { 
  
  for(int i = 0; i < predtr.n_rows; i++) {
    for(int j = 0; j < i; j++) {
      if (predtr(i, k) == predtr(j,k) ) {
        prox(i,j) += 1;
      } 
    }
  }

}
return(prox/m);
}


//maxvote used in tree_pred2
// [[Rcpp::export]]
arma::vec mvote(arma::mat votes){
  arma::vec clval = arma::unique(votes);
  arma::mat counts( clval.size(), votes.n_cols, fill::zeros);
    for(int k = 0; k < clval.size() ; k++) {
     for(int j=0; j< votes.n_cols ;j++){
      for(int i=0; i < votes.n_rows; i++) {
        if( votes(i, j) == k+1 ) {
          counts(k, j) += 1;
        }
        }
      }
    }
    
    arma::vec maxvote(votes.n_cols, fill::zeros);
    for(int i = 0; i < votes.n_cols; i++) {
      maxvote(i)= clval( index_max( counts.col(i) ));
    }
    return(maxvote);
    }
    
    
    
   // order bootstrap samples
  // [[Rcpp::export]]
  NumericMatrix oobindex(List datab, int m){
    int n = as<NumericVector>(datab[0]).size();
    NumericMatrix index(m, n); 
 for(int i = 0; i < m; i++) {
   index(i,_) = as<NumericVector>(datab[i]);
  }
 return(index);
  }
  //index <- as.matrix(plyr::ldply(data.b, function(x) c(pp=(x )))[,-1])
  
  
  
  
  
    // Filter oob observations 
    // [[Rcpp::export]]
    arma::mat oobobs(arma::mat index){ 
    int n = index.n_cols;
    int m = index.n_rows;
    arma::mat oobobs(m, n, fill::zeros);

      for(int i = 0; i < m; i++) {
        for(int j = 0; j< n ;j++){
       arma::uvec oobid = find( index.row(i) == j );
          if(oobid.size()==0){
              oobobs(i,j) = 1;
          }
      }
    }
    return(oobobs);
    }





//maxvote oob oob.pred here votes in pred.tr[[1]] I need to filtrate the oob observations
//to get the oob prediction
// [[Rcpp::export]]
arma::mat mvoteoob(arma::mat votes, arma::mat oobobs){
  
  arma::vec clval = arma::unique(votes);
  arma::mat counts( votes.n_cols, clval.size(), fill::zeros);
  for(int j=0; j< votes.n_cols ;j++){ //obs
    for(int i=0; i < votes.n_rows; i++) {//trees
      for(int k = 0; k < clval.size() ; k++) {
        
        if( (oobobs(i,j) > 0) && (votes(i, j) == k+1) ) { 
          counts(j, k) += 1;
        }
       }
      }
    }
  
  arma::vec maxvote(votes.n_cols, fill::zeros);
  for(int j = 0; j < votes.n_cols; j++) {
    maxvote(j)= clval( index_max( counts.row(j) ));
  }
  
  counts.insert_cols(counts.n_cols, maxvote);
  return counts;
  
  // return Rcpp::List::create(Rcpp::Named("Vote") = counts,
  //                           Rcpp::Named("maxvote") = maxvote);
  }

// [[Rcpp::export]]
arma::vec ooberrortree(arma::mat votes, arma::mat oobobs, arma::vec classe, int m){
  arma::vec err(m);
  for(int i=0; i < votes.n_rows; i++) {//trees
    double n = 0.0;
    double tot = 0.0;
    for(int j=0; j< votes.n_cols ;j++){ //obs
    if((oobobs(i, j) > 0) ) {
      n +=1;
      if(votes(i, j)!= classe(j)){
      tot +=1;
      }
    }
  }
    //Rcout << tot;
    err(i) = tot/n;
    
  }
  return(err);
}


// [[Rcpp::export]]
  List  PPclassification(arma::mat Treestruct, arma::mat testclassindex,
                                  arma::vec IOindex, arma::vec testclass, int id, int rep) {
  
    arma::vec IOindexL(IOindex.size(), fill::zeros);
    arma::vec IOindexR(IOindex.size(), fill::zeros);
    
      if(Treestruct(id, 3) == 0) {
        arma::vec iclass = testclass;
        
        
        
        iclass.elem( find(iclass > 0) ).ones();
        
       
        iclass = 1 - iclass;
        for(int i = 0; i< iclass.size(); i++){
        testclass(i) = testclass(i) + IOindex(i)*iclass(i)*Treestruct(id,2);
        }
        return Rcpp::List::create(Rcpp::Named("testclass") = testclass,
                                  Rcpp::Named("rep") = rep);
        
      }else{
        
        
        
        for(int i=0; i< IOindex.size(); i++){
       
    
          IOindexL(i) = IOindex(i)*testclassindex(rep,i);
          IOindexR(i) = IOindex(i)*(1 - testclassindex(rep,i));
        
         
          // Rcout<< "\n";
          // Rcout<< IOindexR;
        }
       
        rep = rep + 1;
        List a;
        a = PPclassification(Treestruct, testclassindex,
                               IOindexL, testclass, Treestruct(id, 1) -1, rep);

        testclass = as<vec>(a["testclass"]);

        rep = as<int>(a["rep"]);
        a = PPclassification(Treestruct, testclassindex,
                               IOindexR, testclass, Treestruct(id, 2) - 1, rep);
        testclass = as<vec>(a["testclass"]);
        rep = as<int>(a["rep"]);
      }

      return Rcpp::List::create(Rcpp::Named("testclass") = testclass , Rcpp::Named("rep")=rep);
    }


 // [[Rcpp::export]]
List PPclassindex(arma::vec classtemp,arma::mat testclassindex,
                                arma::mat testdata, arma::mat Treestruct, arma::mat AlphaKeep,
                                arma::mat CKeep, int id, int Rule){
      int n = classtemp.size();
      int tn;
      
      arma::uvec tindex(n, fill::zeros);
      arma::vec sele(n, fill::zeros);
     
      
      arma::vec tclass(n, fill::zeros);
      
       if(Treestruct(id, 1) == 0) {
         return Rcpp::List::create(Rcpp::Named("testclassindex") = testclassindex,
                                   Rcpp::Named("classtemp") = classtemp);
       }else{
         
         tclass = classtemp;
         
         
         if(sum( tclass == 0) == 0){
           tn = 0;    
         }else{
           arma::uvec tnaux = find(tclass==0);
           tn = tnaux.size();
         }
         tindex = arma::sort_index(tclass); 
         
         
         if(tn >  0){
            arma::uvec condidx = find(tclass==1);
           tindex = sort(tindex(condidx));
           
         }
           arma::mat tdata = testdata.rows(tindex);
          int idproj = Treestruct(id, 3) -1;
         arma::vec projtest = testdata*AlphaKeep.row(idproj).t();
         
         
         double ckcond = CKeep(idproj, Rule-1);
        
        
         for(int i=0; i< classtemp.size(); i++){
           if(projtest(i)<ckcond){
             classtemp(i) = true;
           }else{
             classtemp(i)  = false;
           }

         }

          testclassindex.insert_rows(testclassindex.n_rows, classtemp.t());
       
           
           List a;
           a = PPclassindex(classtemp, testclassindex,
                               testdata, Treestruct, AlphaKeep, CKeep, Treestruct(id, 1)-1, Rule);

           testclassindex = as<mat>(a["testclassindex"]);

           a = PPclassindex(1 - classtemp,testclassindex,
                               testdata, Treestruct, AlphaKeep, CKeep, Treestruct(id, 2)-1, Rule);
           testclassindex = as<mat>(a["testclassindex"]);
       }
     return Rcpp::List::create(Rcpp::Named("testclassindex") = testclassindex, Rcpp::Named("classtemp") = classtemp);
    
     }
       
       
#include <Rcpp.h>
       using namespace Rcpp;
       
       double innerProduct(NumericVector x, NumericVector y) {
         return std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
       }
       
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
       
       
       // return the prediction using permuted variables in the forest, used to compute the permuted inportance
       //variable measure like in randomForest
       
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
       
       
       

