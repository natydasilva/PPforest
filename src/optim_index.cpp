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
      gsum(k,(origclass(i)-1) )+=projdata(i,k);
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
  //double dw=as<double>(det(wrap(W)));
  double dw=0.0;
   dw =as<double>(det(wrap(WB)));
  index=1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
  return index;
}




// [[Rcpp::export]]
arma::vec LDAindex2(arma::vec origclass, arma::mat origdata, 
                 arma::mat proj, bool weight=true){
  double index=0.0;
  int n=origdata.n_rows, p=origdata.n_cols; 
  //arma::mat proj= proj.zeros(n,p);
  int q=proj.n_cols, p1=proj.n_rows;
  // Environment base("package:base");
  // Function table=base["table"];
  
   
   arma::vec clval = arma::unique(origclass);
   //----
   arma::vec newclass(n, fill::zeros);
   for(int i=0;i<n;i++){
   for (int k=0; k < clval.size(); k++) {
       if (origclass(i) == clval(k) ){
         newclass(i) = k +1;
       }
     }

   }
   origclass=newclass;
   arma::vec gn=tableC(origclass);
   int g=gn.size(); 

   //----
   
   
   if(p1!=p){
     q=p;
     }
    arma::vec allmean(q);
   arma::mat W(q,q),WB(q,q),gsum(q,g);        
   arma::mat projdata(n,q, fill::zeros);
   
  if(p1!=p||p1==1){
    projdata = origdata;
  }else{
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
  double temp=0.0;
  for(int i=0;i<n;i++){
    int l=origclass(i)-1;
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
        temp=((projdata(i,j1)-gsum(j1,l)/gn(l))*
                     (projdata(i,j2)-gsum(j2,l)/gn(l)) +
                     (gsum(j1,l)/gn(l)-allmean(j1))*
                     (gsum(j2,l)/gn(l)-allmean(j2)))/gn(l)*gn1;
        WB(j1,j2)+=temp;
        WB(j2,j1)=WB(j1,j2);
      }
    }
  }
//   Environment base("package:base");
//  Function det=base["det"];
//  double dw = as<double>(det(wrap(W)));
// index = 1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));



  index = 1.0-det(W)/det(WB);
  
   arma::vec out(q, fill::zeros);
   out(0)  =index; 
    return out;
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
      mean_all[i] = mean( origdata.col(i) ) ;// variable mean
      for (int k=0; k < g; k++) { // in groups
        double tot = 0.0;
        for (int j=0; j<n; j++) {//in observations
          if (origclass(j) == clval(k) ) tot += origdata(j,i) ; //total by variable and group 
        }
        mean_g(k, i) = tot/ng(k) ; //mean by variable and group
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
double PDAindex2(arma::vec origclass, arma::mat origdata,
                arma::mat proj,bool weight=true,
                double lambda=0.1){
  double index;
  int n=origdata.n_rows, p=origdata.n_cols; 
  //arma::mat proj= proj.zeros(n,p);
  int q=proj.n_cols, p1=proj.n_rows;
  // Environment base("package:base");
  // Function table=base["table"];
  arma::vec gn=tableC(origclass);
  int g = gn.size();  
  
 
 
 arma::vec clval = arma::unique(origclass);
 //----
 arma::vec newclass(n, fill::zeros);
 for(int i=0;i<n;i++){
   for (int k=0; k < clval.size(); k++) {
     if (origclass(i) == clval(k) ){
       newclass(i) = k +1;
     }
   }
   
 }
 origclass=newclass;
 
  arma::mat W(p,p,fill::zeros),WB(p,p,fill::zeros),gsum(p,g,fill::zeros);
  arma::vec allmean(p);
  
  if(p1!=p)
    q=p;
  for(int i=0;i<n;i++){
    for(int k=0;k<p;k++){
      allmean(k)+=origdata(i,k)/n;
      gsum(k,(origclass(i)-1))+=origdata(i,k);
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
  arma::mat Wt(q,p),WBt(q,p);    
  arma::mat Wtt(q,q),WBtt(q,q); 
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
  //Function det=base["det"];
  //index=1.0-as<double>(det(wrap(Wtt)))/as<double>(det(wrap(WBtt)));   

  index = 1.0-det(Wtt)/det(WBtt);
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
 
  if (PPmethod.compare("LDA")==0){
    //Rcout << "LDA\n";
    a1 = LDAopt(origclass,origdata,1, "LDA",true);
   //double indexbest = LDAindex(origclass,origdata);
   } else {
    //Rcout << "PDA\n";
     a1 = PDAopt(origclass, origdata, 1,"PDA",true, lambda);
  //double indexbest = PDAindex(origclass,origdata)
  }

  
  int  index = arma::index_max(arma::abs(a1));
  double sign = signC(a1(index));
  arma::vec classe = split_rel(origclass, origdata,  origdata*a1); 

  if (g > 2) {
        if (PPmethod.compare("LDA")==0){
         // Rcout << "still LDA\n";
          
      a2 = LDAopt(classe,origdata,1, "LDA",true);
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
 arma::vec ng = tableC(classe); 
 int g = ng.size();
 arma::vec clval = arma::unique(classe);

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
List findprojwrap(arma::vec origclass,arma::mat origdata, std::string PPmethod="LDA",
                  double sizep=1, double lambda =.1){

  int pp = origdata.n_cols;
  List dataspl = datanode(origdata, sizep );
  origdata = as<mat>(dataspl["data"]);
  arma::uvec vrnd = as<uvec>(dataspl["varselected"]);

  List oneDproj = findproj(origclass, origdata, PPmethod, lambda);
  arma::vec projdata = as<vec>(oneDproj["projdata"]);
  arma::vec classe = split_rel(origclass, origdata, projdata);
  //arma::vec classe =  as<vec>(oneDproj["class"]);
  //Rcout << classe;
  arma::mat projbest = as<mat>(oneDproj["projbest"]);

       arma::vec C = nodestr(classe, projdata);

       arma::mat Alpha = zeros<mat>(pp);
       for(int i=0; i<vrnd.size(); i++){
            int v = vrnd(i);
           Alpha(v) = projbest(i,0);

       }

  arma::vec indexbest(1);
  if(PPmethod.compare("LDA")==0){
 indexbest = LDAindex2(classe,origdata,projbest);
  }else{
    indexbest = PDAindex2(classe,origdata,projbest);
  }

       int n = projdata.n_rows;

         // group totals, group means and overall mean
         arma::vec ng = tableC(classe);
        
         int g = ng.size();
         arma::vec clval = arma::unique(classe);

         arma::colvec mean_g(g,fill::zeros);
         
         for (int k=0; k < g; k++) {
           double tot = 0.0;
           for (int j = 0; j<n; j++) {
             if (classe(j) == clval(k) ) {
               tot += projdata(j) ;

               }
           }
          
          mean_g(k) = tot/ng(k);
           //Rcout<<mean_g;
          // if (ng(k) > 0 ) {
          //  mean_g(k) = tot/ng(k) ;
          // } else {
          //   mean_g(k) = 0.0;
          // }
          
          
         }
         
           arma::colvec mLR = sort(mean_g);
           arma::uvec sortLR  = sort_index(mean_g);
           
         arma::vec IOindexL(classe.size());
         arma::vec IOindexR(classe.size());

         
         //-----
         // int ncl=classe.size();
         //  arma::vec newclass(ncl, fill::zeros);
         // for(int i=0;i<ncl;i++){
         //   for (int k=0; k < clval.size(); k++) {
         //     if (classe(i) == clval(k) ){
         //       newclass(i) = k ;
         //     }
         //   }
         //   
         // }
         // classe=newclass;
         // //----
         
         
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



//return indexbest;

return Rcpp::List::create(Rcpp::Named("Index") = indexbest,Rcpp::Named("Alpha") = Alpha, Rcpp::Named("C") = C,
                          Rcpp::Named("IOindexL")=IOindexL, Rcpp::Named("IOindexR") = IOindexR);

}

//Rcpp::Named("indexbest") = indexbest,
// [[Rcpp::export]]
List treeconstruct(arma::vec origclass, arma::mat origdata,arma::mat Treestruct, int id, int rep, int rep1, int rep2, arma::mat projbestnode, arma::mat  splitCutoffnode,
                   std::string PPmethod = "LDA", double lambda = 0.1) {

  //int n = origdata.n_rows;

  arma::vec g = tableC(origclass);
  int cl = g.n_rows;
  arma::vec clid = arma::linspace<vec>(1, cl, cl);

  //matrix with number of observations by class
  arma::mat gg;
  gg.insert_cols(gg.n_cols,clid);
  gg.insert_cols(gg.n_cols, g);

  int G = g.size();
  int ts =sum(Treestruct.row(0));
  int st = sum(splitCutoffnode.row(0));
  int GS = 2*G-1;
  arma::vec part = arma::linspace<vec>(1, GS, GS);
  List a;
  arma::vec C;
  //arma::vec  z = zeros<vec>(GS);
  if(ts == 0){
    Treestruct.insert_cols(0,part);
  }

  if(G == 1){
    //check as.numeric group names
    Treestruct(id, 2) = gg(0, 0);

    List::create( Treestruct, projbestnode,
                  splitCutoffnode, rep, rep1, rep2);
  }else{
    Treestruct(id, 1) = rep1;
    rep1 = rep1 + 1;
    Treestruct(id, 2) = rep1;
    rep1 = rep1 + 1;
    Treestruct(id, 3) = rep2;
    rep2 = rep2 + 1;
    a = findproj(origclass, origdata, PPmethod, lambda);

    // classe =as<vec>(a["class"]);
    C = nodestr(as<vec>(a["class"]),as<vec>(a["projdata"]));

    //splitCutoffnode = join_rows(splitCutoffnode, C);
    if(st == 0){
      splitCutoffnode.row(0) = C.t();
    }else{
      splitCutoffnode.insert_rows( splitCutoffnode.n_rows,C.t());
    }

    // Treestruct(id, 4) = as<int>(a["Index"]);
  }
  //     arma::mat Alpha = as<mat>(a["Alpha"]);
  //     projbestnode = join_rows(projbestnode, Alpha);
  //      arma::vec tclass = origclass;
  //       arma::mat tdata = origdata;
  //       arma::mat IOindexL = as<mat>(a["IOindexL"]);
  //       tclass = tclass*IOindexL;
  //
  //       arma::vec tnaux = tclass(tclass == 0);
  //         int tn = tnaux.size();
  //       arma::uvec tindex = sort_index(tclass);
  //       tindex = sort(tindex(-(arma::linspace<uvec>(0, tn-1, tn))));
  //       tclass = tclass(tindex);
  //       tdata = origdata(tindex);
  //
  //       List  b = treeconstruct(tclass, tdata, Treestruct, Treestruct(id, 2), rep, rep1, rep2, projbestnode,
  //                          splitCutoffnode, PPmethod, lambda);
  //
  //      Treestruct = as<mat>(b["Treestruct"]);
  //       projbestnode = as<vec>(b["projbestnode"]);
  //       splitCutoffnode = as<mat>(b["splitCutoffnode"]);
  //       rep = as<int>(b["rep"]);
  //       rep1 =  as<int>(b["rep1"]);
  //       rep2 =  as<int>(b["rep2"]);
  //       tclass = origclass;
  //       tdata =  origdata;
  //       arma::mat IOindexR = as<mat>(a["IOindexR"]);
  //       tclass = tclass*IOindexR;
  //
  //
  //       tnaux = tclass(tclass == 0);
  //       tn = tnaux.size();
  //       tindex = sort_index(tclass);
  //       tindex = sort(tindex(-(arma::linspace<uvec>(0, tn-1,tn))));
  //       tclass = tclass(tindex);
  //       tdata = origdata(tindex);
  //
  //
  //       n =  tdata.n_rows;
  //       arma::vec g = tableC(tclass);
  //       G = g.size();
  //       b = treeconstruct(tclass, tdata, Treestruct,
  //                           Treestruct(id, 3), rep, rep1, rep2, projbestnode,
  //                           splitCutoffnode, PPmethod, lambda);
  //       arma::mat Treestruct = as<mat>(b["Treestruct"]);
  //       projbestnode = as<vec>(b["projbestnode"]);
  //       splitCutoffnode = as<mat>(b["splitCutoffnode"]);
  //       rep = as<int>(b["rep"]);
  //       rep1 =  as<int>(b["rep1"]);
  //       rep2 =  as<int>(b["rep2"]);

  //}

  return Rcpp::List::create(Rcpp::Named("GS") = GS,
                            Rcpp::Named("gg") = gg, Rcpp::Named("Treestruct")=Treestruct,
                            Rcpp::Named("part")=part, Rcpp::Named("sp")=splitCutoffnode,
                            Rcpp::Named("a")=a,  Rcpp::Named("C")=C);



  // return Rcpp::List::create(Rcpp::Named("Treestruct") = Treestruct,
  //                           Rcpp::Named("projbestnode") = projbestnode,
  //                           Rcpp::Named("splitCutoffnode")= splitCutoffnode,
  //                           Rcpp::Named("rep")= rep,
  //                           Rcpp::Named("rep1")= rep1,
  //                           Rcpp::Named("rep2")= rep2);

}




// List treeconstruct(arma::vec origclass, arma::mat origdata, arma::mat Treestruct, int id, int rep, int rep1, int rep2, arma::mat projbestnode, arma::mat  splitCutoffnode,
//                         std::string PPmethod = "LDA", double lambda = 0.1 ) {
  //   origclass should be an integer
  
  
// 
//   splitCutoffnode = NA_REAL;
//   //Rcpp::Nullable<arma::mat> splitCutoffnode = R_NilValue;
//   projbestnode = NA_REAL;
//   Treestruct = NA_REAL;
//   id = 1;
//   rep1 = 2;
//   rep2 = 1;
//   rep = 1;
//   List Treefinal = treeconstruct(origclass, origdata, Treestruct,
//                                id, rep, rep1, rep2, projbestnode, splitCutoffnode,
//                                PPmethod, lambda);
// 
//     Treestruct = as<mat>(Treefinal["Treestruct"]);
// 
//     // colnames(Tree.Struct) <- c("id", "L.node.ID", "R.F.node.ID",
//     //          "Coef.ID", "Index")
//     projbestnode = as<mat>(Treefinal["projbestnode"]);
//     splitCutoffnode = as<mat>(Treefinal["splitCutoffnode"]);
// //colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
//    // List treeobj =  List::create list(Treestruct, projbestnode,
//    //                  splitCutoffnode, origclass,origdata )
// 
//     //class(treeobj) <- append(class(treeobj), "PPtreeclass")
//     return Rcpp::List::create(Rcpp::Named("Treestruct") = Treestruct,
//                               Rcpp::Named("projbestnode") = projbestnode,
//                               Rcpp::Named("splitCutoffnode")= splitCutoffnode,
//                               Rcpp::Named("origclass")= origclass,
//                               Rcpp::Named("origdata")= origdata);
//   
//     }



  
  
  
///////////
  
  
//
// Tree.construct <- function(origclass, origdata, Tree.Struct,
//                            id, rep, rep1, rep2, projbest.node, splitCutoff.node,
//                            PPmethod, r = NULL, lambda = NULL, ...) {
//   origclass <- as.integer(origclass)
//   n <- nrow(origdata)
//   g <- table(origclass)
//   G <- length(g)
//   if (length(Tree.Struct) == 0) {
//     Tree.Struct <- matrix(1:(2 * G - 1), ncol = 1)
//     Tree.Struct <- cbind(Tree.Struct, 0, 0, 0, 0)
//   }
//   if (G == 1) {
//     Tree.Struct[id, 3] <- as.numeric(names(g))
//     list(Tree.Struct = Tree.Struct, projbest.node = projbest.node,
//          splitCutoff.node = splitCutoff.node, rep = rep,
//          rep1 = rep1, rep2 = rep2)
//   }
//   else {
//     Tree.Struct[id, 2] <- rep1
//     rep1 <- rep1 + 1
//     Tree.Struct[id, 3] <- rep1
//     rep1 <- rep1 + 1
//     Tree.Struct[id, 4] <- rep2
//     rep2 <- rep2 + 1
//     a <- Find.proj(origclass, origdata, PPmethod, weight,
//                    r, lambda, ...)
//     splitCutoff.node <- rbind(splitCutoff.node, a$C)
//     Tree.Struct[id, 5] <- a$Index
//     projbest.node <- rbind(projbest.node, a$Alpha)
//     t.class <- origclass
//     t.data <- origdata
//     t.class <- t.class * a$IOindexL
//     t.n <- length(t.class[t.class == 0])
//     t.index <- sort.list(t.class)
//     t.index <- sort(t.index[-(1:t.n)])
//     t.class <- t.class[t.index]
//     t.data <- origdata[t.index, ]
//     b <- Tree.construct(t.class, t.data, Tree.Struct,
//                         Tree.Struct[id, 2], rep, rep1, rep2, projbest.node,
//                         splitCutoff.node, PPmethod, r, lambda,
//                                    ...)
//                         Tree.Struct <- b$Tree.Struct
//     projbest.node <- b$projbest.node
//     splitCutoff.node <- b$splitCutoff.node
//     rep <- b$rep
//     rep1 <- b$rep1
//     rep2 <- b$rep2
//     t.class <- origclass
//     t.data <- origdata
//     t.class <- (t.class * a$IOindexR)
//     t.n <- length(t.class[t.class == 0])
//     t.index <- sort.list(t.class)
//     t.index <- sort(t.index[-(1:t.n)])
//     t.class <- t.class[t.index]
//     t.data <- origdata[t.index, ]
//     n <- nrow(t.data)
//     G <- length(table(t.class))
//     b <- Tree.construct(t.class, t.data, Tree.Struct,
//                         Tree.Struct[id, 3], rep, rep1, rep2, projbest.node,
//                         splitCutoff.node, PPmethod, r, lambda,
//                                    ...)
//                         Tree.Struct <- b$Tree.Struct
//     projbest.node <- b$projbest.node
//     splitCutoff.node <- b$splitCutoff.node
//     rep <- b$rep
//     rep1 <- b$rep1
//     rep2 <- b$rep2
//   }
//   list(Tree.Struct = Tree.Struct, projbest.node = projbest.node,
//        splitCutoff.node = splitCutoff.node, rep = rep, rep1 = rep1,
//        rep2 = rep2)
// }
//
// splitCutoff.node <- NULL
// projbest.node <- NULL
// Tree.Struct <- NULL
// id <- 1
// rep1 <- 2
// rep2 <- 1
// rep <- 1
// Tree.final <- Tree.construct(origclass, origdata, Tree.Struct,
//                              id, rep, rep1, rep2, projbest.node, splitCutoff.node,
//                              PPmethod, r, lambda, TOL, ...)
//   Tree.Struct <- Tree.final$Tree.Struct
//   colnames(Tree.Struct) <- c("id", "L.node.ID", "R.F.node.ID",
//            "Coef.ID", "Index")
//   projbest.node <- Tree.final$projbest.node
//   splitCutoff.node <- Tree.final$splitCutoff.node
//   colnames(splitCutoff.node) <- paste("Rule", 1:8, sep = "")
//   treeobj <- list(Tree.Struct = Tree.Struct, projbest.node = projbest.node,
//                   splitCutoff.node = splitCutoff.node, origclass = origclass,
//                   origdata = origdata)
//   class(treeobj) <- append(class(treeobj), "PPtreeclass")
//   return(treeobj)
//   }


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


