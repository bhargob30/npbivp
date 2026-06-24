#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


// [[Rcpp::export]]
Rcpp::List pred_cov_cpp(Rcpp::List data_types_train, Rcpp::List data_types_test, arma::vec y, arma::mat W_train, arma::mat W_test, arma::mat samples, bool scale_data=true ){
  //samples is a matrix of size M_reduced * param_reduced
  // param_reduced is the params reqd for computing the covariance
  // M_reduced is the no of MCMC samples over which prediction is done
  //W means covariates (confounders) of size n by c
  //scale_data is bool indicating whether we want to scale data to (0,1) range while predicting. Else by default it will consider original (current) or scaled versions of input matrix. 
  
  //append covariates with data_train/test as another platform
  
  int m_init=data_types_train.size();
  Rcpp::List L_dat(1+m_init);
  for(int j=0; j<m_init; j++){
    L_dat[j]=data_types_train[j];
    
  }
  L_dat[m_init]=W_train;
  
  data_types_train=L_dat;
  
  Rcpp::List L_dat_1(1+m_init);
  for(int j=0; j<m_init; j++){
    L_dat_1[j]=data_types_test[j];
  }
  L_dat_1[m_init]=W_test;
  data_types_test=L_dat_1;
  
  
  arma::mat n_find=Rcpp::as<arma::mat>(data_types_train[0])  ; 
  int n = n_find.n_rows;  // no of observation in train set
  
  n_find=Rcpp::as<arma::mat>(data_types_test[0])  ;
  int n1 = n_find.n_rows;  // no of observation in test set
  
  int m=data_types_train.size();  // no. of data types
  
  arma::field<arma::mat> A_stars_var(m); // distance for each train data matrix
  arma::field<arma::mat> A_stars_cov(m); // distance between train and test
  arma::field<arma::mat> A_stars_var_test(m); // distance for each test data matrix
  arma::uvec p(m);
  arma::mat chol_var(n,n);
  arma::mat U(n,n1);
  arma::vec v;
  arma::mat temp1;
  arma::mat temp2;
  arma::mat temp3;
  arma::mat temp4;
  for(int i=0;i<m;i++){
    temp1=Rcpp::as<arma::mat>(data_types_train[i]);
    temp2=Rcpp::as<arma::mat>(data_types_test[i]);
    
    if(scale_data==true){
      //scale data to lie in (0,1)
      temp1.each_col( [](arma::vec& a){
        a=(a-arma::min(a))/(arma::max(a)-arma::min(a));
        
      } );
      
      temp2.each_col( [](arma::vec& a){
        a=(a-arma::min(a))/(arma::max(a)-arma::min(a));
        
      } );
      
      
    }
    
    
    
    A_stars_var(i)=A_unique_sym(temp1,n);
    A_stars_cov(i)=A_unique_asym(temp2,temp1,n1,n);
    A_stars_var_test(i)=A_unique_sym(temp2,n1);
    p(i)= temp1.n_cols;
  }
  
  arma::uvec q=arma::join_cols(arma::zeros<arma::uvec>(1),p);
  
  arma::mat prediction(n1, samples.n_rows,arma::fill::zeros); // prediction matrix. Each column is a prediction (i.e. conditional mean vector)[E(Y_test|X_test,Y,param)] corresponding to a particular MCMC sample
  arma::mat cond_vars(n1, samples.n_rows,arma::fill::zeros); // matrix of conditional variances. Each column is a vector of conditional variances [diag(Var(Y_test|X_test,Y,param))] corresponding to a particular MCMC sample
  arma::vec rhos;
  arma::vec lambdas;
  arma::vec w_2s;
  double r;
  for(int t=0; t<samples.n_rows; t++){
    rhos=samples.row(t).subvec(0, sum(p)-1).t();
    lambdas=samples.row(t).subvec(sum(p),sum(p)+m).t();
    w_2s=samples.row(t).subvec(sum(p)+m+1,sum(p)+m+m).t();
    r=samples.row(t)(sum(p)+m+m+1);
    
    arma::mat var=(arma::eye(n,n)/r)+(arma::ones(n,n)/lambdas(0));
    arma::mat var_test=(arma::eye(n1,n1)/r)+(arma::ones(n1,n1)/lambdas(0));
    arma::mat cov=(arma::ones(n1,n)/lambdas(0));
    int a;
    int b;
    for(int i=0;i<m;i++){
      a=sum(q.subvec(0,i));
      b=sum(q.subvec(0,i+1))-1;
      var=var+w_2s(i)*generate_C_sym(A_stars_var(i),rhos.subvec(a,b),lambdas(1+i),n);
      cov=cov+w_2s(i)*generate_C_asym(A_stars_cov(i),rhos.subvec(a,b),lambdas(1+i),n1,n);
      var_test=var_test+w_2s(i)*generate_C_sym(A_stars_var_test(i),rhos.subvec(a,b),lambdas(1+i),n1);
    }
    //cholesky decomp of var
    
    chol_var=arma::chol(var,"lower");
    U=arma::solve(arma::trimatl(chol_var), cov.t(), arma::solve_opts::fast);
    v=arma::solve(arma::trimatl(chol_var), y, arma::solve_opts::fast);
    //temp3=arma::inv_sympd(var,arma::inv_opts::allow_approx);
    prediction.col(t)=U.t() * v;
    temp4=var_test - (U.t() * U);
    cond_vars.col(t)=temp4.diag();
    
    
  }
  
  Rcpp::List L(2);
  L[0]=prediction;
  L[1]=cond_vars;
  return L;
  
}
