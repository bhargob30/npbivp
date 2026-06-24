#include "helpers.h"

arma::vec conDMVN(arma::vec mean, arma::mat sigma, arma::uvec dependent_ind, arma::uvec given_ind, arma::vec X_given) {
  arma:: mat B=sigma.submat(dependent_ind, dependent_ind);
  arma:: mat C=sigma.submat(dependent_ind, given_ind);
  arma:: mat D=sigma.submat(given_ind, given_ind);
  arma:: mat L=arma::chol(D,"lower");
  arma:: mat u=arma::solve(arma::trimatl(L),X_given -mean.elem(given_ind) , arma::solve_opts::fast);
  arma:: mat U=arma::solve(arma::trimatl(L),C.t() , arma::solve_opts::fast);
  //arma:: mat CDinv= C * arma::inv_sympd(D, arma::inv_opts::allow_approx);
  arma::vec cMu= mean.elem(dependent_ind)+ U.t() * u;
  arma:: mat cVar=B- U.t() * U;
  arma::vec out={cMu(0),cVar(0,0)};
  return out;
}

arma::mat A_unique(arma::mat X, double n){
  double unique_elem=(pow(n,2)+n)/2;
  int var =(int)X.n_cols;
  arma::mat A_star( unique_elem, var   );
  arma::rowvec temp(var);
  double count=0;
  for(int i=0; i<n; i++){
    for(int j=i;j<n;j++){
      if(i==j){
        A_star.row(count)=arma::zeros(var).t();
      } else{

        temp=X.row(i)-X.row(j);
        A_star.row(count)= temp.transform([](double val) { return (pow(val,2));} );
      }
      count=count+1;
    }

  }

  return A_star;
}

arma::mat generate_C(arma::mat A_star, arma::vec rho,double lambda,double n){ //need to change if we take two lambda params
  int a;
  int b;
  arma::uvec zero_ind=arma::find(rho==1.0);
  rho.shed_rows(zero_ind);
  A_star.shed_cols(zero_ind);
  arma::vec neg_Gvec=A_star *arma::log(rho);
  neg_Gvec=(1.0/lambda)*arma::exp(neg_Gvec); // it is actually Cvec,    change here if there are two lambdas
  arma::mat output_C(n,n);
  output_C.row(0)=neg_Gvec.subvec(0,n-1).t();

  for(int i=0; i<n; i++){
    a=(i*n)- (((i-1)*i)/2);
    b=((i+1)*n) - ((i*(i+1))/2) -1;
    output_C.row(i).subvec(i,n-1)=neg_Gvec.subvec(a,b).t();
  }
  return arma::symmatu(output_C);

}

double neg_log_lik(arma::mat cov,arma:: vec y){
  //cholesky decompo
  arma::mat L=arma::chol(cov, "lower");
  arma::vec u= arma::solve(arma::trimatl(L), y, arma::solve_opts::fast);
  return sum(log(L.diag()))+ 0.5*( arma::as_scalar(u.t() * u));
}

double log_gamma_pdf(double x, double alpha, double beta){
  return  (alpha* log(beta))- lgamma(alpha)+((alpha-1)*log(x))-(beta*x);
}


arma::vec normalise_y(arma::vec y){
  double y0=arma::max(y);
  arma::vec y_new=arma::join_cols(y, arma::zeros<arma::vec>(1));
  y_new=y_new-y0*arma::ones<arma::vec>(y.n_elem+1);
  arma::vec temp_y= arma::exp(y_new);
  return temp_y/arma::accu(temp_y);

}


// Function to calculate A_star for two data matrices X1 and X2. nrow(X1)=n1 <= nrow(X2)=n2


arma::mat A_unique_asym(arma::mat X1, arma::mat X2, int n1, int n2){

  double unique_elem=n1*n2;
  int var =(int)X1.n_cols; //X1 and X2 have same no. of cols
  arma::mat A_star( unique_elem, var   );
  arma::rowvec temp(var);
  double count=0;
  for(int i=0; i<n1; i++){
    for(int j=0;j<n2;j++){
      temp=X1.row(i)-X2.row(j);
      A_star.row(count)= temp.transform([](double val) { return (pow(val,2));} );

      count=count+1;
    }

  }

  return A_star;
}

arma::mat A_unique_sym(arma::mat X, double n){
  double unique_elem=(pow(n,2)+n)/2;
  int var =(int)X.n_cols;
  arma::mat A_star( unique_elem, var   );
  arma::rowvec temp(var);
  double count=0;
  for(int i=0; i<n; i++){
    for(int j=i;j<n;j++){
      if(i==j){
        A_star.row(count)=arma::zeros(var).t();
      } else{

        temp=X.row(i)-X.row(j);
        A_star.row(count)= temp.transform([](double val) { return (pow(val,2));} );
      }
      count=count+1;
    }

  }

  return A_star;
}




//function to generate K(X1,X2) for each dataset


arma::mat generate_C_asym(arma::mat A_star, arma::vec rho,double lambda,double n1, double n2){
  int a;
  int b;
  arma::uvec zero_ind=arma::find(rho==1.0);
  rho.shed_rows(zero_ind);
  A_star.shed_cols(zero_ind);
  arma::vec neg_Gvec=A_star *arma::log(rho);
  neg_Gvec=(1.0/lambda)*arma::exp(neg_Gvec); // it is actually Cvec,    change here if there are two lambdas
  arma::mat output_C(n1,n2);


  for(int i=0; i<n1; i++){
    a=(i*n2);
    b=((i+1)*n2) -1;
    output_C.row(i)=neg_Gvec.subvec(a,b).t();

  }
  return output_C;

}


arma::mat generate_C_sym(arma::mat A_star, arma::vec rho,double lambda,double n){ //need to change if we take two lambda params
  int a;
  int b;
  arma::uvec zero_ind=arma::find(rho==1.0);
  rho.shed_rows(zero_ind);
  A_star.shed_cols(zero_ind);
  arma::vec neg_Gvec=A_star *arma::log(rho);
  neg_Gvec=(1.0/lambda)*arma::exp(neg_Gvec); // it is actually Cvec,    change here if there are two lambdas
  arma::mat output_C(n,n);
  output_C.row(0)=neg_Gvec.subvec(0,n-1).t();

  for(int i=0; i<n; i++){
    a=(i*n)- (((i-1)*i)/2);
    b=((i+1)*n) - ((i*(i+1))/2) -1;
    output_C.row(i).subvec(i,n-1)=neg_Gvec.subvec(a,b).t();
  }
  return arma::symmatu(output_C);

}


arma::vec forward_subst(arma::mat L, arma::vec b){
  arma::vec y(L.n_rows);
  double alpha;
  for(int i=0; i<L.n_rows;i++){
    alpha=b(i);
    for(int j=0; j<i;j++){
      alpha=alpha-L(i,j)*y(j);
    }
    y(i)=alpha*(1.0/L(i,i));
  }
  return y;
}




