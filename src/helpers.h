#ifndef HELPERS_H
#define HELPERS_H

#include <RcppArmadillo.h>
#include <RcppDist.h>

arma::vec conDMVN(
    arma::vec mean,
    arma::mat sigma,
    arma::uvec dependent_ind,
    arma::uvec given_ind,
    arma::vec X_given
);

arma::mat A_unique(arma::mat X, double n);

arma::mat generate_C(
    arma::mat A_star,
    arma::vec rho,
    double lambda,
    double n
);

double neg_log_lik(arma::mat cov, arma::vec y);

double log_gamma_pdf(double x, double alpha, double beta);

arma::vec normalise_y(arma::vec y);

arma::mat A_unique_asym(arma::mat X1, arma::mat X2, int n1, int n2);

arma::mat A_unique_sym(arma::mat X, double n);

arma::mat generate_C_asym(arma::mat A_star, arma::vec rho,double lambda,double n1, double n2);

arma::mat generate_C_sym(arma::mat A_star, arma::vec rho,double lambda,double n);

arma::vec forward_subst(arma::mat L, arma::vec b);

#endif
