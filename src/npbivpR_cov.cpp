#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
Rcpp::List npbivpr_cov_cpp(Rcpp::List data_types, arma::vec y, arma::mat W, arma::vec obs_indicator, int M, int burn_in, arma::mat P_way, arma::mat Cand, arma::vec nu_R, arma::vec lambdas_prior_param, arma::vec r_prior_param, arma::vec r_muL_prior_param, arma::vec r_ev_prior_param,arma::vec tau_prior_param, double beta_prior_tau, double R_prior_q, double sigma_prior, arma::vec sigma_prop, arma::vec c_lambdas, double c_prec, arma::vec init, int prop_scheme, bool weight_upd=true, bool icept_update=true, bool same_lambda=true, bool RN=true, bool scale_data=true){ //sigma_prop--proposal sd for transformed weights //c_lambdas--tuning for lambdas, //c_prec---tuning for precision

  //icept_update is just a bool indicating whether we want to update the lambda_0 or not, if not we have to center the response y
  //scale_data is bool indicating whether we want to scale data to (0,1) range while fitting outcome GP model. Else by default it will consider original (current) or scaled versions of input matrix.
  //Cand is the candidate gene mirna matrix, P_way is Pathway info

  //sigma_prop--proposal sd for transformed weights //c_lambdas--tuning for lambdas, //c_prec---tuning for precision
  //c_lambdas is of length m+1
  //lambdas_prior param is to be decided --we can give all lambdas (m+1 in total) different params
  //update:- so we are taking lambdas_prior_param=c(a1,b1,a2,b2) where the first two correspond to intercept and the last two are common for scaling factors for each data type

  //prop schme is a double '1/2' etc. meanign different proposal scheme we are trying for explicit/implicit weighting scheme. Currently 1 means usual Gamma proposal, i.e. c=1
  /// 2 means keep on updating c in Gamma proposal untill target acceptance is reached. So, c changes with iteration
  // In all these cases, transformed wrights are proposed with proposal adpated
  /// according to usual Normal proposal with continuous tuning like scheme '2' described above


  //in the arguments above W means covariates (confounders)


  //we have to make several changes in the way inputs are provided before making it a package since we will want it to make more user friendly. For example, we should probably just give one option to get weights:- either explitcit,i.e. w2s themselvles or via updating individual lambdas.
  // Otherwise it will make user more confused. Also note that current apporaches for showing posterior consistency would probably relie on working with individual lambdas, afterall this is just one-to-one transformation and we show identifiability in lambdas which automatically indicate identifiability of weights since weights are just reweighted version of those lambdas
  // similarly for NpBIVP same function will be able to be used for continuous as well as survival outcome. We would need some indicator for that as well, so that when 'survival' user has to provide obs_indicator as well. But if 'continuous/binary' then obs_indicator will just be a vector of 1s. Since in C++ we can't give any condition in the function argument itself, it will be better to
  // wrap this cpp function under a R function to get those functionalities. Also need to update dim and other related params as well if we just learn lambdas

  //Also, note weight_upd=F and same_lambda=T does not make any sense. So, maybe later we have to wrap everything in R and there we have to give an arguemnt 'parametrization' with two options,
  // (1) default which will be updating individual lambdas and implicitly learning weights (i.e. both weight_upd=F, same_lambda=F). Or (2) Updating weights on their own.
  //Also have to give a character vector for data transformaton :- 'No', 'min_max', 'mean_sd'


  //to ensure minimal changes in arg input from earlier code/script, so that we can just take the covariate
  // as just another platform (so m=2 changes to m+1=3) with last element of data_types is
  // just the covariates. We just don't have to select any covariates by ensuring gamma=1
  // and sampling rho. Note that m is read later below as no of elements in list data_types.
  // We will just change the input data_types to append W at last, and remove earlier lineat kernel part for W

  // append covariate as last data types. Note many of this scripts's code uses
  // full specifics that we are having two data platforms for Np-BIVPR, and many more.
  //for the general script for Np-BIVP we have to be careful

  Rcpp::List L_dat(3);
  L_dat[0]=data_types[0];
  L_dat[1]=data_types[1];
  L_dat[2]=W;

  data_types = L_dat;


  // obs_ind=observartion indicator.
  arma::vec y_original=y;
  //note that y will be updated once we sample censored values
  // so we save the original value of y to have access to the
  // lower truncation point of censored observations
  arma::uvec cens_ind=arma::find(obs_indicator==0.0); //censored indices
  int total_cens=cens_ind.n_elem; //total no of censored obs
  //arma::uvec obs_ind=arma::find(obs_indicator==1.0); //obs indices (updates:- Prolly not needed)


  arma::mat X_1=Rcpp::as<arma::mat>(data_types[0])  ; // first data matrix--or gene which will be taken as response in RN

  //scale  X_1
  X_1.each_col( [](arma::vec& a){
    a=(a-arma::mean(a))/(arma::stddev( a ));

  } );


  int n = X_1.n_rows;  // no of observation
  int m=data_types.size();  // no. of data types (this will be 3, two platforms+covariates)



  //New comment:- most of these codes are taken directly from the general m case, so while we trivially know there are 2 data types many of things are that's wht being performed as if for general m>=2
  // only some added parts will directly use the fact that we have m=2.

  arma::mat X_2=Rcpp::as<arma::mat>(data_types[1]);  //miRNA or the other data type

  //scale  X_2
  X_2.each_col( [](arma::vec& a){
    a=(a-arma::mean(a))/(arma::stddev( a ));

  } );


  int a;
  int b;

  arma::uvec all_ind=arma::linspace<arma::uvec>(0,n-1,n); // vector of indices 0,1,....n-1
  arma::field<arma::mat> A_stars(m);
  arma::uvec p(m);
  for(int i=0;i<m;i++){
    arma::mat temp=Rcpp::as<arma::mat>(data_types[i]);

    if(scale_data==true){
      //scale data to lie in (0,1)
      temp.each_col( [](arma::vec& a){
        a=(a-arma::min(a))/(arma::max(a)-arma::min(a));

      } );

    }

    A_stars(i)=A_unique(temp,n);
    p(i)= temp.n_cols;
  }




  arma::uvec q=arma::join_cols(arma::zeros<arma::uvec>(1),p);

  //New comment:-under chnages for RN, I think there will be changes. 'dim' might be used in some
  // other later parts, so not changing. dimension of acceptance have to be changed if we decide to keep
  //track of acceptance rate of other params. So, at present not changing anything --02/23/2025

  int dim=2*sum(p)+1+m+m +1 + total_cens; // total parameters in our model, here we are taking the general case, i.e. m+1 lambdas, and also we will save actual weights than the transformed ones
  arma:: mat samples(dim,M-burn_in, arma::fill::zeros); //output of Rcpp function (mcmc samples)
  arma:: mat acceptance(1,1+m+m+1,arma::fill::zeros); //output of Rcpp function (acceptance rates)
  arma::vec nllik(M-burn_in);  //although this is probably just likelihood from outcome model

  int time_count; //output of Rcpp function (total time count of M samples)

  //setting of prior param--each of the following prior can be part args in Rcpp function

  arma:: vec param_r=r_prior_param; // prior param of r, (a,b) in notation
  arma:: vec param_lambda= lambdas_prior_param; // prior param of lambda_j's:-(a1,b1,a2,b2)




  //initialization       #20/1000 initially
  double non_zero_counts = std::floor((20.0/1000.0)*sum(p)); //no of variables to be randomly set as active
  arma::uvec non_zero_indices = arma::randi<arma::uvec>(non_zero_counts, arma::distr_param(0,sum(p)-1));
  arma::vec gammas_old=arma::zeros(sum(p));
  gammas_old.elem(  non_zero_indices).ones();

  arma::vec rhos_old=arma::randu<arma::vec>(sum(p));
  rhos_old.elem( find(gammas_old==0.0)   ).ones();        //rho

  //re-initialize parameters of covariate platform according to what we want, i.e. all gammma=1 and rho~U(0,1)
  int j=m-1;
  for(int k=sum(q.subvec(0,j));k<sum(q.subvec(0,j+1));k++){

    gammas_old(k)= 1.0;
    rhos_old(k)=arma::randu();
  }



  arma::vec lambdas_old =arma::ones(1+m); //lambdas=(lambda_0,lambda_1,...lambda_m)

  lambdas_old(0)=arma::randu();

  arma::vec lambdas_current(1+m) ;

  if(icept_update==false){
    lambdas_old(0)= 1e16;
    lambdas_current(0)=1e16;
  }



  if(same_lambda==true){
    lambdas_old(1)=init(0); //currently looks like init only handles witrh same lambda case.
    for(arma::uword i=2;i<m+1;i++){
      lambdas_old(i)=lambdas_old(1);
      lambdas_current(i)=lambdas_current(1);
    }
  }



  arma::vec y_old;
  arma::vec w_2s_old;
  arma::vec y_current;
  arma::vec lambda_inv(m);
  arma::vec w_2s_old_1(m); // the implicit weights if we just learn lambdas


  y_old=init.subvec(1,m-1); //transformed weights
  w_2s_old= normalise_y(y_old);
  y_current = arma::zeros(m-1);

  if(weight_upd==false){
    lambda_inv=arma::pow(lambdas_old.subvec(1,m), -1.0);
    w_2s_old_1=lambda_inv/arma::accu(lambda_inv);
    w_2s_old=arma::ones(m);
  }




  double r_old=init(m); //r

  //arma::vec cens_obs=y.elem(cens_ind); //sample of censored observations


  //New Comment:- new params initialization
  int L=P_way.n_cols; //no. of total pathways
  arma::mat R_old=arma::zeros<arma::mat>(p(0),p(1));
  //for each gene randomly set 10 percent of the candidate miRNAs to be in RN of that gene
  //Note:- earlier this was 50 percent. we recently changed it to 10%
  // Have to check if this works for real GBM data, earlier GBM analyses were under 50%
  arma::vec zero_row(p(1));
  for(arma::uword i=0;i<p(0);i++){
    arma::uvec non_zero_ind= arma::find(Cand.row(i).t()==1.0);
    non_zero_ind=arma::shuffle(non_zero_ind);
    zero_row=arma::zeros(p(1));
    zero_row.elem(non_zero_ind.head(std::floor(0.1*non_zero_ind.n_elem))).ones();
    R_old.row(i)=zero_row.t();
  }

  arma::mat R_current(p(0),p(1));
  arma::vec r_v_current(p(1)); //vth row of R

  arma::vec R_tau_samples(2,arma::fill::randu);

  arma::vec mu_old(L,arma::fill::randn);
  arma::vec mu_current(L);
  double mu_L_current; //L th elem of mu

  arma::vec r_mu=arma::randg<arma::vec>(L,arma::distr_param(1.0,1.0));
  arma::vec z_1(p(0),arma::fill::randn);
  arma::vec z_2(p(1),arma::fill::randn);
  arma::vec r_e=arma::randg<arma::vec>(p(0),arma::distr_param(1.0,1.0));

  arma::field<arma::vec> R_sums(2); // row and col sums of R, 0 element being row sum


  arma::vec gammas_current(sum(p));
  arma::vec rhos_current(sum(p));
  arma::vec w_2s_current(m);
  double r_current;
  arma::uvec dependent_ind(1);
  arma::uvec given_ind(n-1);
  arma::uvec temp;
  //arma::vec cond_mean_var;
  arma::mat cov_inv;
  arma::mat F11_inv;
  arma::vec u2;
  double d;
  arma::vec cMu;
  arma:: mat cVar;
  arma:: mat B;
  arma:: mat C;
  double log_lik_mu_old; //log lik part needed in metrop of muL's
  double log_lik_mu_current;
  arma::vec rv_temp(p(1));

  //samples.col(0)=arma::join_cols(arma::join_cols(gammas_old,rhos_old), arma::join_cols(lambdas_old,w_2s_old),r_old*arma::ones(1),y.elem(cens_ind)); //set 0th iteration

  arma::field<arma::mat> K_js_old(m); //K_j from gaussian process kernel for all the datasets


  for(int i=0;i<m;i++){
    a=sum(q.subvec(0,i));
    b=sum(q.subvec(0,i+1))-1;
    K_js_old(i) =generate_C(A_stars(i),rhos_old.subvec(a,b),lambdas_old(1+i),n);
  }

  arma::field<arma::mat> K_js_current(m);




  //full covariance matrix under old params

  arma::mat cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));

  for(int i=0;i<m;i++){
    cov=cov+(w_2s_old(i)*K_js_old(i));
  }


  arma::field<arma::vec> A_prods(m); //save the vector A_star * log(rho) for each data type
  for(int i=0;i<m;i++){
    a=sum(q.subvec(0,i));
    b=sum(q.subvec(0,i+1))-1;
    A_prods(i)=A_stars(i) * arma::log(rhos_old.subvec(a,b));
  }



  arma::vec neg_del_Gvec;
  arma::vec neg_del_Gvec_temp;

  double neg_log_lik_old=neg_log_lik(cov,y);
  //nllik(0)=neg_log_lik_old;
  double neg_log_lik_current;



  double gamma_k_prime;
  double rho_k_prime;
  double lambda_prime;
  double g_prime;
  double log_r;
  double alpha;
  double u;
  double gammas_temp;
  int k_star;
  double log_norm_sum;
  double trun_sample_y;

  arma::wall_clock timer;
  timer.tic();

  double probit_prob;
  //for each l=1,..L, get set of genes (or indices) such that p_vl=1
  arma::field<arma::uvec> p_gene_ind(L);
  arma::vec temp_2(p(0));
  for(arma::uword l=0;l<L;l++){
    temp_2=P_way.col(l);
    p_gene_ind(l)=arma::find(temp_2); //find no zero i.e. those which are 1
  }

  double log_det_temp;
  arma::vec mean_mu_temp(n);
  arma::mat temp_inv_mu(n,n);
  arma::uword ind_R;
  arma::uword ind_0R;
  arma::uword ind_1R;
  double log_lik_R_old;
  double log_lik_R_current;

  arma::vec ak_0_s(p(0));
  arma::vec bk_0_s(p(0));
  arma::vec ak_1_s(p(1));
  arma::vec bk_1_s(p(1));

  double A_tau;
  double B_tau;


  arma::mat samples_mu(L,M-burn_in);
  //arma::mat samples_r_mu(L,M);
  //arma::mat samples_z(p(0)+p(1),M);
  arma::mat samples_tau(2, M-burn_in);
  //arma::mat samples_r_e(p(0), M);
  //arma::cube samples_R(p(0),p(1),M);

  //samples_mu.col(0)=mu_old;
  //samples_r_mu.col(0)=r_mu;
  //samples_r_e.col(0)=r_e;
  //samples_z.col(0)=arma::join_cols(z_1,z_2);
  //samples_tau.col(0)=R_tau_samples;
  //samples_R.slice(0)=R_old;
  arma::mat cv_inv(n*n,p(0));

  //save summary of R rather than R for each MCMC iteration
  arma::mat R_sum_save=arma::zeros(p(0),p(1));


  //calculate and save cv_inverse's and update it while updating R:-


  for(arma::uword v=0;v<p(0);v++){

    rv_temp=R_old.row(v).t();
    arma::uvec find_temp_one=arma::find(rv_temp==1.0); //note that there may be no miRNAs for this gene which is related to it
    // so this indices may be empty

    //do the following step iff the above indices vector is non empty
    if(find_temp_one.n_elem >0){
      arma::mat X2_temp=X_2;
      X2_temp.shed_cols(arma::find(rv_temp==0.0)); //X_2_rv
      temp_inv_mu=arma::eye(n,n) - beta_prior_tau*X2_temp*arma::inv_sympd(arma::eye(X2_temp.n_cols,X2_temp.n_cols)+beta_prior_tau*X2_temp.t()*X2_temp)*X2_temp.t();
    } else {
      temp_inv_mu=arma::eye(n,n);
    }

    cv_inv.col(v)=arma::vectorise(temp_inv_mu);

  }

  double Al;
  double Bl;
  arma::vec cv_inv_ones(n);

  //Rcpp::Rcout << "Starting main loop..." << std::endl;

  for (int t=1; t< M; t++){

    //std::cout << "Sample..." << t << std::endl;
    //std::cout << "acceptance" << acceptance/t <<std::endl;


    //calculate row_col sums of R
    R_sums(0)=arma::sum(R_old,1);
    R_sums(1)=arma::sum(R_old,0).t();



    for(int j=0;j<(m-1);j++){ //this loop will be done for actual datatypes excluding covariates

      for(int k=sum(q.subvec(0,j));k<sum(q.subvec(0,j+1));k++){

        //Rcpp::Rcout << "k is:" << k << std::endl;

        if( gammas_old(k) ==1.0  ){
          gamma_k_prime=0.0;
          rho_k_prime=1.0;
        } else{
          gamma_k_prime=1.0;
          rho_k_prime=arma::randu();
        }
        k_star=k-sum(q.subvec(0,j));

        neg_del_Gvec_temp= A_prods(j) + ( A_stars(j).col(k_star)*std::log(rho_k_prime)-A_stars(j).col(k_star)*std::log(rhos_old(k)) );
        neg_del_Gvec=(1.0/lambdas_old(1+j))*arma::exp(neg_del_Gvec_temp);
        arma::mat delta(n,n);
        delta.row(0)=neg_del_Gvec.subvec(0,n-1).t();
        for(int i=0; i<n; i++){
          a=(i*n)- (((i-1)*i)/2);
          b=((i+1)*n) - ((i*(i+1))/2) -1;
          delta.row(i).subvec(i,n-1)=neg_del_Gvec.subvec(a,b).t();
        }

        delta=arma::symmatu(delta);
        K_js_current=K_js_old;
        K_js_current(j)=delta;
        gammas_current=gammas_old;
        rhos_current=rhos_old;
        gammas_current(k)=gamma_k_prime;
        rhos_current(k)=rho_k_prime;

        //printing for error debugging. Delete later
        //a=sum(q.subvec(0,j));
        //b=sum(q.subvec(0,j+1))-1;

        //Rcpp::Rcout << "max difference:" << arma::max(arma::max(arma::abs(K_js_current(j)-generate_C(A_stars(j),rhos_current.subvec(a,b),lambdas_old(1+j),n)))) << std::endl;

        // delete the above part later



        cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));

        for(int i=0;i<m;i++){
          cov=cov+(w_2s_old(i)*K_js_current(i));
        }


        neg_log_lik_current=neg_log_lik(cov,y);

        //Rcpp::Rcout << "CheckPoint 1..." << std::endl;

        //change alpha_k_j to probit probability here if RN is true
        arma::vec temp1=R_sums(j);

        if(RN==true){

          probit_prob=arma::normcdf(nu_R(j)+R_tau_samples(j)*temp1(k_star));

        } else{

          probit_prob=arma::normcdf(nu_R(j));

        }



        log_r=(gamma_k_prime-gammas_old(k))*log( probit_prob/(1.0-probit_prob ) ) + neg_log_lik_old  -neg_log_lik_current;
        gammas_temp= gammas_old(k);

        alpha =std::min(1.0,exp(log_r));
        u = arma::randu();
        if(u<=alpha){
          gammas_old=gammas_current;
          rhos_old=rhos_current;
          K_js_old=K_js_current;
          A_prods(j)=neg_del_Gvec_temp;
          neg_log_lik_old=neg_log_lik_current;

        }

        //within model

        if(gammas_old(k)==1.0){

          gamma_k_prime=1.0;  //actually double prime
          rho_k_prime=arma::randu();
          neg_del_Gvec_temp= A_prods(j) + ( A_stars(j).col(k_star)*std::log(rho_k_prime)-A_stars(j).col(k_star)*std::log(rhos_old(k)) );
          neg_del_Gvec=(1.0/lambdas_old(1+j))*arma::exp(neg_del_Gvec_temp);
          arma::mat delta(n,n);
          delta.row(0)=neg_del_Gvec.subvec(0,n-1).t();
          for(int i=0; i<n; i++){
            a=(i*n)- (((i-1)*i)/2);
            b=((i+1)*n) - ((i*(i+1))/2) -1;
            delta.row(i).subvec(i,n-1)=neg_del_Gvec.subvec(a,b).t();
          }
          delta=arma::symmatu(delta);
          K_js_current=K_js_old;
          K_js_current(j)=delta;
          gammas_current=gammas_old;
          rhos_current=rhos_old;
          gammas_current(k)=gamma_k_prime;
          rhos_current(k)=rho_k_prime;

          //printing for error debugging. Delete later
          //a=sum(q.subvec(0,j));
          //b=sum(q.subvec(0,j+1))-1;

          //Rcpp::Rcout << "max difference:" << arma::max(arma::max(arma::abs(K_js_current(j)-generate_C(A_stars(j),rhos_current.subvec(a,b),lambdas_old(1+j),n)))) << std::endl;

          // delete the above part later


          cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));
          for(int i=0;i<m;i++){
            cov=cov+(w_2s_old(i)*K_js_current(i));
          }


          neg_log_lik_current=neg_log_lik(cov,y);

          //Rcpp::Rcout << "CheckPoint 2..." << std::endl;

          log_r=(gamma_k_prime-gammas_old(k))*log( probit_prob/(1-probit_prob ) ) + neg_log_lik_old  -neg_log_lik_current;
          alpha =std::min(1.0,exp(log_r));
          u = arma::randu();
          if(u<=alpha){
            gammas_old=gammas_current;
            rhos_old=rhos_current;
            K_js_old=K_js_current;
            A_prods(j)=neg_del_Gvec_temp;
            neg_log_lik_old=neg_log_lik_current;
          }

        }



      }




    }

    //now update the parameters for last platform i.e. covariates

    j=m-1;
    for(int k=sum(q.subvec(0,j));k<sum(q.subvec(0,j+1));k++){

      k_star=k-sum(q.subvec(0,j));
      gamma_k_prime=1.0;  //actually double prime
      rho_k_prime=arma::randu();
      neg_del_Gvec_temp= A_prods(j) + ( A_stars(j).col(k_star)*std::log(rho_k_prime)-A_stars(j).col(k_star)*std::log(rhos_old(k)) );
      neg_del_Gvec=(1.0/lambdas_old(1+j))*arma::exp(neg_del_Gvec_temp);
      arma::mat delta(n,n);
      delta.row(0)=neg_del_Gvec.subvec(0,n-1).t();
      for(int i=0; i<n; i++){
        a=(i*n)- (((i-1)*i)/2);
        b=((i+1)*n) - ((i*(i+1))/2) -1;
        delta.row(i).subvec(i,n-1)=neg_del_Gvec.subvec(a,b).t();
      }
      delta=arma::symmatu(delta);
      K_js_current=K_js_old;
      K_js_current(j)=delta;
      gammas_current=gammas_old;
      rhos_current=rhos_old;
      gammas_current(k)=gamma_k_prime;
      rhos_current(k)=rho_k_prime;



      cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));
      for(int i=0;i<m;i++){
        cov=cov+(w_2s_old(i)*K_js_current(i));
      }


      neg_log_lik_current=neg_log_lik(cov,y);


      log_r= neg_log_lik_old  -neg_log_lik_current;
      alpha =std::min(1.0,exp(log_r));
      u = arma::randu();
      if(u<=alpha){
        gammas_old=gammas_current;
        rhos_old=rhos_current;
        K_js_old=K_js_current;
        A_prods(j)=neg_del_Gvec_temp;
        neg_log_lik_old=neg_log_lik_current;
      }



    }


    //update lambdas
    //lambda_0
    if(icept_update==true){
      if(prop_scheme==1){
        c_lambdas(0)=1.0;
        lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(0)*1.0, lambdas_old(0)/c_lambdas(0)   ) );
      }
      if(prop_scheme==2){
        if( (t>200)    ){

          if((t % 50) == 0){
            if(acceptance(0,0)/t< 0.30){
              c_lambdas(0)=c_lambdas(0)+1;

            }
            if(acceptance(0,0)/t> 0.40){
              c_lambdas(0)=std::max(c_lambdas(0)-1,1.0);

            }

          }


          lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(0)*1.0, lambdas_old(0)/c_lambdas(0)   ) );

        } else{
          c_lambdas(0)=1.0;
          lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(0)*1.0, lambdas_old(0)/c_lambdas(0)   ) );
        }

      }

      lambdas_current=lambdas_old;
      lambdas_current(0)=lambda_prime;
      cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_current(0));
      for(int i=0;i<m;i++){
        cov=cov+(w_2s_old(i)*K_js_old(i));
      }
      neg_log_lik_current=neg_log_lik(cov,y);


      log_r=log_gamma_pdf(lambdas_old(0),c_lambdas(0), c_lambdas(0)/lambda_prime  )-log_gamma_pdf(lambda_prime, c_lambdas(0), c_lambdas(0)/lambdas_old(0))+ (param_lambda(0)-1)*log(lambda_prime/lambdas_old(0))+param_lambda(1)*(lambdas_old(0)-lambda_prime) + neg_log_lik_old-neg_log_lik_current;



      alpha=std::min(1.0,exp(log_r));
      u=arma::randu();
      if(u<= alpha){
        lambdas_old=lambdas_current;
        neg_log_lik_old=neg_log_lik_current;
        acceptance(0,0)=acceptance(0,0)+1;

      }


    }


    //update lambda_1,....lambda_m or just lambda_1

    if( same_lambda == false   ){

      for(int j=0; j<m;j++){

        if(prop_scheme==1){
          c_lambdas(1+j)=1.0;
          lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(1+j)*1.0, lambdas_old(1+j)/c_lambdas(1+j)   ) );
        }

        if(prop_scheme==2){
          if( (t>200)    ){

            if((t % 50) == 0){
              if(acceptance(0,1+j)/t< 0.20){
                c_lambdas(1+j)=c_lambdas(1+j)+1;

              }
              if(acceptance(0,1+j)/t> 0.30){
                c_lambdas(1+j)=std::max(c_lambdas(1+j)-1,1.0);

              }

            }


            lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(1+j)*1.0, lambdas_old(1+j)/c_lambdas(1+j)   ) );

          } else{
            c_lambdas(1+j)=1.0;
            lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(1+j)*1.0, lambdas_old(1+j)/c_lambdas(1+j)   ) );
          }

        }


        lambdas_current=lambdas_old;
        lambdas_current(1+j)=lambda_prime;
        K_js_current=K_js_old;

        cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_current(0));

        K_js_current(j)=(lambdas_old(1+j)/lambda_prime)*(K_js_old(j));

        for(int i=0;i<m;i++){
          cov=cov+(w_2s_old(i)*K_js_current(i));
        }

        neg_log_lik_current=neg_log_lik(cov,y);
        log_r=log_gamma_pdf(lambdas_old(1+j),c_lambdas(1+j), c_lambdas(1+j)/lambda_prime  )-log_gamma_pdf(lambda_prime, c_lambdas(1+j), c_lambdas(1+j)/lambdas_old(1+j))+ (param_lambda(2)-1)*log(lambda_prime/lambdas_old(1+j))+param_lambda(3)*(lambdas_old(1+j)-lambda_prime) + neg_log_lik_old-neg_log_lik_current;
        alpha=std::min(1.0,exp(log_r));

        u=arma::randu();
        if(u<= alpha){
          lambdas_old=lambdas_current;
          K_js_old=K_js_current;
          neg_log_lik_old=neg_log_lik_current;
          acceptance(0,1+j)=acceptance(0,1+j)+1;

        }

      }


    } else{

      int j=0;

      if(prop_scheme==1){
        c_lambdas(1+j)=1.0;
        lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(1+j)*1.0, lambdas_old(1+j)/c_lambdas(1+j)   ) );
      }
      if(prop_scheme==2){
        if( (t>200)    ){
          if((t % 50) == 0){

            if(acceptance(0,1)/t< 0.20){
              c_lambdas(1)=c_lambdas(1)+1;

            }
            if(acceptance(0,1)/t> 0.30){
              c_lambdas(1)=std::max(c_lambdas(1)-1,1.0);

            }


          }

          lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(1+j)*1.0, lambdas_old(1+j)/c_lambdas(1+j)   ) );

        } else{
          c_lambdas(1+j)=1.0;
          lambda_prime=arma::randg<double>( arma::distr_param( c_lambdas(1+j)*1.0, lambdas_old(1+j)/c_lambdas(1+j)   ) );


        }



      }



      lambdas_current=lambdas_old;
      for(arma::uword i=0;i<m;i++){
        lambdas_current(1+i)=lambda_prime;  //note under same_lambda=true, lambda_1=lambda_2,lambda_3 etc., so the proposal

      }

      K_js_current=K_js_old;

      cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_current(0));

      for(arma::uword i=0;i<m;i++){
        K_js_current(i)=(lambdas_old(1+i)/lambda_prime)*(K_js_old(i));
        //one lambda, so need to update all the parts of K_js
      }

      for(int i=0;i<m;i++){
        cov=cov+(w_2s_old(i)*K_js_current(i));
      }

      neg_log_lik_current=neg_log_lik(cov,y);


      log_r=log_gamma_pdf(lambdas_old(1+j),c_lambdas(1+j), c_lambdas(1+j)/lambda_prime  )-log_gamma_pdf(lambda_prime, c_lambdas(1+j), c_lambdas(1+j)/lambdas_old(1+j))+ (param_lambda(2)-1)*log(lambda_prime/lambdas_old(1+j))+param_lambda(3)*(lambdas_old(1+j)-lambda_prime) + neg_log_lik_old-neg_log_lik_current;





      alpha=std::min(1.0,exp(log_r));

      u=arma::randu();
      if(u<= alpha){
        lambdas_old=lambdas_current;
        K_js_old=K_js_current;
        neg_log_lik_old=neg_log_lik_current;
        acceptance(0,1+j)=acceptance(0,1+j)+1;

      }


    }

    // update weights

    if(weight_upd==true){
      //update y (here y is transformed weights, don't be confused with data y)
      for(int i=0; i<m-1; i++){
        y_current=y_old;

        if(t>200){
          if((t % 50) == 0){
            if(acceptance(0,1+m+i)/t< 0.30){
              sigma_prop(i)=std::max(1e-6,0.9*sigma_prop(i));

            }
            if(acceptance(0,1+m+i)/t> 0.40){
              sigma_prop(i)=std::min(100.0,1.1*sigma_prop(i));

            }


          }

          y_current(i)=arma::randn<double>(arma::distr_param(y_old(i),sigma_prop(i)));

        } else{
          sigma_prop(i)=1.0;
          y_current(i)=arma::randn<double>(arma::distr_param(y_old(i),sigma_prop(i)));

        }



        w_2s_current = normalise_y(y_current);
        cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));
        for(int ii=0;ii<m;ii++){
          cov=cov+(w_2s_current(ii)*K_js_old(ii));
        }

        neg_log_lik_current=neg_log_lik(cov,y);

        //Rcpp::Rcout << "CheckPoint 5..." << std::endl;


        log_r=arma::log_normpdf(y_current(i),0.0,sigma_prior)-arma::log_normpdf(y_old(i),0.0,sigma_prior)+neg_log_lik_old-neg_log_lik_current;
        alpha=std::min(1.0,exp(log_r));
        u=arma::randu();
        if(u<= alpha){
          w_2s_old=w_2s_current;
          y_old=y_current;
          neg_log_lik_old=neg_log_lik_current;

          acceptance(0,1+m+i)=acceptance(0,1+m+i)+1;

        }
      }

    }


    if( same_lambda == false   ){
      lambda_inv=arma::pow(lambdas_old.subvec(1,m), -1.0);
      w_2s_old_1=lambda_inv/arma::accu(lambda_inv);

    } else{
      w_2s_old_1=w_2s_old;
    }


    //update r

    if(prop_scheme==1){
      c_prec=1.0;
      r_current=arma::randg<double>( arma::distr_param( c_prec*1.0, r_old/c_prec    ) );

    }
    if( prop_scheme==2  ){
      if(t>200){
        if((t % 50) == 0){
          if(acceptance(0,1+m+m)/t< 0.20){
            c_prec=c_prec+1;

          }
          if(acceptance(0,1+m+m)/t> 0.40){
            c_prec=std::max(c_prec-1,1.0);

          }





        }


        r_current=arma::randg<double>( arma::distr_param( c_prec*1.0, r_old/c_prec    ) );


      } else{
        c_prec=1.0;
        r_current=arma::randg<double>( arma::distr_param( c_prec*1.0, r_old/c_prec    ) );

      }






    }



    cov=(arma::eye(n,n)/r_current)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));
    for(int i=0;i<m;i++){
      cov=cov+(w_2s_old(i)*K_js_old(i));
    }
    neg_log_lik_current=neg_log_lik(cov,y);

    log_r=log_gamma_pdf(r_old, c_prec,c_prec/r_current)-log_gamma_pdf(r_current,c_prec,c_prec/r_old)+ (param_r(0)-1)*log(r_current/r_old)+param_r(1)*(r_old-r_current) + neg_log_lik_old-neg_log_lik_current;
    alpha=std::min(1.0,exp(log_r));
    u=arma::randu();
    if(u<= alpha){
      r_old=r_current;
      neg_log_lik_old=neg_log_lik_current;
      acceptance(0,1+m+m)=acceptance(0,1+m+m)+1;

    }

    //update censored observations

    if(total_cens !=0){
      cov=(arma::eye(n,n)/r_old)+(arma::ones<arma::mat>(n,n)/lambdas_old(0));
      for(int i=0;i<m;i++){
        cov=cov+(w_2s_old(i)*K_js_old(i));
      }
      arma::mat cov_inv=arma::inv_sympd(cov, arma::inv_opts::allow_approx);

      for(int i=0;i<total_cens;i++){
        dependent_ind(0)=cens_ind(i);
        temp=all_ind;
        temp.shed_row(cens_ind(i));
        given_ind =temp;
        //
        B=cov.submat(dependent_ind, dependent_ind);
        C=cov.submat(dependent_ind, given_ind);
        F11_inv=cov_inv.submat(given_ind,given_ind);
        d=cov_inv(cens_ind(i),cens_ind(i));
        u2=-(1.0/d)*cov_inv.submat(given_ind,dependent_ind);
        F11_inv=F11_inv-d*u2*u2.t() ;
        cMu= C* F11_inv * ( y.elem(given_ind)         );
        cVar=B- C*F11_inv * C.t();
        //cond_mean_var = conDMVN(arma::zeros<arma::vec>(n),cov,dependent_ind,given_ind,y.elem(given_ind));
        trun_sample_y=r_truncnorm(cMu(0),std::sqrt(cVar(0,0)),y_original(cens_ind(i)),std::numeric_limits<double>::infinity());

        y(cens_ind(i))=trun_sample_y;
      }




    }

    if(t>burn_in-1){
      samples.col(t-burn_in)=arma::join_cols(arma::join_cols(gammas_old,rhos_old),arma::join_cols(lambdas_old,w_2s_old_1),r_old*arma::ones(1),y.elem(cens_ind));
      nllik(t-burn_in)=neg_log_lik_old;

    }


    //std::cout << "samples" << arma::join_cols(lambdas_old,w_2s_old_1,r_old*arma::ones(1)).t() << std::endl;
    //std::cout << "acceptance" << acceptance/t <<std::endl;


    //now just add sampling of additional new params if RN ==true

    if(RN==true){
      //update mu and r_mu
      //loop for l=1,2,,,,L

      //Rcpp::Rcout << "Starting new params mu..." << std::endl;

      for(arma::uword l=0;l<L;l++){

        Al=0.0;
        Bl=0.0;
        //loop for v
        for(arma::uword v : p_gene_ind(l) ){
          temp_inv_mu=arma::reshape(cv_inv.col(v), n, n);
          cv_inv_ones=temp_inv_mu * arma::ones(n);
          Al=Al+arma::as_scalar(r_e(v)*(arma::ones(n)).t() * cv_inv_ones);

          mean_mu_temp=X_1.col(v)-arma::as_scalar(P_way.row(v)*mu_old)*arma::ones(n)+mu_old(l)*arma::ones(n);

          Bl=Bl+r_e(v)*arma::as_scalar( mean_mu_temp.t()  *  cv_inv_ones         );

        }

        mu_old(l)=arma::randn<double>(arma::distr_param(Bl/(r_mu(l)+Al)  , std::sqrt(1.0/(r_mu(l)+Al))   ));


      }

      //update r_mu given mu_old
      for(arma::uword l=0;l<L;l++){
        r_mu(l)=arma::randg<double>(arma::distr_param(r_muL_prior_param(0)+0.5,1.0/(r_muL_prior_param(1)+0.5*mu_old(l)*mu_old(l)) ));
      }



      //update r_e given all
      for(arma::uword v=0;v<p(0);v++){
        temp_inv_mu=arma::reshape(cv_inv.col(v), n, n);
        mean_mu_temp=X_1.col(v)-arma::as_scalar(P_way.row(v)*mu_old)*arma::ones(n);

        r_e(v)=arma::randg<double>(arma::distr_param(r_ev_prior_param(0)+0.5*n, 1.0/( r_ev_prior_param(1)+0.5*( arma::as_scalar(mean_mu_temp.t() *temp_inv_mu*mean_mu_temp     )        ) ) ));


      }

      //Rcpp::Rcout << "Starting new params Z..." << std::endl;

      //update Z
      //update z_1
      R_sums(0)=arma::sum(R_old,1);
      R_sums(1)=arma::sum(R_old,0).t();

      for(arma::uword k=0;k<p(0);k++){
        if(gammas_old(k)==1.0){
          z_1(k)=r_truncnorm(nu_R(0)+R_tau_samples(0)*(R_sums(0))(k),1.0,0.0,std::numeric_limits<double>::infinity());

        } else {
          z_1(k)=r_truncnorm(nu_R(0)+R_tau_samples(0)*(R_sums(0))(k),1.0,(-1.0)*std::numeric_limits<double>::infinity(),0.0);
        }

      }

      //update z_2
      for(arma::uword k=p(0);k<(p(0)+p(1));k++){
        if(gammas_old(k)==1.0){
          z_2(k-p(0))=r_truncnorm(nu_R(1)+R_tau_samples(1)*(R_sums(1))(k-p(0)),1.0,0.0,std::numeric_limits<double>::infinity());

        } else {
          z_2(k-p(0))=r_truncnorm(nu_R(1)+R_tau_samples(1)*(R_sums(1))(k-p(0)),1.0,(-1.0)*std::numeric_limits<double>::infinity(),0.0);
        }

      }

      //Rcpp::Rcout << "Starting new params R..." << std::endl;

      //update R
      for(arma::uword v=0;v<p(0);v++){
        r_v_current=R_old.row(v).t();

        arma::uvec find_one_temp3=arma::find(r_v_current==1.0);


        arma::uvec find_one_temp4=arma::find(Cand.row(v).t()==1.0);
        //not all genes have candidate miRNAs, so for those not having propose a vector of zeros


        arma::uvec find_one_temp5=arma::find( (r_v_current==0.0) && (Cand.row(v).t()==1.0)   );


        if(find_one_temp4.n_elem==0){

          r_v_current=arma::zeros(p(1));

        } else{

          if( ((find_one_temp3.n_elem >0) && (find_one_temp3.n_elem <p(1))) && (find_one_temp5.n_elem >0 )  ){

            u=arma::randu();
            if(u<= 0.5){  //note that we are taking phi as 0.5 in stochastic search algo
              //allowing only candidate moves


              arma::uvec ind_Cand=arma::shuffle(arma::find(Cand.row(v).t()==1.0));
              ind_R= ind_Cand(0);
              //ind_R=arma::randi<arma::uword>(arma::distr_param(0,p(1)-1));
              r_v_current(ind_R)=1.0-r_v_current(ind_R);


            } else{
              arma::uvec ind_0=arma::shuffle(arma::find( (r_v_current==0.0) && (Cand.row(v).t()==1.0)   ));
              ind_0R=ind_0(0);
              arma::uvec ind_1=arma::shuffle(arma::find( (r_v_current==1.0) && (Cand.row(v).t()==1.0)   )); //the second codition is redundant here since by our design all the 1s will be always candidate 1's. nevertheless!
              ind_1R=ind_1(0);
              r_v_current.swap_rows(ind_0R,ind_1R);


            }
          } else{

            arma::uvec ind_Cand=arma::shuffle(arma::find(Cand.row(v).t()==1.0));
            ind_R= ind_Cand(0);
            r_v_current(ind_R)=1.0-r_v_current(ind_R);

          }





        }



        //metrop
        R_current=R_old;
        R_current.row(v)=r_v_current.t();
        //again differently calculate values at old and current. Not sure if parts of the log_lik can be saved
        // and re-used later, but most of the parts will not be saved as they depend on v. So, while going from v1 to v2 those parts will change

        //old

        rv_temp=R_old.row(v).t();
        //calculate row_col sums of R
        R_sums(0)=arma::sum(R_old,1);
        R_sums(1)=arma::sum(R_old,0).t();



        temp_inv_mu=arma::reshape(cv_inv.col(v), n, n);




        mean_mu_temp=X_1.col(v)-arma::as_scalar(P_way.row(v)*mu_old)*arma::ones(n);


        log_lik_R_old=0.5*arma::log_det_sympd(temp_inv_mu)-0.5* r_e(v)*(arma::as_scalar( mean_mu_temp.t()*temp_inv_mu*mean_mu_temp )  )+ arma::accu(rv_temp)*std::log(R_prior_q/(1.0-R_prior_q))+ arma::accu( arma::log_normpdf(z_1, nu_R(0)*arma::ones(p(0))+R_tau_samples(0)*R_sums(0),arma::ones(p(0))   )            ) +  arma::accu( arma::log_normpdf(z_2, nu_R(1)*arma::ones(p(1))+R_tau_samples(1)*R_sums(1),arma::ones(p(1))   )            )        ;



        //current

        rv_temp=R_current.row(v).t();
        //calculate row_col sums of R
        R_sums(0)=arma::sum(R_current,1);
        R_sums(1)=arma::sum(R_current,0).t();

        arma::uvec find_one_temp2=arma::find(rv_temp==1.0);
        if(find_one_temp2.n_elem >0 ){
          arma::mat X2_temp1=X_2;
          X2_temp1.shed_cols(arma::find(rv_temp==0.0)); //X_2_rv



          temp_inv_mu=arma::eye(n,n) - beta_prior_tau*X2_temp1*arma::inv_sympd(arma::eye(X2_temp1.n_cols,X2_temp1.n_cols)+beta_prior_tau*X2_temp1.t()*X2_temp1)*X2_temp1.t();

        } else{
          temp_inv_mu=arma::eye(n,n);
        }





        //mean_mu_temp=X_1.col(v)-arma::as_scalar(P_way.row(v)*mu_old)*arma::ones(n);


        log_lik_R_current=0.5*arma::log_det_sympd(temp_inv_mu)-0.5*r_e(v)*(arma::as_scalar( mean_mu_temp.t()*temp_inv_mu*mean_mu_temp )  )+ arma::accu(rv_temp)*std::log(R_prior_q/(1.0-R_prior_q))+ arma::accu( arma::log_normpdf(z_1, nu_R(0)*arma::ones(p(0))+R_tau_samples(0)*R_sums(0),arma::ones(p(0))   )            ) +  arma::accu( arma::log_normpdf(z_2, nu_R(1)*arma::ones(p(1))+R_tau_samples(1)*R_sums(1),arma::ones(p(1))   )            )        ;


        //metrop step
        //now metrop step
        log_r=log_lik_R_current-log_lik_R_old;
        alpha=std::min(1.0,exp(log_r));
        u=arma::randu();
        if(u<= alpha){
          R_old=R_current;
          cv_inv.col(v)=arma::vectorise(temp_inv_mu);

          //log_lik_R_old=log_lik_R_current; //to check is computation can vbe reduced by saving earlier loglik

        }


      }

      //Rcpp::Rcout << "Starting new params tau..." << std::endl;

      // update tau
      //update tau_1

      ak_0_s=z_1-nu_R(0)*arma::ones(p(0));
      bk_0_s=arma::sum(R_old,1);
      B_tau=arma::as_scalar(ak_0_s.t()*bk_0_s);
      A_tau=arma::accu(   bk_0_s.transform( [](double val) { return (val*val); } )        );
      double mm1=1.0/(A_tau+(1.0/tau_prior_param(0)));
      R_tau_samples(0)=r_truncnorm(B_tau*mm1,std::sqrt(mm1),0.0,std::numeric_limits<double>::infinity());


      //update tau_2

      ak_1_s=z_2-nu_R(1)*arma::ones(p(1));
      bk_1_s=arma::sum(R_old,0).t();
      B_tau=arma::as_scalar(ak_1_s.t()*bk_1_s);
      A_tau=arma::accu(   bk_1_s.transform( [](double val) { return (val*val); } )        );
      mm1=1.0/(A_tau+(1.0/tau_prior_param(1)));
      R_tau_samples(1)=r_truncnorm(B_tau*mm1,std::sqrt(mm1),0.0,std::numeric_limits<double>::infinity());

      //Rcpp::Rcout << "Checkpoint 2..." << std::endl;

      if(t>burn_in-1){
        samples_mu.col(t-burn_in)=mu_old;
        //samples_r_mu.col(t)=r_mu;
        //samples_z.col(t)=arma::join_cols(z_1,z_2);
        samples_tau.col(t-burn_in)=R_tau_samples;
        //samples_r_e.col(t)=r_e;
        //samples_R.slice(t)=R_old;
        R_sum_save=R_sum_save+R_old;

      }



    }


  }

  double n1 =timer.toc();

  arma::mat B_1=samples.t();
  arma::mat B_2=samples_mu.t();
  //arma::mat B_3=samples_r_mu.t();
  //arma::mat B_4=samples_z.t();
  arma::mat B_5=samples_tau.t();
  //arma::mat B_6=samples_r_e.t();



  if(RN==true){
    Rcpp::List LL(7);
    LL[0]=B_1;  //samples from outcome model
    LL[1]=acceptance/(M+1); //acceptance rates of params from outcome model
    LL[2]=n1;        //runtime
    LL[3]=nllik;  // nllik of outcome model
    LL[4]=B_2;    //samples of mu_l's
    LL[5]=B_5;   //samples of tau
    LL[6]=R_sum_save/(M-burn_in);  //posterior prob of RN

    return LL;

  } else{
    Rcpp::List LL(4);
    LL[0]=B_1;  //samples from outcome model
    LL[1]=acceptance/(M+1); //acceptance rates of params from outcome model
    LL[2]=n1;        //runtime
    LL[3]=nllik;  // nllik of outcome model
    return LL;

  }






}
