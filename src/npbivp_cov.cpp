#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "helpers.h"

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


// [[Rcpp::export]]
Rcpp::List npbivp_cov_cpp(Rcpp::List data_types, arma::vec y, arma::mat W, arma::vec obs_indicator, int M, int burn_in, arma::vec alpha_prior_param, arma::vec lambdas_prior_param, arma::vec r_prior_param, double sigma_prior, arma::vec sigma_prop, arma::vec c_lambdas, double c_prec, arma::vec init, int prop_scheme, bool weight_upd=true, bool icept_update=true, bool same_lambda=true, bool scale_data=true){

  //icept_update is just a bool indicating whether we want to update the lambda_0 or not, if not we have to center the response y
  //scale_data is bool indicating whether we want to scale data to (0,1) range while fitting outcome GP model. Else it will take the original data (current) or scale the data by mean and std deviation
  //sigma_prop--proposal sd for transformed weights //c_lambdas--tuning for lambdas, //c_prec---tuning for precision
  //c_lambdas is of length m+1
  //lambdas_prior param is to be decided --we can give all lambdas (m+1 in total) different params
  //update:- so we are taking lambdas_prior_param=c(a1,b1,a2,b2) where the first two correspond to intercept and the last two are common for scaling factors for each data type
  // init is  a vector of initial values for common lambda, transformed weights and r (so we are using explicit weight update case)
  //prop schme is a double '1/2' etc. meanign different proposal scheme we are trying for explicit/implicit weighting scheme. Currently 1 means usual Gamma proposal, i.e. c=1
  /// 2 means keep on updating c in Gamma proposal untill target acceptance is reached. So, c changes with iteration
  // In all these cases, transformed wrights are proposed with proposal adpated
  /// according to usual Normal proposal with continuous tuning like scheme '2' described above


  //in the arguments above W means covariates (confounders)


  //we have to make several changes in the way inputs are provided before making it a package since we will want it to make more user friendly. For example for one data type, we probably won't need many parameter like sigma_prior, sigma_prop etc. Similarly we should probably just give one option to get weights:- either explitcit,i.e. w2s themselvles or via updating individual lambdas.
  // Otherwise it will make user more confused. Also note that current apporaches for showing posterior consistency would probably relie on working with individual lambdas, afterall this is just one-to-one transformation and we show identifiability in lambdas which automatically indicate identifiability of weights since weights are just reweighted version of those lambdas
  // similarly for NpBIVP same function will be able to be used for continuous as well as survival outcome. We would need some indicator for that as well, so that when 'survival' user has to provide obs_indicator as well. But if 'continuous/binary' then obs_indicator will just be a vector of 1s. Since in C++ we can't give any condition in the function argument itself, it will be better to
  // wrap this cpp function under a R function to get those functionalities. Also need to update dim and other related params as well if we just learn lambdas

  //Also, note weight_upd=F and same_lambda=T does not make any sense. So, maybe later we have to wrap everything in R and there we have to give an arguemnt 'parametrization' with two options,
  // (1) default which will be updating individual lambdas and implicitly learning weights (i.e. both weight_upd=F, same_lambda=F). Or (2) Updating weights on their own.
  //Also have to give a character vector for data transformaton :- 'No', 'min_max', 'mean_sd'



  //to ensure minimal changes in arg input from earlier code/script, so that we can just take the covariate
  // as just another platform with last element of data_types is
  // just the covariates. We just don't have to select any covariates by ensuring gamma=1
  // and sampling rho. Note that m is read later below as no of elements in list data_types.
  // We will just change the input data_types to append W at last, and remove earlier linear kernel part for W

  // append covariate as last data types. In m=1 setup, we end up getting two platforms with covariates.
  // So, all problematic stuffs with m=1 case vanishes.

  // currently we dont have any option here to omit covariate

  int m_init=data_types.size();
  Rcpp::List L_dat(1+m_init);
  for(int j=0; j<m_init; j++){
    L_dat[j]=data_types[j];

  }
  L_dat[m_init]=W;

  data_types=L_dat;




  // y means W in this case. obs_ind=observartion indicator.
  arma::vec y_original=y;
  //note that y will be updated once we sample censored values
  // so we save the original value of y to have access to the
  // lower truncation point of censored observations
  arma::uvec cens_ind=arma::find(obs_indicator==0.0); //censored indices
  int total_cens=cens_ind.n_elem; //total no of censored obs
  //arma::uvec obs_ind=arma::find(obs_indicator==1.0); //obs indices (updates:- Prolly not needed)


  arma::mat n_find=Rcpp::as<arma::mat>(data_types[0])  ;
  int n = n_find.n_rows;  // no of observation
  int m=data_types.size();  // no. of data types
  // when making Rcpparmadillo function, make it a function of a Rcpp list object (i.e. argument= Rcpp::List A) and data y(+ other params), define lists' length as m (A.size()). Each entry of the list will be data
  // type X_j, so read n from number of rows of any element of the list. Now to treat an element from the list as armadillo matrix use Rcpp::as<arma::mat>(A[i]).
  //After sourcing, supply with list with the data matrices.


  //override weight and lambda bool if m==1
  if(m==1){
    weight_upd=false;
    same_lambda=false;
  }



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

    } //else{
    //standardize by mean and variance
    //temp.each_col( [](arma::vec& a){
    //a=(a-arma::mean(a))/(arma::stddev( a ));

    // } );

    //}


    A_stars(i)=A_unique(temp,n);
    p(i)= temp.n_cols;
  }




  arma::uvec q=arma::join_cols(arma::zeros<arma::uvec>(1),p);
  int dim=2*sum(p)+1+m+m + 1+ total_cens; // total parameters in our model, here we are taking the general case, i.e. m+1 lambdas, and also we will save actual weights than the transformed ones
  arma:: mat samples(dim,M-burn_in, arma::fill::zeros); //output of Rcpp function (mcmc samples)
  arma:: mat acceptance(1,1+m+m+1,arma::fill::zeros); //output of Rcpp function (acceptance rates)
  arma::vec nllik(M-burn_in);
  int time_count; //output of Rcpp function (total time count of M samples)

  //setting of prior param--each of the following parior can be part args in Rcpp function

  arma:: vec param_r=r_prior_param; // prior param of r, (a,b) in notation
  arma:: vec param_lambda= lambdas_prior_param; // prior param of lambda_j's:-(a1,b1,a2,b2)
  arma:: vec param_bernoulli(sum(p)); //prior for gamma_k_j's
  for(int i=0;i<m;i++){
    a=sum(q.subvec(0,i));
    b=sum(q.subvec(0,i+1))-1;
    param_bernoulli.subvec(a,b)=alpha_prior_param(i)*arma::ones(p(i));
  }

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
    lambdas_old(1)=init(0);

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

  if(m>1){
    y_old=init.subvec(1,m-1); //transformed weights
    w_2s_old= normalise_y(y_old);
    y_current = arma::zeros(m-1);

    if(weight_upd==false){
      lambda_inv=arma::pow(lambdas_old.subvec(1,m), -1.0);
      w_2s_old_1=lambda_inv/arma::accu(lambda_inv);
      w_2s_old=arma::ones(m);
    }


  } else {
    w_2s_old={1.0};
  }


  double r_old=init(m); //r

  //arma::vec cens_obs=y.elem(cens_ind); //sample of censored observations

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

  //if(m>1){

  //samples.col(0)=arma::join_cols(arma::join_cols(gammas_old,rhos_old), arma::join_cols(lambdas_old,y_old),r_old*arma::ones(1),y.elem(cens_ind)); //set 0th iteration
  //} else{
  //samples.col(0)=arma::join_cols(arma::join_cols(gammas_old,rhos_old), arma::join_cols(lambdas_old,r_old*arma::ones(1)),y.elem(cens_ind)); //set 0th iteration

  //}

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


  for (int t=1; t<M; t++){
    for(int j=0;j<(m-1);j++){ //this loop will be done for actual datatypes excluding covariates
      for(int k=sum(q.subvec(0,j));k<sum(q.subvec(0,j+1));k++){
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

        log_r=(gamma_k_prime-gammas_old(k))*log( param_bernoulli(k)/(1-param_bernoulli(k) ) ) + neg_log_lik_old  -neg_log_lik_current;
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

          log_r=(gamma_k_prime-gammas_old(k))*log( param_bernoulli(k)/(1-param_bernoulli(k) ) ) + neg_log_lik_old  -neg_log_lik_current;
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
      neg_log_lik_old = neg_log_lik(cov, y);

    }

    if(t>burn_in-1){
      samples.col(t-burn_in)=arma::join_cols(arma::join_cols(gammas_old,rhos_old),arma::join_cols(lambdas_old,w_2s_old_1),r_old*arma::ones(1),y.elem(cens_ind));
      nllik(t-burn_in)=neg_log_lik_old;

    }

    //std::cout << "samples" << arma::join_cols(lambdas_old,w_2s_old_1,r_old*arma::ones(1)).t() << std::endl;
    //std::cout << "acceptance" << acceptance/t <<std::endl;





  }

  double n1 =timer.toc();

  arma::mat B_1=samples.t();
  Rcpp::List L(4);
  L[0]=B_1;
  L[1]=acceptance/(M+1);
  L[2]=n1;
  L[3]=nllik;

  return L;



}
