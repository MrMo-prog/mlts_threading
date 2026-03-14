functions{
  #include "functions/partial_sum.stan"
  #include "functions/function_calculate_b.stan"
}

data {
  int<lower=0> N;                 // number of observational units
  int<lower=1> G;                 // groups in total
  int<lower=1> D; 	              // number of time-varying constructs
  int<lower=1> D_cen;             // number of constructs to be mean-centered
  int<lower=1, upper=3> maxLag;   // maximum lag
  int<lower=1> N_obs; 	          // observations in total: N * TP
  int<lower=1> n_pars;            // number of parameters
  int<lower=1> n_random;          // number of random effects
  int n_fixed;
  array[1, n_fixed] int is_fixed;
  array[n_random] int is_random;  // which parameters to model person-specific
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  array[D] vector[N_obs] y; 	    // array of observations

  // handling of missing values
  int n_miss;                      // total number of missings across D
  array[D] int n_miss_D;           // missings per D
  array[D,max(n_miss_D)] int pos_miss_D; // array of missings' positions

  //censoring
  real censL_val;
  int n_censL;                     // total number of obs at LB across D
  array[D] int n_censL_D;          // obs at LB per D
  array[D,max(n_censL_D)] int pos_censL_D; // array of obs at LBs' positions
  real censR_val;
  int n_censR;                      // total number of obs at LB across D
  array[D] int n_censR_D;           // obs at LB per D
  array[D,max(n_censR_D)] int pos_censR_D; // array of obs at LBs' positions

  // model adaptions based on user inputs:
  array[D_cen] int<lower=0, upper=1> innos_rand; // 1=person specific (random), 0=fixed
  int n_innos_fix;
  array[D_cen] int innos_pos;
  array[D_cen] int innos_fix_pos;


  // - dynamic model specification per D
  array[D] int<lower=0> N_pred;   // number of predictors per dimension
  array[D, max(N_pred)] int<lower=0> D_pred;   // matrix to index predictors to use per dimension
  array[D, max(N_pred)] int<lower=0> Lag_pred; // matrix to index lag of used predictors
  array[D] int Dpos1;  // index positions of dynamic effect parameters
  array[D] int Dpos2;
  array[D,max(N_pred)] int D_pred2;    // matrix to index predictors to use per dimension
  array[D,max(N_pred)] int Lag_pred2;  // matrix to index lag of used predictors

  // - time-invariant variables:
  // covariates as predictors of random effects
  int<lower=1> n_cov;           // number of covariates - minimum of 1 for intercepts

  // group specific
  array[G] int N_G;  // number of clusters (persons) by group
  array[N] int g_id; // group index per cluster
  array[G, max(N_G)] int g_id_pos; // cluster (person) indexes by group

  array[D] int<lower=0,upper=1> is_wcen;   // parameter should be within centered = 1; should not = 0
  array[D] int<lower=0,upper=D> D_cen_pos; // pos of parameters that should be centered
  int grainsize;
}

transformed data{
  // creating pos and pos_cov for partial_sum
  array[N] int pos_start;
  array[N] int pos_end;
  array[N] int seq_N;
  int pos = 1;
  int obs_id_temp;

  for (n in 1:N){
    seq_N[n] = n;
    obs_id_temp = (N_obs_id[n]);

    pos_start[n] = pos;
    pos_end[n] = pos + obs_id_temp -1;

    pos = pos + obs_id_temp;
  }

  // creating positioning for missings and censoring in partial sum
  array[N, D] int pos_start_miss = rep_array(0, N, D);            //starting position for each person and each dimension in y_impute
  array[N, D] int pos_end_miss = rep_array(0, N, D);
  array[N, D] int seq_N_miss = rep_array(0, N, D);

  int impute_pos = 0;

  for (d in 1:D){
    for (n in 1:N){
      int counter = 0;
      int end_pos = 0;
      for (x in 1: n_miss_D[d]){
        if (pos_miss_D[d, x] >= pos_start[n] && pos_miss_D[d, x] <= pos_end[n]){
          if (counter == 0){
            pos_start_miss[n, d] = impute_pos + x;
          }
          counter = counter + 1;
          end_pos = x;
        }
      }
      if (counter > 0){
        pos_end_miss[n, d] = impute_pos + end_pos;
      }
      seq_N_miss[n, d] = counter;
    }
    impute_pos = impute_pos + n_miss_D[d];
  }
}

parameters {
  array[N] vector[n_random] b_free;      // person-specific parameters
  array[G] vector[n_fixed] b_fix;        // fixed parameters
  array[G] vector<lower=0>[n_random] sd_R;        // random effect SD
  array[G] vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  array[G] cholesky_factor_corr[n_random] L;      // cholesky factor of random effects correlation matrix
  array[G] row_vector[n_random] gammas;           // fixed effect (intercepts)
  vector[n_miss] y_impute;                        // vector to store imputed values
  vector<upper=censL_val>[n_censL] y_impute_censL;
  vector<upper=censR_val>[n_censR] y_impute_censR;
}

transformed parameters{
  array[D_cen] vector[N] sd_noise;
  matrix[N, n_pars] b;

  // transformation of log-innovation variances if modeled as cluster-specific
  b = calculate_b(N, n_pars, n_random, is_random, b_free, n_fixed, G, g_id_pos,
  N_G, b_fix, is_fixed);

  for(d in 1:D_cen){
      if (innos_rand[d] == 0){
        for(g in 1:G){
          sd_noise[d, g_id_pos[g, 1:N_G[g]]] = rep_vector(sigma[g, innos_fix_pos[d]], N_G[g]);
        }
      }
      else{
        sd_noise[d] = sqrt(exp(b[,innos_pos[d]]));  // random effect transformed from log(var) to sd for each person
      }
  }
}


model {

  array[G] matrix[n_random, n_random] SIGMA;

  for(g in 1:G){
    SIGMA[g] = diag_pre_multiply(sd_R[g], L[g]); // covariance matrix of parameters by group
  }

target += reduce_sum(
    partial_sum_log_lik,
    seq_N,
    grainsize,
    N_obs_id, g_id, b_free, gammas, SIGMA, D_cen, maxLag, D,
    is_wcen, y, pos_start, pos_end, b, D_cen_pos, N_pred,
    Lag_pred, D_pred, D_pred2, Lag_pred2, Dpos1, Dpos2, sd_noise, n_miss, n_miss_D,
    pos_miss_D, y_impute, pos_start_miss, pos_end_miss, seq_N_miss, n_censL, n_censL_D, pos_censL_D, y_impute_censL, n_censR,
    n_censR_D, pos_censR_D, y_impute_censR
  );

  for (g in 1:G){
    target += normal_lpdf(gammas[g] | 0, 10);
    target += cauchy_lpdf(sd_R[g] | 0, 2);
    target += lkj_corr_cholesky_lpdf(L[g] | 1);
  }

}

generated quantities{
  array[G] matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
    for(g in 1:G){
        bcorr[g] = multiply_lower_tri_self_transpose(L[g]);
      }
}

