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
}

parameters {
  array[N] vector[n_random] b_free;      // person-specific parameters
  array[G] vector[n_fixed] b_fix;        // fixed parameters
  array[G] vector<lower=0>[n_random] sd_R;        // random effect SD
  array[G] vector<lower=0>[n_innos_fix] sigma;    // SDs of fixed innovation variances
  array[G] cholesky_factor_corr[n_random] L;      // cholesky factor of random effects correlation matrix
  array[G] row_vector[n_random] gammas;           // fixed effect (intercepts)
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
    Lag_pred, D_pred, D_pred2, Lag_pred2, Dpos1, Dpos2, sd_noise
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

