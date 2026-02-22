functions{
  real partial_sum_log_lik(array[] int slice_N, int start, int end,
    
    array[] int starts,
    array[] int N_obs_id,
    vector y,
    array[] vector b_free,
    vector sd_noise,
    vector gammas,
    matrix SIGMA
    
    ) {
    
    real lp = 0;
    
    for(k in 1:size(slice_N)) {
      int pp = slice_N[k];
      
      int st  = starts[pp];
      int n_obs = N_obs_id[pp];

      vector[n_obs] y_cen;
      y_cen = y[st:(st+n_obs-1)] - b_free[pp,1];            
      
      vector[n_obs-1] mus;
      mus = b_free[pp,2] * y_cen[1:(n_obs-1)];
      lp += normal_lpdf(y_cen[2:n_obs] | mus, sd_noise[pp]);
      lp += multi_normal_cholesky_lpdf(b_free[pp,] | gammas, SIGMA);
      }
      
      return(lp);
      
    }
}

data {
  int<lower=0> N;
  int<lower=1> D; 	// number of time-varying constructs
  int<lower=1> N_obs; 	          // observations in total: N * TP
  array[N] int<lower=1> N_obs_id; // number of observations for each unit
  array[D] vector[N_obs] y; 	    // array of observations
  array[N] int<lower=1, upper=N_obs> starts;
  array[N] int<lower=1, upper=N> seq_N;
  int grainsize;
}

transformed data{
  int n_random;
  n_random = 3;
}

parameters {
  array[N] vector[3] b_free;      // person-specific parameter
  vector<lower=0>[3] sd_R;        // random effect SD
  cholesky_factor_corr[3] L;      // cholesky factor of random effects correlation matrix
  vector[3] gammas;           // fixed effect (intercepts)
}

transformed parameters{
  vector[N] sd_noise;
  matrix[n_random, n_random] SIGMA;

  SIGMA = diag_pre_multiply(sd_R, L);
  // transformation of log-innovation variances if modeled as cluster-specific
  sd_noise = sqrt(exp(to_vector(b_free[1:N,3])));
}


model {

  target += reduce_sum(
    partial_sum_log_lik, 
    seq_N,
    grainsize,
    starts,
    N_obs_id,
    y[1,],
    b_free,
    sd_noise,
    gammas,
    SIGMA
    );
    
  target += normal_lpdf(gammas | 0, 10);
  target += cauchy_lpdf(sd_R   | 0, 2);
  target += lkj_corr_cholesky_lpdf(L | 1);

}

generated quantities{
matrix[n_random,n_random] bcorr; // random coefficients correlation matrix
bcorr = multiply_lower_tri_self_transpose(L);
}

