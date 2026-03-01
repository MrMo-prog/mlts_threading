functions{
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
  for(i in 1:n_random){
    b[,is_random[i]] = to_vector(b_free[,i]);
  }
  if(n_fixed>0){
    for(i in 1:n_fixed){
      for(g in 1:G){
        b[g_id_pos[g,1:N_G[g]],is_fixed[1,i]] = rep_vector(b_fix[g,i],N_G[g]);
      }
    }
  }

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

  for(pp in 1:N){
      int obs_id = N_obs_id[pp]; // observations per person
      int gg_p = g_id[pp];       // group

      // level 2 prediction
      target += multi_normal_cholesky_lpdf(b_free[pp] | to_vector(gammas[gg_p]), SIGMA[gg_p]);
      {
      // array of predicted values
      array[D_cen] vector[obs_id-maxLag] mus;
      // create latent mean centered versions of observations
      array[D] vector[N_obs_id[pp]] y_cen;

      // calculating y_cen --> array vector of within centered observations
      for(d in 1:D){
        if(is_wcen[d] == 1){
          y_cen[d] = y[d, pos_start[pp]: pos_end[pp]] - b[pp, D_cen_pos[d]];
        } else {
          y_cen[d] = y[d, pos_start[pp]: pos_end[pp]];
        }
      }

      for(d in 1:D){

        if(is_wcen[d] == 1){
          // build prediction matrix for specific dimensions
          int n_cols; // matrix dimensions
          n_cols = N_pred[d];
          {
            matrix[obs_id - maxLag, n_cols] b_mat; // dimension specific prediction matrix: time points * predictors
            vector[n_cols] b_use; //
            for(nd in 1:N_pred[d]){ // AR effect and CL effects
              int lag_use = Lag_pred[d, nd];
              if(D_pred2[d, nd] == -99){
                b_mat[,nd] = y_cen[D_pred[d, nd], (1+maxLag-lag_use):(obs_id-lag_use)];
              } else { // interactions between two ds
                int lag_use2 = Lag_pred2[d, nd];
                b_mat[,nd] = y_cen[D_pred[d, nd], (1+maxLag-lag_use):(obs_id-lag_use)] .*
                y_cen[D_pred2[d, nd], (1+maxLag-lag_use2):(obs_id-lag_use2)];
              }
            }
            b_use[1:N_pred[d]] = to_vector(b[pp, Dpos1[d]:Dpos2[d]]);
            mus[D_cen_pos[d]] = b_mat * b_use;
          }
          target += normal_lpdf(y_cen[d, (1+maxLag):obs_id] | mus[D_cen_pos[d]], sd_noise[D_cen_pos[d],pp]);
        }
      } // end of loop over dimensions

    } // end of local calculations

  } // end of loop over subjects

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
