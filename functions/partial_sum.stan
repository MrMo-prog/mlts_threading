real partial_sum_log_lik(array[] int slice_N, int start, int end,
    array[] int N_obs_id,
    array[] int g_id,
    array[] vector b_free,
    array[] row_vector gammas,
    array[] matrix SIGMA,
    int D_cen,
    int maxLag,
    int D,
    array[] int is_wcen,
    array[] vector y,
    array[] int pos_start,
    array[] int pos_end,
    matrix b,
    array[] int D_cen_pos,
    array[] int N_pred,
    array[,] int Lag_pred,
    array[,] int D_pred,
    array[,] int D_pred2,
    array[,] int Lag_pred2,
    array[] int Dpos1,
    array[] int Dpos2,
    array[] vector sd_noise

    ) {

      real lp = 0;

      for(i in 1:size(slice_N)){
        int pp = slice_N[i];

        int obs_id = N_obs_id[pp]; // observations per person
        int gg_p = g_id[pp];       // group

        // level 2 prediction
        lp += multi_normal_cholesky_lpdf(b_free[pp] | to_vector(gammas[gg_p]), SIGMA[gg_p]);
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
            lp += normal_lpdf(y_cen[d, (1+maxLag):obs_id] | mus[D_cen_pos[d]], sd_noise[D_cen_pos[d],pp]);
          }
        } // end of loop over dimensions

      } // end of local calculations

    } // end of loop over subjects
  return(lp);

    }
