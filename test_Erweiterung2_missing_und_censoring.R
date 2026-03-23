library(mlts)
library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE, threads_per_chain = 2)

N = 100
TP =100
Q = 1
maxlag = 1
iterations = 3000
missings = 0.20

mod = mlts_model(q = Q, fix_inno_covs = F, inno_covs_dir = "pos", max_lag = maxlag)
simData = mlts_sim(mod, N = N, TP = TP, default = T, seed = 123)

# add some missings:

set.seed(999)
total_obs <- nrow(simData$data)
number_missings <- floor(missings*total_obs)
missing_rows <- sample(1:total_obs, size = number_missings)
simData$data$Y1[missing_rows] <- NA

# censoring
#censL_val <- quantile(simData$data$Y1, 0.15, na.rm = TRUE)
#censR_val <- quantile(simData$data$Y1, 0.85, na.rm = TRUE)

#censL_rows <- which(simData$data$Y1 <= censL_val)
#censR_rows <- which(simData$data$Y1 >= censR_val)

#simData$data$Y1[censL_rows] <- censL_val
#simData$data$Y1[censR_rows] <- censR_val


noFit = mlts_fit(mod, data = simData$data, id = "ID", ts = paste0("Y",1:Q),
                 fit_model = F)
stan_data = noFit$standata

# old model without threading:
fitted_old = stan("AR_Erweiterung2_missings_und_censoring.stan", data = stan_data,
                  pars = c("gammas", "sd_R", "bcorr"),
                  iter = iterations, chains = 2, cores = 2, seed = 1015)

## same model with treading function
stan_data$starts = array(unlist(lapply(1:N, function(x){
  cumsum(stan_data$N_obs_id)[x] - stan_data$N_obs_id[x] +1
})),dim = c(N))
# stan_data$seq_N = 1:N
stan_data$grainsize = 1

fitted = stan("AR_Erweiterung2.1_threading_missings_censoring.stan", data = stan_data,
              pars = c("gammas", "sd_R", "bcorr"),
              iter = iterations, chains = 2, cores = 2, seed = 1015)




# check elapsed times:
## by chain
apply(rstan::get_elapsed_time(fitted),1,sum)
apply(rstan::get_elapsed_time(fitted_old),1,sum)

## max time across chains
max(apply(rstan::get_elapsed_time(fitted),1,sum)) / 60
max(apply(rstan::get_elapsed_time(fitted_old),1,sum)) / 60

## relative time gain
1 - ( max(apply(rstan::get_elapsed_time(fitted),1,sum)) /
        max(apply(rstan::get_elapsed_time(fitted_old),1,sum)) )


# same parameter estimates?
summary(fitted)$summary == summary(fitted_old)$summary
# no longer the same when sampling based on y^w instead of y

# percentage of equal values by column after rounding
apply(round(summary(fitted)$summary,1) == round(summary(fitted_old)$summary,1),2,
      function(x){sum(x,na.rm = T)/sum(!is.na(x))})

# check estimates and convergence criteria
sums     = monitor(fitted,     digits_summary = 3, print = F)
sums_old = monitor(fitted_old, digits_summary = 3, print = F)

# compare model convergence indices
## Rhat
sum(sums$Rhat < sums_old$Rhat) / length(sums$Rhat)
max(sums$Rhat)
max(sums_old$Rhat)
plot(sums$Rhat, sums_old$Rhat)
abline(a = 0, b = 1)
abline(v = 1.01, h = 1.01, lty = 2)

## ESS
plot(sums$Bulk_ESS, sums_old$Bulk_ESS,
     xlim = c(0, max(c(sums$Bulk_ESS, sums_old$Bulk_ESS))),
     ylim = c(0, max(c(sums$Bulk_ESS, sums_old$Bulk_ESS))))
abline(a = 0, b = 1)
abline(v = 400, h = 400, lty = 2)

plot(sums$Tail_ESS, sums_old$Tail_ESS,
     xlim = c(0, max(c(sums$Tail_ESS, sums_old$Tail_ESS))),
     ylim = c(0, max(c(sums$Tail_ESS, sums_old$Tail_ESS))))
abline(a = 0, b = 1)
abline(v = 400, h = 400, lty = 2)

plot(sums$n_eff, sums_old$n_eff)
abline(a = 0, b = 1)
abline(v = 100, h = 100, lty = 2)


# check individual estimates based on their posterior means
means     = get_posterior_mean(fitted)
means_old = get_posterior_mean(fitted_old)

# select random mean levels
mus = means[rownames(means) %in% paste0("b_free[",1:N,",1]"),3]
mus_old = means_old[rownames(means_old) %in% paste0("b_free[",1:N,",1]"),3]


plot(mus, mus_old)
abline(a = 0, b = 1)
cor(mus, mus_old)
summary(mus)
summary(mus_old)



