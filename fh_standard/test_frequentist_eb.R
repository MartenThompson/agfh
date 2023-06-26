rm(list=ls())

source('fh_standard/frequentist_eb.R')
source('data_gen.R')

set.seed(2023)
M <- 100
p <- 3
D <- rep(0.1, M)
dat <- gamma_err_gen(M, p, D, 5, shape = 1/2, rate = 10)#null_gen(M=M, p=p, D=D, lambda=5)
RM_theta_var_moment_est(dat$X, dat$Y, dat$D)
adjpfl_thvar <- adj_profile_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
plot(seq(1e-5, 20, length.out=100), sapply(seq(1e-5, 20, length.out=100), adjpfl_thvar))
#YL_theta_var_proflik_est_opt(pfl_thvar)
theta_var_est_grid(adjpfl_thvar)

adjreml_thvar <- adj_resid_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
plot(seq(1e-5, 20, length.out=100), sapply(seq(1e-5, 20, length.out=100), adjreml_thvar))
#YL_theta_var_proflik_est_opt(pfl_thvar)
theta_var_est_grid(adjreml_thvar)


reml_thvar <- resid_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
plot(seq(1e-5, 20, length.out=100), sapply(seq(1e-5, 20, length.out=100), reml_thvar))
#YL_theta_var_proflik_est_opt(pfl_thvar)
theta_var_est_grid(reml_thvar)

ests <- data.frame(
  lamb = rep(NA, 400),
  yl_est = rep(NA, 400)
)
itr <- 1
for (i in 1:1e2) {
  # confirm no error outs
  for (lamb in c(0.5, 5)) {
    for (D.val in c(0.1, 0.01)) {
      D <- rep(D.val, M)
      dat <- null_gen(M=M, p=p, D=D, lambda=lamb)
      pfl_thvar <- adj_profile_log_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
      est <- YL_theta_var_proflik_est(pfl_thvar)
      ests[itr,] <- c(lamb, est)
      itr <- itr+1
      
      dat <- gamma_err_gen(M, p, D, lamb, shape = 1/2, rate = 10)
      pfl_thvar <- adj_profile_log_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
      est <- YL_theta_var_proflik_est(pfl_thvar)
      ests[itr,] <- c(lamb, est)
      itr <- itr+1
      
      dat <- beta_err_gen(M, p, D, lamb, 1/8, 1/8)
      pfl_thvar <- adj_profile_log_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
      est <- YL_theta_var_proflik_est(pfl_thvar)
      ests[itr,] <- c(lamb, est)
      itr <- itr+1
      
      dat <- beta_err_gen(M, p, D, lamb, 1/12, 1/6)
      pfl_thvar <- adj_profile_log_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
      est <- YL_theta_var_proflik_est(pfl_thvar)
      ests[itr,] <- c(lamb, est)
      itr <- itr+1
    }
  }
}

summary(ests)
