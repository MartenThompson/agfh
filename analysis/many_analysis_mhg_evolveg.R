source('data_gen.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')



run_it <- function(seed, slug.data.name, dat, kern.a0, kern.a1, kern.fuzz, mhg.iter, thin) {
  save_slug <- paste0('analysis/output/', slug.data.name, '_s', seed,
                      '_a0', kern.a0, 'a1', kern.a1, 'fz', kern.fuzz,
                      'mhgiter', mhg.iter, 'thin', thin, '/')
  dir.create(file.path(save_slug))
  dir.create(file.path(save_slug, 'robjects/', fsep=''))
  dir.create(file.path(save_slug, 'rscripts/', fsep=''))
  file.copy('./analysis/many_analysis_mhg_evolveg.R', paste0(save_slug, 'rscripts/many_analysis_mhg_evolveg.R'))
  file.copy('./mcmc/mh_within_gibbs_evolveg.R', paste0(save_slug, 'rscripts/mh_within_gibbs_evolveg.R'))
 
  S.init <- seq(-10,10,length.out=1e2)
  mh.scale.thetas <- 1
  mh.scale.alphas <- list(a0=0, a1=0)
  saveRDS(mh.scale.thetas, paste0(save_slug, 'robjects/', 'mh_scale_thetas.RData'))
  saveRDS(mh.scale.alphas, paste0(save_slug, 'robjects/', 'mh_scale_alphas.RData'))

  thet.var.prior.a <- 1
  thet.var.prior.b <- 1
  
  mhg_sampler <- make_mhingibbs_sampler(X=dat$X, Y=dat$Y, D=dat$D,
                                        var_gamma_a=thet.var.prior.a, var_gamma_b=thet.var.prior.b,
                                        S=S.init,
                                        kern.a0=kern.a0,
                                        kern.a1=kern.a1,
                                        kern.fuzz=kern.fuzz)
  init <- list(beta=rep(0,p),
               theta= predict(lm(dat$Y ~ dat$X), data.frame(dat$X)), # rep(0,M),
               theta.var=1 ,
               gamma = rep(0, length(S.init)))


  set.seed(seed)
  mhg.out <- mhg_sampler(init, mhg.iter, thin, mh.scale.thetas, mh.scale.alphas)
  saveRDS(mhg.out, paste0(save_slug, 'mhg_out_goodthetastarts.RData'))
  
  
  
#  init <- list(beta=rep(0,p),
#               theta= rep(0,M),
#               theta.var=1 ,
#               gamma = rep(0, length(S.init))) 
#
#  set.seed(seed)
#  mhg.out <- mhg_sampler(init, mhg.iter, thin, mh.scale.thetas, mh.scale.alphas)
#  saveRDS(mhg.out, paste0(save_slug, 'mhg_out_0thetastarts.RData'))

  # Hier Bayes FH Sampling
  set.seed(seed)
  params.init <- list(beta=rep(0,p), 
                      theta=predict(lm(dat$Y ~ dat$X), data.frame(dat$X)),
                      theta.var=1) 
  sampler <- make_gibbs_sampler(dat$X, dat$Y, dat$D, thet.var.prior.a, thet.var.prior.b)
  gibbs.out <- sampler(params.init, mhg.iter, thin) 
  saveRDS(sampler, paste0(save_slug, 'robjects/', 'gibbs_sampler.RData'))
  saveRDS(gibbs.out, paste0(save_slug, 'robjects/', 'gibbs_out.RData'))
  
  # Linear Model 
  m <- lm(dat$Y ~ dat$X)
  saveRDS(m, paste0(save_slug, 'lm.RData'))
  
  
  # Additional Saving
  saveRDS(dat, paste0(save_slug, 'data.RData'))
  saveRDS(mhg_sampler, paste0(save_slug, 'robjects/', 'mhg_sampler.RData'))
  saveRDS(init, paste0(save_slug, 'robjects/', 'mhg_param_init.RData'))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rng.seed <- 2022
M <- 100
p <- 3
kernel.a0 <- 0.1
kernel.a1 <- 0.1
kernel.fuzz <- 1e-2
mhg.iter <- 1e3
mhg.thin <- 100


set.seed(rng.seed)
D <- rep(1/100, M)
dat <- null_gen(M=M, p=p, D=D, lambda=5)
run.name <- 'no_recthv_gkde_nullgen_D001_lamb5'
run_it(rng.seed, run.name, dat, kernel.a0, kernel.a1, kernel.fuzz, mhg.iter, mhg.thin)


set.seed(rng.seed)
D <- rep(1/100, M)
dat <- gamma_err_gen(M, p, D, 5, shape = 1/2, rate = 10)
run.name <- 'no_recthv_gkde_gamma_D001_lamb5'
run_it(rng.seed, run.name, dat, kernel.a0, kernel.a1, kernel.fuzz, mhg.iter, mhg.thin)


set.seed(rng.seed)
D <- rep(1/100, M)
dat <- beta_err_gen(M, p, D, 5, 1/8, 1/8)
run.name <- 'no_recthv_gkde_beta_D001_lamb5'
run_it(rng.seed, run.name, dat, kernel.a0, kernel.a1, kernel.fuzz, mhg.iter, mhg.thin)


set.seed(rng.seed)
D <- rep(1/100, M)
dat <- beta_err_gen(M, p, D, 5, 1/12, 1/6)
run.name <- 'no_recthv_gkde_betaassym_D001_lamb5'
run_it(rng.seed, run.name, dat, kernel.a0, kernel.a1, kernel.fuzz, mhg.iter, mhg.thin)
