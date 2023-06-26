source('data_gen.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')



run_it <- function(dat_maker, reps, slug.data.name, dat, kern.a0, kern.a1, kern.fuzz, mhg.iter, thin) {
  save_slug <- paste0('analysis/output/', slug.data.name, '_r', reps,
                      '_a0', kern.a0, 'a1', kern.a1, 'fz', kern.fuzz,
                      'mhgiter', mhg.iter, 'thin', thin, '/')
  dir.create(file.path(save_slug))
  dir.create(file.path(save_slug, 'robjects/', fsep=''))
  dir.create(file.path(save_slug, 'rscripts/', fsep=''))
  saveRDS(data.frame(reps=reps), paste0(save_slug, 'robjects/', 'reps.RData'))
  file.copy('./analysis/rep_analysis_mhg_evolveg.R', paste0(save_slug, 'rscripts/rep_analysis_mhg_evolveg.R'))
  file.copy('./mcmc/mh_within_gibbs_evolveg.R', paste0(save_slug, 'rscripts/mh_within_gibbs_evolveg.R'))
  
  S.init <- seq(-10,10,length.out=1e2)
  mh.scale.thetas <- 1
  mh.scale.alphas <- list(a0=0, a1=0)
  saveRDS(mh.scale.thetas, paste0(save_slug, 'robjects/', 'mh_scale_thetas.RData'))
  saveRDS(mh.scale.alphas, paste0(save_slug, 'robjects/', 'mh_scale_alphas.RData'))
  
  for (r in 1:reps) {
    dat <- dat_maker()
    mhg_sampler <- make_mhingibbs_sampler(X=dat$X, Y=dat$Y, D=dat$D,
                                          var_gamma_a=0.01, var_gamma_b=0.01,
                                          S=S.init,
                                          kern.a0=kern.a0,
                                          kern.a1=kern.a1,
                                          kern.fuzz=kern.fuzz)
    init <- list(beta=rep(0,p),
                 theta= predict(lm(dat$Y ~ dat$X), data.frame(dat$X)), # rep(0,M),
                 theta.var=1 ,
                 gamma = rep(0, length(S.init)))
    
    mhg.out <- mhg_sampler(init, mhg.iter, thin, mh.scale.thetas, mh.scale.alphas)
    saveRDS(mhg.out, paste0(save_slug, 'mhg_out_goodthetastarts',r,'.RData'))
    
    # Hier Bayes FH Sampling
    params.init <- list(beta=rep(0,p), 
                        theta=predict(lm(dat$Y ~ dat$X), data.frame(dat$X)),
                        theta.var=1) 
    sampler <- make_gibbs_sampler(dat$X, dat$Y, dat$D, 1, 1/5)
    gibbs.out <- sampler(params.init, 1e3, 50)
    saveRDS(sampler, paste0(save_slug, 'robjects/', 'gibbs_sampler',r,'.RData'))
    saveRDS(gibbs.out, paste0(save_slug, 'robjects/', 'gibbs_out',r,'.RData'))
    
    # Linear Model 
    #m <- lm(dat$Y ~ dat$X)
    #saveRDS(m, paste0(save_slug, 'lm',r,'.RData'))
    
    
    # Additional Saving
    saveRDS(dat, paste0(save_slug, 'data',r,'.RData'))
    saveRDS(mhg_sampler, paste0(save_slug, 'robjects/', 'mhg_sampler', r, '.RData'))
    saveRDS(init, paste0(save_slug, 'robjects/', 'mhg_param_init', r, '.RData'))  
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rng.seed <- 2023
M <- 100
p <- 3
kernel.a0 <- 0.1
kernel.a1 <- 0.1
kernel.fuzz <- 1e-2
mhg.iter <- 1e3
mhg.thin <- 10
reps <- 30


set.seed(rng.seed)
D <- rep(1/10, M)
data_maker <- function() {
#  null_gen(M=M, p=p, D=D, lambda=0.5)
#  gamma_err_gen(M, p, D, 0.5, shape = 1/2, rate = 10)
#  beta_err_gen(M, p, D, 0.5, 1/8, 1/8)
  beta_err_gen(M, p, D, 0.5, 1/12, 1/6)
}
run.name <- 'rep_betaassym_D01_lamb05'
run_it(data_maker, reps, run.name, dat, kernel.a0, kernel.a1, kernel.fuzz, mhg.iter, mhg.thin)

