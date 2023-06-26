source('data_gen.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')

run_it <- function(dat_maker, M.tr, reps, slug.data.name, kern.a0, kern.a1, kern.fuzz, mhg.iter, thin) {
  save_slug <- paste0('analysis/output/', slug.data.name, '_r', reps,
                      '_a0', kern.a0, 'a1', kern.a1, 'fz', kern.fuzz,
                      'mhgiter', mhg.iter, 'thin', thin, '/')
  dir.create(file.path(save_slug))
  dir.create(file.path(save_slug, 'robjects/', fsep=''))
  dir.create(file.path(save_slug, 'rscripts/', fsep=''))
  saveRDS(data.frame(reps=reps), paste0(save_slug, 'robjects/', 'reps.RData'))
  file.copy('./analysis/rep_oospred_betaassym.R', paste0(save_slug, 'rscripts/rep_oospred_betaassym.R'))
  file.copy('./mcmc/mh_within_gibbs_evolveg.R', paste0(save_slug, 'rscripts/mh_within_gibbs_evolveg.R'))
  
  S.init <- seq(-10,10,length.out=M.tr)
  mh.scale.thetas <- 1
  
  saveRDS(mh.scale.thetas, paste0(save_slug, 'robjects/', 'mh_scale_thetas.RData'))
  
  thet.var.prior.a <- 1
  thet.var.prior.b <- 1
  
  for (r in 1:reps) {

    dat.all <- dat_maker()
    dat.train <- dat.all$train
    dat.test <- dat.all$test

    saveRDS(dat.train, paste0(save_slug, 'robjects/', 'data_train', r, '.RData'))
    saveRDS(dat.test, paste0(save_slug, 'robjects/','data_test', r, '.RData'))
    
    
    S.init <- seq(-10,10,length.out=M.train)
    mhg_sampler <- make_mhingibbs_sampler(X=dat.train$X, Y=dat.train$Y, D=dat.train$D, 
                                          var_gamma_a=thet.var.prior.a, var_gamma_b=thet.var.prior.b,
                                          S=S.init, 
                                          kern.a0=kern.a0, 
                                          kern.a1=kern.a1, 
                                          kern.fuzz=kern.fuzz) 
    init <- list(beta=rep(0,p), 
                 theta= predict(lm(dat.train$Y ~ dat.train$X), data.frame(dat.train$X)), # rep(0,M), #dat$theta,
                 theta.var=1 , 
                 gamma = rep(0, length(S.init)))
    
    mh.scale.thetas <- 1
    
    mhg.out <- mhg_sampler(init, mhg.iter, thin, mh.scale.thetas)
    saveRDS(mhg.out, paste0(save_slug, 'robjects/', 'mhg_out_goodthetastarts', r, '.RData'))
  
    
    params.init <- list(beta=rep(0,p), 
                        theta=predict(lm(dat.train$Y ~ dat.train$X), data.frame(dat.train$X)),#rep(0,M), 
                        theta.var=1) 
    sampler <- make_gibbs_sampler(dat.train$X, dat.train$Y, dat.train$D, thet.var.prior.a, thet.var.prior.b)
    gibbs.out <- sampler(params.init, mhg.iter, thin)
    saveRDS(sampler, paste0(save_slug, 'robjects/', 'gibbs_sampler', r, '.RData'))
    saveRDS(gibbs.out, paste0(save_slug, 'robjects/', 'gibbs_out', r, '.RData'))
    
    saveRDS(mhg_sampler, paste0(save_slug, 'robjects/', 'mhg_sampler', r, '.RData'))
    saveRDS(init, paste0(save_slug, 'robjects/', 'mhg_param_init', r, '.RData'))
  }
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rng.seed <- 2023
M.train <- 50
M.test <- 50
p <- 3
kernel.a0 <- 0.1
kernel.a1 <- 0.1
kernel.fuzz <- 1e-2
mhg.iter <- 1e3
mhg.thin <- 10
reps <- 30


set.seed(rng.seed)
D <- rep(1/100, M.train+M.test)
data_maker <- function() {
  split_traintest(beta_err_gen(M.train+M.test, p, D, 0.5, 1/12, 1/6), M.train, M.test)
}
run.name <- 'oospred_rep_betaassym_D001_lamb05'
run_it(data_maker, M.train, reps, run.name, kernel.a0, kernel.a1, kernel.fuzz, mhg.iter, mhg.thin)

