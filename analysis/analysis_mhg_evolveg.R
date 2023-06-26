source('data_gen.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')

save_slug <- paste0('analysis/output/', 'gkde_nullgen_D01_lamb5_s2023_a00.1a10.1fz0.01mhgiter1000thin10/')
dir.create(file.path(save_slug))
dir.create(file.path(save_slug, 'robjects/', fsep=''))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Data Gen ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set.seed(2023)
M <- 50
p <- 3
D <- rep(0.1, M)
dat <- null_gen(M=M, p=p, D=D, lambda=5)
#dat <- gamma_err_gen(M, p, D, 1/2, shape = 1/2, rate = 10)
#dat <- beta_err_gen(M, p, D, 5, 1/8, 1/8)
#dat <- beta_err_gen(M, p, D, 1/2, 1/12, 1/6)

thet.var.prior.a <- 1
thet.var.prior.b <- 1

saveRDS(dat, paste0(save_slug, 'data.RData'))

# dev.off(dev.list()["RStudioGD"])
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### MH within Gibbs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#]
n.mcmc <- 1e3
n.thin <- 10
S.init <- seq(-10,10,length.out=M)
mhg_sampler <- make_mhingibbs_sampler(X=dat$X, Y=dat$Y, D=dat$D, 
                                      var_gamma_a=thet.var.prior.a, var_gamma_b=thet.var.prior.b,
                                      S=S.init, 
                                      kern.a0=0.1, 
                                      kern.a1=0.1, 
                                      kern.fuzz=1e-2) # higher overall fuzz -> less sensitive (pulled less) to any 1 point
init <- list(beta=rep(0,p), 
             theta= predict(lm(dat$Y ~ dat$X), data.frame(dat$X)), # rep(0,M), #dat$theta,
             theta.var=1 , 
             gamma = rep(0, length(S.init)))

mh.scale.thetas <- 1

set.seed(2023)
mhg.out <- mhg_sampler(init, n.mcmc, n.thin, mh.scale.thetas)
saveRDS(mhg.out, paste0(save_slug, 'mhg_out_1e3_10thin_goodthetastarts.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Hier Bayes FH Sampling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
params.init <- list(beta=rep(0,p), 
                    theta=predict(lm(dat$Y ~ dat$X), data.frame(dat$X)),#rep(0,M), 
                    theta.var=1) 
sampler <- make_gibbs_sampler(dat$X, dat$Y, dat$D, thet.var.prior.a, thet.var.prior.b)
gibbs.out <- sampler(params.init, n.mcmc, n.thin)
saveRDS(sampler, paste0(save_slug, 'robjects/', 'gibbs_sampler.RData'))
saveRDS(gibbs.out, paste0(save_slug, 'robjects/', 'gibbs_out.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Linear Model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
m <- lm(dat$Y ~ dat$X)
saveRDS(m, paste0(save_slug, 'lm.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Additional Saving ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
saveRDS(mhg_sampler, paste0(save_slug, 'robjects/', 'mhg_sampler.RData'))
saveRDS(init, paste0(save_slug, 'robjects/', 'mhg_param_init.RData'))

