source('data_gen.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs.R')

save_slug <- 'analysis/output/mhg_bnorm251_M100p3_0thetastarts_1/'
dir.create(file.path(save_slug))
dir.create(file.path(save_slug, 'robjects/', fsep=''))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Data Gen ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set.seed(2022)
M <- 100
p <- 3
D <- rep(1, M)
dat <- bimodal_err_gen(M=M, p=p, D=D, psi=exp(1))
cat('True params.\nBeta: ', dat$Beta, '\nPsi: ', dat$psi, '\n')
saveRDS(dat, paste0(save_slug, 'data.RData'))
saveRDS(null_gen, paste0(save_slug, 'robjects/', 'null_data_gen.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Hier Bayes FH Sampling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
params.init <- list(beta=rep(0,p), theta=rep(0,M), theta.var=1) # note that we sample theta directly here.
sampler <- make_gibbs_sampler(dat$X, dat$Y, dat$D, 1, 1/5)
gibbs.out <- sampler(params.init, 1e4, 50)
saveRDS(sampler, paste0(save_slug, 'robjects/', 'gibbs_sampler.RData'))
saveRDS(gibbs.out, paste0(save_slug, 'robjects/', 'gibbs_out.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### MH within Gibbs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
mhg_sampler <- make_mhingibbs_sampler(dat$X, dat$Y, dat$D, 0.01, 0.01)
init <- list(beta=rep(0,p), 
             theta= rep(0,M), #predict(lm(dat$Y ~ dat$X), data.frame(dat$X)),
             theta.var=1, 
             alpha.0 = 0.5, alpha.1 = 100)

mh.scale.thetas <- 0.01
mh.scale.alphas <- c(0.01,10)
saveRDS(mh.scale.thetas, paste0(save_slug, 'mh_scale_thetas.RData'))
saveRDS(mh.scale.alphas, paste0(save_slug, 'mh_scale_alphas.RData'))

set.seed(2022)
mhg_out <- mhg_sampler(init, 1e2, 100, mh.scale.thetas, mh.scale.alphas)
saveRDS(mhg_out, paste0(save_slug, 'mhg_out_1e2_100thin.RData'))


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

