#library(lineprof)
library(mcmc)
source('data_gen.R')
source('mcmc/mcmc_hybrid_fapprox.R')
source('fh_standard/hier_bayes.R')

save_slug <- 'analysis/output/noscalelik_M100p3_6/'
dir.create(file.path(save_slug))
dir.create(file.path(save_slug, 'robjects/', fsep=''))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Data Gen ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set.seed(2022)
M <- 100
p <- 3
dat <- null_gen(M=M, p=p)
cat('True params.\nBeta: ', dat$Beta, '\nPsi: ', dat$psi, '\n')
saveRDS(dat, paste0(save_slug, 'data.RData'))
saveRDS(null_gen, paste0(save_slug, 'robjects/', 'null_data_gen.RData'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Hier Bayes FH Sampling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
params.init <- list(beta=rep(0,p), theta=rep(0,M), theta.var=1) # note that we sample theta directly here.
#sampler <- make_gibbs_sampler(dat$X, dat$Y, dat$D, 1, 1/5)
#gibbs.out <- sampler(params.init, 1e4, 50)
#saveRDS(sampler, paste0(save_slug, 'robjects/', 'gibbs_sampler.RData'))
#saveRDS(gibbs.out, paste0(save_slug, 'robjects/', 'gibbs_out.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### AGFH Sampling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
S <- seq(-50, 50, length.out=50)
lup <- lu_post_maker(dat$X, dat$Y, dat$D, S)

lup.no.g <- lu_post_maker(dat$X, dat$Y, dat$D, S, include_g = FALSE)
params.init <- c(0,0,0,1,.5,100) #c(0,0,0,0,.5,100) 
params.init <- c(3.16,1.47,-3.8, 1, 0.5, 100) # beta 1, 2, 3, psi, alpha.0, alpha.1
#lup(params.init)  
#lup.no.g(params.init[1:4])

met.scale <- 0.1*c(.1,.1,.1,1, .1, 10)
cat('Scale: ', met.scale, '\n')
set.seed(2022)
mcmc.out <- metrop(lup.no.g, params.init, nbatch=1e2, nspac = 100, 
                   scale = met.scale)

mcmc.out.2 <- metrop(mcmc.out, nbatch=1e2, nspac = 100, 
                   scale = met.scale)

for (i in 1:6) {
  plot(ts(mcmc.out.2$batch[,i]))
}

#saveRDS(mcmc.out, paste0(save_slug, 'mcmc_out_1e3-100_beta0init.RData'))

mcmc.nog.out <- metrop(lup.no.g, params.init[1:4], nbatch=1e3, nspac = 100, 
                   scale = met.scale[1:4])
saveRDS(mcmc.nog.out, paste0(save_slug, 'mcmc_out_nog_1e3-100_beta0init.RData'))
saveRDS(params.init, paste0(save_slug, 'mcmc_params_init.RData'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Linear Model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
m <- lm(dat$Y ~ dat$X)
saveRDS(m, paste0(save_slug, 'lm.RData'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Additional Saving ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
saveRDS(S, paste0(save_slug, 'robjects/', 'S.RData'))
saveRDS(met.scale, paste0(save_slug, 'robjects/', 'metrop_scale.RData'))
saveRDS(lup, paste0(save_slug, 'robjects/', 'log_unnormalized_posterior.RData'))
saveRDS(lu_post_maker, paste0(save_slug, 'robjects/', 'log_unnormalized_posterior_maker.RData'))
saveRDS(lup.no.g, paste0(save_slug, 'robjects/', 'log_unnormalized_posterior_nog.RData'))
saveRDS(g_maker, paste0(save_slug, 'robjects/', 'g_maker.RData'))
#notes <- paste0('M=', M, ', p=', p, '\n1e4 batch w/ nspac=50')
#write.csv(notes, paste0(save_slug, 'notes.txt'), quote=FALSE, row.names=FALSE)
