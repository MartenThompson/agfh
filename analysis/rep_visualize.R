rm(list=ls())

setwd('~/Git/agfh/')
source('fh_standard/frequentist_eb.R')
source('mcmc/mcmc_hybrid_fapprox.R')
source('mcmc/mh_within_gibbs_evolveg.R')

slug <- 'rep_gamma_D001_lamb05_r30_a00.1a10.1fz0.01mhgiter1000thin10'
reps <- readRDS(paste0('analysis/output/', slug, '/robjects/reps.RData'))$reps

M <- 100
n.sample.map.est <- 100

agfh.mses <- rep(NA, reps)
hb.mses <- rep(NA, reps)
lm.mses <- rep(NA, reps)
eb.mses <- rep(NA, reps)
adj.prof.mses <- rep(NA, reps) # theta var est w/ adj prof lik
adj.reml.mses <- rep(NA, reps)# theta var est w/ adj REML
reml.mses <- rep(NA, reps)  # theta var est w/ REML

fig.path <- paste0('analysis/output/',slug,'/figs/')
dir.create(fig.path)

for (r in 1:reps) {
  cat(r, '/', reps, '\n')
  #if (7 == r) next
  dat <- readRDS(paste0('analysis/output/', slug, '/data',r,'.RData'))
  mhg.out <- readRDS(paste0('analysis/output/', slug, '/mhg_out_goodthetastarts',r,'.RData'))
  linmod <- lm(dat$Y ~ dat$X)
  gibbs.out <- readRDS(paste0('analysis/output/', slug, '/robjects/gibbs_out',r,'.RData'))
  theta.ests.lm <- predict(linmod, data.frame(dat$X))
  theta.ests.freqeb <- RM_theta_eblup(dat$X, dat$Y, dat$D)
  
  adj_pfl_thvar <- adj_profile_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
  theta.var.est.adj.prof <- theta_var_est_grid(adj_pfl_thvar)
  adj_reml_thvar <- adj_resid_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
  theta.var.est.adj.reml <- theta_var_est_grid(adj_reml_thvar)
  reml_thvar <- resid_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
  theta.var.est.reml <- theta_var_est_grid(reml_thvar)
  
  theta.ests.adj.prof <- RM_theta_eblup_varother(dat$X, dat$Y, dat$D, theta.var.est.adj.prof)
  theta.ests.adj.reml <- RM_theta_eblup_varother(dat$X, dat$Y, dat$D, theta.var.est.adj.reml)
  theta.ests.reml <- RM_theta_eblup_varother(dat$X, dat$Y, dat$D, theta.var.est.reml)
  
  n.mcmc <- dim(mhg.out$param.samples.list$theta)[1]
  
  theta.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta[(n.mcmc-n.sample.map.est):n.mcmc,m])})
  #theta.recon.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta.recon[(n.mcmc-n.sample.map.est):n.mcmc,m])})
  theta.maps.gibbs <- sapply(1:M, function(m) {map_from_density(gibbs.out$param.samples.list$theta[,m])})
  
  agfh.mses[r] <- mean((t(dat$theta)-theta.maps.agfh)^2)
  hb.mses[r] <- mean((t(dat$theta)-theta.maps.gibbs)^2)
  lm.mses[r] <- mean((t(dat$theta)-theta.ests.lm)^2)
  eb.mses[r] <- mean((dat$theta-theta.ests.freqeb)^2)
  adj.prof.mses[r] <- mean((dat$theta-theta.ests.adj.prof)^2)
  adj.reml.mses[r] <- mean((dat$theta-theta.ests.adj.reml)^2)
  reml.mses[r] <- mean((dat$theta-theta.ests.reml)^2)
}

# agfh.mses <- na.omit(agfh.mses)
# hb.mses <- na.omit(hb.mses)
# lm.mses <- na.omit(lm.mses)
# eb.mses <- na.omit(eb.mses)

output <- data.frame(agfh.mse.mean = mean(na.omit(agfh.mses)),
                     agfh.mse.sd = sd(agfh.mses),
                     hb.mse.mean = mean(na.omit(hb.mses)),
                     hb.mse.sd = sd(hb.mses),
                     lm.mse.mean = mean(na.omit(lm.mses)),
                     lm.mse.sd = sd(lm.mses),
                     eb.mse.mean = mean(na.omit(eb.mses)),
                     eb.mse.sd = sd(eb.mses),
                     adj.prof.mse.mean = mean(na.omit(adj.prof.mses)),
                     adj.prof.mse.sd = sd(adj.prof.mses),
                     adj.reml.mse.mean = mean(na.omit(adj.reml.mses)),
                     adj.reml.mse.sd = sd(adj.reml.mses),
                     reml.mse.mean = mean(na.omit(reml.mses)),
                     reml.mse.sd = sd(reml.mses))

round(output,2)


#png(paste0(fig.path, 'hb_agfh_mse_comp.png'), height=6, width=6, units='in', res = 100)
#plot(hb.mses, agfh.mses, xlim = c(min(hb.mses, agfh.mses), max(hb.mses, agfh.mses)), 
#     ylim = c(min(hb.mses, agfh.mses), max(hb.mses, agfh.mses)),
#     pch=16, ylab = 'AGFH MSE', xlab = 'HB MSE')
#lines(c(-1,100), c(-1,100))
#dev.off()

png(paste0(fig.path, 'all_agfh_mse_comp', slug, '.png'), height=6, width=6, units='in', res = 100)
plot(agfh.mses, adj.reml.mses, xlim = c(min(eb.mses, agfh.mses), max(eb.mses, agfh.mses)), 
     ylim = c(min(eb.mses, agfh.mses), max(eb.mses, agfh.mses)),
     pch=16, xlab = 'AGFH MSE', ylab = 'Other Method MSE', col='magenta')
points(agfh.mses, lm.mses, pch=16, col='blue')
points(agfh.mses, hb.mses, pch=16, col='grey')
lines(c(-1,100), c(-1,100))
legend('bottomright', legend=c('FH HB MAP', 'LM', 'Adj REML'), col=c('gray', 'blue', 'magenta'), pch=16)
dev.off()

write.csv(agfh.mses, paste0(fig.path, 'agfh_mses.csv'), row.names = FALSE)
write.csv(hb.mses, paste0(fig.path, 'hb_mses.csv'), row.names = FALSE)
write.csv(lm.mses, paste0(fig.path, 'lm_mses.csv'), row.names = FALSE)
write.csv(eb.mses, paste0(fig.path, 'eb_mses.csv'), row.names = FALSE)
write.csv(adj.prof.mses, paste0(fig.path, 'adj_prof_mses.csv'), row.names = FALSE)
write.csv(adj.reml.mses, paste0(fig.path, 'adj_reml_mses.csv'), row.names = FALSE)
write.csv(reml.mses, paste0(fig.path, 'reml_mses.csv'), row.names = FALSE)
write.csv(output, paste0(fig.path, 'mses.csv'), row.names = FALSE)
write.csv(round(output,2), paste0(fig.path, 'mses_rounded.csv'), row.names = FALSE)


# runtimes 
#runtimes <- readRDS(paste0('analysis/output/', slug, '/robjects/runtimes.RData'))
#summary(runtimes)
