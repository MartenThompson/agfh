rm(list=ls())

setwd('~/Git/agfh/')
source('fh_standard/frequentist_eb.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')
library(ggplot2)

slug <- 'oospred_nullgen_D01_lamb5_s2023_a00.1a10.1fz0.01mhgiter1000thin10'

dat.train <- readRDS(paste0('analysis/output/', slug, '/data_train.RData'))
dat.test <- readRDS(paste0('analysis/output/', slug, '/data_test.RData'))
dat.train$Beta
dat.test$Beta
mhg.out <- readRDS(paste0('analysis/output/', slug, '/mhg_out_goodthetastarts.RData'))
linmod <- readRDS(paste0('analysis/output/', slug, '/lm.RData'))
gibbs.out <- readRDS(paste0('analysis/output/', slug, '/robjects/gibbs_out.RData'))

fig.path <- paste0('analysis/output/',slug,'/figs/')
dir.create(fig.path)

mhg.out$theta.avg.accept


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Training Error ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

theta.train.ests.lm <- linmod$fitted.values
theta.train.ests.freqeb <- RM_theta_eblup(dat.train$X, dat.train$Y, dat.train$D)


n.mcmc <- dim(mhg.out$param.samples.list$theta)[1]
n.sample.map.est <- 100

M <- dim(dat.train$X)[1]
theta.train.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta[(n.mcmc-n.sample.map.est):n.mcmc,m])})
theta.train.means.agfh <- apply(mhg.out$param.samples.list$theta, 2, mean)
theta.train.medians.agfh <- apply(mhg.out$param.samples.list$theta, 2, median)
theta.train.maps.gibbs <- sapply(1:M, function(m) {map_from_density(gibbs.out$param.samples.list$theta[,m])})


png(paste0(fig.path, 'theta_train_ests.png'), height=6, width=6, units='in', res = 100)
plot(dat.train$theta, theta.train.ests.lm, pch=16, col='red',
     #xlim=c(0,11), ylim=c(0,11),
     xlab=expression(theta[true]), ylab=expression(theta[est]),
     main='Theta Estimates')
lines(c(-100,100), c(-100,100))
points(dat.train$theta, theta.train.ests.freqeb, pch=16, col='orange')
points(dat.train$theta, theta.train.maps.gibbs, pch=16, col='gray')
points(dat.train$theta, theta.train.maps.agfh, pch=16, col='green')
legend('topleft', legend=c('AGFH MAP', 'FH HB MAP', 'LM', 'Freq/EB'), col=c('green', 'gray', 'red', 'orange'), pch=16)
dev.off()


mse <- function(x,y) {
  mean((x-y)^2)
}

mses.train <- data.frame(
  AGFH_mean   =mse(dat.train$theta, theta.train.means.agfh),
  AGFH_median   =mse(dat.train$theta, theta.train.medians.agfh),
  AGFH_map   =mse(dat.train$theta, theta.train.maps.agfh),
  HB   =mse(dat.train$theta, theta.train.maps.gibbs),
  LM     =mse(dat.train$theta, theta.train.ests.lm),
  FreqEB =mse(dat.train$theta, theta.train.ests.freqeb)
)
xtable::xtable(mses.train)
# compare to post variance exp psi/(D + exp psi)

write.csv(mses.train, paste0(fig.path, 'mses_train.csv'), row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Testing Error ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
theta.test.ests.lm <- matrix(linmod$coefficients, nrow=1)%*%t(dat.test$X)
beta.est.eb <- RM_beta_eblue(dat.train$X, dat.train$Y, dat.train$D, 
                             theta_var_est = RM_theta_var_moment_est(dat.train$X, dat.train$Y, dat.train$D))
theta.test.ests.freqeb <- RM_theta_new_pred(dat.test$X, beta.est.eb)


n.mcmc <- dim(mhg.out$param.samples.list$theta)[1]
n.sample.map.est <- 100

M <- dim(dat.test$X)[1]

theta.test.samp.agfh <- apply(dat.test$X, 1, function(x) {agfh_theta_new_pred(x, t(mhg.out$param.samples.list$beta), mhg.out$param.samples.list$theta.var)})
theta.test.samp.gibbs <- apply(dat.test$X, 1, function(x) {hb_theta_new_pred(x, t(mhg.out$param.samples.list$beta), mhg.out$param.samples.list$theta.var)})

saveRDS(theta.test.samp.agfh, paste0(fig.path, 'theta_test_samp_agfh.RData'))
saveRDS(theta.test.samp.gibbs, paste0(fig.path, 'theta_test_samp_gibbs.RData'))

plot.data <- data.frame(theta = c(theta.test.samp.agfh[,10], theta.test.samp.gibbs[,10]), 
                        method = rep(c('AGFH', 'HB'), each=n.mcmc))

ggplot() + 
  geom_histogram(data=plot.data, aes(theta, fill=method), position='dodge')

theta.test.maps.agfh <- sapply(1:M, function(m) {map_from_density(theta.test.samp.agfh[(n.mcmc-n.sample.map.est):n.mcmc,m])})
theta.test.means.agfh <- apply(theta.test.samp.agfh, 2, mean)
theta.test.medians.agfh <- apply(theta.test.samp.agfh, 2, median)
theta.test.maps.gibbs <- sapply(1:M, function(m) {map_from_density(theta.test.samp.gibbs[(n.mcmc-n.sample.map.est):n.mcmc,m])})


png(paste0(fig.path, 'theta_test_ests.png'), height=6, width=6, units='in', res = 100)
plot(dat.test$theta, theta.test.ests.lm, pch=16, col='red',
     #xlim=c(0,11), ylim=c(0,11),
     xlab=expression(theta[true]), ylab=expression(theta[est]),
     main='Theta Estimates')
lines(c(-100,100), c(-100,100))
points(dat.test$theta, theta.test.ests.freqeb, pch=16, col='orange')
points(dat.test$theta, theta.test.maps.gibbs, pch=16, col='gray')
points(dat.test$theta, theta.test.maps.agfh, pch=16, col='green')
legend('topleft', legend=c('AGFH MAP', 'FH HB MAP', 'LM', 'Freq/EB'), col=c('green', 'gray', 'red', 'orange'), pch=16)
dev.off()



mses.test <- data.frame(
  AGFH_mean   =mse(dat.test$theta, theta.test.means.agfh),
  AGFH_median   =mse(dat.test$theta, theta.test.medians.agfh),
  AGFH_map   =mse(dat.test$theta, theta.test.maps.agfh),
  HB   =mse(dat.test$theta, theta.test.maps.gibbs),
  LM     =mse(dat.test$theta, theta.test.ests.lm),
  FreqEB =mse(dat.test$theta, theta.test.ests.freqeb)
)
xtable::xtable(mses.test)
# compare to post variance exp psi/(D + exp psi)

write.csv(mses.test, paste0(fig.path, 'mses_test.csv'), row.names = FALSE)

