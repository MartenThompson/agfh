setwd('~/Git/agfh/')
source('fh_standard/frequentist_eb.R')
source('mcmc/mcmc_hybrid_fapprox.R')

slug <- 'noscalelik_M100p3_6'

dat <- readRDS(paste0('analysis/output/', slug, '/data.RData'))
mcmc.out <- readRDS(paste0('analysis/output/', slug, '/mcmc_out_1e3-100_beta0init.RData'))
mcmc.out.nog <- readRDS(paste0('analysis/output/', slug, '/mcmc_out_nog_1e3-100_beta0init.RData'))
linmod <- readRDS(paste0('analysis/output/', slug, '/lm.RData'))
gibbs.out <- readRDS(paste0('analysis/output/', slug, '/robjects/gibbs_out.RData'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Traceplots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## AGFH ##
par(mfrow=c(2,3))
titles <- c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(psi), expression(alpha[0]), expression(alpha[1]))
for (i in 1:ncol(mcmc.out$batch)) {
  plot(ts(mcmc.out$batch[,i]), main=paste0(titles[i], ' Traceplot'), ylab='')
}

## AGFH No G ##
par(mfrow=c(2,2))
for (i in 1:ncol(mcmc.out.nog$batch)) {
  plot(ts(mcmc.out.nog$batch[,i]), main=paste0(titles[i], ' Traceplot'), ylab='')
}

## HB ##
for (i in 1:ncol(gibbs.out$param.samples.list$beta)) {
  plot(ts(gibbs.out$param.samples.list$beta[,i]), main=paste0('HB ', titles[i], ' Traceplot'), ylab='')
}

plot(ts(gibbs.out$param.samples.list$theta.var), )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### ACF Plots #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
acf(mcmc.out$batch[,1:3])
acf(mcmc.out$batch[,4:6])

acf(mcmc.out.nog$batch)

acf(gibbs.out$param.samples.list$beta)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Beta Densities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(1,1))
dens_plot <- function(mcmc.ts, mcmc.ts.nog, gibbs.ts, lm.est, eb.est, true.val, title, ymax=0.2) {
  plot(density(mcmc.ts), col='green', main=title, xlab='param', lwd=2, ylim=c(0,ymax))
  lines(density(mcmc.ts.nog), col='blue', lwd=2)
  lines(density(gibbs.ts), col='gray', lwd=2)
  lines(c(eb.est, eb.est), c(0,100), col='orange', lwd=4)
  lines(c(true.val, true.val), c(0,100), col='black', lwd=2)
  
  if (!is.na(lm.est)) {
    lines(c(lm.est, lm.est), c(0,100), col='red', lwd=2)
    legend(min(mcmc.ts), ymax, 
         legend=c('AGFH w/ G', 'AGFH w/o G', 'HB FH', 'LM Est', 'Freq/EB', 'True Val'),
         col=c('green', 'blue', 'gray', 'red', 'orange', 'black'), lty=1, lwd=2)
  } else {
    legend(min(mcmc.ts), ymax, 
           legend=c('AGFH w/ G', 'AGFH w/o G', 'HB FH', 'Freq/EB', 'True Val'),
           col=c('green', 'blue', 'gray', 'orange', 'black'), lty=1, lwd=2)
  }
}



beta.est.eb <- RM_beta_eblue(dat$X, dat$Y, dat$D, theta_var_est = RM_theta_var_moment_est(dat$X, dat$Y, dat$D))
for (i in 1:3) {
  gibbs.burn.in <- 4
  n.gibbs <- nrow(gibbs.out$param.samples.list$beta)
  
  #png(paste0('analysis/output/M100p3_1/images/', 'beta', i, '_dens.png'), height=6, width=8, units='in', res = 100)
  dens_plot(mcmc.out$batch[,i], mcmc.out.nog$batch[,i], 
            gibbs.out$param.samples.list$beta[4:n.gibbs,i],
            linmod$coefficients[i+1], beta.est.eb[i],
            dat$Beta[i], titles[i], 3)
  #dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Psi Densities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# NOTE: remember to exponentiate agfh ests.

gibbs.burn.in <- 4
#png(paste0('analysis/output/M100p3_1/images/', 'psi_dens.png'), height=6, width=8, units='in', res = 100)
n.gibbs <- nrow(gibbs.out$param.samples.list$beta)
theta.var.est.eb <- RM_theta_var_moment_est(dat$X, dat$Y, dat$D)
dens_plot(c(0, mcmc.out$batch[,4]), mcmc.out.nog$batch[,4], 
          gibbs.out$param.samples.list$theta.var[4:n.gibbs],
          NA, theta.var.est.eb, exp(dat$psi), 'Theta Variance (exp(psi))', 2)
#dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Theta Ests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
theta.ests.lm <- predict(linmod, data.frame(dat$X))
theta.ests.freqeb <- RM_theta_eblup(dat$X, dat$Y, dat$D)

map_from_density <- function(param.ts, plot=FALSE) {
  d <- density(param.ts) 
  map <- d$x[which.max(d$y)]
  
  if (plot) {
    plot(d)
    lines(c(map,map), c(0,100))
  }
  
  return(map)
}

M <- dim(dat$X)[1]
theta.samples.agfh <- monte_carlo_posterior_theta(1e3, dat$X, mcmc.out$batch, FALSE, 2022)
theta.maps.agfh <- sapply(1:M, function(m) {map_from_density(theta.samples.agfh[,m])})

# 'quasi' FH in the sense that we used AGFH w/ g \equiv 0
theta.samples.agfhnog <- monte_carlo_posterior_theta(1e3, dat$X, mcmc.out.nog$batch, FALSE, 2022)
theta.maps.agfhnog <- sapply(1:M, function(m) {map_from_density(theta.samples.agfhnog[,m])})

theta.maps.gibbs <- sapply(1:M, function(m) {map_from_density(gibbs.out$param.samples.list$theta[,m])})


#png(paste0('analysis/output/M100p3_1/images/', 'theta_ests.png'), height=6, width=8, units='in', res = 100)

plot(dat$theta, theta.ests.lm, ylim=c(-35,35), xlim=c(-35,35), pch=16, col='red',
     xlab=expression(theta[true]), ylab=expression(theta[est]))
lines(c(-100,100), c(-100,100))
points(dat$theta, theta.ests.freqeb, pch=16, col='orange')
points(dat$theta, theta.maps.agfhnog, pch=16, col='blue')
points(dat$theta, theta.maps.gibbs, pch=16, col='gray')
points(dat$theta, theta.maps.agfh, pch=16, col='green')
legend(-35, 35, legend=c('AGFH MAP', 'AGFH w/o G MAP', 'FH HB MAP', 'LM', 'Freq/EB'), col=c('green', 'blue', 'gray', 'red', 'orange'), pch=16)

#dev.off()

# HB and EB nearly identical
mean((theta.maps.gibbs- theta.ests.freqeb)^2)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Tables ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(xtable)
mse <- function(x,y) {
  mean((x-y)^2)
}

mses <- data.frame(
  agfh   =mse(dat$theta, theta.maps.agfh),
  agfhnog=mse(dat$theta, theta.maps.agfhnog),
  hb     =mse(dat$theta, theta.maps.gibbs),
  lm     =mse(dat$theta, theta.ests.lm),
  freqeb =mse(dat$theta, theta.ests.freqeb)
)
xtable(mses)
# compare to post variance exp psi/(D + exp psi)
