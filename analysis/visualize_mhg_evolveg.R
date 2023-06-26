rm(list=ls())

setwd('~/Git/agfh/')
source('fh_standard/frequentist_eb.R')
source('mcmc/mcmc_hybrid_fapprox.R')
source('mcmc/mh_within_gibbs_evolveg.R')

slug <- 'rep_nullgen_D001_p1b1_lamb05_r30_a00.1a10.1fz0.01mhgiter1000thin10'
rep <- 1

dat <- readRDS(paste0('analysis/output/', slug, '/data', rep, '.RData'))
dat$Beta
n.beta <- length(dat$Beta)
mhg.out <- readRDS(paste0('analysis/output/', slug, '/mhg_out_goodthetastarts', rep, '.RData'))
#linmod <- readRDS(paste0('analysis/output/', slug, '/lm.RData'))
linmod <- lm(dat$Y ~ dat$X)
gibbs.out <- readRDS(paste0('analysis/output/', slug, '/robjects/gibbs_out', rep, '.RData'))

fig.path <- paste0('analysis/output/',slug,'/figs/')
dir.create(fig.path)

mhg.out$theta.avg.accept

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Traceplots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## AGFH ##
png(paste0(fig.path, 'traceplots.png'), height=4, width=6, units='in', res = 100)
par(mfrow=c(1,2))
titles <- c(expression(beta[1]), expression(beta[2]), expression(beta[3]), expression(psi), expression(alpha[0]), expression(alpha[1]))
for (i in 1:n.beta) {
  plot(ts(mhg.out$param.samples.list$beta[,i]), main=paste0(titles[i], ' Traceplot'), ylab='')
  lines(c(0,1e5), c(dat$Beta[i],dat$Beta[i]), col='red', lwd=2)
}
plot(ts(mhg.out$param.samples.list$theta.var), main='Theta Var Traceplot')
lines(c(0,1e5), rep(dat$lambda, 2), col='red', lwd=2)
dev.off()


plot(ts(mhg.out$param.samples.list$kern.a0))
plot(ts(mhg.out$param.samples.list$kern.a1))
plot(mhg.out$g.final$S,mhg.out$param.samples.list$gamma[499,])
plot(mhg.out$g.final$S,mhg.out$param.samples.list$gamma[dim(mhg.out$param.samples.list$gamma)[1],])


dim(mhg.out$param.samples.list$gamma)

#plot(ts(log(mhg.out$param.samples.list$e.eta)))

par(mfrow=c(1,1))
g <- mhg.out$g.final
(mass.g <- integrate(function(u) {(1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))}, -100,100, subdivisions = 1e5))
(mean.g <- integrate(function(u) {u*(1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))}, -100,100, subdivisions = 1e5))
(var.g <- integrate(function(u) {(u^2)*(1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))}, -100,100, subdivisions = 1e5))
plot(g$S, g$gamma, xlim=c(-5,5), ylim=c(min(-3, min(g$gamma)),max(3, max(g$gamma))), pch=16)
points(mhg.out$S, sapply(mhg.out$S, function(x){g$g(x)}), col='green')
points(seq(-5,5,length.out=50), sapply(seq(-5,5,length.out=50), function(x){g$g(x)}), col='red')

g.summary <- data.frame(mass = mass.g$value, mean=mean.g$value, variance=var.g$value)
write.csv(g.summary, paste0(fig.path, 'g_summary.csv'), row.names = F)


u_dens <- function(u) {
  (1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))
}

if(grepl('nullgen', slug)) {
  ud.xlim <- c(-4,4)
  ud.ylim <- c(0,1)
} else if(grepl('gamma', slug)) {
  ud.xlim <- c(-2,4)
  ud.ylim <- c(0,1.4)
} else if(grepl('betaassym', slug)) {
  ud.xlim <- c(-2,2)
  ud.ylim <- c(0,4)
} else if(grepl('beta', slug)) {
  ud.xlim <- c(-2,2)
  ud.ylim <- c(0,3)
} else {
  stop('Uknown file name type')
}


png(paste0(fig.path, 'u_dens.png'), height=6, width=6, units='in', res = 100)
x.dense <- seq(-12,12,length.out=500)
u.dens <- sapply(x.dense, u_dens)
hist(dat$err, freq = FALSE, xlim=ud.xlim, ylim=ud.ylim, xlab='U', main='U Density', breaks = 12)
legend(-.7,ud.ylim[2],legend=c('Empirical Density', 'Estimated Density'), fill=c('gray', 'green'))
lines(x.dense, u.dens, lwd=2, col='green')#, main='forcing theta.var=1, bimod data, it learns!')
#lines(x.dense, dnorm(x.dense, 0, sd=1), col='red')
dev.off()

#par(mfrow=c(1,3))
#plot(dat$theta, dat$Y)
#plot(x.dense, u.dens, ylim=c(0,0.4))
#lines(x.dense, dnorm(x.dense, 0, sd=1), col='red')
#plot(ts(mhg.out$param.samples.list$theta.var))


# ## AGFH No G ##
# par(mfrow=c(2,2))
# for (i in 1:ncol(mcmc.out.nog$batch)) {
#   plot(ts(mcmc.out.nog$batch[,i]), main=paste0(titles[i], ' Traceplot'), ylab='')
# }

## HB ##
#for (i in 1:ncol(gibbs.out$param.samples.list$beta)) {
#  plot(ts(gibbs.out$param.samples.list$beta[,i]), main=paste0('HB ', titles[i], ' Traceplot'), ylab='')
#}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### ACF Plots #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#acf(mcmc.out$batch[,1:3])
#acf(mcmc.out$batch[,4:6])

#acf(mcmc.out.nog$batch)

#acf(gibbs.out$param.samples.list$beta)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Beta Densities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(1,1))
dens_plot <- function(mcmc.ts, gibbs.ts, lm.est, eb.est, true.val, title, leg.pos=NA) {
  ylim <- c(0, max(density(mcmc.ts)$y, density(gibbs.ts)$y))
  xlim <- c(min(mcmc.ts, gibbs.ts, lm.est, eb.est, true.val), max(mcmc.ts, gibbs.ts, lm.est, eb.est, true.val))
  plot(density(mcmc.ts), col='green', main=title, xlab=NA, lwd=2, xlim=xlim, ylim=ylim)
  lines(density(gibbs.ts), col='gray', lwd=2)
  lines(c(eb.est, eb.est), c(0,100), col='orange', lwd=2)
  lines(c(true.val, true.val), c(0,100), col='black', lwd=2)
  
  if (!is.na(lm.est)) {
    lines(c(lm.est, lm.est), c(0,100), col='blue', lwd=2)
    if (!is.na(leg.pos)) {
      legend(leg.pos[1], leg.pos[2], 
             legend=c('AGFH', 'HB FH', 'LM Est', 'Freq/EB', 'True Val'),
             col=c('green', 'gray', 'blue', 'orange', 'black'), lty=1, lwd=3)  
    } else {
      legend('topright', 
             legend=c('AGFH', 'HB FH', 'LM Est', 'Freq/EB', 'True Val'),
             col=c('green', 'gray', 'blue', 'orange', 'black'), lty=1, lwd=3)
    }
    
  } else {
    if (!is.na(leg.pos)) {
      legend(leg.pos[1], leg.pos[2], 
             legend=c('AGFH', 'HB FH', 'Freq/EB', 'True Val'),
             col=c('green', 'gray', 'orange', 'black'), lty=1, lwd=3)  
    } else {
      legend('topright', 
             legend=c('AGFH', 'HB FH', 'Freq/EB', 'True Val'),
             col=c('green', 'gray', 'orange', 'black'), lty=1, lwd=3)
    }
    
  }
}



beta.est.eb <- RM_beta_eblue(dat$X, dat$Y, dat$D, theta_var_est = RM_theta_var_moment_est(dat$X, dat$Y, dat$D))
#xlims <- matrix(c(4,7,-4,0,2,5), ncol=2, byrow=T)
#ylims <- matrix(c(0,3, 0,3, 0,3), ncol=2, byrow = T)
for (i in 1:n.beta) {
  gibbs.burn.in <- 0
  n.gibbs <- nrow(gibbs.out$param.samples.list$theta)
  
  png(paste0(fig.path, 'beta', i, '_dens.png'), height=5, width=5, units='in', res = 100)
  dens_plot(mhg.out$param.samples.list$beta[,i], 
            gibbs.out$param.samples.list$beta[gibbs.burn.in:n.gibbs],
            linmod$coefficients[i+1], beta.est.eb[i],
            dat$Beta[i], titles[i], 'topleft')
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Psi Densities ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
theta.var.est.eb <- RM_theta_var_moment_est(dat$X, dat$Y, dat$D)
adj_pfl_thvar <- adj_profile_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
theta.var.est.adj.prof <- theta_var_est_grid(adj_pfl_thvar)
adj_reml_thvar <- adj_resid_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
theta.var.est.adj.reml <- theta_var_est_grid(adj_reml_thvar)
reml_thvar <- resid_likelihood_theta_var_maker(dat$X, dat$Y, dat$D)
theta.var.est.reml <- theta_var_est_grid(reml_thvar)

png(paste0(fig.path, 'lambda_dens.png'), height=5, width=5, units='in', res = 100)
plot(density(mhg.out$param.samples.list$theta.var), col='green', main='Theta Variance', 
     ylim=c(0, max(density(mhg.out$param.samples.list$theta.var)$y, density(gibbs.out$param.samples.list$theta.var)$y)), 
     xlim=c(-0.05, 15),#max(gibbs.out$param.samples.list$theta.var)), 
     xlab=NA, lwd=2)
lines(density(gibbs.out$param.samples.list$theta.var), col='gray', lwd=2)
lines(rep(theta.var.est.eb, 2), c(0,100), col='orange', lwd=2)
lines(rep(theta.var.est.adj.prof, 2), c(0,100), col='light blue', lwd=2)
lines(rep(theta.var.est.adj.reml, 2), c(0,100), col='magenta', lwd=2)
lines(rep(dat$lambda, 2), c(0,100), col='black', lwd=2)
legend('topright', 
           legend=c('AGFH', 'HB FH', 'Freq/EB', 'Adj Prof', 'Adj REML', 'True Val'),
           col=c('green', 'gray', 'orange', 'light blue', 'magenta', 'black'), lty=1, lwd=3)  
    
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Theta Ests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
theta.ests.lm <- predict(linmod, data.frame(dat$X))
theta.ests.freqeb <- RM_theta_eblup(dat$X, dat$Y, dat$D)
theta.ests.adj.reml <- RM_theta_eblup_varother(dat$X, dat$Y, dat$D, theta.var.est.adj.reml)

map_from_density <- function(param.ts, plot=FALSE) {
  d <- density(param.ts) 
  map <- d$x[which.max(d$y)]
  
  if (plot) {
    plot(d)
    lines(c(map,map), c(0,100))
  }
  
  return(map)
}

n.mcmc <- dim(mhg.out$param.samples.list$theta)[1]
n.sample.map.est <- 100

M <- dim(dat$X)[1]
theta.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta[(n.mcmc-n.sample.map.est):n.mcmc,m])})
#theta.recon.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta.recon[(n.mcmc-n.sample.map.est):n.mcmc,m])})
theta.maps.gibbs <- sapply(1:M, function(m) {map_from_density(gibbs.out$param.samples.list$theta[,m])})

lim.low <- min(dat$theta, theta.ests.lm, theta.ests.freqeb, theta.ests.adj.reml, theta.maps.gibbs, theta.maps.agfh)
lim.high <- max(dat$theta, theta.ests.lm, theta.ests.freqeb, theta.ests.adj.reml, theta.maps.gibbs, theta.maps.agfh)

png(paste0(fig.path, 'theta_ests.png'), height=5, width=5, units='in', res = 100)
plot(dat$theta, theta.ests.lm, pch=16, col='blue',
     xlim=c(lim.low,lim.high), ylim=c(lim.low,lim.high),
     xlab=expression(theta[true]), ylab=expression(theta[est]),
     main='Theta Estimates')
lines(c(-100,100), c(-100,100))
points(dat$theta, theta.ests.freqeb, pch=16, col='orange')
points(dat$theta, theta.ests.adj.reml, pch=16, col='magenta')
points(dat$theta, theta.maps.gibbs, pch=16, col='gray')
points(dat$theta, theta.maps.agfh, pch=16, col='green')
legend('topleft', legend=c('AGFH MAP', 'FH HB MAP', 'LM', 'Freq/EB', 'Adj REML'), col=c('green', 'gray', 'blue', 'orange', 'magenta'), pch=16)
dev.off()

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
  AGFH   =mse(dat$theta, theta.maps.agfh),
  #AGFHRecon = mse(dat$theta, theta.recon.maps.agfh),
  HB   =mse(dat$theta, theta.maps.gibbs),
  LM     =mse(dat$theta, theta.ests.lm),
  FreqEB =mse(dat$theta, theta.ests.freqeb)
)
xtable(mses)
# compare to post variance exp psi/(D + exp psi)

write.csv(mses, paste0(fig.path, 'mses.csv'), row.names = FALSE)

