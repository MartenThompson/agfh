rm(list=ls())

setwd('~/Git/agfh/')

source('~/Git/agfh/package/agfh/R/agnostic_fh.R')

run.name <- 'mh1e+08nmcmc1000thin100/'
base.path <- paste0('analysis/output/saipe_state/', run.name)
fig.path <- paste0(base.path, 'figs/')
dir.create(file.path(fig.path))
mhg.out <- readRDS(paste0(base.path, 'robjects/', 'mhg_out_thin_goodthetastarts.RData'))
mhg.out$theta.avg.accept
gibbs.out <- readRDS(paste0(base.path, 'robjects/', 'gibbs_out.RData'))
linmod <- readRDS('./analysis/output/saipe_state/lm.RData')
preds.all <- read.csv(paste0(base.path, 'robjects/', 'dat_preds.csv'))
beta.ests.eblue <- read.csv(paste0(base.path, 'robjects/rm_beta_eblue_est.csv'))

p <- dim(mhg.out$param.samples.list$beta)[2]

#png(paste0(fig.path, 'traceplots.png'), height=6, width=6, units='in', res = 100)
#par(mfrow=c(2,2))
for(i in 1:p) {
  plot(ts(mhg.out$param.samples.list$beta[,i]))
  lines(c(0,1e6), rep(beta.ests.eblue[i,],2), lwd=2, col='red')
}
#dev.off()
par(mfrow=c(1,1))

g <- mhg.out$g.final
(mass.g <- integrate(function(u) {(1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))}, -100,100, subdivisions = 1e5))
(mean.g <- integrate(function(u) {u*(1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))}, -100,100, subdivisions = 1e5))
(var.g <- integrate(function(u) {(u^2)*(1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))}, -100,100, subdivisions = 1e5))
plot(g$S, g$gamma, xlim=c(-5,5), ylim=c(min(-3, min(g$gamma)),max(3, 10)))#max(g$gamma))), pch=16)
points(mhg.out$S, sapply(mhg.out$S, function(x){g$g(x)}), col='green')
points(seq(-5,5,length.out=50), sapply(seq(-5,5,length.out=50), function(x){g$g(x)}), col='red')

g.summary <- data.frame(mass = mass.g$value, mean=mean.g$value, variance=var.g$value)
write.csv(g.summary, paste0(fig.path, 'g_summary.csv'), row.names = F)


u_dens <- function(u) {
  (1/sqrt(2*pi))*exp(-u^2/2 + g$g(u))
}

#png(paste0(fig.path, 'u_dens.png'), height=6, width=6, units='in', res = 100)
x.dense <- seq(-5,5,length.out=500)
u.dens <- sapply(x.dense, u_dens)
plot(NA, NA, xlim=c(-5,5), ylim=c(0,max(u.dens)), xlab='u', ylab='')
lines(x.dense, u.dens, lwd=2, col='green')#, main='forcing theta.var=1, bimod data, it learns!')
lines(x.dense, dnorm(x.dense, 0, sd=1), lwd=2, col='red')
legend('topleft', c('Est. U Density', 'N(0,1) Density'), col=c('green', 'red'), lty=c(1,1), lwd=3)
#dev.off()

theta.ests.eblup <- preds.all$eblup
theta.maps.agfh <- preds.all$agfh

plot(theta.ests.eblup, theta.maps.agfh)
points(theta.ests.eblup, theta.maps.agfh)

# mass slightly shifted positive -> more positive sampling errors 
# -> theta smaller than other methods' estimates
png(paste0(fig.path, 'eblup_agfh_hist.png'), height = 4, width=4, units='in', res=100)
hist(theta.ests.eblup-theta.maps.agfh, breaks=10, main='EBLUP - AGFH Differences', xlab='People') 
dev.off()
write.csv(data.frame(as.matrix(summary(theta.ests.eblup-theta.maps.agfh))), paste0(fig.path, 'eblup_agfh_summary.csv'))

plot(NA,NA, xlim=c(0,100), ylim=c(5e4, 5e6))
for (i in 1:20) {
  lines(ts(mhg.out$param.samples.list$theta[,i]))  
}



png(paste0(fig.path, 'box.png'), height=4, width=5, units='in', res=100)
boxplot(preds.all[,c('lm', 'eblup', 'hb', 'agfh')],
        names=c('LM', 'EBLUP', 'HB', 'AGFH'),
        ylab='People in Poverty', 
        width=rep(0.5, 4))
dev.off()


eblup.thvar <- read.csv(paste0(base.path, 'robjects/', 'rm_thvar_est.csv'))

png(paste0(fig.path, 'lambda_dens.png'), height=5, width=5, units='in', res=100)
plot(density(mhg.out$param.samples.list$theta.var), xlab=expression(lambda), main='Theta Variance', lwd=2, col='green')
lines(density(gibbs.out$param.samples.list$theta.var), lwd=2, col='grey')
lines(rep(eblup.thvar, 2), c(0, 1e14), lwd=2, col='orange') 
legend('topleft', legend=c('AGFH', 'HB', 'Freq/EB'), lwd=2, col=c('green', 'grey', 'orange'))
dev.off()

png(paste0(fig.path, 'beta_dens.png'), height=5, width=5, units='in', res=100)
plot(density(mhg.out$param.samples.list$beta), lwd=2, main=expression(beta), col='green')
lines(density(gibbs.out$param.samples.list$beta), lwd=2, col='grey')
lines(rep(beta.ests.eblue,2), c(0, 10), lwd=2, col='orange')
lines(rep(linmod$coefficients[2],2), c(0,10), lwd=2, col='red')
legend('topleft', legend=c('AGFH', 'HB', 'Freq/EB', 'LM'), lwd=2, col=c('green', 'grey', 'orange', 'red'))
dev.off()
