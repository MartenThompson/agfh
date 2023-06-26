rm(list=ls())

setwd('~/Git/agfh/')
source('fh_standard/frequentist_eb.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')
library(ggplot2)

slug <- 'oospred_rep_beta_D01_lamb05_r30_a00.1a10.1fz0.01mhgiter1000thin10' # UPDATE BELOW TOO
reps <- readRDS(paste0('analysis/output/', slug, '/robjects/reps.RData'))$reps

fig.path <- paste0('analysis/output/',slug,'/figs/')
dir.create(fig.path)

M <- 100
n.sample.mcmc.est <- 100


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Training Error ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

agfh.means.mses <- rep(NA, reps)
agfh.medians.mses <- rep(NA, reps)
agfh.maps.mses <- rep(NA, reps)
hb.maps.mses <- rep(NA, reps)
lm.mses <- rep(NA, reps)
eb.mses <- rep(NA, reps)

for (r in 1:reps) {
  dat.train <- readRDS(paste0('analysis/output/', slug, '/robjects/data_train', r, '.RData'))
  mhg.out <- readRDS(paste0('analysis/output/', slug, '/robjects/mhg_out_goodthetastarts', r, '.RData'))
  linmod <- lm(dat.train$Y ~ dat.train$X -1)
  gibbs.out <- readRDS(paste0('analysis/output/', slug, '/robjects/gibbs_out', r, '.RData'))
  
  theta.train.ests.lm <- linmod$fitted.values
  theta.train.ests.freqeb <- RM_theta_eblup(dat.train$X, dat.train$Y, dat.train$D)
  
  
  n.mcmc <- dim(mhg.out$param.samples.list$theta)[1]
  
  M <- dim(dat.train$X)[1]
  theta.train.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta[(n.mcmc-n.sample.mcmc.est):n.mcmc,m])})
  theta.train.means.agfh <- apply(mhg.out$param.samples.list$theta, 2, mean)
  theta.train.medians.agfh <- apply(mhg.out$param.samples.list$theta, 2, median)
  theta.train.maps.gibbs <- sapply(1:M, function(m) {map_from_density(gibbs.out$param.samples.list$theta[,m])})

  
  agfh.means.mses[r] <- mse(dat.train$theta, theta.train.means.agfh)
  agfh.medians.mses[r] <- mse(dat.train$theta, theta.train.medians.agfh)
  agfh.maps.mses[r] <- mse(dat.train$theta, theta.train.maps.agfh)
  hb.maps.mses[r] <- mse(dat.train$theta, theta.train.maps.gibbs)
  lm.mses[r] <- mse(dat.train$theta, theta.train.ests.lm)
  eb.mses[r] <- mse(dat.train$theta, theta.train.ests.freqeb)
}

output <- data.frame(
  agfh_means = agfh.means.mses,
  agfh_medians = agfh.medians.mses,
  agfh_maps = agfh.maps.mses, 
  hb = hb.maps.mses,
  lm = lm.mses,
  eb = eb.mses
)
summary(output)
write.csv(output, paste0(fig.path, 'mse_train.csv'), row.names = FALSE)


plot.data <- data.frame(
  MSE = c(agfh.means.mses, agfh.medians.mses, agfh.maps.mses, hb.maps.mses, lm.mses, eb.mses),
  Method = rep(c('AGFH Mean', 'AGFH Median', 'AGFH MAP', 'HB MAP', 'LM', 'Freq/EB'), each=reps)
)

png(paste0(fig.path, 'mse_comp_train.png'), height=5, width=5, units='in', res = 100)
ggplot(data=plot.data[plot.data$Method != 'AGFH Median',]) +
  geom_density(aes(x=MSE, fill=Method), alpha=0.5, color=NA) + 
  scale_fill_manual(values=c('red', 'orange', 'green', 'blue', 'purple')) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlim(min(plot.data$MSE)-0.1*diff(range(plot.data$MSE)), max(plot.data$MSE)+0.1*diff(range(plot.data$MSE))) +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 2))
dev.off()


png(paste0(fig.path, 'mse_boxplot_train.png'), height=5, width=5, units='in', res = 100)
ggplot(data=plot.data) +
  geom_boxplot(aes(x=Method, y=MSE)) +
  xlab(element_blank()) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, vjust=0.5)) 
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Test Error ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
rm(list=ls())

setwd('~/Git/agfh/')
source('fh_standard/frequentist_eb.R')
source('fh_standard/hier_bayes.R')
source('mcmc/mh_within_gibbs_evolveg.R')
library(ggplot2)

slug <- 'oospred_rep_beta_D01_lamb05_r30_a00.1a10.1fz0.01mhgiter1000thin10'
reps <- readRDS(paste0('analysis/output/', slug, '/robjects/reps.RData'))$reps

fig.path <- paste0('analysis/output/',slug,'/figs/')

M <- 100
n.sample.mcmc.est <- 100

agfh.means.mses <- rep(NA, reps)
agfh.medians.mses <- rep(NA, reps)
agfh.maps.mses <- rep(NA, reps)
hb.maps.mses <- rep(NA, reps)
lm.mses <- rep(NA, reps)
eb.mses <- rep(NA, reps)

for (r in 1:reps) {
  dat.train <- readRDS(paste0('analysis/output/', slug, '/robjects/data_train', r, '.RData'))
  dat.test <- readRDS(paste0('analysis/output/', slug, '/robjects/data_test', r, '.RData'))
  mhg.out <- readRDS(paste0('analysis/output/', slug, '/robjects/mhg_out_goodthetastarts', r, '.RData'))
  linmod <- lm(dat.train$Y ~ dat.train$X -1)
  gibbs.out <- readRDS(paste0('analysis/output/', slug, '/robjects/gibbs_out', r, '.RData'))
  
  theta.ests.lm <- matrix(linmod$coefficients, nrow=1)%*%t(dat.test$X)
  beta.est.eb <- RM_beta_eblue(dat.train$X, dat.train$Y, dat.train$D, 
                               theta_var_est = RM_theta_var_moment_est(dat.train$X, dat.train$Y, dat.train$D))
  #cat('RM: ', beta.est.eb, '. ', RM_theta_var_moment_est(dat.train$X, dat.train$Y, dat.train$D), '. LM:', linmod$coefficients, '\n' )
  theta.ests.freqeb <- RM_theta_new_pred(dat.test$X, beta.est.eb)
  
  
  n.mcmc <- dim(mhg.out$param.samples.list$theta)[1]
  
  M <- dim(dat.train$X)[1]
  
  theta.samp.agfh <- apply(dat.test$X, 1, function(x) {agfh_theta_new_pred(x, t(mhg.out$param.samples.list$beta), mhg.out$param.samples.list$theta.var)})
  theta.samp.gibbs <- apply(dat.test$X, 1, function(x) {hb_theta_new_pred(x, t(mhg.out$param.samples.list$beta), mhg.out$param.samples.list$theta.var)})
  
  theta.maps.agfh <- sapply(1:M, function(m) {map_from_density(theta.samp.agfh[(n.mcmc-n.sample.mcmc.est):n.mcmc,m])})
  theta.means.agfh <- apply(theta.samp.agfh, 2, mean)
  theta.medians.agfh <- apply(theta.samp.agfh, 2, median)
  theta.maps.gibbs <- sapply(1:M, function(m) {map_from_density(theta.samp.gibbs[(n.mcmc-n.sample.mcmc.est):n.mcmc,m])})
  
  agfh.means.mses[r] <- mse(dat.test$theta, theta.means.agfh)
  agfh.medians.mses[r] <- mse(dat.test$theta, theta.medians.agfh)
  agfh.maps.mses[r] <- mse(dat.test$theta, theta.maps.agfh)
  hb.maps.mses[r] <- mse(dat.test$theta, theta.maps.gibbs)
  lm.mses[r] <- mse(dat.test$theta, theta.ests.lm)
  eb.mses[r] <- mse(dat.test$theta, theta.ests.freqeb)
}

output <- data.frame(
  agfh_means = agfh.means.mses,
  agfh_medians = agfh.medians.mses,
  agfh_maps = agfh.maps.mses, 
  hb = hb.maps.mses,
  lm = lm.mses,
  eb = eb.mses
)
summary(output)
write.csv(output, paste0(fig.path, 'mse_test.csv'), row.names = FALSE)


plot.data <- data.frame(
  MSE = c(agfh.means.mses, agfh.medians.mses, agfh.maps.mses, hb.maps.mses, lm.mses, eb.mses),
  Method = rep(c('AGFH Mean', 'AGFH Median', 'AGFH MAP', 'HB MAP', 'LM', 'Freq/EB'), each=reps)
)

png(paste0(fig.path, 'mse_comp_test.png'), height=5, width=5, units='in', res = 100)
ggplot(data=plot.data[plot.data$Method != 'AGFH Median',]) +
  geom_density(aes(x=MSE, fill=Method), alpha=0.5, color=NA) + 
  scale_fill_manual(values=c('red', 'orange', 'green', 'blue', 'purple')) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlim(min(plot.data$MSE)-0.1*diff(range(plot.data$MSE)), max(plot.data$MSE)+0.1*diff(range(plot.data$MSE))) +
  theme_bw() + 
  theme(panel.grid=element_blank(),
        legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 2))
dev.off() # LM and EB same


png(paste0(fig.path, 'mse_boxplot_test.png'), height=5, width=5, units='in', res = 100)
ggplot(data=plot.data) +
  geom_boxplot(aes(x=Method, y=MSE)) +
  xlab(element_blank()) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, vjust=0.5)) 
dev.off()


