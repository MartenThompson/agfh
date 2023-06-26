rm(list=ls())

mh.scale.thetas <- 1e8
n.mcmc <- 1e3
n.thin <- 1e2

setwd('~/Git/agfh/')
source('package/agfh/R/frequentist_eb.R')
source('package/agfh/R/hier_bayes.R')
source('package/agfh/R/agnostic_fh.R')
source('package/agfh/R/u_normal.R')

save.slug <- paste0('analysis/output/', 'saipe_state/', 'mh',mh.scale.thetas,'nmcmc', n.mcmc, 'thin', n.thin, '/')
dir.create(file.path(save.slug))
dir.create(file.path(save.slug, 'robjects/', fsep=''))
dir.create(file.path(save.slug, 'rscripts/', fsep=''))

file.copy('./analysis/analysis_saipe_state.R', paste0(save.slug, 'rscripts/', 'analysis_saipe_state.R'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Data Read ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat.raw <- readxl::read_xls('./data/est19us.xls', sheet = 1)
n.rows <- length(dat.raw$...2)
state.names <- dat.raw$...3[5:n.rows]
poverty.90l <- as.numeric(dat.raw$...5[5:n.rows])
poverty.50 <- as.numeric(dat.raw$...4[5:n.rows])
poverty.90h <- as.numeric(dat.raw$...6[5:n.rows])
summary(poverty.50)


#income.90l <- as.numeric(dat.raw$...23[5:n.rows])
income.50 <- as.numeric(dat.raw$...22[5:n.rows])
#income.90h <- as.numeric(dat.raw$...24[5:n.rows])


Y <- poverty.50
# 90% CI so 0.05 - 0.95 percent
D <- ((2*1.645)/(poverty.90h - poverty.90l))^2
X <- as.matrix(income.50, ncol=1)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### EBLUP Fit ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theta.var.est <- RM_theta_var_moment_est(X,Y,D)
cat('Theta var:', theta.var.est, '\n')
beta.ests.eblue <- RM_beta_eblue(X,Y,D,theta.var.est)
cat('Param   :', colnames(X),'\n')
cat('Beta est:', beta.ests.eblue)

write.csv(theta.var.est, paste0(save.slug, 'robjects/', 'rm_thvar_est.csv'), row.names = F)
write.csv(beta.ests.eblue, paste0(save.slug, 'robjects/', 'rm_beta_eblue_est.csv'), row.names = F)

# this is what we care about
theta.ests.eblup <- RM_theta_eblup(X, Y, D)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### AGFH Fit ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
thet.var.prior.a <- 1
thet.var.prior.b <- 1
M <- dim(X)[1]
p <- dim(X)[2]
S.init <- seq(-10,10,length.out=M)
mhg_sampler <- make_agfh_sampler(X=X, Y=Y, D=D, 
                                      var_gamma_a=thet.var.prior.a, var_gamma_b=thet.var.prior.b,
                                      S=S.init, 
                                      kern.a0=0.1, 
                                      kern.a1=0.1, 
                                      kern.fuzz=1e-2) # higher overall fuzz -> less sensative (pullled less) to any 1 point
init <- list(beta=rep(0,p), 
             theta= predict(lm(Y ~ as.numeric(X)), data.frame(X)), # rep(0,M), #dat$theta,
             theta.var=1 , 
             gamma = rep(0, length(S.init)))

set.seed(2022)
mhg.out <- mhg_sampler(init, n.mcmc, n.thin, mh.scale.thetas)
saveRDS(mhg.out, paste0(save.slug, 'robjects/', 'mhg_out_thin_goodthetastarts.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Hier Bayes FH Sampling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
params.init <- list(beta=rep(0,p), 
                    theta=predict(lm(Y ~ as.numeric(X)), data.frame(X)),#rep(0,M), 
                    theta.var=1) 
sampler <- make_gibbs_sampler(X, Y, D, thet.var.prior.b, thet.var.prior.b)
gibbs.out <- sampler(params.init, n.mcmc, n.thin)
saveRDS(sampler, paste0(save.slug, 'robjects/', 'gibbs_sampler.RData'))
saveRDS(gibbs.out, paste0(save.slug, 'robjects/', 'gibbs_out.RData'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Linear Model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
linmod <- lm(Y ~ as.numeric(X))
saveRDS(linmod, paste0('analysis/output/saipe_state/lm.RData'))
theta.ests.lm <- predict(linmod, data.frame(X))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Saving Preds ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

map_from_density <- function(param.ts, plot=FALSE) {
  d <- density(param.ts) 
  map <- d$x[which.max(d$y)]
  
  if (plot) {
    plot(d)
    lines(c(map,map), c(0,100))
  }
  
  return(map)
}

burn.in <- 10
M <- dim(X)[1]
theta.maps.agfh <- sapply(1:M, function(m) {map_from_density(mhg.out$param.samples.list$theta[burn.in:n.mcmc,m])})
theta.maps.gibbs <- sapply(1:M, function(m) {map_from_density(gibbs.out$param.samples.list$theta[burn.in:n.mcmc,m])})


dat <- data.frame(state=state.names, Y=Y, D=D, 
                  agfh=theta.maps.agfh, 
                  hb = theta.maps.gibbs,
                  eblup=theta.ests.eblup,
                  lm = theta.ests.lm)
dat <- cbind(dat, X)
write.csv(dat, paste0(save.slug, 'robjects/', 'dat_preds.csv'), row.names = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Testing for Normality ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
source('~/Git/agfh/package/agfh/R/u_normal.R')
#mhg.out <- readRDS(paste0(save.slug, '/robjects/mhg_out_1e2_1e2thin_goodthetastarts.RData'))
g <- mhg.out$g.final

hyptest.out <- data.frame(
  SW_p = c(NA),
  SW_stat = c(NA),
  KS_p = c(NA),
  KS_stat = c(NA),
  CM_p = c(NA),
  CM_stat = c(NA),
  AD_p = c(NA),
  AD_stat = c(NA)
)

test <- test_u_normal(g$S, 'SW')
hyptest.out$SW_p[1] <- test$p.value
hyptest.out$SW_stat[1] <- test$statistic
test <- test_u_normal(g$S, 'KS')
hyptest.out$KS_p[1] <- test$p.value
hyptest.out$KS_stat[1] <- test$statistic
test <- test_u_normal(g$S, 'CM')
hyptest.out$CM_p[1] <- test$p.value
hyptest.out$CM_stat[1] <- test$statistic
test <- test_u_normal(g$S, 'AD')
hyptest.out$AD_p[1] <- test$p.value
hyptest.out$AD_stat[1] <- test$statistic
hyptest.out
write.csv(hyptest.out, paste0(save.slug, '/robjects/hyptest.csv'), row.names=FALSE)
