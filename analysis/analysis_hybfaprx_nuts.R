source('./nuts.R')
source('../data_gen.R')
source('./mcmc_hybrid_fapprox.R')


pd_func <- function(theta, log) {
  dnorm(1, theta, 1, log = log)
}
nuts_results <- sampler_nuts(pd_func, start = 0, iterations = 1e3, epsilon=0.1)
plot(ts(nuts_results))



set.seed(2022)
dat <- null_gen(n=100, p=3)
cat('True params.\n', 'Beta: ', dat$Beta, '\nPsi: ', dat$psi)

S <- seq(-20, 20, length.out=10)
lup <- lu_post_maker(dat$X, dat$Y, dat$D, S)
params.init <- c(0,0,0, 1, .1,1) #rep(0.1, 3+1+2) # beta 1, 2, 3, psi, alpha.0, alpha.1

lup(params.init)  
#lup(c(3.2, 1.4, -3.8,1, .1,.1))

lup_wrapped <- function(params, log_unused) {
  return(lup(params))
}

set.seed(2022)
mcmc.out <- sampler_nuts(lup_wrapped, start = params.init, iterations = 5e2, epsilon=1e-2)

for (i in 1:ncol(mcmc.out)) {
  plot(ts(mcmc.out[,i]))
}
