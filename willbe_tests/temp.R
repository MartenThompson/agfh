# from data_gen.R

#M <- 1e5
#a <- null_gen(M, 2, rep(1/2,M), exp(1))
#mean(a$Y - a$theta)
#var(a$Y - a$theta) # should be 2 according to D's specification.
#var(a$err) # should always be 1
#summary(a$Y)
#par(mfrow=c(1,1))
#hist(a$theta)
#hist(a$Y)

#M <- 1e3
#a <- beta_err_gen(M, 1, rep(1/19,M),lambda = 1/2, 1/12, 1/8)
#plot(a$X, a$Y, col=rgb(0,0,0,0.5),pch=16)
#mean(a$Y-a$theta)
#var(a$Y-a$theta) # should be 19
#var(a$err) # should be 1

# M <- 1e3
# a <- gamma_err_gen(M, 2, rep(2,M), exp(1), 2, 10)
# plot(a$theta, a$Y)
# mean(a$Y - a$theta)
# var(a$Y - a$theta) # should be 1/2
# var(a$err) # should always be 1
# hist(a$err)




# from hier_bayes.R
#source('data_gen.R')
#set.seed(2022)
#M <- 100
#p <- 3
#dat <- null_gen(M=M, p=p)
#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#beta_marginal(dat$X, theta.ests, 1)

#theta_marginal(dat$X, dat$Y, dat$D, c(3,1,-1), 1)

#th_var_marg <- theta_var_marginal_maker_hb(100, 100, dat$X)
#th_var_marg(dat$X, c(3,1,-1), theta.ests)

#params.init <- list(beta=rep(0,p), theta=rep(0,M), theta.var=1)
#sampler <- make_gibbs_sampler(dat$X, dat$Y, dat$D, 0.001, 0.001)
#gibbs.out <- sampler(params.init, 1e3, 10)

#plot(ts(gibbs.out$param.samples.list$beta[,3]))
#plot(ts(gibbs.out$param.samples.list$theta.var))
#acf(gibbs.out$param.samples.list$beta)

# x <- matrix(1, 1, 3)
# bsamp <- matrix(c(0,0,0,0,
#                   0,1,10,-10,
#                   0,2,10,-10), nrow=3, byrow=T)
# tvsamp <- c(1e-3, 1e1, 1e3, 1)
# hb_theta_new_pred(x, bsamp, tvsamp)







# from mcmc/mh_within_gibs_evolveg aka anostic_fh.R

#set.seed(2022)
#M <- 100
#p <- 3
#dat <- null_gen(M=M, p=p)


#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#theta_scalar_llik(theta.ests[1], dat$X[1,], dat$Y[1], dat$D[1], c(3,1,-3), 1, 0.5, 1)
#theta_scalar_llik(-theta.ests[1], dat$X[1,], dat$Y[1], dat$D[1], c(3,1,-3), 1, 0.5, 1)


#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#plot(theta.ests, theta_mh(dat$X, dat$Y, dat$D, c(3,1,-3), theta.ests, 1, 0.5, 1, 1))

#S <- seq(-10,10,length.out=100)
#plot(S, mvtnorm::rmvnorm(1, c(rep(0,20),rep(2,10),rep(-1,70)), rbf_kernel(S,0.001,1,1e-10)))



# my_kernel <- function(X, alpha.0, alpha.1, phi) {
#   #K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
#   #K <- (1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2
#
#   n <- length(X) # X is Uhat and always 1 dim
#   K <- matrix(0, n,n)
#   for (i in 1:n) {
#     for (j in 1:n) {
#       K[i,j] <- alpha.0*exp((1/alpha.1)*(-(X[i]  - X[j])^2 - 0.1*X[i]^4 - 0.1*X[j]^4))
#     }
#   }
#
#
#   diag(K) <- diag(K) + phi
#   image(K)
#   return(K)
# }

#S <- seq(-10,10,length.out=100)
#K <- my_kernel(S,10,100,1e-10)
#plot(S, mvtnorm::rmvnorm(1, rep(0, length(S)), K))


#S <- seq(-10, 10, length.out=1e2)
#g <- init_g(c(rep(1,10),rep(0,50),rep(2,20), rep(-1,20)), S, 0.01, 1, 1e-6)
# integrate(function(u) {exp(-u^2/2 + g$g(u))}, -100,100)
#plot(S, g$gamma, xlim=c(-25,25), ylim=c(-3,6))
#points(S, sapply(S, function(x){g$g(x)}), col='green')
#points(seq(-25,25,length.out=500), sapply(seq(-25,25,length.out=500), function(x){g$g(x)}), col='red')


#plot(S, g$gamma, xlim=c(-55,-45))
#points(seq(-55,-45,length.out=1e3), sapply(seq(-55,-45,length.out=1e3), function(x){g$g(x)}), col='red')
# weird tail behavior caused by us just saying g(.) = 0 outside support


# rMVNormP <- function(n, mu, Prec){
#   p <- length(mu)
#   Z <- matrix(rnorm(p*n), p, n)
#   U <- chol(Prec) # By default R's chol fxn returns upper cholesky factor
#   X <- backsolve(U, Z) # more efficient and stable than acctually inverting
#   X <- t(sweep(X, 1, mu, FUN=`+`))
#   return(X)
# }


# g.curr/g.prop a bit of a misnomer: we let gamma update freely and
# really curr/prop on alpha values.
# sample_new_g <- function(g.current, a0.mh.scale=0.01, a1.mh.scale=0.01) {
#   a0.curr <- g.current$a0
#   a1.curr <- g.current$a1
#
#   a0.prop <- rnorm(1, a0.curr, sd=sqrt(a0.mh.scale))
#   a0.prop <- ifelse(a0.prop > 0, a0.prop, a0.curr)
#   a1.prop <- rnorm(1, a1.curr, sd=sqrt(a1.mh.scale))
#   a1.prop <- ifelse(a1.prop > 0, a1.prop, a1.curr)
#   phi.prop <- rnorm(1, g.current$phi, sd=sqrt(1))
#   phi.prop <- ifelse(phi.prop >0, phi.prop, g.current$phi)
#   g.prop <- init_g(g.current$gamma, g.current$S, a0.prop, a1.prop, phi.prop)
# }

# profiling
profile_mhg <- function() {
  source('data_gen.R')
  set.seed(2022)
  M <- 50
  p <- 1
  D <- rep(0.1, M)
  dat <- beta_err_gen(M, p, D, 5, 1/8, 1/8)
  thet.var.prior.a <- 1
  thet.var.prior.b <- 1
  S.init <- seq(-10,10,length.out=M)
  mhg_sampler <- make_mhingibbs_sampler(X=dat$X, Y=dat$Y, D=dat$D,
                                        var_gamma_a=thet.var.prior.a, var_gamma_b=thet.var.prior.b,
                                        S=S.init,
                                        kern.a0=0.1,
                                        kern.a1=0.1,
                                        kern.fuzz=1e-2) # higher overall fuzz -> less sensative (pullled less) to any 1 point
  init <- list(beta=rep(0,p),
               theta= predict(lm(dat$Y ~ dat$X), data.frame(dat$X)), # rep(0,M),
               theta.var=1 ,
               gamma = rep(0, length(S.init)))

  mh.scale.thetas <- 1

  set.seed(2022)
  mhg.out <- mhg_sampler(init, 3, 1, mh.scale.thetas)
}
