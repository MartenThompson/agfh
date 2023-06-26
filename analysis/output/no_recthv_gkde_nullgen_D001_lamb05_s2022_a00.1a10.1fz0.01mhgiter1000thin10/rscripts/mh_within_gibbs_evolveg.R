# Here model used exp(psi) := lambda

#source('data_gen.R')
#source('mcmc/mcmc_hybrid_fapprox.R')

#set.seed(2022)
#M <- 100
#p <- 3
#dat <- null_gen(M=M, p=p)


beta_marginal <- function(X, theta.est, theta.var.est) {
  M <- dim(X)[1]
  theta.est <- matrix(theta.est, nrow=M, ncol=1)
  
  mu <- solve(t(X)%*%X)%*%t(X)%*%theta.est # TODO: this is correct? not Y
  sigma <- solve(t(X)%*%X)*theta.var.est # TODO: had theta.var.est in solve(...)
  mvtnorm::rmvnorm(1, mean=mu, sigma=sigma) # 1xp vector
}

#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#beta_marginal(dat$X, theta.ests, 1)

# blows up when (theta.est - X%*%beta.est)^2 blows up.
# ie bad theta/beta alignment.
theta_var_marginal_maker <- function(prior.gamma.shape.a, prior.gamma.rate.b, X) {
  M <- dim(X)[1]
  a <- prior.gamma.shape.a
  b <- prior.gamma.rate.b
  
  theta_var_marginal <- function(X, Y, D, beta.est, theta.est) {
#    u.hat<- sqrt(D)*(Y - theta.est)
#    u.reconstructed <- (u.hat - mean(u.hat))/sd(u.hat)
#    theta.est <- Y-(1/sqrt(D))*u.reconstructed
    
    beta.est <- matrix(beta.est, ncol=1)
    theta.est <- matrix(theta.est, nrow=M, ncol=1)
    
    shape <- a + M/2
    rate <- b + 0.5*sum((theta.est - X%*%beta.est)^2)
    prec <- rgamma(1, shape=shape, rate=rate)
    
    return(list(theta_var=1/prec,
                theta_recon=theta.est))
  }
  
  return(theta_var_marginal)
}

#th_var_marg <- theta_var_marginal_maker(.1, .1, dat$X)
#th_var_marg(dat$X, c(3,1,-1), theta.ests)

mh_scalar_proposal <- function(current, scale) {
  return(rnorm(1, mean=current, sd=sqrt(scale)))
}

theta_scalar_llik <- function(theta, x, y, d, beta, theta.var, g) {
  x < matrix(x, ncol=1)

  llik <- -0.5*d*(y-theta)^2 + g$g(sqrt(d)*(y-theta)) - 0.5*(1/theta.var)*(theta - beta%*%x)^2
  return(llik)
}

#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#theta_scalar_llik(theta.ests[1], dat$X[1,], dat$Y[1], dat$D[1], c(3,1,-3), 1, 0.5, 1)
#theta_scalar_llik(-theta.ests[1], dat$X[1,], dat$Y[1], dat$D[1], c(3,1,-3), 1, 0.5, 1)

theta_mh <- function(X, Y, D, beta, theta.current, theta.var, g, mh.scale.th) {
  M <- length(theta.current)
  
  for (m in 1:M) {
    current <- theta.current[m]
    llik.current <- theta_scalar_llik(current, X[m,], Y[m], D[m], beta, theta.var, g)
    proposed <- mh_scalar_proposal(current, mh.scale.th)
    llik.proposed <- theta_scalar_llik(proposed, X[m,], Y[m], D[m], beta, theta.var, g)
    
    mh_rat <- max(0, exp(llik.proposed - llik.current))
    r <- runif(1)
    
    if (mh_rat > r) {
      theta.current[m] <- proposed
    } 
  }
  
  return(theta.current)
}

#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#plot(theta.ests, theta_mh(dat$X, dat$Y, dat$D, c(3,1,-3), theta.ests, 1, 0.5, 1, 1))


rbf_kernel <- function(X, alpha.0, alpha.1, phi) {
  K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
  diag(K) <- diag(K) + phi
  #image(K)
  return(K)
}

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


init_g <- function(gamma, S, a0, a1, phi) {
  Sigma.virtual <- rbf_kernel(S, a0, a1, phi)
  Sigma.virtual.inv <- solve(Sigma.virtual)
  n <- length(S)
  
  # g may get scalars (e.g. theta_llik, or vectors (e.g. integrate U dens, sample new g))
  g <- function(x) {
    Sigma.new.virtual.full <- rbf_kernel(c(x, S), a0, a1, phi)
    m <- length(x)
    n <- length(S)
    Sigma.new.virtual <- Sigma.new.virtual.full[1:m,(m+1):(m+n)]
    post.mean <- Sigma.new.virtual %*% Sigma.virtual.inv %*% gamma
    return(post.mean)
  }
  
  mass_u <- tryCatch({
   mass_u <- integrate(function(u){exp(-u^2/2 + g(u))}, -100, 100)
  }, error=function(cond) {
    cat('Bad Normalizing Integral for g.\n')
    mass_u <- list(value=NA)
  })

  #plot(S, gamma, xlim=c(-5,5), ylim = c(-1,1))
  #lines(seq(-5,5,length.out=1e2), sapply(seq(-5,5,length.out=1e2), function(x){g(x)}), col='green')
  
  return(list(
    g = g,#_scaled,
    gamma = gamma,
    K = Sigma.virtual, 
    K.inv = Sigma.virtual.inv,
    S = S,
    a0 = a0,
    a1 = a1,
    phi = phi,
    e.eta = as.numeric(mass_u$value)
  ))
}

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


g_update <- function(Y, D, theta, g.current) {
  u.hat<- sqrt(D)*(Y - theta)
  u.reconstructed <- (u.hat - mean(u.hat))/sd(u.hat)
  
  kde.u <- ks::kde(u.reconstructed)
  kde.est <- predict(kde.u, x = u.reconstructed)
  g.hat <- log(sqrt(2*pi)*kde.est) + (u.reconstructed^2)/2
  g.curr <- init_g(g.hat, u.reconstructed, g.current$a0, g.current$a1, g.current$phi)
  return(g.curr)
  
  # g.prop <- sample_new_g(g.curr, mh.scales$a0, mh.scales$a1)
  # mh_rat <- exp(sum(g.prop$g(sqrt(D)*(Y-theta)) - g.current$g(sqrt(D)*(Y-theta))))#*(g_prior(g.prop)/g_prior(g.current))
  # 
  # r <- runif(1)
  # if (mh_rat < r) {
  #   return(g.current)
  # } else {
  #   return(g.prop)
  # }
}


make_mhingibbs_sampler <- function(X, Y, D, var_gamma_a=0.01, var_gamma_b=0.01, 
                                   S=seq(-50,50,length.out=50), kern.a0=0.5, kern.a1=1, kern.fuzz=0.1) {
  
  theta_var_marginal <- theta_var_marginal_maker(var_gamma_a, var_gamma_b, X)
  
  # n.thin like batch length e.g. 10 means take every tenth
  mh_in_gibbs_sampler <- function(params.init, n.samples, n.thin, mh.scale.thetas, mh.scale.alphas) {
    beta.samples <- matrix(NA, n.samples, length(params.init$beta))
    theta.samples <- matrix(NA, n.samples, length(params.init$theta))
    theta.reconst.samples <- matrix(NA, n.samples, length(params.init$theta))
    theta.var.samples <- matrix(NA, n.samples, 1)
    gamma.samples <- matrix(NA, n.samples, length(params.init$gamma))
    e.eta.samples <- rep(NA, n.samples)
    kern.a0.samples <- rep(NA, n.samples)
    kern.a1.samples <- rep(NA, n.samples)
    
    beta.samples[1,] <- params.init$beta
    theta.samples[1,] <- params.init$theta
    theta.var.samples[1,] <- params.init$theta.var
    gamma.samples[1,] <- params.init$gamma
    kern.a0.samples[1] <- kern.a0
    kern.a1.samples[1] <- kern.a1
    
    if(is.null(params.init$theta) || is.null(params.init$theta.var) || is.null(params.init$gamma)) {
      stop('some param inits NULL.')      
    }
    
    theta.init <- params.init$theta
    theta.var.init <- params.init$theta.var
    g.init <- init_g(params.init$gamma, S, kern.a0, kern.a1, kern.fuzz) 
    e.eta.samples[1] <- g.init$e.eta
    
    beta.new <- beta_marginal(X, theta.init, theta.var.init)
    theta.new <- theta_mh(X, Y, D, beta.new, theta.init, theta.var.init, g.init, mh.scale.thetas)
    temp <- theta_var_marginal(X, Y, D, beta.new, theta.new)
    theta.var.new <- temp$theta_var
    theta.recon.new <- temp$theta_recon
    g.new <- g_update(Y, D, theta.new, g.init)
    
    beta.samples[2,] <- beta.new
    theta.samples[2,] <- theta.new
    theta.reconst.samples[2,] <- theta.recon.new
    theta.var.samples[2,] <- theta.var.new
    gamma.samples[2,] <- g.new$gamma
    kern.a0.samples[2] <- g.new$a0
    kern.a1.samples[2] <- g.new$a1
    e.eta.samples[2] <- g.new$e.eta
    
    theta.accpt <- rep(NA, n.thin*(n.samples-2))
    alphas.accpt <- matrix(NA, n.thin*(n.samples-2), 2)
    overall.iter <- 0
    
    for(i in 3:n.samples) {
      cat('sample: ', i, '\n')
      for (j in 1:n.thin) {
        theta.hold <- theta.new # just for accept/reject ratio
        a0.hold <- g.new$a0
        a1.hold <- g.new$a1
        overall.iter <- overall.iter + 1
        
        beta.new <- beta_marginal(X, theta.new, theta.var.new)
        theta.new <- theta_mh(X, Y, D, beta.new, theta.new, theta.var.new, g.new, mh.scale.thetas)
        temp <- theta_var_marginal(X, Y, D, beta.new, theta.new)  
        theta.var.new <- temp$theta_var
        theta.recon.new <- temp$theta_recon
        g.new <- g_update(Y, D, theta.new, g.new)
        
        theta.accpt[overall.iter] <- mean(theta.hold != theta.new)
        alphas.accpt[overall.iter,] <- c(a0.hold != g.new$a0, a1.hold != g.new$a1)
      }
      
      beta.samples[i,] <- beta.new
      theta.samples[i,] <- theta.new
      theta.reconst.samples[i,] <- theta.recon.new
      theta.var.samples[i,] <- theta.var.new
      gamma.samples[i,] <- g.new$gamma
      
      kern.a0.samples[i] <- g.new$a0
      kern.a1.samples[i] <- g.new$a1
      e.eta.samples[i] <- g.new$e.eta
    }
    
    param.samples.list <- list(
      beta = beta.samples,
      theta = theta.samples,
      theta.recon = theta.reconst.samples,
      theta.var = theta.var.samples,
      gamma = gamma.samples,
      kern.a0 = kern.a0.samples,
      kern.a1 = kern.a1.samples,
      e.eta = e.eta.samples
    ) 
    
    #theta.accept = mean(theta.accpt)

    output = list(
      param.init = params.init,
      S.init = S,
      mh.scale.thetas = mh.scale.thetas,
      mh.scale.alphas = mh.scale.alphas,
      n.samples = n.samples,
      n.thin = n.thin,
      #param.samples = param.samples,
      param.samples.list = param.samples.list,
      kern.a0=kern.a0, 
      kern.a1=kern.a1, 
      kern.fuzz=kern.fuzz,
      g.final =g.new,
      theta.avg.accept = mean(theta.accpt),
      alphas.accept = list(a0=mean(alphas.accpt[,1]),
                           a1=mean(alphas.accpt[,2]))
    )
    
    return(output)
  }
  
  
  return(mh_in_gibbs_sampler)  
}
