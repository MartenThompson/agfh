# Here model used exp(psi) := lambda

#source('data_gen.R')
source('mcmc/mcmc_hybrid_fapprox.R')

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


theta_var_marginal_maker <- function(prior.gamma.shape.a, prior.gamma.rate.b, X) {
  M <- dim(X)[1]
  a <- prior.gamma.shape.a
  b <- prior.gamma.rate.b
  
  theta_var_marginal <- function(X, beta.est, theta.est) {
    beta.est <- matrix(beta.est, ncol=1)
    theta.est <- matrix(theta.est, nrow=M, ncol=1)
    
    shape <- a + M/2
    rate <- b + 0.5*sum((theta.est - X%*%beta.est)^2)
    prec <- rgamma(1, shape=shape, rate=rate)
    #cat(prec, ', shape:', shape, ', rate:', rate, '\n')
    return(1/prec)
  }
  
  return(theta_var_marginal)
}

#th_var_marg <- theta_var_marginal_maker(.1, .1, dat$X)
#th_var_marg(dat$X, c(3,1,-1), theta.ests)

mh_scalar_proposal <- function(current, scale) {
  return(rnorm(1, mean=current, sd=sqrt(scale)))
}

theta_scalar_llik <- function(theta, x, y, d, beta, theta.var, alpha.0, alpha.1) {
  x < matrix(x, ncol=1)
  S <- seq(-50, 50, length.out=50) # TODO pass in 
  g <- g_maker(c(alpha.0, alpha.1), S)

  llik <- -0.5*d*(y-theta)^2 + g$g(sqrt(d)*(y-theta)) - 0.5*(1/theta.var)*(theta - beta%*%x)^2
  return(llik)
}

#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#theta_scalar_llik(theta.ests[1], dat$X[1,], dat$Y[1], dat$D[1], c(3,1,-3), 1, 0.5, 1)
#theta_scalar_llik(-theta.ests[1], dat$X[1,], dat$Y[1], dat$D[1], c(3,1,-3), 1, 0.5, 1)

theta_mh <- function(X, Y, D, beta, theta.current, theta.var, alpha.0, alpha.1, mh.scale.th) {
  M <- length(theta.current)
  
  for (m in 1:M) {
    current <- theta.current[m]
    llik.current <- theta_scalar_llik(current, X[m,], Y[m], D[m], beta, theta.var, alpha.0, alpha.1)
    proposed <- mh_scalar_proposal(current, mh.scale.th)
    llik.proposed <- theta_scalar_llik(proposed, X[m,], Y[m], D[m], beta, theta.var, alpha.0, alpha.1)
    
    mh_rat <- max(0, exp(llik.proposed - llik.current))
    r <- runif(1)
    
    #cat('lik prop:', llik.proposed, ', curr:', llik.current, ', mhratio:', mh_rat, ', r:', r, '\n')
    
    if (mh_rat > r) {
      theta.current[m] <- proposed
      #cat('accepted \n')
    } #else {
      #cat('rejected \n')
    #}
  }
  
  return(theta.current)
}

#theta.ests <- predict(lm(dat$Y ~ dat$X), data.frame(dat$X))
#plot(theta.ests, theta_mh(dat$X, dat$Y, dat$D, c(3,1,-3), theta.ests, 1, 0.5, 1, 1))

g_maker_zero <- function() {
  g <- function(x) {
    return(rep(0, length(x)))
  }
  
  return(list(g=g))
}

alphas_llik <- function(Y, D, theta, alpha.0, alpha.1) {
  alpha0_prior <- function(alpha) {
    dbeta(alpha, 1, 3)
  }
  alpha1_prior <- function(alpha) {
    dgamma(alpha, shape=2, rate=1/50)
  }
  
  if(alpha.0 <= 0 || alpha.1 <= 0) {
    return(-Inf)
  }
  
  S <- seq(-50, 50, length.out=50) # TODO pass in 
  g <- g_maker(c(alpha.0, alpha.1), S)
  #g <- g_maker_zero()
  M <- length(Y)
  
  lliks <- sapply(1:M, function(m) {
    -0.5*D[m]*(Y[m] - theta[m])^2 + g$g(sqrt(D[m])*(Y[m]-theta[m]))
  })
  
  #plot(S, g$gamma)
  #points(sqrt(D)*(Y-theta), g$g(sqrt(D)*(Y-theta)), col='red')
  
  return(sum(lliks) + log(alpha0_prior(alpha.0)) + log(alpha1_prior(alpha.1)))
}

#alphas_llik(dat$Y, dat$D, theta.ests, -0.5, 10)
#alphas_llik(dat$Y, dat$D, theta.ests, 0.5, 100)

alphas_mh <- function(Y, D, theta, alpha.0.curr, alpha.1.curr, mh.scale.al) {

  llik.current <- alphas_llik(dat$Y, dat$D, theta, alpha.0.curr, alpha.1.curr)
  alpha.0.prop <- mh_scalar_proposal(alpha.0.curr, mh.scale.al[1])
  alpha.1.prop <- mh_scalar_proposal(alpha.1.curr, mh.scale.al[2])
  llik.proposed <- alphas_llik(dat$Y, dat$D, theta, alpha.0.prop, alpha.1.prop)
  
  mh_rat <- max(0, exp(llik.proposed - llik.current))
  #cat('curr: ', llik.current, ', prop: ', llik.proposed, ', a0.c:', alpha.0.curr, 
  #    ', a1.c:', alpha.1.curr, ', a0.p:', alpha.0.prop, ', a1.p:', alpha.1.prop, '\n')
  r <- runif(1)
  if (-Inf == llik.proposed || mh_rat < r) {
    return(list(alpha.0 = alpha.0.curr,
                alpha.1 = alpha.1.curr)) 
  } else {
    return(list(alpha.0 = alpha.0.prop,
                alpha.1 = alpha.1.prop))
  }
}


# this will walk over alphas (maybe make one that walks over gammas as a proxy for g?)
make_mhingibbs_sampler <- function(X, Y, D, var_gamma_a=0.01, var_gamma_b=0.01) {
  
  
  theta_var_marginal <- theta_var_marginal_maker(var_gamma_a, var_gamma_b, X)
  
  # n.thin like batch length e.g. 10 means take every tenth
  mh_in_gibbs_sampler <- function(params.init, n.samples, n.thin, mh.scale.thetas, mh.scale.alphas) {
    params.init.flat <- c(params.init$beta, params.init$theta, params.init$theta.var, 
                          params.init$alpha.0, params.init$alpha.1)
    n.params <- length(params.init.flat)
    
    param.samples <- matrix(NA, n.samples, n.params)
    param.samples[1,] <- params.init.flat
    p <- length(params.init$theta)
    if (length(mh.scale.alphas) != 2) {
      stop('mh.scale should be ', 2)
    }
    
    theta.init <- params.init$theta
    theta.var.init <- params.init$theta.var
    alpha.0.init <- params.init$alpha.0
    alpha.1.init <- params.init$alpha.1
    
    if(is.null(theta.init) || is.null(theta.var.init) || is.null(alpha.0.init) || is.null(alpha.1.init)) {
      stop('some param inits NULL.')      
    }
    
    beta.new <- beta_marginal(X, theta.init, theta.var.init)
    theta.new <- theta_mh(X, Y, D, beta.new, theta.init, theta.var.init, alpha.0.init, alpha.1.init, mh.scale.thetas)
    theta.var.new <- theta_var_marginal(X, beta.new, theta.new)
    alphas.new <- alphas_mh(Y, D, theta.new, alpha.0.init, alpha.1.init, mh.scale.alphas)
    alpha.0.new <- alphas.new$alpha.0
    alpha.1.new <- alphas.new$alpha.1
    
    param.samples[2,] <- c(beta.new, theta.new, theta.var.new, alpha.0.new, alpha.1.new)
    
    theta.accpt <- rep(NA, n.thin*(n.samples-3))
    alpha.0.accpt <- rep(NA, n.thin*(n.samples-3))
    alpha.1.accpt <- rep(NA, n.thin*(n.samples-3))
    overall.iter <- 0
    
    for(i in 3:n.samples) {
      cat('sample: ', i, '\n')
      for (j in 1:n.thin) {
        theta.hold <- theta.new # just for accept/reject ratio
        alpha.0.hold <- alphas.new$alpha.0
        alpha.1.hold <- alphas.new$alpha.1
        overall.iter <- overall.iter + 1
        
        beta.new <- beta_marginal(X, theta.new, theta.var.new)
        theta.new <- theta_mh(X, Y, D, beta.new, theta.new, theta.var.new, alpha.0.new, alpha.1.new, mh.scale.thetas)
        theta.var.new <- theta_var_marginal(X, beta.new, theta.new)  
        alphas.new <- alphas_mh(Y, D, theta.new, alpha.0.new, alpha.1.new, mh.scale.alphas)
        alpha.0.new <- alphas.new$alpha.0
        alpha.1.new <- alphas.new$alpha.1
        
        theta.accpt[overall.iter] <- mean(theta.hold == theta.new)
        alpha.0.accpt[overall.iter] <- alpha.0.hold == alpha.0.new
        alpha.1.accpt[overall.iter] <- alpha.1.hold == alpha.1.new
      }
      param.samples[i,] <- c(beta.new, theta.new, theta.var.new, alpha.0.new, alpha.1.new)
    }
    
    n.beta <- length(params.init$beta)
    n.theta <- length(params.init$theta)
    param.samples.list <- list(
      beta = param.samples[,1:n.beta],
      theta = param.samples[,(n.beta+1):(n.beta+n.theta)],
      theta.var = param.samples[,(n.beta+n.theta+1)],
      alpha.0 = param.samples[,(n.beta+n.theta+2)],
      alpha.1 = param.samples[,(n.beta+n.theta+3)]
    )
    
    theta.accept = mean(theta.accpt)
    alpha.0.accept = mean(alpha.0.accpt)
    alpha.1.accept = mean(alpha.1.accpt)
    
    output = list(
      param.init = params.init,
      n.samples = n.samples,
      n.thin = n.thin,
      #param.samples = param.samples,
      param.samples.list = param.samples.list,
      theta.accept = theta.accept,
      alpha.0.accept = alpha.0.accept,
      alpha.1.accept = alpha.1.accept
    )
    
    return(output)
  }
  
  
  return(mh_in_gibbs_sampler)  
}
