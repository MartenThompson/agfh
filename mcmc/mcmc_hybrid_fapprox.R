library(mvtnorm)
source('data_gen.R')

g_maker <- function (params, S) {
  alpha.0 <- params[1]
  alpha.1 <- params[2]
  
  if (alpha.0 < 0 || alpha.1 < 0) {
    stop('kern params less than zero: ', alpha.0, ', ', alpha.1)
  }
  
  rbf_kernel <- function(X) {
    K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
    return(K)
  }
  
  Sigma.virtual <- rbf_kernel(S)
  diag(Sigma.virtual) <- diag(Sigma.virtual) + 1e-4
  Sigma.virtual.inv <- solve(Sigma.virtual)
  gamma <- mvtnorm::rmvnorm(1, mean=rep(0, length(S)), sigma=Sigma.virtual)  # virtual observations
  
  #lhs <- mean(gamma[S < 0])
  #rhs <- mean(gamma[S >= 0])
  #cat(lhs, rhs, '\n')
  #gamma[S < 0] <- gamma[S < 0]/lhs
  #gamma[S >=0] <- gamma[S >=0]/rhs

  g <- function(x) {
    Sigma.new.virtual.full <- rbf_kernel(c(x, S))
    m <- length(x)
    n <- length(S)
    Sigma.new.virtual <- Sigma.new.virtual.full[1:m,(m+1):(m+n)]
    post.mean <- Sigma.new.virtual %*% Sigma.virtual.inv %*% t(gamma)
    return(post.mean)
  }
  
  mass_g <- integrate(function(u) {exp(-u^2/2 + g(u))}, -100,100)
  
  # so exp(-u^2/2 + g(u)) integrates to one
  g_scaled <- function(x) {
    g(x) + log(1/as.numeric(mass_g$value))
  }
  
  #lhs <- integrate(function(u) {u*exp(-u^2/2 + g_scaled(u))}, -100,0)
  #rhs <- integrate(function(u) {u*exp(-u^2/2 + g_scaled(u))}, 0,100)
  #cat(lhs$value, rhs$value, '\n')
  
  return(list(g=g_scaled,
              g_raw=g,
              gamma=gamma,
              sig.virt=Sigma.virtual))
}

#S <- seq(-100, 100, length.out=1e2)
#g <- g_maker(c(1/2,100), S)
#plot(S, g$gamma, col='green', ylim=c(-4,3), xlim=c(-140,140))
#lines(S, g$g_raw(S), col='green')
#lines(S, g$g(S))
#snew <- sort(runif(30, -140, 140))
#points(snew, g$g_raw(snew), col='red', pch=16)
#points(snew, g$g(snew), pch=16)

# mass
#K <- integrate(function(u) {exp(-u^2/2 + g$g_raw(u))}, -100,100)
#integrate(function(u) {exp(-u^2/2 + log(1/K$value) + g$g_raw(u))}, -100,100)
#integrate(function(u) {exp(-u^2/2 + g$g(u))}, -100,100)

# mean
# integrate(function(u) {u*exp(-u^2/2 + g$g_raw(u))}, -100,100)
# integrate(function(u) {u*exp(-u^2/2 + g$g(u))}, -100,100)
# m <- sapply(1:1e2, function(x) {
#  g <- g_maker(c(1/2,10), S)
#  i <- integrate(function(u) {u*exp(-u^2/2 + g$g(u))}, -100,100)
#  return(i$value)})
#hist(m)
#summary(m)


# var 
# without any further modifications, this is pretty close to 1.
# slightly worse for larger values of alpha.1.
#v <- sapply(1:5e2, function(x) {
#  g <- g_maker(c(1/2,100), S)
#  i <- integrate(function(u) {u^2*exp(-u^2/2 + g$g(u))}, -100,100)
#  return(i$value)})
#hist(v)
#summary(v)

lu_post_maker <- function(X, Y, D, S, include_g=TRUE) {
  p <- ncol(X)
  M <- length(D)
  beta_prior <- function(beta) {
    mvtnorm::dmvnorm(beta, mean=rep(0,p), sigma = diag(10, p))
  }
  psi_prior <- function(psi) {
    dnorm(psi, mean=0, sd=sqrt(.01)) 
    # log normal has mean exp(mu + sigma^2/2)
  }
  alpha0_prior <- function(alpha) {
    dbeta(alpha, 1, 3)
  }
  alpha1_prior <- function(alpha) {
    dgamma(alpha, shape=2, rate=1/50)
  }
  
  log_unnormalized_post <- function(params) {
    #cat('params: ', params, '\n')
    Beta <- params[1:p]
    psi <- params[p+1]
    
    if (include_g) {
      alpha.0 <- params[p+2] # exp(params[p+2])
      alpha.1 <- params[p+3] # exp(params[p+3])
      if(alpha.0 <= 0 || alpha.1 <= 0) {
        return(-Inf)
      }
      
      g.mkr <- g_maker(c(alpha.0, alpha.1), S)
      g <- g.mkr$g
    } else {
      # this will just at a const to lik
      alpha.0 <- 0.5
      alpha.1 <- 1
      g <- function(x) {
        return(rep(0, length(x)))
      }
    }
    
    Utilde <- sqrt(D)*(Y - X%*%Beta)
    
    if (length(D) != length(Utilde)) {
      stop('length mismatch:', length(D), length(Utilde))
    }
    
    # log Likelihood is sum_M sum_n_m f_mi(U_mi|g,params).
    # We set n_m \equiv 1 so llog is sum_M f_mi(U_mi|g,params).
    # f_mi(u) is an integral w/ integrand exp{fh(u-s) + g(u-s)} over s.
    # We approximate this integral with a histogram, evaluated at points S.
    # 
    
    Delta.S <- diff(S)
    fapprox <- rep(NA, M)
    fapprox_nog <- rep(NA, M)
    for (m in 1:M) {
      g_evald <- g(Utilde[m]-S)

      fh <- sapply(S, function(s) {
        -0.5*(Utilde[m]-s)^2 - (s^2)/(2*D[m]*exp(psi))
      })
      
      fapprox[m] <- sum(exp(fh + g_evald)*c(Delta.S, 0))
      fapprox_nog[m] <- sum(exp(fh)*c(Delta.S,0))
    }
    #cat(sum(log(fapprox[fapprox > 1e-5])), '\n')
    #plot(Utilde, fapprox_nog, main='f(u)', col='red', ylim=c(0,4))
    points(Utilde, fapprox, col='green', pch=16)
    
    if (include_g) {
      log_lik <- sum(log(fapprox))   
    } else {
      log_lik <- sum(log(fapprox_nog)) 
    }
    
    log_beta_pr <- log(beta_prior(Beta))
    log_psi_pr <- log(psi_prior(psi)) 
    log_alpha0_pr <- log(alpha0_prior(alpha.0)) 
    log_alpha1_pr <- log(alpha1_prior(alpha.1))
    cat(log_lik, ',', log_beta_pr, ',', log_psi_pr, ',', log_alpha0_pr, ',', log_alpha1_pr, '\n')
    lup <- log_lik + log_beta_pr + log_psi_pr + log_alpha0_pr + log_alpha1_pr
    
    return(lup)
  }
}


# given samples of beta, theta.var, return samples of theta
monte_carlo_posterior_theta <- function(n.sample, X, params.ts, include.g=FALSE, seed=NA) {
  M <- nrow(X)
  p <- ncol(X)
  theta.samples <- matrix(NA, nrow=n.sample, ncol=M)
  n.mcmc <- nrow(params.ts)
  idx <- sample(1:n.mcmc, n.sample, replace=TRUE)
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  if (include.g) {
    
  } else {
    # theta_m = x_m Beta + N(0, e^psi)
    for (i in 1:n.sample) {
      params <- params.ts[idx[i],]
      beta <- params[1:p]
      psi <- params[p+1]
      sigma <- diag(rep(exp(psi), M))
      theta.samples[i,] <- mvtnorm::rmvnorm(1, X%*%beta, sigma)
      #cat(dim(X%*%beta), '\n')
    }
  }
  
  return(theta.samples)
}


profile_me <- function() {
  set.seed(2022)
  dat <- null_gen(n=100, p=3)
  cat('True params.\n', 'Beta: ', dat$Beta, '\nPsi: ', dat$psi)
  
  S <- seq(-20, 20, length.out=100)
  lup <- lu_post_maker(dat$X, dat$Y, dat$D, S)
  params.init <- rep(0.1, 3+1+2) # beta 1, 2, 3, psi, alpha.0, alpha.1
  
  lup(params.init)  
}



#plot(seq(0,100,length.out=1000), sapply(seq(0,100,length.out=1000), function(x){dgamma(x, 2,scale=10)}))

# log gamma density
#temp <- function(psi) {
#  a <- 1
#  b <- 1
#  return(exp(b*psi)*exp(-exp(psi)/a)/(a^b * gamma(b)))
#}
#plot(sapply(seq(-10,10,length.out=100), temp))
#plot(exp(sapply(seq(-10,10,length.out=100), temp)))

