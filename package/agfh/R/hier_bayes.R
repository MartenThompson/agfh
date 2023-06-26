# Gibbs sampler described by Rao&Molina 2015 Sec 10.3.3; Eqs 10.3.13-15


beta_marginal_hb <- function(X, theta.est, theta.var.est) {
  M <- dim(X)[1]
  theta.est <- matrix(theta.est, nrow=M, ncol=1)

  mu <- solve(t(X)%*%X)%*%t(X)%*%theta.est # TODO: this is correct? not Y
  sigma <- solve(t(X)%*%X)*theta.var.est
  mvtnorm::rmvnorm(1, mean=mu, sigma=sigma) # 1xp vector
}


theta_marginal <- function(X, Y, D, beta.est, theta.var.est) {
  beta.est <- matrix(beta.est, ncol=1)
  gamma <- theta.var.est/(1/D + theta.var.est)
  mu <- gamma*Y + (1-gamma)*X%*%beta.est
  sigma <- gamma/D
  mvtnorm::rmvnorm(1, mean = mu, sigma=diag(sigma)) # 1 x M vector
}


theta_var_marginal_maker_hb <- function(prior.gamma.shape.a, prior.gamma.rate.b, X) {
  M <- dim(X)[1]
  a <- prior.gamma.shape.a
  b <- prior.gamma.rate.b

  theta_var_marginal <- function(X, beta.est, theta.est) {
    beta.est <- matrix(beta.est, ncol=1)
    theta.est <- matrix(theta.est, nrow=M, ncol=1)

    shape <- a + M/2
    rate <- b + 0.5*sum((theta.est - X%*%beta.est)^2)
    prec <- rgamma(1, shape=shape, rate=rate)
    return(1/prec)
  }

  return(theta_var_marginal)
}


# beta has flat prior, and 1/theta_var term has Gamma(a,b) prior.
make_gibbs_sampler <- function(X, Y, D, var_gamma_a=1, var_gamma_b=1) {
  theta_var_marginal <- theta_var_marginal_maker_hb(var_gamma_a, var_gamma_b, X)

  # n.thin like batch length e.g. 10 means take every tenth
  RM_gibbs_sampler <- function(params.init, n.samples, n.thin) {
    params.init.flat <- c(params.init$beta, params.init$theta, params.init$theta.var)
    n.params <- length(params.init.flat)

    param.samples <- matrix(NA, n.samples, n.params)
    param.samples[1,] <- params.init.flat

    theta.init <- params.init$theta
    theta.var.init <- params.init$theta.var

    beta.new <- beta_marginal_hb(X, theta.init, theta.var.init)
    theta.new <- theta_marginal(X, Y, D, beta.new, theta.var.init)
    theta.var.new <- theta_var_marginal(X, beta.new, theta.new)
    param.samples[2,] <- c(beta.new, theta.new, theta.var.new)

    for(i in 3:n.samples) {
      for (j in 1:n.thin) {
        beta.new <- beta_marginal_hb(X, theta.new, theta.var.new)
        theta.new <- theta_marginal(X, Y, D, beta.new, theta.var.new)
        theta.var.new <- theta_var_marginal(X, beta.new, theta.new)
      }
      param.samples[i,] <- c(beta.new, theta.new, theta.var.new)
    }

    n.beta <- length(params.init$beta)
    n.theta <- length(params.init$theta)
    param.samples.list <- list(
      beta = param.samples[,1:n.beta],
      theta = param.samples[,(n.beta+1):(n.beta+n.theta)],
      theta.var = param.samples[,(n.beta+n.theta+1)]
    )

    output = list(
      params.init = params.init,
      n.samples = n.samples,
      n.thin = n.thin,
      #param.samples = param.samples,
      param.samples.list = param.samples.list
    )

    return(output)
  }


  return(RM_gibbs_sampler)
}


# X_new is 1 by p.
# beta is p by k, theta_var is length k.
# RM b/t Eq 10.13.13 & 10.13.14
hb_theta_new_pred <- function(X_new, beta_samples, theta_var_samples) {
  if (ncol(beta_samples) != length(theta_var_samples)) {
    stop(ncol(beta_samples), ' samples of beta, but ',  length(theta_var_samples), ' variance samples.')
  }
  means <- X_new%*%beta_samples
  return(rnorm(length(theta_var_samples), means, sqrt(theta_var_samples)))
}

