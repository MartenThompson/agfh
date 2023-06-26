mh_scalar_proposal <- function(current, scale) {
  return(rnorm(1, mean=current, sd=sqrt(scale)))
}

theta_scalar_llik <- function(theta, x, y, d, beta, theta.var, g) {
  x < matrix(x, ncol=1)

  llik <- -0.5*d*(y-theta)^2 + g$g(sqrt(d)*(y-theta)) - 0.5*(1/theta.var)*(theta - beta%*%x)^2
  return(llik)
}


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


rbf_kernel <- function(X, alpha.0, alpha.1, phi) {
  K <- alpha.0*exp(-(1/alpha.1)*as.matrix(dist(X, diag=TRUE, upper=TRUE))^2)
  diag(K) <- diag(K) + phi
  #image(K)
  return(K)
}


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

  # mass_u <- tryCatch({
  #  mass_u <- integrate(function(u){exp(-u^2/2 + g(u))}, -100, 100)
  # }, error=function(cond) {
  #   cat('Bad Normalizing Integral for g.\n')
  #   mass_u <- list(value=NA)
  # })
  mass_u <- list(value=NA)

  return(list(
    g = g,
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


g_update <- function(Y, D, theta, g.current) {
  u.hat<- sqrt(D)*(Y - theta)
  u.reconstructed <- (u.hat - mean(u.hat))/sd(u.hat)

  kde.u <- ks::kde(u.reconstructed)
  kde.est <- predict(kde.u, x = u.reconstructed)
  g.hat <- log(sqrt(2*pi)*kde.est) + (u.reconstructed^2)/2
  g.curr <- init_g(g.hat, u.reconstructed, g.current$a0, g.current$a1, g.current$phi)
  return(g.curr)
}


make_agfh_sampler <- function(X, Y, D, var_gamma_a=0.01, var_gamma_b=0.01,
                                   S=seq(-50,50,length.out=length(D)), kern.a0=0.5, kern.a1=1, kern.fuzz=0.1) {

  theta_var_marginal <- theta_var_marginal_maker_hb(var_gamma_a, var_gamma_b, X)

  # n.thin like batch length e.g. 10 means take every tenth
  mh_in_gibbs_sampler <- function(params.init, n.samples, n.thin, mh.scale.thetas, trace=TRUE) {
    beta.samples <- matrix(NA, n.samples, length(params.init$beta))
    theta.samples <- matrix(NA, n.samples, length(params.init$theta))
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

    beta.new <- beta_marginal_hb(X, theta.init, theta.var.init)
    theta.new <- theta_mh(X, Y, D, beta.new, theta.init, theta.var.init, g.init, mh.scale.thetas)
    theta.var.new <- theta_var_marginal(X, beta.new, theta.new)
    g.new <- g_update(Y, D, theta.new, g.init)

    beta.samples[2,] <- beta.new
    theta.samples[2,] <- theta.new
    theta.var.samples[2,] <- theta.var.new
    gamma.samples[2,] <- g.new$gamma
    kern.a0.samples[2] <- g.new$a0
    kern.a1.samples[2] <- g.new$a1
    e.eta.samples[2] <- g.new$e.eta

    theta.accpt <- rep(NA, n.thin*(n.samples-2))
    overall.iter <- 0

    for(i in 3:n.samples) {
      if(trace) {
        cat('sample: ', i, '\n')
      }

      for (j in 1:n.thin) {
        theta.hold <- theta.new # just for accept/reject ratio
        overall.iter <- overall.iter + 1

        beta.new <- beta_marginal_hb(X, theta.new, theta.var.new)
        theta.new <- theta_mh(X, Y, D, beta.new, theta.new, theta.var.new, g.new, mh.scale.thetas)
        theta.var.new <- theta_var_marginal(X, beta.new, theta.new)
        g.new <- g_update(Y, D, theta.new, g.new)

        theta.accpt[overall.iter] <- mean(theta.hold != theta.new)
      }

      beta.samples[i,] <- beta.new
      theta.samples[i,] <- theta.new
      theta.var.samples[i,] <- theta.var.new
      gamma.samples[i,] <- g.new$gamma

      kern.a0.samples[i] <- g.new$a0
      kern.a1.samples[i] <- g.new$a1
      e.eta.samples[i] <- g.new$e.eta
    }

    param.samples.list <- list(
      beta = beta.samples,
      theta = theta.samples,
      theta.var = theta.var.samples,
      gamma = gamma.samples,
      kern.a0 = kern.a0.samples,
      kern.a1 = kern.a1.samples,
      e.eta = e.eta.samples
    )


    output = list(
      param.init = params.init,
      S.init = S,
      mh.scale.thetas = mh.scale.thetas,
      n.samples = n.samples,
      n.thin = n.thin,
      #param.samples = param.samples,
      param.samples.list = param.samples.list,
      kern.a0=kern.a0,
      kern.a1=kern.a1,
      kern.fuzz=kern.fuzz,
      g.final =g.new,
      theta.avg.accept = mean(theta.accpt)
    )

    return(output)
  }


  return(mh_in_gibbs_sampler)
}

map_from_density <- function(param.ts, plot=FALSE) {
  d <- density(param.ts)
  map <- d$x[which.max(d$y)]

  if (plot) {
    plot(d)
    lines(c(map,map), c(0,100))
  }

  return(map)
}


mse <- function(x,y) {
  mean((x-y)^2)
}

# X_new is p by 1
# beta is p by k, theta_var is length k.
# identical to hb_theta_new_pred
agfh_theta_new_pred <- function(X_new, beta_samples, theta_var_samples) {
  if (ncol(beta_samples) != length(theta_var_samples)) {
    stop(ncol(beta_samples), ' samples of beta, but ',  length(theta_var_samples), ' variance samples.')
  }
  means <- X_new%*%beta_samples
  return(rnorm(length(theta_var_samples), means, sqrt(theta_var_samples)))
}


