# Frequentist EBLUP are the same as EB estimators (Rao & Molina, 2015).
# RM for Rao-Molina


# Rao Molina Eq 6.1.5 for Eq 6.1.12
# theta_var_est: estimate of theta latent variable's variance e.g. exp(psi) or sigma_nu^2 (latter is RM)
# assumes each small area has one obs i.e. Y_i bar is just Y_i
# D precisions
RM_beta_eblue <- function(X, Y, D, theta_var_est) {
  X.weighted <- X/(1/D + theta_var_est)
  Y.weighted <- Y/(1/D + theta_var_est)

  beta.eblue <- solve(t(X)%*%X.weighted)%*%t(X)%*%Y.weighted
  beta.eblue
}


# RM Eq 6.1.15
RM_theta_var_moment_est <- function(X, Y, D) {
  M <- dim(X)[1]
  p <- dim(X)[2]

  beta.hat <- solve(t(X)%*%X)%*%t(X)%*%Y
  H.diag <- diag(X%*%solve(t(X)%*%X)%*%t(X))
  var.est <- 1/(M-p)*(sum((Y-X%*%beta.hat)^2) - sum((1/D*(1-H.diag)))) # TODO unclear abt this second term

  var.est <- max(0, var.est)
  var.est
}
#RM_theta_var_moment_est(dat$X, dat$Y, dat$D)

# Yoshimori & Lahiri (2014) adjusted (profile) likelihood of theta var
# D precisions
adj_profile_likelihood_theta_var_maker <- function(X, Y, D) {
  adj_profile_likelihood_theta_var <- function(theta_var) {
    # Adjustment Factor
    B <- (1/D)/(theta_var + 1/D)
    M <- length(Y)
    trace.IB <- sum(rep(1, M)-B)
    h.YL <- (atan(trace.IB))^(1/M)

    # Profile Likelihood
    V <- diag(theta_var + 1/D)
    V.inv <- solve(V)
    P <- V.inv - V.inv%*%X%*%solve(t(X)%*%V.inv%*%X)%*%t(X)%*%V.inv

    Log_Lik_prf <- (-1/2)*log(det(V)) -(1/2)*t(Y)%*%P%*%Y
    l <- log(h.YL) + Log_Lik_prf
    return(h.YL * (det(V)^(-1/2))*exp(-(1/2)*t(Y)%*%P%*%Y))
  }

  return(adj_profile_likelihood_theta_var)
}

# Yoshimori & Lahiri (2014) adjusted (residual) likelihood of theta var
# D precisions
adj_resid_likelihood_theta_var_maker <- function(X, Y, D) {
  adj_resid_likelihood_theta_var <- function(theta_var) {
    # Adjustment Factor
    B <- (1/D)/(theta_var + 1/D)
    M <- length(Y)
    trace.IB <- sum(rep(1, M)-B)
    h.YL <- (atan(trace.IB))^(1/M)

    # Residual Likelihood
    V <- diag(theta_var + 1/D)
    V.inv <- solve(V)
    P <- V.inv - V.inv%*%X%*%solve(t(X)%*%V.inv%*%X)%*%t(X)%*%V.inv

    lik.resid <- (det(t(X)%*%V%*%X)^(-1/2))*(det(V)^(-1/2))*exp(-(1/2)*t(Y)%*%P%*%Y)
    return(h.YL * lik.resid)
  }

  return(adj_resid_likelihood_theta_var)
}


# Yoshimori & Lahiri (2014) non-adjusted (residual) likelihood of theta var
# D precisions
resid_likelihood_theta_var_maker <- function(X, Y, D) {
  resid_likelihood_theta_var <- function(theta_var) {
    V <- diag(theta_var + 1/D)
    V.inv <- solve(V)
    P <- V.inv - V.inv%*%X%*%solve(t(X)%*%V.inv%*%X)%*%t(X)%*%V.inv

    lik.resid <- (det(t(X)%*%V%*%X)^(-1/2))*(det(V)^(-1/2))*exp(-(1/2)*t(Y)%*%P%*%Y)
    return(lik.resid)
  }

  return(resid_likelihood_theta_var)
}


# expects you to work small -> large over subsequent calls.
grid_opt <- function(f, x.min, x.max, n.x) {
  grid <- seq(x.min, x.max, length.out=n.x)
  max.id <- which.max(sapply(grid, f))
  opt.val <- grid[max.id]

  if (x.max == opt.val) {
    return(list(
      increase_range = TRUE,
      opt_val = opt.val
    ))
  } else {
    return(list(
      increase_range = FALSE,
      opt_val = opt.val
    ))
  }
}

theta_var_est_grid <- function(likelihood_theta_var) {
  var.min <- 1e-6
  var.max <- 20
  n.eval <- 1e3
  opt.out <- grid_opt(likelihood_theta_var, var.min, var.max, n.eval)

  if (opt.out$increase_range) {
    var.min <- var.max
    var.max <- 100
    n.eval <- 1e3
    opt.out <- grid_opt(likelihood_theta_var, var.min, var.max, n.eval)
    return(opt.out$opt_val)
  } else {
    return(opt.out$opt_val)
  }
}


# Rao Molina Eq 6.1.2, 6.1.12
# D precisions
RM_theta_eblup <- function(X, Y, D, theta.var.est=NA) {
  M <- dim(X)[1]
  p <- dim(X)[2]

  if (p > M) {
    stop('M should be > p. Actual values, M=', M, ', p=', p)
  }

  if (is.na(theta.var.est)) {
    # if an specific flavor of theta.var.est is supplied, use what is passed in;
    # otherwise, use the basic moment estimator.
    theta.var.est <- RM_theta_var_moment_est(X, Y, D)
  }

  gamma <- theta.var.est/(1/D + theta.var.est)
  beta.eblue <- RM_beta_eblue(X, Y, D, theta.var.est)
  theta.eblup <- gamma*Y + (1-gamma)*X%*%beta.eblue
  theta.eblup
}



# X_new is n by p. beta is p by 1.
# Roa Molina Eq 6.1.13
RM_theta_new_pred <- function(X.new, beta.est) {
  return(X.new%*%beta.est)
}
