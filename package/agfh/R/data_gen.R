check_inputs_datagen <- function(M, D, lamb) {
  if (M != length(D)) {
    stop('D should have length M=', M)
  }
  if (lamb <= 0) {
    stop('Theta variance (lambda) should be > 0.')
  }
}

beta_gen <- function(p) {
  return(runif(p, -5,5))
  #return(rep(1,p))
}

X_gen <- function(M,p) {
  matrix(runif(M*p, 0,10), nrow = M, ncol=p)
}

# Shared by all methods
theta_gen <- function(M, p, X, Beta, lambda) {
  theta <- X%*%Beta + rnorm(M, 0, sd=sqrt(lambda))
}


# Y = theta + N(0, D^-1) = theta + D^-1/2 N(0,1)
# theta = x Beta + N(0, lambda)
null_gen <- function(M, p, D, lambda) {
  check_inputs_datagen(M, D, lambda)
  Beta <- beta_gen(p)
  X <- X_gen(M, p)
  theta <- theta_gen(M, p, X, Beta, lambda)

  err <- rnorm(M,0,1) # mean 0, var 1
  Y <- theta + (1/sqrt(D))*err

  ret <- list(
    D = D,
    Beta = Beta,
    lambda = lambda,
    X = X,
    theta = theta,
    Y = Y,
    err=err,
    name='null gen (normal normal)'
  )
  return(ret)
}


# U shape: 1/8,1/8 or 1/2, 1/2
# Assym U: 1/2, 7/8
beta_err_gen <- function(M, p, D, lambda, a, b) {
  check_inputs_datagen(M, D, lambda)
  Beta <- beta_gen(p)
  X <- X_gen(M, p)
  theta <- theta_gen(M, p, X, Beta, lambda)

  em <- a/(a+b)
  ev <- (a*b)/(((a+b)^2)*(a+b+1))
  err <- (1/sqrt(ev))*(rbeta(M,a,b) - em) # mean 0, var 1
  Y <- theta + (1/sqrt(D))*err

  ret <- list(
    D = D,
    Beta = Beta,
    lambda = lambda,
    X = X,
    theta = theta,
    Y = Y,
    err=err,
    a=a,
    b=b,
    name=paste0('Beta gen (err [Beta(1/2,1/2)-1]*sqrt(8)). a:',a,' b:', b)
  )
  return(ret)
}


# Y = theta + D^-1/2 Gamma(s,r)*
# * massaged so mean zero, unit variance
# theta = x Beta + N(0, lambda)
gamma_err_gen <- function(M, p, D, lambda, shape, rate) {
  check_inputs_datagen(M, D, lambda)
  Beta <- beta_gen(p)
  X <- X_gen(M, p)
  theta <- theta_gen(M, p, X, Beta, lambda)

  em <- shape/rate
  ev <- shape/rate^2
  err <- (1/sqrt(ev))*(rgamma(M, shape=shape, rate=rate) -em) # mean 0 var 1
  Y <- theta + (1/sqrt(D))*err

  ret <- list(
    D = D,
    Beta = Beta,
    lambda=lambda,
    X = X,
    theta = theta,
    Y = Y,
    err=err,
    shape=shape,
    rate=rate,
    name=paste0('gamma errors. shape: ', shape, ' rate: ', rate)
  )
  return(ret)
}


split_traintest <- function(data.obj, M.tr, M.te) {
  return(list(train=
                list(
                  D=data.obj$D[1:M.tr],
                  Beta=data.obj$Beta,
                  lambda=data.obj$lambda,
                  X = data.obj$X[1:M.tr,],
                  theta=data.obj$theta[1:M.tr],
                  Y = data.obj$Y[1:M.tr],
                  err= data.obj$err[1:M.tr],
                  name = paste0(data.obj$name, ' train')
                ),
              test=list(
                D=data.obj$D[(1+M.tr):(M.tr+M.te)],
                Beta=data.obj$Beta,
                lambda=data.obj$lambda,
                X = data.obj$X[(1+M.tr):(M.tr+M.te),],
                theta=data.obj$theta[(1+M.tr):(M.tr+M.te)],
                Y = data.obj$Y[(1+M.tr):(M.tr+M.te)],
                err= data.obj$err[(1+M.tr):(M.tr+M.te)],
                name = paste0(data.obj$name, ' test')
              )))
}

