library(mcmc)
setwd('~/Git/agfh/mcmc/')
source('../data_gen.R')


utilde_dens <- function(v, D, psi, g) {
  if (1 != length(v)) {
    stop('Utilde dense is for scalars, recieved shape ', length(v))
  }

  f <- function(s) {
    if (1 == length(s)) {
      exp(-0.5*(v-s)^2 - (s^2)/(2*D*exp(psi)) + g(v-s))
    } else {
      sapply(s, function(x){exp(-0.5*(v-x)^2 - (x^2)/(2*D*exp(psi)) + g(v-x))})
    }
  } 
  
  ans <- integrate(f, -1e2,1e2,subdivisions = 1e3, stop.on.error = FALSE)
  return(ans$value)
}

g_maker <- function (mu) {
  K <- 1
  if (length(mu) != 9) {
    stop('should be 9, is: ', length(mu))
  }
  xi <- seq(-15,15,length.out=10)
  mu.final <- log(K/(sqrt(2*pi)*exp(sum(mu))))
  mu <- c(mu, mu.final)
  
  g <- function(x) {
    if (1 != length(x)) {
      stop('x should be scalar, is length ', length(x))
    }
    alpha.0 <- 1 
    alpha.1 <- 1
    rbf_kernel <- function(X) {
      X <- matrix(X, ncol=1)
      n <- nrow(X)
      K <- matrix(NA,n,n)
      for (i in 1:n) {
        for (j in i:n) {
          K[i,j] <- alpha.0*exp(-alpha.1*abs(X[i,1] - X[j,1]))
        }
      }
      
      K[lower.tri(K)] = t(K)[lower.tri(K)]
      return(K)
    }
    
    Sigma.new.virtual.full <- rbf_kernel(c(xi, x))
    m <- ncol(Sigma.new.virtual.full)
    Sigma.virtual <- Sigma.new.virtual.full[1:(m-1), 1:(m-1)]
    Sigma.new.virtual <- Sigma.new.virtual.full[m, 1:(m-1)]
    
    post.mean <- Sigma.new.virtual %*% solve(Sigma.virtual) %*% mu
    return(post.mean)
  }
  
  return(g)
}

lu_post_maker <- function(X, Y, D) {
  # TODO currently uninformed prior
  log_unnormalized_post <- function(params) {
    p <- dim(X)[2]
    Beta <- params[1:p]
    psi <- params[p+1]
    mu <- params[(p+2):(p+2+9-1)]
    
    Utilde <- sqrt(D)*(Y - X%*%Beta)
    g <- g_maker(mu)
    
    if (length(D) != length(Utilde)) {
      stop('length mismatch:', length(D), length(Utilde))
    }
    
    return(
      sum(
        mapply(
          function(um, Dm) {utilde_dens(um, Dm, psi, g)}, 
          Utilde, D
               )
        )
      )
  }
}




dat <- null_gen(n=100, p=3)
cat('True params.\n', 'Beta: ', dat$Beta, '\nPsi: ', dat$psi)

lup <- lu_post_maker(dat$X, dat$Y, dat$D)
params.init <- rep(0, 3+1+9) # beta 1, 2, 3, psi, 9 virtual obs
lup(params.init)

set.seed(2022)
mcmc.out <- metrop(lup, params.init, 1e3, scale = 0.5)
mcmc.out$accept
saveRDS(mcmc.out, 'mcmc_out.RData')

mcmc.out <- readRDS('mcmc_out.RData')


t <- ts(mcmc.out$batch)



virt.ten <- sapply(apply(t[,5:13], 1, sum), function(row){log(1/(sqrt(2*pi)*exp(sum(row))))})
plot(virt.ten)

for (i in 1:ncol(t)) {
  plot(t[,i], main=paste0('Trace ', i))
}

par(mfrow=c(2,2))
plot(density(t[,1]))
plot(density(t[,2]))
plot(density(t[,3]))
plot(density(t[,4]))
par(mfrow = c(1,1))








# ========
xi.garbo <- seq(-1,1,length.out=10)
mu.garbo <- sin(2*xi.garbo[1:9])
mu.g.final <- log(1/(sqrt(2*pi)*exp(sum(mu.garbo))))

mu.garbo <- c(mu.garbo, mu.g.final)

g.garbo <- function(x) {
  # x should be scalar
  alpha.0 <- 1 
  alpha.1 <- 1
  rbf_kernel <- function(X) {
    X <- matrix(X, ncol=1)
    n <- nrow(X)
    K <- matrix(NA,n,n)
    for (i in 1:n) {
      for (j in i:n) {
        K[i,j] <- alpha.0*exp(-alpha.1*abs(X[i,1] - X[j,1]))
      }
    }
    
    K[lower.tri(K)] = t(K)[lower.tri(K)]
    return(K)
  }
  
  Sigma.new.virtual.full <- rbf_kernel(c(xi.garbo, x))
  m <- ncol(Sigma.new.virtual.full)
  Sigma.virtual <- Sigma.new.virtual.full[1:(m-1), 1:(m-1)]
  Sigma.new.virtual <- Sigma.new.virtual.full[m, 1:(m-1)]
  
  #return(Sigma.new.virtual.full)
  post.mean <- Sigma.new.virtual %*% solve(Sigma.virtual) %*% mu.garbo
  return(post.mean)
}

x.garbo <- seq(-2,2,length.out=20)
plot(x.garbo, sapply(x.garbo, g.garbo), ylim=c(-5,5))
points(xi.garbo, mu.garbo, col='red')





