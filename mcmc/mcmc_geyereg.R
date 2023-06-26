library(mcmc)
data(logit)
out <- glm(y ~ x1 + x2 + x3 + x4, data = logit, family = binomial, x = TRUE)



lupost_factory <- function(x,y) function(beta) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)), -log1p(exp(-eta)))
  logq <- ifelse(eta < 0, -log1p(exp(eta)), -eta - log1p(exp(-eta)))
  logl <- sum(logp[y==1]) + sum(logq[y==0])
  return(logl - sum(beta^2)/8)
}

lupost <- lupost_factory(out$x, out$y)

set.seed(42)
beta.init <- as.numeric(coefficients(out))
mcmc.out <- metrop(lupost, beta.init, 1e3)


plot(ts(mcmc.out$batch))
acf(mcmc.out$batch)
