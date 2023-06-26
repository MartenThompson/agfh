library(agfh)

set.seed(2023)
n <- 10
p <- 1
X <- matrix(1:n, ncol=p)
Y <- 2*X + rnorm(n, sd=1.1)
D <- rep(1, n)
ag <- make_agfh_sampler(X, Y, D)

params.init <- list(
  beta=1,
  theta=rep(0,n),
  theta.var=1,
  gamma=rep(0,n)
)
n.mcmc <- 5
ag.out <- ag(params.init, n.mcmc, 1, 0.1)


test_that('init params save correctly', {
  expect_equal(c(ag.out$param.init$beta, ag.out$param.init$theta, ag.out$param.init$theta.var, ag.out$param.init$gamma),
               c(params.init$beta, params.init$theta, params.init$theta.var, params.init$gamma))
})

test_that('correct number of samples are returned', {
  expect_equal(dim(ag.out$param.samples.list$beta), c(n.mcmc, p))
})

