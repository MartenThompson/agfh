#https://rdrr.io/github/lucas-castillo/SampleR/src/R/samplers.R

#' Markov Chain Monte Carlo Sampler
#'
#'
#' This sampler navigates the proposal distribution following a random walk. At each step, it generates a new proposal from a proposal distribution (in this case a Gaussian centered at the current position) and chooses to accept it or reject it following the Metropolis-Hastings rule: it accepts it if the density of the posterior distribution at the proposed point is higher than at the current point. If the current position is denser, it still may accept the proposal with probability `proposal_density / current_density`.
#'
#' As mentioned, the proposal distribution is a Normal distribution. Its mean is the current position, and its variance is equal to the `sigma_prop` parameter, which defaults to the identity matrix if not specified.
#'
#' @param pdf Probability Density Function of the posterior distribution. Takes a vector as input
#' @param start Vector. Starting point for the sampler
#' @param iterations Numeric. Number of times the sampler runs
#' @param sigma_prop Variance of the univariate proposal distribution. For multivariate proposals, covariance matrix of the proposal.
#'
#' @return A list containing
#' 1. the history of visited places (a n x d matrix, n = iterations; d = dimensions)
#' 1. acceptance ratio - the proportions of proposals that were accepted (numeric)
#' @export
#'
#' @examples
#' # generate a density function that only takes a vector as input
#' pdf_func <- function(x){return(stats::dnorm(x, 0, 1))}
#' # Not giving a sigma_prop issues a warning, but the sampler runs anyway with a default value
#' chain <- sampler_mcmc(pdf_func, start = 0, iterations = 20)
#'
#' # pdf functions can be easily created with the make_*_pdf helpers
#' density_function <- make_distr_pdf(distr::Beta(.5,.5))
#' chain <- sampler_mcmc(
#'  density_function,
#'  iterations = 20,
#'  start = 0,
#'  sigma_prop = .5
#'  )

sampler_mcmc <- function(pdf, start, iterations=1024, sigma_prop=NULL){
  # Initialize variables ---------------------------------
  acceptances <- 0
  n_dim <- length(start)
  chain <- matrix(0, nrow=iterations, ncol = n_dim)
  chain[1,] <- start
  
  # Create a proposal function -----------------------------
  
  # If no variance for the proposal distribution was given...
  if (is.null(sigma_prop)){
    # ... give it the arbitrary value of 1 in the univariate case and issue a warning
    if (n_dim == 1){
      sigma_prop <- 1
      warning("The variance of the proposal distribution was not given and defaulted to 1")
      # ... give it the arbitrary value of the Identity Matrix in the multivariate case and issue a warning
    } else{
      sigma_prop <- diag(n_dim)
      warning("The variance of the proposal distribution was not given and defaulted to the identity matrix")
    }
  }
  
  # get a proposal function depending on dimensions. The proposal function is a normal distribution centered around the current position and with variance = sigma_prop
  if (n_dim == 1){
    proposal_f <- function(mu){stats::rnorm(1, mean = mu, sd = sqrt(sigma_prop))}
  } else{
    proposal_f <- function(mu){mvtnorm::rmvnorm(1, mean = mu, sigma = sigma_prop)}
  }
  
  
  # Run the sampler ------------------------------------------------
  
  for (i in 2:iterations){
    # the metropolis_step helper function creates a new proposal based on the function above and chooses whether to accept it or not
    step <- metropolis_step(chain[i-1,], pdf, proposal_f)
    chain[i,] <- step[[1]]
    acceptances <- acceptances + step[[2]]
    # print(paste(format((i / iterations)*100, digits = 2),"%", sep=""))
  }
  
  return(list(
    chain = chain,
    acceptance_ratio = acceptances/(iterations - 1)))
}

#' Metropolis-coupled MCMC sampler (MC3)
#'
#' This sampler is a variant of MCMC in which multiple parallel chains are run at different temperatures. The chains stochastically swap positions which allows the coldest chain to visit regions far from its starting point (unlike in MCMC). Because of this, an MC3 sampler can explore far-off regions, whereas an MCMC sampler may become stuck in a particular point of high density.
#'
#'
#' @param pdf Probability Density Function of the posterior distribution. Takes a vector as input
#' @param start Vector. Starting point for the sampler
#' @param nChains number of parallel chains to be run.
#' @param iterations Numeric. Number of times the sampler runs
#' @param sigma_prop Variance of the univariate proposal distribution. For multivariate proposals, covariance matrix of the proposal.
#' @param delta_T numeric, >1. Temperature increment parameter. The bigger this number, the steeper the increase in temperature between the cold chain and the next chain
#' @param swap_all Boolean. If true, every iteration attempts floor(nChains / 2) swaps. If false, only one swap per iteration.
#'
#' @return List with:
#' 1. array of iterations x target_dimensions x nChains,
#' 1. acceptance ratio for each chain,
#' 1. history of swaps and
#' 1. swap ratio.
#'
#' If nChains == 1 items 3 and 4 are not returned
#' @export
#'
#' @examples
#' pdf_function <- function(x){return(mvtnorm::dmvnorm(x, mean = c(0,1), sigma = diag(2)))}
#' mc3_results <- sampler_mc3(
#'  pdf_function,
#'  start = c(0,0),
#'  iterations=10,
#'  nChains = 3,
#'  sigma_prop = diag(2) / 8
#' )
#'

sampler_mc3 <-function(pdf, start, nChains = 6, iterations=1024, sigma_prop = NULL, delta_T = 4, swap_all = TRUE){
  
  # Initialize Variables ----------------------
  acceptances <- vector(mode = "numeric", length = nChains)
  swaps <- matrix(ncol = 3, dimnames = list(NULL, c("chain1","chain2", "iteration")))
  swap_attempts <- 0
  swap_accepts <- 0
  n_dim <- length(start)
  chain <- array(0, dim = (c(iterations, n_dim, nChains)))
  chain[1,,] <- start
  beta <- 1 / (1 + delta_T * seq(0,nChains-1,length.out=nChains))
  if (swap_all){
    nSwaps <- floor(nChains/2)
  } else {
    nSwaps <- 1
  }
  
  # Create a proposal function -----------------------------
  
  # If sigma not given, get one depending on n_dimensions
  if (is.null(sigma_prop)){
    # ... give it the arbitrary value of 1 in the univariate case and issue a warning
    if (n_dim == 1){
      sigma_prop <- 1
      warning("The variance of the proposal distribution was not given and defaulted to 1")
      # ... give it the arbitrary value of the Identity Matrix in the multivariate case and issue a warning
    } else{
      sigma_prop <- diag(n_dim)
      warning("The variance of the proposal distribution was not given and defaulted to the identity matrix")
    }
  }
  
  # get a proposal function depending on dimensions. The proposal function is a normal distribution centered around the current position and with variance = sigma_prop
  if (n_dim == 1){
    proposal_f <- function(mu){stats::rnorm(1, mean = mu, sd = sqrt(sigma_prop))}
  } else{
    proposal_f <- function(mu){mvtnorm::rmvnorm(1, mean = mu, sigma = sigma_prop)}
  }
  
  # Run Sampler -----------------------------
  for (iter in 2:iterations){
    for (c in 1:nChains){
      # new step for every chain, carried out by metropolis_step
      step <- metropolis_step(chain[iter-1,,c], pdf = pdf, proposal_f = proposal_f, beta = beta[c])
      chain[iter,,c] <- step[[1]]
      acceptances[c] <- acceptances[c] + step[[2]]
    }
    
    # Swap Chains ---------------------
    if (nChains > 1){
      # arrange chains randomly
      sChain = pracma::randperm(nChains, nChains)
      # swap nSwaps times (depending on swap_all)
      for (k in 1:nSwaps){
        swap_attempts <- swap_attempts + 1
        
        m <- sChain[k*2 - 1]
        n <- sChain[k*2]
        
        # chains are swapped with probability alpha, which is the ratio between:
        #   the product of the density of each chain's location at the temperature of the other chain, and
        #   the product of the density of each chain's location at their own temperatures
        
        top <- pdf(chain[iter,,m]) ^ beta[n] * pdf(chain[iter,,n]) ^ beta[m]
        bottom <- pdf(chain[iter,,m]) ^ beta[m] * pdf(chain[iter,,n]) ^ beta[n]
        
        if ((bottom != 0 & stats::runif(1) <= top/bottom) || (bottom == 0 && top > 0)){
          temp <- chain[iter,,m]
          chain[iter,,m] <- chain[iter,,n]
          chain[iter,,n] <- temp
          swaps <- rbind(swaps, c(m,n, iter))
          swap_accepts <- swap_accepts + 1
        }
      }
    }
  }
  # swap information is not returned if nChains == 1
  if (nChains > 1){
    return (list(chain, acceptances/iterations, swaps[-1,], swap_accepts / swap_attempts))
  } else{
    return (list(chain, acceptances/iterations))
  }
}

#' Metropolis Step
#'
#' Helper function for MCMC-based samplers.Calculates a proposed move from the current position and decides whether the chain should move there based on a Metropolis-Hastings acceptance rule.
#'
#' @param current_x vector. Current position of the chain
#' @param pdf Probability density function of the target distribution
#' @param proposal_f Function from which to draw proposals, given a mean
#' @param beta Temperature of the chain. Acceptance of the proposed state becomes more likely as beta approaches 0.
#'
#' @return List with the position for the next step of the chain and a boolean representing whether the proposal has been accepted or not
#' @export
#'
#' @examples
#' # This function is used internally in MCMC based samplers
#'
#' pdf_function <- function(x){mvtnorm::dmvnorm(x, mean = c(0,1), sigma = diag(2))}
#' proposal_f <- function(mu){mvtnorm::rmvnorm(1, mean = mu, sigma = diag(2))}
#' current_position = c(2,3)
#' metropolis_step(current_position,pdf_function,proposal_f)

metropolis_step <- function(current_x, pdf, proposal_f, beta = 1){
  # generate proposal
  proposal <- proposal_f(current_x)
  # calculate current and proposal probabilities
  prob_curr <- pdf(current_x)
  prob_prop <- pdf(proposal)
  accept <- FALSE
  
  # proposal is accepted with probability prob_prop / prob_curr
  if (prob_curr != 0){
    ratio <- prob_prop / prob_curr
    if (ratio >= 1){
      accept <- TRUE
      # The beta parameter (temperature), beta <= 1, increases the value of the ratio making hotter chains more likely to accept proposals
    } else if (stats::runif(1) < (ratio ^ beta)){
      accept <- TRUE
    }
  } else {
    # in case the current probability is 0 (in which case the ratio cannot be calculated), the step is accepted if the probability of the proposal is not 0
    if (prob_prop > 0) {
      accept <- TRUE
    }
  }
  
  if (accept){
    return(list(proposal, TRUE))
  } else {
    return(list(current_x, FALSE))
  }
}


#' Leapfrog Step
#'
#' Helper Function for Hamiltonian Monte Carlo based samplers.The leapfrog function is a symplectic integrator which deals with error accumulation in the generation of Hamiltonian trajectories.
#' It requires that the distribution of the momentum is independent of position.
#'
#' @param theta A vector of length d with the current position
#' @param momentum A vector of length d with the current momentum
#' @param epsilon A numeric value with the size of the leapfrog step
#' @param log_posterior Function of the log probability density of a value theta in the target distribution.
#'
leapfrog_step <- function(theta, momentum, epsilon, log_posterior, L = 1){
  # start with half step for momentum
  momentum <- momentum + (epsilon/2) * as.vector(rootSolve::gradient(log_posterior, theta))
  # alternate full steps for position and momentum
  for (i in 1:L){
    theta <- theta + epsilon * momentum
    if (i != L){
      # full step for momentum except for the end of trajectory
      momentum <- momentum + epsilon * as.vector(rootSolve::gradient(log_posterior, theta))
    }
  }
  
  # make a half step (instead of a full one) for the momentum at the end
  momentum <- momentum + (epsilon/2) * as.vector(rootSolve::gradient(log_posterior, theta))
  # negate momentum to make the proposal symmetric
  momentum = -momentum
  
  return(list(theta = theta, momentum = momentum))
}

#' Hamiltonian Monte-Carlo Sampler.
#'
#' Hamiltonian Monte-Carlo, also called Hybrid Monte Carlo, is a sampling algorithm that uses Hamiltonian Dynamics to approximate a posterior distribution. Unlike MCMC and MC3, HMC uses not only the current position, but also a sense of momentum, to draw future samples. An introduction to HMC can be read [here](http://arxiv.org/abs/1701.02434)
#'
#'
#' This implementations assumes that the momentum is drawn from a normal distribution with mean 0 and identity covariance matrix (p ~ N (0, I) )
#'
#' @param pdf pdf Probability Density Function of the target distribution.Takes 2 arguments: a vector with a position, and, optionally, a boolean determining whether to return the log density instead.
#' @param start Vector. Starting point for the sampler
#' @param epsilon Size of the leapfrog step
#' @param L Number of leapfrog steps per iteration
#' @param iterations Number of times the sampler runs
#'
#' @return Chain with a history of visited places
#' @export
#'
#' @examples
#' # pdfs can be made easily using the make_*_pdf helpers
#' pd_func <- make_distr6_pdf(distr6::MultivariateNormal$new(mean = c(2,5)))
#' hmc_results <- sampler_hmc(pd_func, start = c(0,0), epsilon = 1, L = 10, iterations=10)
#'
#'

sampler_hmc <- function(pdf, start, epsilon = .5, L = 10, iterations=1024){
  # initialize variables -------------
  log_pdf <- function(x){return(pdf(x, log = TRUE))}
  dim <- length(start)
  chain <- matrix(nrow = iterations, ncol=dim)
  momentums <- matrix(nrow = iterations, ncol=dim)
  acceptances <- 0
  chain[1,] <- start
  
  # Run Sampler --------------
  for (i in 2:iterations){
    # draw a sample of momentum from N(0, M)
    momentum <- as.vector(mvtnorm::rmvnorm(1, rep(0,dim), sigma = diag(dim)))
    # initialize vars
    theta_prime <- chain[i-1,]
    momentum_prime <- momentum
    # Leapfrog Integrator for each L step
    leapInfo <- leapfrog_step(theta_prime, momentum_prime, epsilon, log_pdf, L)
    theta_prime <- leapInfo[[1]]
    momentum_prime <- leapInfo[[2]]
    
    # Metropolis - Hastings Acceptance, using the joint density of position + momentum
    top <- exp(joint_d(theta_prime, momentum_prime, log_pdf))
    bottom <- exp(joint_d(chain[i-1, ], momentum, log_pdf))
    
    alpha <- top/bottom
    
    if (stats::runif(1, 0,1) <= top/bottom){
      chain[i,] <- theta_prime
      momentums[i,] <- momentum_prime
      acceptances <- acceptances + 1
    } else{
      chain[i,] <- chain[i-1,]
      momentums[i,] <- momentums[i-1,]
    }
    
  }
  
  return(list(
    chain = chain,
    momentums = momentums[-1,],
    acceptance_ratio = acceptances / iterations))
  
}


#' Joint Density Calculator
#'
#' Helper function for Hamiltonian Samplers that calculates the log joint density of a position and momentum pair.
#'
#' @param theta Vector. Position
#' @param momentum Vector, momentum.
#' @param logf function. Log probability density function  of the target distribution
#'
#' @return Numeric. Log joint density of the theta/momentum pair
#' @export
#'
#' @examples
#'
#' # This function is internally used
#'
#' target <- distr6::Normal$new()
#' log_func <- function (x){return(target$pdf(log=TRUE, data = matrix(x, nrow=1)))}
#' theta <- target$mean()
#' momentum <- as.vector(mvtnorm::rmvnorm(1, 0))
#' joint_d(theta, momentum, log_func)
#'
joint_d <- function(theta, momentum, logf){
  
  return (logf(theta) - .5 * momentum %*% momentum)
}

#' No U-Turn Sampler.
#'
#' Adapted from Hoffman and Gelman (2014). The No U-Turn Sampler (NUTS) aims to eliminate the need to set a number of steps L that is present in Hamiltonian Monte Carlo, which may lead to undesirable behaviour in HMC if not set correctly.NUTS does so by recursively building a set of candidate points that span the target distribution, and stopping when it starts to double back (hence its name). More information can be found [here](https://arxiv.org/abs/1111.4246)
#'
#'
#' @param pdf Probability Density Function of the target distribution.Takes 2 arguments: a vector with a position, and, optionally, a boolean determining whether to return the log density instead.
#' @param start starting point
#' @param epsilon step size. If left to the default value, a suitable step-size will be estimated
#' @param delta_max Measure of the required accuracy of the simulation. The authors recommend a large value (1000)
#' @param iterations Times the sampler runs
#' @param iterations_adapt burn-in period
#' @param delta desired acceptance rate (default at 0.6)
#'
#' @return matrix of size iterations x dimensions with the points visited by the sampler.
#' @export
#'
#' @examples
#' pd_func <- make_distr_pdf(distr::Norm())
#' nuts_results <- sampler_nuts(pd_func, start = 0, iterations = 20)

# sampler_nuts <- function(pdf, start, epsilon = "estimate", delta_max = 1000, iterations = 1024, iterations_adapt = floor(iterations / 10), delta = .60){
#   # BuildTree Function -----------------
#   build_tree <- function(theta, momentum, u, v, j, epsilon, theta_0, momentum_0,  logf, delta_max){
#     print(paste("Tree Depth:"))
#     if (j == 0){
#       # base case - take one leapfrog step in direction v
#       lf <- leapfrog_step(theta, momentum, epsilon = v * epsilon, logf)
#       theta_prime <- lf[[1]]
#       momentum_prime <- lf[[2]]
#       # is the new point in the slice?
#       n_prime <- as.numeric(u <= exp(joint_d(theta_prime, momentum_prime, logf)))
#       # is the simulation wildly inaccurate?
#       s_prime <- as.numeric(u < exp(delta_max + joint_d(theta_prime, momentum_prime, logf)))
#       # acceptance probability
#       alpha <- min(1, exp(joint_d(theta_prime, momentum_prime, logf) - joint_d(theta_0,momentum_0, logf)))
#       return(list(theta_prime,momentum_prime, theta_prime, momentum_prime, theta_prime, n_prime, s_prime, alpha, 1))
#     } else {
#
#       # recursion - implicitly build the left and right subtrees
#       bt <- build_tree(theta, momentum, u, v, j - 1, epsilon, theta_0, momentum_0, logf, delta_max)
#       theta_minus <- bt[[1]]
#       momentum_minus <- bt[[2]]
#       theta_plus <- bt[[3]]
#       momentum_plus <- bt[[4]]
#       theta_prime <- bt[[5]]
#       n_prime <- bt[[6]]
#       s_prime <- bt[[7]]
#       alpha_prime <- bt[[8]]
#       n_prime_alpha <- bt[[9]]
#       # was the stopping criteria met in the first subtree?
#       if (s_prime == 1){
#
#         if (v == -1){
#           bt <- build_tree(theta_minus, momentum_minus, u, v, j - 1, epsilon, theta_0, momentum_0, logf, delta_max)
#           theta_minus <- bt[[1]]
#           momentum_minus <- bt[[2]]
#           theta_prime_2 <- bt[[5]]
#           n_prime_2 <- bt[[6]]
#           s_prime_2 <- bt[[7]]
#           alpha_prime_2 <- bt[[8]]
#           n_prime_alpha_2 <- bt[[9]]
#         } else{
#           bt <- build_tree(theta_plus, momentum_plus, u, v, j - 1, epsilon, theta_0, momentum_0, logf, delta_max)
#           theta_plus <- bt[[3]]
#           momentum_plus <- bt[[4]]
#           theta_prime_2 <- bt[[5]]
#           n_prime_2 <- bt[[6]]
#           s_prime_2 <- bt[[7]]
#           alpha_prime_2 <- bt[[8]]
#           n_prime_alpha_2 <- bt[[9]]
#         }
#         # which subtree to propagate a sample from
#         if (stats::runif(1,0,1) <= (n_prime_2 / max(1, (n_prime + n_prime_2)))) {
#           theta_prime <- theta_prime_2
#         }
#         # update acceptance probability``
#         alpha_prime <- alpha_prime + alpha_prime_2
#         n_prime_alpha <- n_prime_alpha + n_prime_alpha_2
#         # update stopping criterion
#         diff <- theta_plus - theta_minus
#         s_prime <- s_prime_2 *
#           as.numeric(diff %*% momentum_minus >= 0) *
#           as.numeric(diff %*% momentum_plus >= 0)
#         # Update number of valid points
#         n_prime <- n_prime + n_prime_2
#       }
#       return(list(theta_minus,momentum_minus, theta_plus, momentum_plus, theta_prime, n_prime, s_prime, alpha_prime, n_prime_alpha))
#     }
#   }
#
#   # Initialize Variables --------------------------------
#   # use the log pdf henceforth
#   logf <- function(x) {
#     return(pdf(x, log = TRUE))
#     }
#
#   # decide whether epsilon needs to be estimated
#   if (!is.numeric(epsilon)){
#     epsilon <- estimate_epsilon(start, logf)
#   }
#
#   mu <- log(10 * epsilon)
#   epsilonbar <- 1
#   Hbar <- 0
#   gamma <- .05
#   t_0 <- 10
#   kappa <- .75
#   dim <- length(start)
#   chain <- matrix( nrow=iterations, ncol = dim)
#   chain[1, ] <- start
#   epsilons <- vector(mode = "numeric", length = iterations)
#   epsilons[1] <- epsilon
#
#   # Run Sampler
#   for (iter in 2:iterations){
#     print(paste("Iteration:", iter))
#     momentum_0 <- as.vector(mvtnorm::rmvnorm(1, rep(0,dim), diag(dim))) # resample momentum
#
#     u <- stats::runif(1, 0, exp(logf(chain[iter-1,]) - .5* momentum_0 %*% momentum_0))
#     # Initialize Variables
#     theta_minus <- chain[iter-1,]
#     theta_plus <- chain[iter-1,]
#     momentum_minus <- momentum_0
#     momentum_plus <- momentum_0
#     j <- 0
#     chain[iter, ] <- chain[iter-1, ]
#     n <- 1
#     s <- 1
#     while (s == 1) {
#       v_j <- sample(c(-1,1), 1) # random direction
#       if (v_j == -1){
#         bt <- build_tree(theta_minus, momentum_minus, u, v_j, j, epsilon, chain[iter-1,], momentum_0, logf, delta_max)
#         theta_minus <- bt[[1]]
#         momentum_minus <- bt[[2]]
#         theta_prime <- bt[[5]]
#         n_prime <- bt[[6]]
#         s_prime <- bt[[7]]
#         alpha <- bt[[8]]
#         n_alpha <- bt[[9]]
#       } else {
#         bt <- build_tree(theta_plus, momentum_plus, u, v_j, j, epsilon, chain[iter-1,], momentum_0, logf, delta_max)
#         theta_plus <- bt[[3]]
#         momentum_plus <- bt[[4]]
#         theta_prime <- bt[[5]]
#         n_prime <- bt[[6]]
#         s_prime <- bt[[7]]
#         alpha <- bt[[8]]
#         n_alpha <- bt[[9]]
#       }
#
#       if (s_prime == 1){
#         if (stats::runif(1,0,1) <= (n_prime / n)){
#           chain[iter,] <- theta_prime
#         }
#       }
#       n <- n + n_prime
#       d_theta <- theta_plus - theta_minus
#       s <- s_prime *
#         as.numeric((d_theta %*% momentum_minus) >= 0) *
#         as.numeric((d_theta %*% momentum_plus) >= 0)
#       j <- j + 1
#     }
#     if (iter <= iterations_adapt){
#       m <- iter - 1
#       Hbar <- (1 - (1 / (m + t_0))) * Hbar + (1/ (m + t_0)) * (delta  - alpha / n_alpha)
#       epsilon = exp(mu - sqrt(m) / gamma * Hbar)
#       eta = m ** -kappa
#       epsilonbar = exp((1. - eta) * log(epsilonbar) + eta * log(epsilon))
#
#     } else{
#       epsilon <- epsilonbar
#     }
#   }
#
#
#   return(chain)
# }
#
sampler_nuts <- function(pdf, start, epsilon, iterations = 1024, delta_max = 1000){
  # BuildTree Function -----------------
  build_tree <- function(theta, momentum, u, v, j, epsilon){
    if (j == 0){
      # base case - take one leapfrog step in direction v
      lf <- leapfrog_step(theta, momentum, epsilon = v * epsilon, logf)
      theta_prime <- lf[[1]]
      momentum_prime <- lf[[2]]
      
      # is the new point in the slice?
      n_prime <- as.numeric(u <= exp(joint_d(theta_prime, momentum_prime, logf)))
      # is the simulation wildly inaccurate?
      s_prime <- as.numeric(log(u, base = 10) - delta_max < joint_d(theta_prime, momentum_prime, logf))
      
      return(list(theta_prime,momentum_prime, theta_prime, momentum_prime, theta_prime, n_prime, s_prime))
      
    } else {
      # recursion - implicitly build the left and right subtrees
      bt <- build_tree(theta, momentum, u, v, j - 1, epsilon)
      theta_minus <- bt[[1]]
      momentum_minus <- bt[[2]]
      theta_plus <- bt[[3]]
      momentum_plus <- bt[[4]]
      theta_prime <- bt[[5]]
      n_prime <- bt[[6]]
      s_prime <- bt[[7]]
      # was the stopping criteria met in the first subtree?
      if (s_prime == 1){
        if (v == -1){
          bt <- build_tree(theta_minus, momentum_minus, u, v, j - 1, epsilon)
          theta_minus <- bt[[1]]
          momentum_minus <- bt[[2]]
          theta_prime_2 <- bt[[5]]
          n_prime_2 <- bt[[6]]
          s_prime_2 <- bt[[7]]
        } else{
          bt <- build_tree(theta_plus, momentum_plus, u, v, j - 1, epsilon)
          theta_plus <- bt[[3]]
          momentum_plus <- bt[[4]]
          theta_prime_2 <- bt[[5]]
          n_prime_2 <- bt[[6]]
          s_prime_2 <- bt[[7]]
        }
        # which subtree to propagate a sample from
        if (stats::runif(1,0,1) <= (n_prime_2 / (n_prime + n_prime_2))) {
          theta_prime <- theta_prime_2
        }
        
        # update stopping criterion
        diff <- theta_plus - theta_minus
        s_prime <- s_prime_2 *
          as.numeric(diff %*% momentum_minus >= 0) *
          as.numeric(diff %*% momentum_plus >= 0)
        # Update number of valid points
        n_prime <- n_prime + n_prime_2
      }
      return(list(theta_minus,momentum_minus, theta_plus, momentum_plus, theta_prime, n_prime, s_prime))
    }
  }
  
  # Initialize Variables --------------------------------
  
  # use the log pdf henceforth
  logf <- function(x) {
    return(pdf(x, log = TRUE))
  }
  
  dim <- length(start)
  chain <- matrix( nrow=iterations, ncol = dim)
  chain[1, ] <- start
  
  # Run Sampler -------------------
  for (iter in 2:iterations){
    # resample momentum
    momentum_0 <- as.vector(mvtnorm::rmvnorm(1, rep(0,dim), diag(dim)))
    # get a slice u
    u <- stats::runif(1, 0, exp(logf(chain[iter-1,]) - .5* momentum_0 %*% momentum_0))
    
    # Initialize Variables
    theta_minus <- chain[iter-1,]
    theta_plus <- chain[iter-1,]
    momentum_minus <- momentum_0
    momentum_plus <- momentum_0
    j <- 0
    chain[iter, ] <- chain[iter-1, ]
    n <- 1
    s <- 1
    
    while (s == 1) {
      # random direction
      v_j <- sample(c(-1,1), 1)
      
      if (v_j == -1){
        bt <- build_tree(theta_minus, momentum_minus, u, v_j, j, epsilon)
        theta_minus <- bt[[1]]
        momentum_minus <- bt[[2]]
        theta_prime <- bt[[5]]
        n_prime <- bt[[6]]
        s_prime <- bt[[7]]
        # alpha <- bt[[8]]
        # n_alpha <- bt[[9]]
      } else {
        bt <- build_tree(theta_plus, momentum_plus, u, v_j, j, epsilon)
        theta_plus <- bt[[3]]
        momentum_plus <- bt[[4]]
        theta_prime <- bt[[5]]
        n_prime <- bt[[6]]
        s_prime <- bt[[7]]
        # alpha <- bt[[8]]
        # n_alpha <- bt[[9]]
      }
      
      if (s_prime == 1){
        if (stats::runif(1,0,1) <= (n_prime / n)){
          chain[iter,] <- theta_prime
        }
      }
      n <- n + n_prime
      d_theta <- theta_plus - theta_minus
      s <- s_prime *
        as.numeric((d_theta %*% momentum_minus) >= 0) *
        as.numeric((d_theta %*% momentum_plus) >= 0)
      j <- j + 1
    }
  }
  
  return(chain)
}


#' Epsilon estimator for Hamiltonian-based samplers.
#'
#' Adapted from Hoffman and Gelman (2014). This function estimates step size (epsilon) for a Hamiltonian-based sampler (e.g. HMC or NUTS).
#'
#' This function uses the following heuristic: Starting with an epsilon of 1, do a leapfrog step with a randomly sampled momentum p ~ N(0, I). If the joint density of the new position and momentum pair is not at least half the joint density of the starting position and momentum pair, double epsilon and try again.-
#'
#' @param theta Vector. Start position of the sampler
#' @param logf Function - given theta, returns the its log probability
#'
#' @return Epsilon to use in a Hamiltonian sampler
#' @export
#'
#' @examples
#' target <- distr6::Normal$new()
#' log_func <- function (x){return(target$pdf(log=TRUE, data = matrix(x, nrow=1)))}
#' estimate_epsilon(0, log_func)
estimate_epsilon <-function(theta, logf){
  # initialize vars
  epsilon <- 1
  momentum <- as.vector(mvtnorm::rmvnorm(1, rep(0,length(theta)), diag(length(theta))))
  lf <- leapfrog_step(theta, momentum, epsilon, logf)
  theta_prime <- lf[[1]]
  momentum_prime <- lf[[2]]
  
  k <- 1
  while (is.infinite(logf(theta_prime)) | any(is.infinite(rootSolve::gradient(logf, theta)))){
    k <- .5 * k
    lf <- leapfrog_step(theta, momentum, epsilon * k, logf)
    theta_prime <- lf[[1]]
    momentum_prime <- lf[[2]]
  }
  
  top <- exp(joint_d(theta_prime, momentum_prime, logf))
  bottom <- exp(joint_d(theta, momentum, logf))
  x <- top/bottom
  alpha <- 2 * as.numeric(x > .5) - 1
  
  while (x ** alpha > 2 ** (-alpha)){
    epsilon <- epsilon * (2 ** alpha)
    lf <- leapfrog_step(theta, momentum, epsilon, logf)
    theta_prime <- lf[[1]]
    momentum_prime <- lf[[2]]
    
    top <- exp(joint_d(theta_prime, momentum_prime, logf))
    bottom <- exp(joint_d(theta, momentum, logf))
    x <- top/bottom
  }
  return(epsilon)
}