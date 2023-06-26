require(goftest)

test_u_normal <- function(samples, test='JB') {
  check_inputs(samples, test)

  if ('JB' == test) { # Jarque-Bera
    return(jarque_bera(samples))
  } else if ('SW' == test) { # Shapiro-Wilk
    return( shapiro_wilk(samples))
  } else if ('KS' == test) { # Kolmogorov-Smirnov
    return(kolmogorov_smirnov(samples))
  } else if ('CM' == test) { # Cramer von Mises
    return(cramer_vonmises(samples))
  } else if ('AD' == test) { # Anderson Darling
    return(anderson_darling(samples))
  }
}

check_inputs <- function(samples, test) {
  if (!is.vector(samples) && ncol(samples) > 1) {
    stop('Samples should be vector, not matrix. Actual dim:', dim(samples))
  }
  if (!(test %in% c('JB', 'SW', 'KS', 'CM', 'AD'))) {
    stop('Uknown test type:', test, '. Should be one of (JB, SW, KS, CM, AD).')
  }
}

jarque_bera <- function(samples) {
  n <- length(samples)
  m1 <- mean(samples)
  m2 <- sum((samples-m1)^2)/n
  m3 <- sum((samples-m1)^3)/n
  m4 <- sum((samples-m1)^4)/n
  S <- m3/(m2^(3/2))
  K <- m4/(m2^2)
  stat <- (n/6)*(S^2 + (1/4)*(K-3)^2)
  pval <- 1 - pchisq(stat, df = 2)

  name <- 'Jarque Bera'
  return(list(statistic = stat,
              p.value = pval,
              name = name)
         )

}

shapiro_wilk <- function(samples) {
  out <- shapiro.test(samples)
  return(list(statistic = out$statistic,
              p.value = out$p.value,
              name = 'Shapiro Wilk'))
}

kolmogorov_smirnov <- function(samples) {
  out <- ks.test(samples, 'pnorm')
  return(list(statistic = out$statistic,
              p.value = out$p.value,
              name = 'Kolmogorov Smirnov'))
}

cramer_vonmises <- function(samples) {
  out <- cvm.test(samples, 'pnorm', mean=0, sd=1, estimated=FALSE)
  return(list(statistic = out$statistic,
              p.value = out$p.value,
              name = 'Cramer Von Mises'))
}

anderson_darling <- function(samples) {
  out <- ad.test(samples, 'pnorm', mean=0, sd=1, estimated=FALSE)
  return(list(statistic = out$statistic,
              p.value = out$p.value,
              name = 'Anderson Darling'))
}
