rm(list=ls())

setwd('~/Git/agfh/')
source('hyptest/u_normal.R')

slug <- 'analysis/output/'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Reading 'Rep' Artifacts ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

output.files <- list.files(slug)
n <- length(output.files)

concat_hyptest_out <- function(hyptest.output, new.vals) {
  if(is.na(hyptest.output)) {
    hyptest.output <- data.frame(
      file=c(new.vals$file),
      err_gen=c(new.vals$err_gen),
      D=c(new.vals$D),
      theta_var=c(new.vals$theta_var),
      JB_p=c(new.vals$JB_p),
      JB_stat=c(new.vals$JB_stat),
      SW_p=c(new.vals$SW_p),
      SW_stat=c(new.vals$SW_stat),
      KS_p=c(new.vals$KS_p),
      KS_stat=c(new.vals$KS_stat),
      CM_p=c(new.vals$CM_p),
      CM_stat=c(new.vals$CM_stat),
      AD_p=c(new.vals$AD_p),
      AD_stat=c(new.vals$AD_stat)
    )
  } else {
    temp <- data.frame(
      file=c(new.vals$file),
      err_gen=c(new.vals$err_gen),
      D=c(new.vals$D),
      theta_var=c(new.vals$theta_var),
      JB_p=c(new.vals$JB_p),
      JB_stat=c(new.vals$JB_stat),
      SW_p=c(new.vals$SW_p),
      SW_stat=c(new.vals$SW_stat),
      KS_p=c(new.vals$KS_p),
      KS_stat=c(new.vals$KS_stat),
      CM_p=c(new.vals$CM_p),
      CM_stat=c(new.vals$CM_stat),
      AD_p=c(new.vals$AD_p),
      AD_stat=c(new.vals$AD_stat)
    )
    hyptest.output <- rbind(hyptest.output, temp)
  }
  return(hyptest.output)
}

hyptest.out <- NA
for (i in 1:n) {
  filename <- strsplit(output.files[i], '_')[[1]]
  
  if ('rep' == filename[1] & grepl('iter1000', output.files[i]) &
      grepl('r30', output.files[i]) ) {
    cat(output.files[i], '\n')
    
    err.name <- filename[2]
    D <- substr(filename[3], 2, nchar(filename[3]))
    
    if ('p1b1' == filename[4]) {
      theta.var <- substr(filename[5], 5, nchar(filename[5]))
    } else {
      theta.var <- substr(filename[4], 5, nchar(filename[4]))
    }
    
    for(r in 1:30) {
      mhg.out <- readRDS(paste0('analysis/output/', output.files[i], '/mhg_out_goodthetastarts', r,'.RData'))
      g <- mhg.out$g.final

      
      test <- test_u_normal(g$S, 'JB')
      JB_p <- test$p.value
      JB_stat <- test$statistic
      test <- test_u_normal(g$S, 'SW')
      SW_p <- test$p.value
      SW_stat <- test$statistic
      test <- test_u_normal(g$S, 'KS')
      KS_p <- test$p.value
      KS_stat <- test$statistic
      test <- test_u_normal(g$S, 'CM')
      CM_p <- test$p.value
      CM_stat <- test$statistic
      test <- test_u_normal(g$S, 'AD')
      AD_p <- test$p.value
      AD_stat <- test$statistic
      
      new.vals <- list(file = output.files[i],
                       err_gen = err.name,
                       D = D,
                       theta_var = theta.var,
                       JB_p=JB_p,
                       JB_stat=JB_stat,
                       SW_p=SW_p,
                       SW_stat=SW_stat,
                       KS_p=KS_p,
                       KS_stat=KS_stat,
                       CM_p=CM_p,
                       CM_stat=CM_stat,
                       AD_p=AD_p,
                       AD_stat=AD_stat)
      
      hyptest.out <- concat_hyptest_out(hyptest.out, new.vals)
    }
    
  }
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Nominal Coverage Rates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_nominal_rates <- function(hyptest.out, err.gen, alpha) {
  relevant.data <- hyptest.out[hyptest.out$err_gen == err.gen,]  
  
  return(data.frame(
    JB = mean(relevant.data$JB_p < alpha),
    SW = mean(relevant.data$SW_p <  alpha),
    KS = mean(relevant.data$KS_p <  alpha),
    CM = mean(relevant.data$CM_p <  alpha),
    AD = mean(relevant.data$AD_p <  alpha)
  ))
}

# n.alpha <- 10
# nominal.rates <- data.frame(
#   err_name <- rep(NA, n.alpha),
#   alpha <- rep(NA, n.alpha),
#   JB_rate <- rep(NA, n.alpha),
#   SW_rate <- rep(NA, n.alpha),
#   KS_rate <- rep(NA, n.alpha),
#   CM_rate <- rep(NA, n.alpha),
#   AD_rate <- rep(NA, n.alpha)
# )
# alphas <- seq(1e-3, 1, length.out=n.alpha)
# 
# for (i in 1:n.alpha) {
#   alpha <- alphas[i]
#   get_nominal_rates(hyptest.out, 'beta', alpha)
#   get_nominal_rates(hyptest.out, 'betaassym', alpha)
#   get_nominal_rates(hyptest.out, 'gamma', alpha)
#   get_nominal_rates(hyptest.out, 'nullgen', alpha)
#   
#   nominal.rates
# }

alpha <- 0.05
nom.rates <- data.frame()
nom.rates <- rbind(nom.rates, c('Beta', get_nominal_rates(hyptest.out, 'beta', alpha)))
colnames(nom.rates)[1] <- 'err_gen'
nom.rates <- rbind(nom.rates, c(err_gen='Beta Assym', get_nominal_rates(hyptest.out, 'betaassym', alpha)))
nom.rates <- rbind(nom.rates, c(err_gen='Gamma', get_nominal_rates(hyptest.out, 'gamma', alpha)))
nom.rates <- rbind(nom.rates, c(err_gen='Normal', get_nominal_rates(hyptest.out, 'nullgen', alpha)))
nom.rates

xtable::xtable(nom.rates, digits = 2, align = c('l', rep('c',6)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### P Value Histograms ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


hyptest.long <- data.frame(tidyr::pivot_longer(hyptest.out, cols=c('JB_p', 'SW_p', 'KS_p', 'CM_p', 'AD_p'), names_to = c('test')))


library(ggplot2)

facet.labels <- c('beta'='Symmetric Beta',
                  'betaassym'='Asymmetric Beta',
                  'gamma'='Gamma',
                  'nullgen'='Normal')

g <- ggplot(dat=hyptest.long) + 
  geom_histogram(aes(value, fill=test), bins=10, alpha=0.75, position='dodge') + 
  facet_wrap(.~err_gen, labeller = as_labeller(facet.labels)) +
  geom_line(data=data.frame(x=rep(0.05,2), y=c(0,250)), aes(x=x,y=y), color='red', lty=2) + 
  labs(x='p value', y='frequency') + 
  scale_y_continuous(expand=c(0,0))+
  scale_fill_discrete(name='',
                      breaks=c('JB_p', 'AD_p', 'CM_p', 'KS_p', 'SW_p'),
                      labels=c('Jarque\nBera', 'Anderson\nDarling','Cramer\nvon Mises','Kolmogorov\nSmirnov','Shapiro\nWilks')) + 
  theme_bw() + 
  theme(legend.position='bottom', panel.grid = element_blank())

g

png(paste0(slug, 'sampling_err/', 'rep30exantruns.png'),  height=6, width=6, units='in', res = 100)
g
dev.off()


