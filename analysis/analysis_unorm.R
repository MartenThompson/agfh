rm(list=ls())

setwd('~/Git/agfh/')
source('hyptest/u_normal.R')

slug <- 'analysis/output/'

output.files <- list.files(slug)
n <- length(output.files)
hyptest.out <- data.frame(
  file=rep(NA, n),
  err_gen=rep(NA,n),
  D=rep(NA,n),
  theta_var=rep(NA,n),
  SW_p=rep(NA,n),
  SW_stat=rep(NA,n),
  KS_p=rep(NA,n),
  KS_stat=rep(NA,n),
  CM_p=rep(NA,n),
  CM_stat=rep(NA,n),
  AD_p=rep(NA,n),
  AD_stat=rep(NA,n)
)


for (i in 1:n) {
  filename <- strsplit(output.files[i], '_')[[1]]
  err.name <- filename[2]
  D <- substr(filename[3], 2, nchar(filename[3]))
  theta.var <- substr(filename[4], 5, nchar(filename[4]))
  
  hyptest.out$file[i] <- output.files[i]
  hyptest.out$err_gen[i] <- err.name
  hyptest.out$D[i] <- D
  hyptest.out$theta_var[i] <- theta.var
  
  mhg.out <- readRDS(paste0('analysis/output/', output.files[i], '/mhg_out_goodthetastarts.RData'))
  g <- mhg.out$g.final

  test <- test_u_normal(g$S, 'SW')
  hyptest.out$SW_p[i] <- test$p.value
  hyptest.out$SW_stat[i] <- test$statistic
  test <- test_u_normal(g$S, 'KS')
  hyptest.out$KS_p[i] <- test$p.value
  hyptest.out$KS_stat[i] <- test$statistic
  test <- test_u_normal(g$S, 'CM')
  hyptest.out$CM_p[i] <- test$p.value
  hyptest.out$CM_stat[i] <- test$statistic
  test <- test_u_normal(g$S, 'AD')
  hyptest.out$AD_p[i] <- test$p.value
  hyptest.out$AD_stat[i] <- test$statistic
}


hyptest.long <- data.frame(tidyr::pivot_longer(hyptest.out, cols=c('SW_p', 'KS_p', 'CM_p', 'AD_p'), names_to = c('test')))


library(ggplot2)


facet.labels <- c('beta'='Symmetric Beta',
                     'betaassym'='Asymmetric Beta',
                     'gamma'='Gamma',
                     'nullgen'='Normal')

ggplot(dat=hyptest.long) + 
  geom_histogram(aes(value, fill=test), bins=10, alpha=0.75, position='dodge') + 
  facet_wrap(.~err_gen, labeller = as_labeller(facet.labels)) +
  geom_line(data=data.frame(x=rep(0.05,2), y=c(0,10)), aes(x=x,y=y), color='red', lty=2) + 
  labs(x='p value', y='frequency') + 
  scale_y_continuous(expand=c(0,0))+
  scale_fill_discrete(name='',
                      breaks=c('AD_p', 'CM_p', 'KS_p', 'SW_p'),
                      labels=c('Anderson\nDarling','Cramer\nvon Mises','Kolmogorov\nSmirnov','Shapiro\nWilks')) + 
  theme_bw() + 
  theme(legend.position='bottom', panel.grid = element_blank())




##### 
mhg.out <- readRDS(paste0('../health_climate/code/output/main_mort_endc_ssp3-70_85_1/robjects/mhg_out_1e2_1e2thin_goodthetastarts.RData'))
g <- mhg.out$g.final
gS <- matrix(g$S, ncol = 1)

test_u_normal(gS, 'JB')
test_u_normal(gS, 'SW')
#hyptest.out$SW_p[i] <- test$p.value
#hyptest.out$SW_stat[i] <- test$statistic
test_u_normal(gS, 'KS')
#hyptest.out$KS_p[i] <- test$p.value
#hyptest.out$KS_stat[i] <- test$statistic
test_u_normal(gS, 'CM')
#hyptest.out$CM_p[i] <- test$p.value
#hyptest.out$CM_stat[i] <- test$statistic
test_u_normal(gS, 'AD')
#hyptest.out$AD_p[i] <- test$p.value
#hyptest.out$AD_stat[i] <- test$statistic
