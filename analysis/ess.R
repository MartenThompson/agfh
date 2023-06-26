rm(list=ls())

require(mcmcse)

setwd('~/Git/agfh/')

string_reverse <- function(st) {
  paste0(rev(strsplit(st,'')[[1]]), collapse = '')
}

ess <- function(samples) {
  if(is.null(dim(samples))) {
    dim <- 1
    n <- length(samples)
  } else {
    dim <- ncol(samples)
    n <- nrow(samples)
  }
  
  samples <- matrix(samples, ncol=dim)
  
  combo.ess <- min(multiESS(samples), n)
  
  if (dim > 1) {
    indiv.ess <- apply(samples, 2, function(col){multiESS(matrix(col,ncol=1))})
  } else {
    indiv.ess <- NA
  }
  
  return(list(combo.ess = combo.ess,
              indiv.ess = indiv.ess))
}

save.slug <- 'analysis/output'
dirs <- list.dirs(save.slug)
dirs <- dirs[2:length(dirs)] # first is self dir

for (i in 1:length(dirs)) {
  run.name <- strsplit(dirs[i], '/')[[1]][3]
  run.name.bk <- string_reverse(run.name)
  
  if ('rep' == substr(run.name, 1,3) & 
      '01niht0001' == substr(run.name.bk, 1, 10)) {
    
    run.path <- paste0(save.slug, '/', run.name, '/')
    run.name.parts <- strsplit(run.name, '_')[[1]]
    data.gen.name <- run.name.parts[2]
    D <- run.name.parts[3]
    temp <- run.name.parts[4]
    if ('p' == substr(temp,1,2)) {
      lambda <- run.name.parts[5]
    } else {
      lambda <- run.name.parts[4]
    }
    
    reps <- readRDS(paste0(run.path, 'robjects/reps.RData'))$reps
    
    if (30 == reps) {
      cat(run.path, '\n')  
      
      dat1 <- readRDS(paste0(run.path, 'data1.RData'))
      dim.theta <- length(dat1$theta)
      dim.beta <- length(dat1$Beta)
      
      gibbs.thvar.ess <- rep(NA, reps)
      gibbs.beta.ess.combo<- rep(NA, reps) 
      gibbs.beta.ess.indiv <- matrix(NA, nrow=reps, ncol=dim.beta)
      gibbs.theta.ess.combo <- rep(NA, reps)
      gibbs.theta.ess.indiv <- matrix(NA, nrow=reps, ncol=dim.theta)
      agfh.thvar.ess <- rep(NA, reps)
      agfh.beta.ess.combo<- rep(NA, reps) 
      agfh.beta.ess.indiv <- matrix(NA, nrow=reps, ncol=dim.beta)
      agfh.theta.ess.combo <- rep(NA, reps)
      agfh.theta.ess.indiv <- matrix(NA, nrow=reps, ncol=dim.theta)
      
      for (r in 1:reps) {
        gibbs.samples <- readRDS(paste0(run.path, 'robjects/gibbs_out', r, '.RData'))
         
        gibbs.thvar.ess[r] <- ess(gibbs.samples$param.samples.list$theta.var)$combo.ess
        ess.beta <- ess(gibbs.samples$param.samples.list$beta)
        gibbs.beta.ess.combo[r] <- ess.beta$combo.ess
        gibbs.beta.ess.indiv[r,] <- ess.beta$indiv.ess
        ess.theta <- ess(gibbs.samples$param.samples.list$theta)
        gibbs.theta.ess.combo[r] <- ess.theta$combo.ess
        gibbs.theta.ess.indiv[r,] <- ess.theta$indiv.ess
        
        agfh.samples <- readRDS(paste0(run.path, '/mhg_out_goodthetastarts', r, '.RData'))
        
        agfh.thvar.ess[r] <- ess(agfh.samples$param.samples.list$theta.var)$combo.ess
        ess.beta <- ess(agfh.samples$param.samples.list$beta)
        agfh.beta.ess.combo[r] <- ess.beta$combo.ess
        agfh.beta.ess.indiv[r,] <- ess.beta$indiv.ess
        ess.theta <- ess(agfh.samples$param.samples.list$theta)
        agfh.theta.ess.combo[r] <- ess.theta$combo.ess
        agfh.theta.ess.indiv[r,] <- ess.theta$indiv.ess
        
      }
      ess.all <- list(gibbs.thvar=gibbs.thvar.ess,
                      gibbs.beta.combo=gibbs.beta.ess.combo,
                      gibbs.beta.indiv=gibbs.beta.ess.indiv,
                      gibbs.theta.combo=gibbs.theta.ess.combo,
                      gibbs.theta.indiv=gibbs.theta.ess.indiv,
                      agfh.thvar=agfh.thvar.ess,
                      agfh.beta.combo=agfh.beta.ess.combo,
                      agfh.beta.indiv=agfh.beta.ess.indiv,
                      agfh.theta.combo=agfh.theta.ess.combo,
                      agfh.theta.indiv=agfh.theta.ess.indiv)
      
      if (is.na(mean(gibbs.beta.ess.indiv))) {
        cat('Gibbs:', gibbs.beta.ess.indiv,
            '\n AGFH:', agfh.beta.ess.indiv)
      }
      
      gibbs.ess.summary <- data.frame(thvar.avg=round(mean(gibbs.thvar.ess)),
                                thvar.sd = round(sd(gibbs.thvar.ess)),
                                beta.combo.avg=round(mean(gibbs.beta.ess.combo)),
                                beta.combo.sd=round(sd(gibbs.beta.ess.combo)),
                                beta.indiv.avg=round(mean(gibbs.beta.ess.indiv)),
                                beta.indiv.sd=round(sd(gibbs.beta.ess.indiv)),
                                theta.combo.avg=round(mean(gibbs.theta.ess.combo)),
                                theta.combo.sd=round(sd(gibbs.theta.ess.combo)),
                                theta.indiv.avg=round(mean(gibbs.theta.ess.indiv)),
                                theta.indiv.sd=round(sd(gibbs.theta.ess.indiv)))
      
      agfh.ess.summary <- data.frame(thvar.avg=round(mean(agfh.thvar.ess)),
                                      thvar.sd = round(sd(agfh.thvar.ess)),
                                      beta.combo.avg=round(mean(agfh.beta.ess.combo)),
                                      beta.combo.sd=round(sd(agfh.beta.ess.combo)),
                                      beta.indiv.avg=round(mean(agfh.beta.ess.indiv)),
                                      beta.indiv.sd=round(sd(agfh.beta.ess.indiv)),
                                      theta.combo.avg=round(mean(agfh.theta.ess.combo)),
                                      theta.combo.sd=round(sd(agfh.theta.ess.combo)),
                                      theta.indiv.avg=round(mean(agfh.theta.ess.indiv)),
                                      theta.indiv.sd=round(sd(agfh.theta.ess.indiv)))
      
      #saveRDS(ess.all, paste0(run.path, '/robjects/ess_all.RData'))
      #write.csv(gibbs.ess.summary, paste0('analysis/output/ess_rep_summary/gibbs_', run.name, '_', D, '_beta', dim.beta, '_', lambda, '.csv'), row.names = FALSE)
      #write.csv(agfh.ess.summary, paste0('analysis/output/ess_rep_summary/agfh_', run.name, '_', D, '_beta', dim.beta, '_', lambda, '.csv'), row.names = FALSE)
    }
    
    #break
  }
}

temp <- 'analysis/output/rep_beta_D01_lamb05_r30_a00.1a10.1fz0.01mhgiter1000thin10'
temp <- readRDS(paste0(temp, '/robjects/ess_all.RData'))

#gibbs.samples <- readRDS(paste0(run.path, 'robjects/gibbs_out', r, '.RData'))





