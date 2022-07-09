# PARALLEL NIMBLE ON CGRB
# SGE_Batch -c "R CMD BATCH bear_nimble_cgrb.R bear_nimble_cgrb_OUT" -r bear_nimble_cgrb -q otter -P 3


rm(list=ls())

run_js_par <- function(seed){
  
  load("bear_nimble_data.RData")
  
  data<-list(y=yaug,X=X,xlim=xlim,ylim=ylim)
  constants<-list(M=M,J=J,K=K)
  monitors<-c('N','sigma','p0','psi')
  inits<-function() list(z=rep(1,M),p0=0.1,psi=0.5,sigma=runif(1,2,3),s=cbind(rep(mean(X[,1]),M),rep(mean(X[,2]),M)))
  
  library(nimble)
  
  bear_nimble_code<-nimbleCode({
    
    psi~dbeta(1,1)
    p0~dbeta(1,1)
    sigma~dunif(0,20)
    sigma2<-sigma^2
    
    for(i in 1:M){
      z[i]~dbern(psi)
      s[i,1]~dunif(xlim[1], xlim[2])
      s[i,2]~dunif(ylim[1], ylim[2])
      
      for(j in 1:J){
        d2[i,j]<-(s[i,1]-X[j,1])^2 + (s[i,2]-X[j,2])^2
        p[i,j]<-p0*exp(-d2[i,j]/2*sigma2)*z[i]
        y[i,j]~dbin(p[i,j],K)
      }
    }
    
    N<-sum(z[1:M])
    
  })
  
  bear_nimble<-nimbleMCMC(code=bear_nimble_code,constants=constants,data=data,monitors=monitors,inits=inits,
                          niter=100000,nburnin=50000,nchains=1,progressBar=T,samplesAsCodaMCMC=T,summary=F,setSeed = seed)  
  
  return(bear_nimble)
}

this_cluster <- parallel::makeCluster(3, type="SOCK")

bear_nimble_parallel <- parallel::parLapply(cl = this_cluster, X = 1:3, 
                                            fun = runNimbleParallel)

parallel::stopCluster(this_cluster)

library(coda)
bear_nimble_results <- as.mcmc.list(bear_nimble_parallel)
summary(bear_nimble_results)
gelman.diag(bear_nimble_results, multivariate = F)

save(bear_nimble_results, file="bear_nimble_results.RData")