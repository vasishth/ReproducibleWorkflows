library(MASS)
gen_sim_norm <- function(dat=NULL,
                           alpha=NULL,beta=NULL,
                           Sigma_u=NULL,Sigma_w=NULL,sigma_e=NULL){
  ## subject random effects:
  nsubj<-length(unique(dat$subj))
  u<-mvrnorm(n=nsubj,
             mu=c(0,0),Sigma=Sigma_u)
  u<-data.frame(subjid=unique(dat$subj),u)
  ## item random effects
  nitem<-length(unique(dat$item))
  w<-mvrnorm(n=nitem,
             mu=c(0,0),Sigma=Sigma_w)
  
  ## generate data row by row:
  N<-dim(dat)[1]
  simbase100<-rep(NA,N)
  for(i in 1:N){
    simbase100[i] <- rnorm(1,alpha +
                      u[dat[i,]$subj,1] +
                      w[dat[i,]$item,1] +
                      (beta+u[dat[i,]$subj,2]+
                         w[dat[i,]$item,2])*dat$zcloze[i],
                    sigma_e)}
  dat$simbase100<-simbase100
  simdat<-dat
  simdat
  }
