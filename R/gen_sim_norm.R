library(MASS)
gen_sim_norm <- function(dat=NULL,
                         alpha=NULL,beta=NULL,
                         Sigma_u=NULL,Sigma_w=NULL,sigma_e=NULL){
  
  nsubj<-length(unique(dat$subj))
  u<-mvrnorm(n=nsubj,
             mu=c(0,0),Sigma=Sigma_u)
  u<-data.frame(subjid=unique(dat$subj),u)
  ## item random effects
  nitem<-length(unique(dat$item))
  w<-mvrnorm(n=nitem,
             mu=c(0,0),Sigma=Sigma_w)
  w<-data.frame(itemid=unique(dat$item),w)
  
  ## generate data row by row:
  N<-dim(dat)[1]
  simbase100<-rep(NA,N)
  for(i in 1:N){
    current_subjid<-dat[i,]$subj
    current_itemid<-dat[i,]$item
    
    simbase100[i] <- rnorm(1,alpha +
                             u[u$subjid==current_subjid,2] +
                             w[w$itemid==current_itemid,2] +
                             (beta+
                                u[u$subjid==current_subjid,3] +
                                w[w$itemid==current_itemid,3])*dat$zcloze[i],
                           sigma_e)}
  dat$simbase100<-simbase100
  simdat<-dat
  simdat
}
