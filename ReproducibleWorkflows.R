## ----setup, include = FALSE-------------------------------------------------------
library("papaja")


## ----analysis-preferences---------------------------------------------------------
# Seed for random number generation
set.seed(42)
knitr::opts_chunk$set(cache.extra = knitr::rand_seed,
                      echo=TRUE)
library(brms)
library(lme4)
library(tidyverse)
library(ggplot2)
library(lattice)


## ----loaddata---------------------------------------------------------------------
articles <-read.delim("NieuwlandEtAl2018/public_article_data.txt", 
                      quote="")
head(articles)


## ----originalcode-----------------------------------------------------------------
# R script for the DeLong replication data.
# base100 is the N400 dependent measure with the pre-registered 100 ms baseline, which will be used throughout this script to keep it short
# base200 has the 200 ms baseline, base500 the 500 ms baseline
# filt01 has the 0.01 hz filtered data, and pre-stim the -500 to -100 ms data 
#read article data: loaded above by SV.
#articles <-read.delim("~/Desktop/DELONGR/public_article_data.txt", quote="")

# make sure data is right format
articles$item <-as.factor(articles$item)
articles$cloze <- as.numeric(as.character(articles$cloze))
## no column called n400:
##articles$base100 <- ##as.numeric(as.character(articles$n400))

# Z-transform cloze and store the mean and SD
articles$zcloze <- scale(articles$cloze, center = TRUE, scale = TRUE)
clozemean <- mean( articles$cloze, na.rm=T )
clozesd <- sd( articles$cloze, na.rm=T )


## ----articleanalysis,cache=TRUE---------------------------------------------------
# article models
model1a <- lmer(base100 ~ lab*zcloze + 
                  ( zcloze | subject) + (zcloze  | item),
                contrasts=list( lab=contr.sum(9) ), data = articles , 
                control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )

model2a <- lmer(base100 ~ lab/zcloze + 
                  ( zcloze | subject) + (zcloze  | item), 
                contrasts=list( lab=contr.sum(9) ), data = articles , 
                control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )

model3a <- lmer(base100 ~  zcloze + ( zcloze | subject) +
                  (zcloze  | item),  
                data = articles,
                control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )

model4a <- lmer(base100 ~         ( zcloze | subject) +
                  (zcloze  | item),  data = articles ,
                control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )

# compare models
anova(model1a,model2a)


## ----nounanalysis,cache=TRUE,eval=FALSE-------------------------------------------
## ########################
## #### NOUNS
## 
## #nouns <-read.delim("~/Desktop/DELONGR/public_noun_data.txt", quote="")
## nouns <-read.delim("NieuwlandEtAl2018/public_noun_data.txt", quote="")
## 
## nouns$item <-as.factor(nouns$item)
## nouns$cloze <- as.numeric(as.character(nouns$cloze))
## nouns$n400 <- as.numeric(as.character(nouns$n400))
## 
## nouns$zcloze <- scale(nouns$cloze, center = TRUE, scale = TRUE)
## clozemean <- mean( nouns$cloze, na.rm=T )
## clozesd <- sd( nouns$cloze, na.rm=T )
## 
## model1n <- lmer(n400 ~  lab*zcloze +
##                   ( zcloze | subject) + ( zcloze  | item), contrasts=list( lab=contr.sum(9) ),
##                 data = nouns ,
##                 control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )
## model2n <- lmer(n400 ~  lab + zcloze +
##                   (zcloze | subject) + (zcloze  | item), contrasts=list( lab=contr.sum(9) ),
##                 data = nouns ,
##                 control=lmerControl(optCtrl=list(maxfun=1e9)), REML = FALSE )
## model3n <- lmer(n400 ~  zcloze +
##                   (zcloze | subject) + (zcloze  | item), data = nouns ,
##                 control=lmerControl(optCtrl=list(maxfun=1e9)), REML = FALSE )
## model4n <- lmer(n400 ~
##                   (zcloze | subject) + (zcloze  | item), data = nouns ,
##                 control=lmerControl(optCtrl=list(maxfun=1e9)), REML = FALSE )
## 
## # compare models
## anova(model1n,model2n)
## 


## ----model3a----------------------------------------------------------------------
model3a <- lmer(base100 ~  zcloze + 
                  ( zcloze | subject) + (zcloze  | item),  
                data = articles, 
                control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )
model3aNULL <- lmer(base100 ~  1 + 
                      ( zcloze | subject) + (zcloze  | item),  
                    data = articles, 
                    control=lmerControl(optCtrl=list(maxfun=1e5)), REML = FALSE )
anova(model3a,model3aNULL)


## ----dotplotmodel3a,fig.height=10,fig.width=8-------------------------------------
print(dotplot(ranef(model3a,condVar=TRUE)))


## ----subsetdata-------------------------------------------------------------------
dat<-dplyr::select(articles,subject,item,lab,zcloze,base100)


## ----renamesubjects---------------------------------------------------------------
## get subject names
subjnames<-unique(dat$subject)
## convert to numerical values:
subjects_numerical<-data.frame(subjname=subjnames,subj=factor(1:length(subjnames)))

## create new column with numerical id's:
dat2<-merge(dat,subjects_numerical,by.x="subject",by.y="subjname")
dat<-dat2
head(dat)


## ----bysubject,fig.width=7,fig.height=7-------------------------------------------
## save lab names:
labnames<-as.character(unique(dat$lab))

## Choose first lab's data:
current_lab<-filter(dat,lab==labnames[1])

xyplot(base100~zcloze|subj,current_lab,
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.abline(lm(y ~ x))
       },main="By subjects data")


## ----byitem,fig.width=7,fig.height=7----------------------------------------------
xyplot(base100~zcloze|item,
       current_lab,
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.abline(lm(y ~ x))
       },,main="By items data")


## ----definexyplotfunction---------------------------------------------------------
gg_xyplot <- function(x, y, formula, shape, size, 
                      xlabel="zcloze", 
                      ylabel="base100 (microvolts)",
                      data=current_lab){
    ggplot(data = data, aes(x = data[,x],
                            y = data[,y]))  +
    facet_wrap(formula) +
    geom_smooth(method="lm")+
    geom_point(color = "blue", 
               shape = shape, size = size) +
    theme(panel.grid.minor = element_blank()) +
    theme_bw() +
    ylab(ylabel) +
    xlab(xlabel)
}


## ----ggxyplotsubj,fig.width=7,fig.height=7----------------------------------------
gg_xyplot(x = "zcloze", y = "base100",  ~ subj,  
          shape = 1, size = 3, 
          data = current_lab)
          


## ----ggxyplotitem,fig.width=7,fig.height=7----------------------------------------
gg_xyplot(x = "zcloze", y = "base100",  ~ item,  
          shape = 1, size = 3, 
          data = current_lab)
          


## ----lmList-----------------------------------------------------------------------
## fit separate linear models for each subject: 
m_lmlist<-lmList(base100~zcloze|subj,dat)
lmlist_coef<-summary(m_lmlist)$coefficients
## extract by subject slopes:
slopes<-lmlist_coef[,,2]
means<-slopes[,1]
lower<-means-2*slopes[,2]
upper<-means+2*slopes[,2]

subjects<-unique(dat$subj)

slopes_summary<-data.frame(subj=subjects,
                       means=means,
                       lower=lower,upper=upper)

## reorder means by magnitude:
slopes_summary<-slopes_summary[order(slopes_summary$means),]

p<-ggplot(slopes_summary,aes(x=subjects, y=means)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                 position=position_dodge(0.05))+theme_bw()
p


## ----extractestimates-------------------------------------------------------------
## extract estimates of fixed-effects parameters:
beta<-summary(model3a)$coefficients[,1]
## extract standard deviation estimate:
sigma_e<-attr(VarCorr(model3a),"sc")
## assemble variance covariance matrix for subjects:
subj_ranefsd<-attr(VarCorr(model3a)$subj,"stddev")
subj_ranefcorr<-attr(VarCorr(model3a)$subj,"corr")
Sigma_u<-round(SIN::sdcor2cov(stddev=subj_ranefsd,
                        corr=subj_ranefcorr),6)

## assemble variance covariance matrix for items:
item_ranefsd<-attr(VarCorr(model3a)$item,"stddev")
item_ranefcorr<-attr(VarCorr(model3a)$item,"corr")
Sigma_w<-SIN::sdcor2cov(stddev=item_ranefsd,
                        corr=item_ranefcorr)


## ----generatefakedata-------------------------------------------------------------
source("R/gen_sim_norm.R")


## ----producesimulatedvalues,eval=FALSE--------------------------------------------
## ## Example simulated data, generated several times:
## nsim<-10
## sim_values<-matrix(rep(NA,nsim*dim(dat)[1]),ncol = nsim)
## for(i in 1:nsim){
## simdat<-gen_sim_norm(dat,
##                alpha=beta[1],beta=beta[2],
##                Sigma_u=Sigma_u,Sigma_w=Sigma_w,
##                sigma_e=sigma_e)
## sim_values[,i]<-simdat$simbase100
## }
## head(sim_values)


## ----comparerealsim,fig.height=5,fig.weight=5-------------------------------------
## Compare fake and simulated data:
## load precomputed simulated values:
load("data/sim_values.Rda")
## the observed data:
hist(dat$base100,freq=FALSE,
     ylim=c(0,0.06),
     main="Real vs simulated data",
     xlab="base100")
## simulated data:
for(i in 1:ncol(sim_values)){
lines(density(sim_values[,i]),lty=2)
}


## ---------------------------------------------------------------------------------
b1<-round(summary(model3a)$coefficients[2,1],3)
b1se<-round(summary(model3a)$coefficients[2,2],3)


## ---------------------------------------------------------------------------------
b1-2*b1se; b1+2*b1se


## ----poweranalysismean,cache=TRUE-------------------------------------------------
## store for power:
nsim<-100
tvals<-c()
for(i in 1:nsim){
  #print(paste("i=",i,sep=""))
simdat<-gen_sim_norm(dat,
               alpha=beta[1],
               beta=b1, ## using the estimated slope
               Sigma_u=Sigma_u,Sigma_w=Sigma_w,
               sigma_e=sigma_e)
m<-lmer(simbase100~zcloze+(1+zcloze|subj)+
          (1+zcloze|item),simdat,
        control=lmerControl(calc.derivs=FALSE))
tvals[i]<-summary(m)$coefficients[2,3]
}
power_mean<-mean(abs(tvals)>2)
power_mean
save(power_mean,file="data/power_mean.Rda")


## ----poweranalysislower,cache=TRUE,eval=FALSE-------------------------------------
## ## store for power:
## nsim<-100
## tvals<-c()
## for(i in 1:nsim){
##   #print(paste("i=",i,sep=""))
## simdat<-gen_sim_norm(dat,
##                alpha=beta[1],
##                beta=b1-2*b1se, ## using the lower bound
##                Sigma_u=Sigma_u,Sigma_w=Sigma_w,
##                sigma_e=sigma_e)
## m<-lmer(simbase100~zcloze+(1+zcloze|subj)+(1+zcloze|item),
##         simdat,
##         control=lmerControl(calc.derivs=FALSE))
## tvals[i]<-summary(m)$coefficients[2,3]
## }
## power_lower<-mean(abs(tvals)>2)
## power_lower


## ----poweranalysisupper,cache=TRUE,eval=FALSE-------------------------------------
## ## store for power:
## nsim<-100
## tvals<-c()
## for(i in 1:nsim){
##   #print(paste("i=",i,sep=""))
## simdat<-gen_sim_norm(dat,
##                alpha=beta[1],
##                beta=b1+2*b1se, ## using the upper bound
##                Sigma_u=Sigma_u,Sigma_w=Sigma_w,
##                sigma_e=sigma_e)
## m<-lmer(simbase100~zcloze+(1+zcloze|subj)+(1+zcloze|item),
##         simdat,
##         control=lmerControl(calc.derivs=FALSE))
## tvals[i]<-summary(m)$coefficients[2,3]
## }
## power_upper<-mean(abs(tvals)>2)
## power_upper


## ----datareplicationexample-------------------------------------------------------
library(data.table)
df <- data.frame(a = c(1,2,3), b = c(1,2,3))
dt <- as.data.table(df)
n <- 3
dt[rep(dt[, .I], n)]


## ----selectbirmdata---------------------------------------------------------------
birm_dat<-dplyr::select(current_lab,
                        subj,item,zcloze)
head(birm_dat)


## ----replicatebirmdata------------------------------------------------------------
birm_dat<-as.data.table(birm_dat)
subjid<-birm_dat$subj
n <- 3
repl_dat<-birm_dat[rep(birm_dat[, .I], n)]
dim(birm_dat)[1]*3
dim(repl_dat)[1]

## create new subject vector extended n times
add_id<-rep(seq(100,100*n,by=100),each=dim(birm_dat)[1])
repl_dat$subj<-rep(as.numeric(as.character(birm_dat$subj)),n)+add_id
head(repl_dat)
length(unique(birm_dat$subj))*3
length(unique(repl_dat$subj))


## ----computesamplesize,cache=TRUE-------------------------------------------------
nsim<-100
tvals<-c()
for(i in 1:nsim){
simdat<-gen_sim_norm(repl_dat,
               alpha=beta[1],
               beta=b1, ## using the estimated slope
               Sigma_u=Sigma_u,Sigma_w=Sigma_w,
               sigma_e=sigma_e)
m<-lmer(simbase100~zcloze+(1+zcloze|subj)+(1+zcloze|item),simdat,
        control=lmerControl(calc.derivs=FALSE))
tvals[i]<-summary(m)$coefficients[2,3]
}
## number of subjects:
length(unique(simdat$subj))
power_mean<-mean(abs(tvals)>2)
power_mean


## ----computepowerfunction---------------------------------------------------------
compute_power<-function(dat=birm_dat,
                        replicates=3,nsims=100){
dat<-as.data.table(dat)
subjid<-dat$subj
nrep <- replicates
repl_dat<-dat[rep(dat[, .I], nrep)]

## create new subject vector extended n times
add_id<-rep(seq(100,100*nrep,by=100),each=dim(dat)[1])
repl_dat$subj<-rep(as.numeric(as.character(birm_dat$subj)),
                   nrep)+add_id
nsim<-nsims
tvals<-c()
for(i in 1:nsim){
simdat<-gen_sim_norm(repl_dat,
               alpha=beta[1],
               beta=b1, ## using the estimated slope
               Sigma_u=Sigma_u,Sigma_w=Sigma_w,
               sigma_e=sigma_e)
m<-lmer(simbase100~zcloze+(1+zcloze|subj)+(1+zcloze|item),
        simdat,
        control=lmerControl(calc.derivs=FALSE))
tvals[i]<-summary(m)$coefficients[2,3]
}
power_mean<-mean(abs(tvals)>2)
## return result with sample size:
res<-c(length(unique(simdat$subj)),
                power_mean)
res
}


## ----computesamplesize2,cache=TRUE------------------------------------------------
repl4<-compute_power(dat=birm_dat,
              replicates = 4,
              nsims=10)
repl4
repl6<-compute_power(dat=birm_dat,
              replicates = 6,
              nsims=10)
repl6


## ----powersamplesizeplot----------------------------------------------------------
results<-rbind(repl4,repl6)
results<-data.frame(results)
colnames(results)<-c("nsubj","power estimate")
results


## ----generateranefs---------------------------------------------------------------
nsubj<-length(unique(dat$subj))
nitem<-length(unique(dat$item))
u<-mvrnorm(n=nsubj, # number of subjects
             mu=c(0,0),Sigma=Sigma_u)
w<-mvrnorm(n=nitem, # number of items
             mu=c(0,0),Sigma=Sigma_w)

## add subject and item id's:
u<-data.frame(subjid=unique(dat$subj),u)
w<-data.frame(itemid=unique(dat$item),w)



## ---------------------------------------------------------------------------------
head(u)
head(w)


## ---------------------------------------------------------------------------------
u[1,2]


## ---------------------------------------------------------------------------------
u[1,3]


## ---------------------------------------------------------------------------------
dat[1,]


## ---------------------------------------------------------------------------------
## intercept adjustment:
u[u$subjid==1,2]
## slope adjustment:
u[u$subjid==1,3]


## ---------------------------------------------------------------------------------
beta[1] + u[u$subjid==1,2] + 
  w[w$itemid==102,2]+ 
  (beta[2] + 
     u[u$subjid==1,3] +
     w[w$itemid==102,3])*dat[1,]$zcloze +
  rnorm(1,0,sigma_e)


## ---------------------------------------------------------------------------------
dat[1,]$subj
dat[1,]$item


## ---------------------------------------------------------------------------------
i<-1
current_subjid<-dat[i,]$subj
current_itemid<-dat[i,]$item

beta[1] + u[u$subjid==current_subjid,2] +
  w[w$itemid==current_itemid,2]+ 
  (beta[2] + u[u$subjid==current_subjid,3] +
      w[w$itemid==current_itemid,3])*dat[i,]$zcloze + rnorm(1,0,sigma_e)


## ---------------------------------------------------------------------------------
## number of rows
N<-dim(dat)[1]
# save vector of simulated data:
simbase100 <- rep(NA,N)
for(i in 1:N){
current_subjid<-dat[i,]$subj
current_itemid<-dat[i,]$item
simbase100[i] <-  beta[1] + u[u$subjid==current_subjid,2] +
  w[w$itemid==current_itemid,2]+ 
  (beta[2] + u[u$subjid==current_subjid,3] +
      w[w$itemid==current_itemid,3])*dat[i,]$zcloze + rnorm(1,0,sigma_e)
}


## ---------------------------------------------------------------------------------
hist(dat$base100)
lines(density(simbase100))


## ----sessioninfo------------------------------------------------------------------
sessionInfo()


## ----create_r-references,echo=FALSE-----------------------------------------------
r_refs(file = "bibliography.bib")

