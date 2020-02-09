
###############################################################################
### Compute Bayes factors after a replication study for correlations 
### to find out whether the original findings are replicated
### Last revised: 30-11-2014
### Authors: Josine Verhagen, josineverhagen@gmail.com, www.josineverhagen.com
### Eric- Jan wagenmakers, ejwagenmakers@gmail.com
### INPUT:  r original, n original, r replication, n replication
### OUTPUT: Replication Bayes factor, JZS Bayes factors one and two sided,
### plots of the prior and posterior distribution for teh single studies and 
### the replication Bayes factor 

### This is the file with supporting functions
###############################################################################

###Required packages
library(hypergeo)


#In this file: 
# Replication test
# CorrelationReplicationBF : This function will calculate the Replication Bayes factor (Verhagen & Wagenmakers, 2014) for correlations
# corBF.beta: JZS Bayes Factor
# corBF.beta.1s : one-sided JZS Bayes Factor

#Supporting functions:
#PosteriorRho : Posterior for the correlation coefficient #Jeffreys (1961), p. 175, but normalized 
#LLrho :  p(rho|r,n)  

# Plotting functions: 
#repposteriorplot: 
# This function will plot the posterior when the posterior of the original study is used as prior distribution
# for the effect size in a replicated study. Also the outcome of a Savage Dickey density ratio 
#(Verdinelli & Wasserman, 1995) is given based on the height of the prior and posterior at H0: rho= 0. 


#posteriorplotBeta11 : 
#This function will plot the posterior when a uniform distribution is used as prior distribution.
#Also the outcome of a Savage Dickey density ratio (Verdinelli & Wasserman, 1995) is given based on the height of
#the prior and posterior at H0: rho= 0. One sided or two-sided. 

#repposteriorplot2 and posteriorplotBeta112 provide possibilities to plot multiple
# plots in a windows, removing axis labels for the inner plots (la = , lg = ) 

# Replication test 
CorrelationReplicationBF = function(r.orig, n.orig, r.rep, n.rep)
{ 
  #This function will calculate the Replication Bayes factor (Verhagen & Wagenmakers, 2014) for correlations
  #Replace n with n-1
  posterior.proportional <- function(rho) { PosteriorRho(rho,(n.orig-1),r.orig) *LLrho(rho,(n.rep-1),r.rep) }    
  marginal.likelihood   <- integrate(posterior.proportional, lower=-1, upper=1)$value
  ReplicationPosteriorRho <- function(rho) {posterior.proportional(rho)/marginal.likelihood} 
  BF01 <- 1/integrate(posterior.proportional, lower=-1, upper=1)$value
  BF10 <- 1/BF01
  
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  return(BFS)
}

# JZE Bayes factor: two sided 

corBF.beta = function(r, n)
{  
  #This function will calculate the Regular Bayes factor for correlations
  #Based on a uniform prior distribution
  #Replace n with n-1
  n<- n-1
  integrand <- function(rho){LLrho(rho,n,r) * (.5* dbeta(abs(rho),1,1)) }
  marginal.likelihood   <- integrate(integrand, lower=-1, upper=1)$value
  BF01 <- 1/marginal.likelihood
  BF10 <- 1/BF01
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  return(BFS)
}

# JZS Bayes factor: one sided
corBF.beta.1s = function(r, n, side)
{  
  #This function will calculate the Regular Bayes factor for correlations
  #Based on a uniform prior distirbution
  #Replace n with n-1
  n<- n-1
  integrand <- function(rho){LLrho(rho,n,r) * dbeta (abs(rho),1,1) }
  if( side == "upper") {
    integrand <- function(rho){LLrho(rho,n,r) * dbeta (rho,1,1) }
    marginal.likelihood   <- integrate(integrand, lower=0, upper=1)$value }
  if( side == "lower") {
    integrand <- function(rho){LLrho(rho,n,r) * dbeta (-rho,1,1) }
    marginal.likelihood   <- integrate(integrand, lower=-1, upper=0)$value }
  BF01 <- 1/marginal.likelihood
  BF10 <- 1/BF01
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  return(BFS)
}


#Function for the posterior of the correlation coefficient
#Creates the density for each value: but now we need to sample from it

PosteriorRho = function(rho, n, r)
{
  #This function will calculate the posterior density for a correlation 
  #coefficient (Jeffreys, 1961, p. 175), normalized. Note that if the means
  #are subtracted, n needs to be replaced by n-1
  #Author: Eric Jan Wagenmakers 
  hypgeo    <- hypergeo(.5, .5, n+.5, .5+.5*(r*rho))
  if (n <= 160){
    prior.times.likelihood  <- (LLrho(rho,n,r)) * (sqrt(pi/2)*gamma(n)/gamma(n+.5)) * hypgeo
    integrand <- function(rho){as.double({as.double((LLrho(rho,n,r))*(sqrt(pi/2)
                                                                      *gamma(n)/gamma(n+.5)) * hypergeo(.5, .5, n+.5, .5+.5*(r*rho)))})}
  }
  # The gamma() function cannot handle large numbers, so we need to work around it 
  # with a beta function. This CAN lead to an error in the integral in case of  
  # high r and high n, however 
  # Author: Josine Verhagen 
  if (n>160){
    prior.times.likelihood  <- (LLrho(rho,n,r))*(sqrt(pi/2)*beta(n, .5)/gamma(.5))*hypgeo
    integrand <- function(rho){as.double({as.double((LLrho(rho,n,r))*(sqrt(pi/2)
                                                                      *beta(n, .5)/gamma(.5)) * hypergeo(.5, .5, n+.5, .5+.5*(r*rho)))})}
  }
  marginal.likelihood   <- integrate(integrand, lower=-1, upper=1)$value
  posterior.density      <- prior.times.likelihood/marginal.likelihood
  return(as.numeric(posterior.density))
}

# Basic likelihood function for correlation

LLrho <- function(rho,n,r) {
  n <- n-1
  (1-rho^2)^(n/2) /(1-rho*r)^(n-.5)
}


# Posterior plot with beta prior

posteriorplotBeta11 <- function(n, r, side = "two", max.y = 0, BFIN = "no", BF01 = 0)
{
  # This function will plot the posterior when a uniform distribution is used as prior distribution.
  # Also the outcome of a Savage Dickey density ratio (Verdinelli & Wasserman, 1995) is given based on the height of
  # the prior and posterior at H0: rho= 0. 
  
  seqrho <- seq(-1,1, by = .01) 
  
  par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
      font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
  
  if(side == "two"){  
    if ( max.y == 0) {
      vmax.y <- max(PosteriorRho(seqrho,(n-1),r))
      max.y <- vmax.y + vmax.y/5
    }
    
    plot ( function(rho) PosteriorRho(rho,(n-1),r), -1, 1 
           , xlim = c(-1,1), ylim = c(0,max.y), lwd=2, lty=1, col= 1,
           ylab="Density", xlab=" ", main = " ") 
    axis(1)
    axis(2)
    mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=2)
    
    par ( new = TRUE)
    
    plot( function(x) .5*dbeta(abs(x),1,1), xlim = c(-1,1), ylim = c(0,max.y), lty =3, lwd = 2 , main = " ", xlab=" ", ylab="Density")
    
    prior0 <- .5
    post0 <- PosteriorRho(0,(n-1),r)
    
    xleg = -1.15
    xleg2 <- 0.5
  }
  
  if( side == "upper") {
    area <- integrate(PosteriorRho, lower=0, upper=1,n=(n-1), r=r)$value    
    if ( max.y == 0) {
      vmax.y <- max(PosteriorRho(seqrho[101:200],(n-1),r))*(1/area)
      max.y <- vmax.y + vmax.y/5
    }  
    
    plot ( function(rho) PosteriorRho(rho,(n-1),r)*(1/area), 0, 1 
           , xlim = c(0,1), ylim = c(0,max.y), lwd=2, lty=1, col= 1,
           ylab="Density", xlab=" ", main = " ") 
    axis(1)
    axis(2)
    mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=2)
    par ( new = TRUE)
    plot( function(x) dbeta(abs(x),1,1), xlim = c(0,1), ylim = c(0,max.y), lty =3, lwd = 2 , main = " ", xlab=" ", ylab="Density")
    
    prior0 <- 1
    post0 <- PosteriorRho(0,(n-1),r) * (1/area)
    
    xleg = -0.075
    xleg2 <- .75
  }
  if( side == "lower") {
    area <- integrate(PosteriorRho, lower=-1, upper=0,n=(n-1), r=r)$value  
    
    if ( max.y == 0) {
      vmax.y <- max(PosteriorRho(seqrho[0:101],(n-1),r))*(1/area) 
      max.y <- vmax.y + vmax.y/5
    }  
    
    plot ( function(rho) PosteriorRho(rho,(n-1),r)*(1/area), -1, 0 
           , xlim = c(-1,0), ylim = c(0,max.y), lwd=2, lty=1, col= 1,
           ylab="Density", xlab=" ", main = " ") 
    axis(1)
    axis(2)
    mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=2)
    par ( new = TRUE)
    plot( function(x) dbeta(abs(x),1,1), xlim = c(-1,0), ylim = c(0,max.y), lty =3, lwd = 2 , main = " ", xlab=" ", ylab="Density")
    
    prior0 <- 1
    post0 <- PosteriorRho(0,(n-1),r) * (1/area)
    
    xleg = -1.075
    xleg2 <- -.25
  }
  
  par ( new = FALSE)
  
  points(0, prior0 , pch=21, cex=2, bg="grey")
  points(0, post0, pch=21, cex=2, bg="grey")
  
  legend(xleg,(max.y + max.y/10), c("Prior","Posterior"), lwd= 2, lty=c(3,1) ,bty = 'n' , cex= 1.5)  
  
  
  if (BFIN == "no") {
    BF01 <- post0/prior0
  }
  BF10 <- 1/BF01
  
  
  if( BF10 >= 1){ 
    text(xleg2, (max.y), labels = substitute(paste(BF[10], " = ", v), list(v=signif(BF10, digits=3))) , cex = 1.3)
  }
  if( BF10 < 1 & BF10 > .000501){ 
    text(xleg2, (max.y) , labels = substitute(paste(BF[10], " = ", v), list(v=round(BF10,digits=2))) , cex = 1.3)
  }
  if( BF10 < .000501){ 
    text(xleg2, (max.y) , labels = substitute(paste(BF[10], " = ", v), list(v= "0.000")) , cex = 1.3)
  }
  
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  
  return(BFS)
}


repposteriorplot <- function(n.orig, r.orig, n.rep, r.rep, BFIN = "no", BF01 = 0)
{
  # This function will plot the posterior when the posterior of the original study is used as prior distribution
  # for the effect size in a replicated study. Also the outcome of a Savage Dickey density ratio 
  #(Verdinelli & Wasserman, 1995) is given based on the height of the prior and posterior at H0: rho= 0. 
  
  posterior.proportional <- function(rho) {PosteriorRho(rho,(n.orig-1),r.orig) *LLrho(rho,(n.rep-1),r.rep) }    
  
  marginal.likelihood   <- integrate(posterior.proportional, lower=-1, upper=1)$value
  ReplicationPosteriorRho <- function(rho) {posterior.proportional(rho)/marginal.likelihood} 
  
  par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
      font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
  
  my <- 10
  mix <- -1
  max <- 1
  
  plot ( function(rho) PosteriorRho(rho,(n.orig-1),r.orig), -1, 1 
         , xlim = c(mix,max), ylim = c(0,my), lwd=2, lty=3, col= 1,
         ylab="Density", xlab=" ", main = " ") 
  axis(1)
  axis(2)
  mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=2)
  par ( new = TRUE)
  plot ( function(rho) ReplicationPosteriorRho(rho), -1, 1, 
         xlim = c(mix,max), ylim = c(0,my), lwd=2, lty=1, col= 1,
         ylab="Density", xlab=" ", main = " ") 
  
  par ( new = FALSE)
  
  post0 <- ReplicationPosteriorRho(0)
  prior0 <- PosteriorRho(0,(n.orig-1),r.orig)
  
  points(0, prior0 , pch=21, cex=2, bg="grey")
  points(0, post0, pch=21, cex=2, bg="grey")
  
  legend(mix , my+1 , c("Prior","Posterior"), lwd= 2, lty=c(3,1) ,bty = 'n' , cex= 1.5)  
  
  if (BFIN == "no") {
    BF01 <- post0/prior0
  }
  
  BF10 <- 1/BF01
  
  if( BF10 >= 1){ 
    text(.5, my , labels = substitute(paste(BF[10], " = ", v), list(v=signif(BF10, digits=3))) , cex = 1.3)
  }
  if( BF10 < 1 & BF10 > .000501){ 
    text(.5, my , labels = substitute(paste(BF[10], " = ", v), list(v=round(BF10,digits=2))) , cex = 1.3)
  }
  if( BF10 < .000501){ 
    text(.5, my , labels = substitute(paste(BF[10], " = ", v), list(v= "0.000")) , cex = 1.3)
  }
  
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  
  
  return(BFS)
}


repposteriorplot2 <- function(n.orig, r.orig, n.rep, r.rep, title = "",lg= 1,la=1, max.y =0 )
{
  posterior.proportional <- function(rho) {PosteriorRho(rho,(n.orig-1),r.orig) *LLrho(rho,(n.rep-1),r.rep) }    
  marginal.likelihood   <- integrate(posterior.proportional, lower=-1, upper=1)$value
  ReplicationPosteriorRho <- function(rho) {posterior.proportional(rho)/marginal.likelihood} 
  
  my <- 12
  if ( max.y == 0 ) { my == max.y}
  
  mix <- -1
  max <- 1
  
  if ( la ==1 ) {
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
        font.lab = 2, cex.axis = 1.3, bty = "n", las=1, mai=c(0.7,0.7,0.3,0.1))
    
    plot ( function(rho) PosteriorRho(rho,(n.orig-1),r.orig), -1, 1 
           , xlim = c(mix,max), ylim = c(0,my), lwd=2, lty=3, col= 1,
           ylab="Density", xlab=" ", main = title) 
    
    axis(2)
  }
  
  if ( la == 0 ) {
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
        font.lab = 2, cex.axis = 1.3, bty = "n", las=1, mai=c(0.7,0.2,0.3,0.1))
    
    
    plot ( function(rho) PosteriorRho(rho,(n.orig-1),r.orig), -1, 1 
           , xlim = c(mix,max), ylim = c(0,my), lwd=2, lty=3, col= 1,
           ylab="", xlab=" ", , main = title, yaxt='n') 
  }
  axis(1)
  mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=1)
  par ( new = TRUE)
  
  plot ( function(rho) ReplicationPosteriorRho(rho), -1, 1, 
         xlim = c(mix,max), ylim = c(0,my), lwd=2, lty=1, col= 1,
         ylab="Density", xlab=" ", main = " ", yaxt='n') 
  
  
  par ( new = FALSE)
  
  post0 <- ReplicationPosteriorRho(0)
  prior0 <- PosteriorRho(0,(n.orig-1),r.orig)
  
  
  points(0, prior0 , pch=21, cex=2, bg="grey")
  points(0, post0, pch=21, cex=2, bg="grey")
  
  if (lg == 1) {
    legend(mix , my+1 , c("Prior","Posterior"), lwd= 2, lty=c(3,1) ,bty = 'n' , cex= 1.5)  
  }
  
  BF01 <- post0/prior0
  BF10 <- 1/BF01
  
  if( BF10 >= 1){ 
    text(.5, my , labels = substitute(paste(BF[10], " = ", v), list(v=signif(BF10, digits=3))) , cex = 1.3)
  }
  if( BF10 < 1 & BF10 > .000501){ 
    text(.5, my , labels = substitute(paste(BF[10], " = ", v), list(v=round(BF10,digits=2))) , cex = 1.3)
  }
  if( BF10 < .000501){ 
    text(.5, my , labels = substitute(paste(BF[10], " = ", v), list(v= "0.000")) , cex = 1.3)
  }
  
  
  
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  
  
  return(BFS)
}


posteriorplotBeta112 <- function(n, r, side = "two", title, lg= 1,la=1, max.y = 0)
{
  # This function will plot the posterior when a uniform distribution is used as prior distribution.
  # Also the outcome of a Savage Dickey density ratio (Verdinelli & Wasserman, 1995) is given based on the height of
  # the prior and posterior at H0: rho= 0. 
  
  seqrho <- seq(-1,1, by = .01) 
  
  if ( la ==1 ) {
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
        font.lab = 2, cex.axis = 1.3, bty = "n", las=1, mai=c(0.7,0.7,0.3,0.1))
    taxis = "Density" 
    vyaxt <- NULL
  }  
  
  if ( la == 0 ) {
    
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
        font.lab = 2, cex.axis = 1.3, bty = "n", las=1, mai=c(0.7,0.2,0.3,0.1))
    taxis = " " 
    vyaxt <- "n"
  }
  
  
  if(side == "two"){  
    if ( max.y == 0) {
      vmax.y <- max(PosteriorRho(seqrho,(n-1),r))
      max.y <- vmax.y + vmax.y/5
    }
    
    plot ( function(rho) PosteriorRho(rho,(n.orig-1),r.orig), -1, 1 
           , xlim = c(mix,max), ylim = c(0,my), lwd=2, lty=3, col= 1,
           ylab="Density", xlab=" ", main = " " ,yaxt = vyaxt) 
    
    if ( la ==1 ) { axis(2)}
    axis(1)
    mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=1)
    
    par ( new = TRUE)
    
    plot( function(x) .5*dbeta(abs(x),1,1), xlim = c(-1,1), ylim = c(0,max.y), lty =3, lwd = 2 , main = " ", xlab=" ", ylab="Density",yaxt = vyaxt)
    
    prior0 <- .5
    post0 <- PosteriorRho(0,(n-1),r)
    
    xleg = -1.15
    xleg2 <- 0.5
  }
  
  if( side == "upper") {
    area <- integrate(PosteriorRho, lower=0, upper=1,n=(n-1), r=r)$value    
    if ( max.y == 0) {
      vmax.y <- max(PosteriorRho(seqrho[101:200],(n-1),r))*(1/area)
      max.y <- vmax.y + vmax.y/5
    }  
    
    plot ( function(rho) PosteriorRho(rho,(n-1),r)*(1/area), 0, 1 
           , xlim = c(0,1), ylim = c(0,max.y), lwd=2, lty=1, col= 1,
           ylab=taxis, xlab=" ", main = title,yaxt = vyaxt) 
    if ( la ==1 ) { axis(2)}
    axis(1)
    mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=1)
    par ( new = TRUE)
    plot( function(x) dbeta(abs(x),1,1), xlim = c(0,1), ylim = c(0,max.y), lty =3, lwd = 2 , main = " ", xlab=" ", ylab="Density",yaxt = vyaxt)
    
    prior0 <- 1
    post0 <- PosteriorRho(0,(n-1),r) * (1/area)
    
    xleg = -0.075
    xleg2 <- .75
  }
  
  
  if( side == "lower") {
    area <- integrate(PosteriorRho, lower=-1, upper=0,n=(n-1), r=r)$value  
    
    if ( max.y == 0) {
      vmax.y <- max(PosteriorRho(seqrho[0:101],(n-1),r))*(1/area) 
      max.y <- vmax.y + vmax.y/5
    }  
    
    plot ( function(rho) PosteriorRho(rho,(n-1),r)*(1/area), -1, 0 
           , xlim = c(-1,0), ylim = c(0,max.y), lwd=2, lty=1, col= 1,
           ylab= taxis, xlab=" ", main = title,yaxt = vyaxt) 
    if ( la ==1 ) { axis(2)}
    axis(1)
    mtext(expression(Correlation ~ rho), side=1, line = 2.8, cex=1)
    par ( new = TRUE)
    plot( function(x) dbeta(abs(x),1,1), xlim = c(-1,0), ylim = c(0,max.y), lty =3, lwd = 2 , main = " ", xlab=" ", ylab="Density",yaxt = vyaxt)
    
    prior0 <- 1
    post0 <- PosteriorRho(0,(n-1),r) * (1/area)
    
    xleg = -1.075
    xleg2 <- -.25
  }
  
  par ( new = FALSE)
  
  points(0, prior0 , pch=21, cex=2, bg="grey")
  points(0, post0, pch=21, cex=2, bg="grey")
  
  if (lg == 1) {
    legend(xleg,(max.y + max.y/10), c("Prior","Posterior"), lwd= 2, lty=c(3,1) ,bty = 'n' , cex= 1.5)  
  }
  
  
  BF01 <- post0/prior0
  BF10 <- 1/BF01
  
  
  if( BF10 >= 1){ 
    text(xleg2, (max.y), labels = substitute(paste(BF[10], " = ", v), list(v=signif(BF10, digits=3))) , cex = 1.3)
  }
  if( BF10 < 1 & BF10 > .000501){ 
    text(xleg2, (max.y) , labels = substitute(paste(BF[10], " = ", v), list(v=round(BF10,digits=2))) , cex = 1.3)
  }
  if( BF10 < .000501){ 
    text(xleg2, (max.y) , labels = substitute(paste(BF[10], " = ", v), list(v= "0.000")) , cex = 1.3)
  }
  
  BFS<- c(BF01,BF10)
  names(BFS) <- c("BF01","BF10")
  
  return(BFS)
}
