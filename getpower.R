################################################################################
## Copyright (C) 2010 Peter Murakami <peter.murakami@gmail.com>
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

#Functions for calculating power, for given sample sizes.

########################################################################
mufunc = function(eta,family,offset=0){
    if(family=="gaussian")      mu= eta
    if(family=="binomial")      mu= exp(eta)/(1+exp(eta))
    #if(family=="quasibinomial")
    if(family=="poisson") mu= exp(eta)*exp(offset)
    #if(family=="quasipoisson")
    return(mu)
}

#Instead of calling glm() below, can calculate the stats directly to make it go faster, but i haven't implemented that yet.

#n is the sample size.
#n.sim is number of iterations to run.
#covariate is the values of the covariate. It can be a numeric or integer vector of length n (e.g., rep(c(0,1),n/2)), or an expression for obtaining a random sample (e.g., expression(rnorm(n,0,2)) or expression(rpois(n,3))).
#b0 is true value of beta_0.  Values that are not extreme make the random samples that get generated more representative. E.g., if we're only taking a sample of 30 and lambda in the rpois is .005, the resulting sample will not be too good.  Also, in logistic regression for example, setting b0 determines what probabiilities the desired OR corresponds to.  E.g., if you want to detect an OR difference of 1.25 that corresponds to p0=.05 and p1=.062, set b1=log(1.25) and b0=-2.995.
#b1 is smallest value of b1 you want to be able detect.
#alpha is the significance level at which you reject the null hypothesis.
#getps tells whether you want the function to return all the n p-values for b1 (getps=TRUE) or just the power (getps=FALSE).
#fam can be "gaussian","binomial", or "poisson"
#Must specify sigma.y if fam="gaussian"
#Offset, if given, must already by logged, since that's what glm takes.
getpower= function(n,n.sim=2000,covariate,b1,b0=0,alpha=.05,getps=FALSE,fam,sigma.y,Offset=rep(0,n),...){
    if(fam=="gaussian" & missing(sigma.y)) stop("need sigma.y for gaussian regression")
    
    if(!inherits(covariate,c("numeric","integer","expression"))) stop("covariate must be numeric or integer vector, or expression.")
    if(inherits(covariate,c("numeric","integer"))){
        stopifnot(n==length(covariate))
        mu = mufunc(eta=b0+b1*covariate,family=fam,offset=Offset)
    }
    if(inherits(covariate,"expression")) {
        if(length(eval(covariate))!=n) stop("covariate incorrectly specified.")
    }
 
    pv = vector("numeric",length=n.sim)
    pb <- txtProgressBar(min=1, max=n.sim, initial=0, style=3)
    for(i in 1:n.sim){
        if(inherits(covariate,"expression")) {
            covariate = eval(covariate)
            mu = mufunc(eta=b0+b1*covariate,family=fam,offset=Offset)
        }
        
        if(fam=="gaussian") y = rnorm(n=n,mean=mu,sd=sigma.y)
        if(fam=="poisson")  y = rpois(n=n,lambda=mu)
        if(fam=="binomial") y = rbinom(n=n,size=1,prob=mu)
        
        if(fam!="poisson") reg= summary(glm(y~covariate,family=fam,...))
        if(fam=="poisson") reg= summary(glm(y~covariate,family=fam,offset=Offset,...))
        #canonical link by default.
        
        pv[i] = reg[[12]]["covariate",4]
        setTxtProgressBar(pb, i)
    }
    close(pb)

    if(getps)  return(pv)
    if(!getps) return(mean(pv<=alpha))
}

##################################################################
##################################################################
#Code from chapter 20 of Gelman/Hill for random intercept and coefficient model.
#See page 450 for the model statements. In this situations there are treatment=0 and treatment=1 subjects, 1/2 of each, followed over time.  Function could be generalized to allow for other distributions of the independent variable, but that's not implemented for now.  Note that random effects are assumed to be normally distributed.
#The call to lmer() below is the bottleneck.

#J is number of subjects
#K is number of measurements/subject
#fam can be "gaussian","binomial", or "poisson". Link functions will be the canonical ones.
#mu.a.true is true average value of random intercept. For gaussian outcomes, this is the true link(outcome) at t=0
#g.0.true is true average slope for the controls
#g.1.true is smallest true average effect of the treatment at each time you want to be able to detect, except time=0 since treatment doesn't affect the slope in this model since this model assumes treatment can have no effect at time=0.
#sigma.y.true is true sd of the outcome. Needed if outcome is gaussian.
#sigma.a.true is true sd of the random intercept.
#sigma.b.true is true sd of the random coefficient.
#rho is the true correlation between intercepts and slopes. 0 (independence, since they're gaussian) by default.
#Offsets must already be logged, as in glm
getdata <- function(J,K,family,mu.a.true,g.0.true,g.1.true,sigma.y.true,sigma.a.true,sigma.b.true,offset=0,Rho=0){
  require(splus2R) # need for rmvnorm()
  time       <- rep(seq(0,1,length=K), J) # K measurements during the year
  person     <- rep(1:J, each=K)          # person ID's
  treatment  <- sample(rep (0:1, J/2))    # J total, half 0, half 1, in random order.
  treatment1 <- treatment[person]

  randoms <- rmvnorm(n=J, mean=cbind(mu.a.true, g.0.true+g.1.true*treatment), sd=cbind(rep(sigma.a.true,J), rep(sigma.b.true,J)), rho=Rho)
  a.true <- randoms[,1]
  b.true <- randoms[,2]
  zeta <- a.true[person] + b.true[person]*time
  #Data:
  if(family=="gaussian"){
      y <- rnorm(J*K, zeta, sigma.y.true)
  }
  if(family=="binomial"){
      p1 <- exp(zeta)/(1+exp(zeta)) #zeta is logit(p1)
      y  <- rbinom(n=J*K, size=1, prob=p1)
  }
  if(family=="poisson"){
      mu <- exp(zeta)*exp(offset) #aka lambda
      y  <- rpois(J*K, mu)
  }
  return(data.frame(y, time, person, treatment1))
}
getpower.mlm <- function (J,K,n.sims=1000,fam="gaussian",mu.a.true,g.0.true,g.1.true,sigma.y.true,sigma.a.true,sigma.b.true,Offset=rep(0,J*K),rho=0,...){
  if(J%%2==1) stop("J should be even in order for group sizes to be equal.")
  require(lme4)
  require(arm) # need for se.fixef()
  signif <- rep (NA, n.sims)
  pb <- txtProgressBar(min=1, max=n.sims, initial=0, style=3)
  for (s in 1:n.sims){
    fake <- getdata(J, K, family=fam, mu.a.true=mu.a.true, g.0.true=g.0.true, g.1.true=g.1.true, sigma.y.true=sigma.y.true, sigma.a.true=sigma.a.true, sigma.b.true=sigma.b.true,offset=Offset,Rho=rho)
    if(fam=="gaussian") {
        lme.power = lmer(y ~ time + time:treatment1 +(1 + time | person), data=fake, ...)
    } else if(fam=="poisson") {
        lme.power = glmer(y ~ time + time:treatment1 +(1 + time | person), data=fake, family=fam, offset=Offset, ...)
    } else {
        lme.power = glmer(y ~ time + time:treatment1 +(1 + time | person), data=fake, family=fam, ...)
    }
    theta.hat <- fixef(lme.power)["time:treatment1"]
    theta.hat.se  <- se.fixef(lme.power)["time:treatment1"]
    signif[s] <- (theta.hat - 2*theta.hat.se) > 0 #there's a question about what the p-values should be in such models (because how many degrees of freedom should be used is in question), so just use t quantile= 2 which should be approximately right, for alpha=.05.
    #can't get t-stat or p-value from lme.power directly, and summary() takes too long..
    setTxtProgressBar(pb, s)
  }
  close(pb)
  power <- mean(signif)
  return(power)
}

########################################################################################
########################################################################################
#Random intercept model:
#This can also be used to get power for gee gaussian model with uniform correlation structures, since a gaussian marginal model with uniform correlation structure is equivalent to a gaussian model with just a random intercept.

#y ~ N( a.true + beta_1*covariate, sigma.y.true)
#where a.true = mu.a.true + gamma_i ~ N(mu.a.true, sigma.a.true)
#covariate should be numeric or integer vector of length K, or an expression for obtaining a random sample of size J*K.
getdata.ri <- function(J,K,family,cvs,b1,mu.a.true,sigma.a.true,sigma.y.true,offset=0){
  if(!inherits(cvs,c("numeric","integer","expression"))) stop("covariate must be numeric or integer vector, or expression.")
  if(inherits(cvs,c("numeric","integer"))) cv = rep(cvs,J)
  if(inherits(cvs,"expression"))           cv = eval(cvs)
  if(length(cv)!=J*K) stop("covariate incorrectly specified.")
  cluster <- rep(1:J, each=K) # person ID's

  a.true <- rnorm(J, mu.a.true, sigma.a.true)
  zeta   <- a.true[cluster] + b1*cv
  #Data:
  if(family=="gaussian"){
      y <- rnorm(J*K, zeta, sigma.y.true)
  }
  if(family=="binomial"){
      p1 <- exp(zeta)/(1+exp(zeta)) #zeta is logit(p1)
      y  <- rbinom(n=J*K, size=1, prob=p1)
  }
  if(family=="poisson"){
      mu <- exp(zeta)*exp(offset) #aka lambda
      y  <- rpois(J*K, mu)
  }
  return(data.frame(y, cluster, cv))
}
getpower.mlm.ri <- function (J,K,n.sims=1000,fam="gaussian",covariate,B1,mu.a.true,sigma.a.true,sigma.y.true,Offset=rep(0,J*K),...){
  require(lme4)
  require(arm) # need for se.fixef()
  signif <- rep (NA, n.sims)
  pb <- txtProgressBar(min=1, max=n.sims, initial=0, style=3)
  for (s in 1:n.sims){
    fake <- getdata.ri(J, K, family=fam, cvs = covariate, b1=B1, mu.a.true=mu.a.true, sigma.a.true=sigma.a.true, sigma.y.true=sigma.y.true, offset=Offset)
    if(fam=="gaussian") {
        lme.power = lmer(y ~ cv + (1 | cluster), data=fake,...)
    } else if(fam=="poisson") {
        lme.power = glmer(y ~ cv + (1 | cluster), data=fake, family=fam, offset=Offset,...)
    } else {
        lme.power = glmer(y ~ cv + (1 | cluster), data=fake, family=fam,...)
    }
    theta.hat <- fixef(lme.power)["cv"]
    theta.hat.se  <- se.fixef(lme.power)["cv"]
    signif[s] <- (theta.hat - 2*theta.hat.se) > 0 #there's a question about what the p-values should be in such models (because how many degrees of freedom should be used is in question), so just use t quantile= 2 which should be approximately right, for alpha=.05.
    #can't get t-stat or p-value from lme.power directly, and summary() takes too long..
    setTxtProgressBar(pb, s)
  }
  close(pb)
  power <- mean(signif)
  return(power)
}
