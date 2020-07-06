#' Returns a list of occurrence lists for all species in a given simulation,
#' either the real or a range of occurrence times.
#'
#' \code{SampleClade} takes times of speciation and extinction, information to create a sampling rate with \code{MakeRate}, a vector of geologic time 
#' intervals and whether one wants the true return times or a range based on \code{IntVec}.
#' 
#' @param \code{S} a species number (i.e. its identity) vector to be sampled. Could be only a subset of the species if the user wishes.
#' @param \code{sim} a simulation, usually the output of \code{BDSim}.
#' @param \code{stages} a vector of time intervals corresponding to geological time ranges. If \code{returnTrue} is false, \code{SampleClade} returns the occurrence time as a range. In this way, we simulate the granularity in real world fossil records. If \code{returnTrue} is true, this is ignored.
#' @param \code{rr} a sampling rate function. May be a constant, a time-dependent function, a function dependent on time and environment, or a vector of rates corresponding to the times in \code{rshifts}. If \code{env_rr} and \code{rshifts} are NULL, it will be either treated as an exponential rate. Must be a constant if \code{dFUN} is not \code{NULL}.
#' @param \code{tmax} the maximum simulation time, used by \code{rexp_var}.
#' @param \code{env_rr} a matrix containing time points and values of anenviromental variable, like temperature, for each time point. This will be used to create a sampling rate, so \code{rr} must be a function of time and said variable if \code{env_rr} is not NULL.
#' @param \code{rshifts} vector of rate shifts. First element must be starting time for simulation (0, or tmax). Vector must have the same length as \code{rr} E.g. \code{rr = c(0.15, 0.1, 0.2)}, \code{pshifts = c(0, 20, 30)} means r will be 0.15 from 0 to 20, 0.1 from 20 to 30, and 0.2 from 30 to \code{tmax}. Note that using this method for step-function rates is currently slower than using \code{ifelse}. Using \code{c(0, 10, tmax)} is the same as \code{c(tmax, tmax - 10, 0)}.
#' @param \code{returnTrue} if set to true, a list of true ocurrence times for each species will be returned. If set to true, we call \code{Grid} to take the list of samples and return a list of upper bounds for occurrences.
#' @param \code{dFUN} A density function representing the age-dependent preservation model. It must be a density function, and consequently (1) integrate in 1 (though this condition is not verified by the function, its the user responsability to check this property), (2) should describe the density of sampling a lineage in a given point \code{t} in geological time (3) should be parameterized in absolute geological time (i.e. do not function in respect to "age", the relative time since speciation, but should be relative to absolute geological time, in Mya), (4) should be limited between \code{s} (i.e. the lineage's speciation/origination geological, absolute, time) and \code{e} (i.e. the lineage's extinction geological, absolute, time), with \code{s} > \code{e}, (5) have the arguments \code{t}, \code{s}, \code{e} and \code{sp}.
#' @param \code{dFUNmax} A function that calculates the maximum (density) value of dFUN. It can also be a number representing the maximum density.
#' @param \code{...} Additional parameters related to dFUN and dFUNmax.
#' @return a list of occurrences for that species, expected to be around \code{(Ts-TE)*rr} occurrences, with their distribution in time given by the function provided by the user.
#' @author written by Matheus Januario and Bruno do Rosario Petrucci.
#' 
#' @examples
#' Note: for more detailed examples on the construction of the age-dependen preservation model, see the detailed examples in the help page of the SAmple()  and SampleADPP() functions
#'
#' first we can try constant sampling
sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived lineage which will obscure the pattern
  sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
}
#'
r <- 100 # high so we can see the pattern
#' the resolution of the fossil dataset:
stages=seq(from=10, to = 0, # from simulation's tmax to present
           by = -.1) # note that we will provide a very high resolution (= 0.1 My) just to test the function:
#'
dt<-SampleClade(1:length(sim$TE), sim, r, tmax=10, stages=stages)
ids<-unique(dt$Species)
mids<-(dt$MaxT-dt$MinT)+dt$MinT
#'
for(i in 1:length(ids)){
  
  sp<-unique(as.numeric(gsub("spp_", "", ids[i])))
  
  hist(mids[dt$Species==ids[i]], main=paste0("spp = ", sp, "; duration ~ ", round(sim$TS[sp]-sim$TE[sp], digits = 2), "my"), 
       xlab="My", breaks=seq(round(sim$TS[i]), round(sim$TE[i]), -1), xlim=c(sim$TS[i], sim$TE[i]))
  
  mid<-par[sp]
  abline(h=r)
}
#'
#' now let us try a linearly increasing r
sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived lineage which will obscure the pattern
  sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
}
#'
r <- function(t) {
  return(200-5*t)
}
#' the resolution of the fossil dataset:
stages=seq(from=10, to = 0, # from simulation's tmax to present
           by = -.1) # note that we will provide a very high resolution (= 0.1 My) just to test the function:
#'
dt<-SampleClade(1:length(sim$TE), sim, r, tmax=10, stages=stages)
ids<-unique(dt$Species)
mids<-(dt$MaxT-dt$MinT)+dt$MinT
#'
for(i in 1:length(ids)){
  sp<-unique(as.numeric(gsub("spp_", "", ids[i])))
  
  hist(mids[dt$Species==ids[i]], main=paste0("spp = ", sp, "; duration ~ ", round(sim$TS[sp]-sim$TE[sp], digits = 2), "my"), 
       xlab="My", breaks=seq(round(sim$TS[i]), round(sim$TE[i]), -1), xlim=c(sim$TS[i], sim$TE[i]))
  mid=par[sp]
  t <- seq(sim$TE[i], sim$TS[i], 0.1)
  lines(t, rev(r(t)))
}
#'
#' sampling could be any function of time, of course, such as a step function
sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived lineage which will obscure the pattern
  sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
}
#'
rlist <- c(100, 50, 120)
rshifts <- c(0, 4, 8) # this could also be c(10, 6, 2)
#' this makes the simulation slower than when using \code{ifelse} to make a stepfunction, but it is an option
#' make it a function so we can plot it
r <- MakeRate(rlist, 10, fshifts=rshifts)
#' the resolution of the fossil dataset:
stages=seq(from=10, to = 0, # from simulation's tmax to present
           by = -.1) # note that we will provide a very high resolution (= 0.1 My) just to test the function:
#'
dt<-SampleClade(1:length(sim$TE), sim, rlist, tmax=10, rshifts=rshifts, stages=stages)
ids<-unique(dt$Species)
mids<-(dt$MaxT-dt$MinT)+dt$MinT
#'
for(i in 1:length(ids)){
  sp<-unique(as.numeric(gsub("spp_", "", ids[i])))
  
  hist(mids[dt$Species==ids[i]], main=paste0("spp = ", sp, "; duration ~ ", round(sim$TS[sp]-sim$TE[sp], digits = 2), "my"), 
       xlab="My", breaks=seq(round(sim$TS[i]), round(sim$TE[i]), -1), xlim=c(sim$TS[i], sim$TE[i]))
  mid=par[sp]
  t <- seq(sim$TE[i], sim$TS[i], 0.1)
  lines(t, rev(r(t)))
}
#'
#' finally, \code{SampleClade} also accepts an environmental variable
library(RPANDA)
data(InfTemp)
#'
sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived lineage which will obscure the pattern
  sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
}
#'
env_r <- InfTemp
#' we can make sampling dependent on the temperature
r <- function(t, env) {
  return(25*env)
}
#' this makes the simulation slower than when using \code{ifelse} to make a stepfunction, but it is an option
#' make it a function so we can plot it
rr <- MakeRate(r, env_f=env_r)
#' let us check r is high enough to see a pattern
plot(1:10, rr(1:10), type='l', main="Sampling rate", xlab="My", ylab="r")
#' the resolution of the fossil dataset:
stages=seq(from=10, to = 0, # from simulation's tmax to present
           by = -.1) # note that we will provide a very high resolution (= 0.1 My) just to test the function:
#'
dt<-SampleClade(1:length(sim$TE), sim, r, tmax=10, env_rr=env_r, stages=stages)
ids<-unique(dt$Species)
mids<-(dt$MaxT-dt$MinT)+dt$MinT
#'
for(i in 1:length(ids)){
  sp<-unique(as.numeric(gsub("spp_", "", ids[i])))
  
  hist(mids[dt$Species==ids[i]], main=paste0("spp = ", sp, "; duration ~ ", round(sim$TS[sp]-sim$TE[sp], digits = 2), "my"), 
       xlab="My", breaks=seq(round(sim$TS[i]), round(sim$TE[i]), -1), xlim=c(sim$TS[i], sim$TE[i]))
  mid=par[sp]
  t <- seq(sim$TE[i], sim$TS[i], 0.1)
  lines(t, rev(rr(t)))
}
#' 
#' we will now do some tests with age-dependent rates. For more details, check \code{SampleADPP}.
#' 
#' let us start with a hat-shaped increase through the duration of a species
sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
while((sim$TS[1]-sim$TE[1])<10){ # in case first simulation has short-lived lineage which will obscure the pattern
  sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
}
#'
# preservation function in respect to age
# here we will use the PERT function. It is described in: 
# Silvestro, D., Schnitzler, J., Liow, L. H., Antonelli, A., & Salamin, N. (2014). Bayesian estimation of speciation and extinction from incomplete fossil occurrence data. Systematic biology, 63(3), 349-367.
#'
dPERT<-function(t,s,e,sp,a=3,b=3, log=F){ #this function have all criteria cited in the help page
  
  if(e>=s){
    message("There is no PERT with e>=s") 
    return(rep(NaN, times=length(t)))
  }
  id1<-which(t<=e | t>=s)
  id2<-which(!(t<=e | t>=s))
  t<-t[id2]
  
  res<-vector()
  if(log){
    res[id1]<--Inf
  }else{
    res[id1]<-0
  }
  
  if(log){
    res[id2]<-log(((s-t)^2)*((-e+t)^2)/((s-e)^5*beta(a,b)))
  } else{
    res[id2]<-((s-t)^2)*((-e+t)^2)/((s-e)^5*beta(a,b))
  }
  return(res)
}
#'
dPERTmax<-function(s,e,sp){
  return(((s-e)/2)+e)
}
#'
#' the resolution of the fossil dataset:
stages<-seq(from=10, to = 0, # from simulation's tmax to present
            by = -.1) # note that we will provide a very high resolution (= 0.1 My) just to test the function:
#'
dt<-SampleClade(S = 1:length(sim$TE), sim, rr=10, tmax = 100, stages = stages, dFUN = dPERT, dFUNmax = dPERTmax)
#'
ids<-unique(dt$Species)
mids<-(dt$MaxT-dt$MinT)+dt$MinT
#'
for(i in 1:length(ids)){
  
  sp<-unique(as.numeric(gsub("spp_", "", ids[i])))
  
  hist(mids[dt$Species==ids[i]],  probability = T,
       main=paste0("spp = ", sp, "; duration ~ ", round(sim$TS[sp]-sim$TE[sp], digits = 2), "my"))
  mid<-par[sp]#+par1[sp]
  curve(dPERT(x, s = sim$TS[sp], e = sim$TE[sp], sp=sp),from = sim$TE[sp], to = sim$TS[sp], add=T, col="red", n = 100) # expected by model
}
# note that the sampling at fossil stages distorts the quantiles of the distribution, even with the high preservation rate and the high resolution of fossil stages.
#'
#' Now, a hat-shaped increase through the duration of a species dependent on two parameters
#'
sim<-BDSim(N0 = 1, pp = 0.1, qq = 0.1, tmax = 10)
while(length(sim$TE)<20){ #generating a reasonable amount of lineages with some (expected) reasonable duration
  sim<-BDSim(N0 = 1, pp = .1, qq = 0.1, tmax = 10)
}
#'
# preservation function in respect to age
# here we will use the triangular distribution as a model of a "hat-shaped" function. In this model, preservation responds to the age of a linegae (in absolute time) and a "mode" of the triangle. Thi mode, in this example, is the result of the interaction betweentwo parameters: par and par1
dTRImod2<-function(t,s,e,sp){
  
  if(e>=s){
    message("There is no TRI with e>=s")
    return(rep(NaN, times=length(t)))
  }
  
  md<-par[sp]+par1[sp]
  
  if(md<e | md>s){
    message("There is no TRI with md outside [s, e] interval") 
    return(rep(NaN, times=length(t)))
  }
  
  id1<-which(t>=e & t<md)
  id2<-which(t==md)
  id3<-which(t>md & t<=s)
  id4<-which( !(1:length(t) %in% c(id1,id2,id3)))
  
  res<-vector()
  
  res[id1]<-(2*(t[id1]-e))/((s-e)*(md-e)) 
  res[id2]<-2/(s-e) 
  res[id3]<-(2*(s-t[id3]))/((s-e)*(s-md)) 
  res[id4]<-0 
  
  return(res) #for more details in this function, see the deatiled examples in the help page of the SampleADPP() function
}
#'
dTRImaxmod2<-function(s,e,sp){ ##note that now we don't have the "md", just like in the dTRImod1 (example 4) function
  return(2/(s-e))
}
par<-runif(n = length(sim$TE), min = sim$TE, max = sim$TS) # a random point inside each lineage's duration
par1<-(((sim$TS-sim$TE)/2)+sim$TE)-par #a distance between "par" and the lineage's duration middle 
#'
#' the resolution of the fossil dataset:
stages<-seq(from=10, to = 0, #from simulation's tmax to present
            by = -.1) #note that we will provide a very high resolution (= 0.1 My) just to test the function:
dt<-SampleClade(S = 1:length(sim$TE), sim, rr=10, tmax = 100, stages = stages, dFUN = dTRImod2, dFUNmax = dTRImaxmod2)
#'
ids<-unique(dt$Species)
mids<-(dt$MaxT-dt$MinT)+dt$MinT
#'
for(i in 1:length(ids)){
  
  sp<-unique(as.numeric(gsub("spp_", "", ids[i])))
  
  hist(mids[dt$Species==ids[i]],  probability = T,
       main=paste0("spp = ", sp, "; duration ~ ", round(sim$TS[sp]-sim$TE[sp], digits = 2), "my"))
  mid<-par[sp]#+par1[sp]
  curve(dTRImod2(x, e=sim$TE[sp], s=sim$TS[sp], sp=sp),from = sim$TE[sp], to = sim$TS[sp], add=T, col="red", n = 100)  #expected by model
}
#' note that the sampling at fossil stages distorts the quantiles of the distribution, even with the high preservation rate and the high resolution of fossil stages.
#'

SampleClade<-function(S, sim, rr,tmax,env_rr=NULL,rshifts=NULL,returnTrue=FALSE,stages=NULL, dFUN=NULL, dFUNmax=NULL,...){
  # get the speciation and extinction times vectors
  TE <- sim$TE[S]
  TS <- sim$TS[S]
  
  # check if it is age-dependent
  if(is.null(dFUN)){
    rr <- MakeRate(rr, tmax, env_rr, rshifts)
  } else{
    if(!is.numeric(rr) | length(rr)>1){
      stop("ADPP cannot be used with time-varing preservation rates")
    }
  }
  
  # adjusting stages
  stages<-sort(stages, decreasing = T)
  
  # sample using Poisson process:
  if(is.null(dFUN)){ # independent of age (i.e. occurrences uniformly distributed through the lineage's age)
    point_estimates<-lapply(S,Sample,TE=TE,TS=TS,rr=rr,tmax=tmax)  
  } else{ #dependent of age (i.e. occurrences distributed through the lineage's age accourding to the function provided by the user)
    point_estimates<-SampleADPP(S, TS=TS, TE=TE, rr=rr, dFUN = dFUN, dFUNmax = dFUNmax, ...)
  }
  
  #wrapping data:
  if(!returnTrue){ # output as fossil occurrence binned within stages/bins
    
    res<-data.frame(matrix(nrow=0, ncol=4))
    colnames(res)<-c("Species", "Extant", "MaxT", "MinT")
    for(i in 1:length(point_estimates)){
      binned_occs<-binner(point_estimates[[i]], bins=stages)
      for(k in 1:(length(stages)-1)){
        if(binned_occs[k]>0){
          # make a row of the data frame
          aux<-data.frame(Species=i,Extant=NA,MaxT=rep(stages[k], times=binned_occs[k]),MinT=stages[k+1])
          res<-rbind(res, aux)
        }
      }
    }
    # make the Extant column
    res$Extant<-FALSE
    res$Extant[res$Species %in% which(sim$EXTANT)]<-TRUE
    
    # and the species column
    res$Species<-paste0("spp_", res$Species) 
  } else{ # output as the "true" times of preservation of each lineage
    # if returnTrue=TRUE, get a a data frame with the real sampling times only
    res<-data.frame(matrix(nrow=length(unlist(point_estimates)), ncol=3))
    colnames(res)<-c("Species", "Extant", "SampT")
    res$Species<-rep(S, times=lapply(point_estimates, length))
    res$Extant<-FALSE
    res$Extant[res$Species %in% which(sim$EXTANT)]<-TRUE
    res$Species<-paste0("spp_", res$Species)
    res$SampT<-unlist(point_estimates)
  }
  
  return(res)
}
