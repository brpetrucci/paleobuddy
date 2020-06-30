#' Returns a list of occurrence lists for all species in a given simulation,
#' either the real or a range of occurrence times.
#'
#' \code{SampleClade} takes times of speciation and extinction, information to
#' create a sampling rate with \code{MakeRate}, a vector of geologic time 
#' intervals and whether one wants the true return times or a range based on 
#' \code{IntVec}.
#' 
#' @param \code{TE} a vector of extinction times, usually an output of 
#' \code{BDSim}.
#' 
#' @param \code{TS} a vector of speciation times.
#'
#' @param \code{IntVec} a vector of time intervals corresponding to 
#' geological time ranges. If \code{returnTrue} is false, \code{SampleClade} 
#' returns the member of \code{IntVec} right after the true occurrence time. 
#' In this way, we simulate the granularity in real world fossil records. If 
#' \code{returnTrue} is true, this is ignored.
#'
#' @param \code{rr} a sampling rate function. May be a constant, a 
#' time-dependent function, a function dependent on time and environment, or a 
#' vector of rates corresponding to the times in \code{rshifts}. If 
#' \code{env_rr} and \code{rshifts} are NULL, it will be either treated as an 
#' exponential rate, if \code{rshape} is 0, or a weibull scale, if 
#' \code{rshape} > 0.
#' 
#' @param \code{tmax} the maximum simulation time, used by 
#' \code{rexp_var}.
#'
#' @param \code{rshape} a shape for weibull distributions. Sampling is not 
#' usually modelled as a Weibull, but we provide this feature while a more 
#' common age-dependent sampling model is not implemented.
#'
#' @param \code{env_rr} a matrix containing time points and values of an
#' enviromental variable, like temperature, for each time point. This will be
#' used to create a sampling rate, so \code{rr} must be a function of time and
#' said variable if \code{env_rr} is not NULL.
#'
#' @param \code{rshifts} vector of rate shifts. First element must be 
#' starting time for simulation (0). Vector must have the same length as 
#' \code{rr} E.g. \code{rr = c(0.15, 0.1, 0.2)}, \code{pshifts = c(0, 20, 30)} 
#' means r will be 0.15 from 0 to 20, 0.1 from 20 to 30, and 0.2 from 30 to 
#' \code{tmax}. Note that using this method for step-function rates is 
#' currently slower than using \code{ifelse}.
#'
#' @param \code{returnTrue} if set to true, a list of true ocurrence times
#' for each species will be returned. If set to true, we call \code{Grid} to 
#' take the list of samples and return a list of upper bounds for occurrences.
#'
#' @return a list of either the true occurrence times for a species, or the
#' upper bound for those times taking into account \code{IntVec} (one can then 
#' use \code{IntVec} to find the exact range). In both cases, a list of lists 
#' (i.e. \code{return[[i]]} is a list of occurrences/bounds on occurrences for
#' species i).
#'
#' @author written by Bruno do Rosario Petrucci.
#'
#' @examples
#'
#' to start we need a simulation, let us generate a simple one
sim <- BDSim(1, 0.1, 0.09, 50)
length(sim$TE) # we need this to have more than 1 species for the tests to work
#' start with a constant rate
r <- 1
#'
nrep<-1000
samples<-lapply(1:nrep, function(x) SampleClade(sim$TE, sim$TS, r, returnTrue=TRUE))
#'
Species<-1:length(sim$TE)
countMean <- c()
countTotal <- matrix(0,length(Species),nrep)
countExpected <- c()
#'
for (s in Species) {
  countNum <- c()
  for (n in 1:nrep) {
    sampled <- samples[[n]][[s]]
    countNum <- c(countNum, length(sampled))
    countTotal[s, n] <- length(sampled)
  }
  countMean <- c(countMean, mean(countNum))
  countExpected <- c(countExpected, integrate(Vectorize(r), ifelse(sim$TS[s] < 0, 0, sim$TS[s]), 
                                              ifelse(sim$TE[s] > tmax, tmax, sim$TE[s]))$value)
}
#' first we create a boxplot to check that the median is where we would expect
boxplot.matrix(countTotal, outline=FALSE, use.cols=FALSE, main="Occurrence count per species boxplot",
               ylab="Occurrences")
lines(Species, countExpected, type='p', pch=19, col='GREEN')
#' then a plot of the mean occurrence number to compare to the expected curve
plot(Species, countMean, type='o', main="Mean occurrence number per species", ylab="Occurrences")
lines(Species, countExpected, type='o', col='GREEN')
#'
#' now let us do the same for a more complicated r
r <- function(t) {
  return(1+0.1*t)
}
samples<-lapply(1:nrep, function(x) SampleClade(sim$TE, sim$TS, r, returnTrue=TRUE))
#'
countMean <- c()
countTotal <- matrix(0,length(Species),nrep)
countExpected <- c()
#'
for (s in Species) {
  countNum <- c()
  for (n in 1:nrep) {
    sampled <- samples[[n]][[s]]
    countNum <- c(countNum, length(sampled))
    countTotal[s, n] <- length(sampled)
  }
  countMean <- c(countMean, mean(countNum))
  countExpected <- c(countExpected, integrate(Vectorize(r), ifelse(sim$TS[s] < 0, 0, sim$TS[s]), 
                                              ifelse(sim$TE[s] > tmax, tmax, sim$TE[s]))$value)
}
boxplot.matrix(countTotal, outline=FALSE, use.cols=FALSE, main="Occurrence count per species boxplot",
               ylab="Occurrences")
lines(Species, countExpected, type='p', pch=19, col='GREEN')
plot(Species, countMean, type='o', main="Mean occurrence number per species", ylab="Occurrences")
lines(Species, countExpected, type='o', col='GREEN')
#' we see that it fits pretty well
#' r can be any time-dependent function, and we could make it a weibull scale
#' that is not common in the literature, however
#' 
#' seen as the \code{returnTrue=TRUE} tests are working, we refer to the successful
#' tests of \code{Grid} to note that the \code{returnTrue=FALSE} tests should also
#' be successful

SampleClade<-function(TE,TS,rr,tmax,rshape=0,env_rr=NULL,rshifts=NULL,returnTrue=FALSE,IntVec=NULL){
  
  # get species number
  NSpecies<-length(TE)
  
  # use MakeRate to create the sampling rate
  r<-MakeRate(rr, env_rr, rshifts)
  
  # make it a function if it is a
  rr<-ifelse(is.numeric(r), function(t) return(r+0*t), r)
  
  # sample using the function above
  sampled<-lapply(1:NSpecies,Sample,TE=TE,TS=TS,rr=rr,tmax=tmax,rshape=rshape)
  
  if (returnTrue==TRUE){
    return(sampled)
  }
  else{
    GridedSample<-lapply(1:NSpecies,function(x){as.numeric(Grid(sampled[[x]], IntVec))})
    return(GridedSample)
  }
}
