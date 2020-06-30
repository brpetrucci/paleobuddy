#' Returns a list of occurrences for a species \code{S} during their lifespan,
#' following a Poisson Point process
#'
#' \code{Sample} takes a species number, a vector of speciation and extinction
#' times, a sampling rate with possible shape and a maximum time for simulation.
#'
#' @parameter \code{S} the species number to be sampled. Since \code{Sample} 
#' will be called by a wrapper using \code{lapply}, it is through \code{S} 
#' that we apply this function.
#' 
#' @parameter \code{TE} a vector of extinction times, usually an output of
#' \code{BDSim}.
#' 
#' @parameter \code{TS} a vector of speciation times.
#' 
#' @parameter \code{rr} a sampling rate function, usually created by 
#' \code{MakeRate}. Can be an exponential rate or a weibull scale, if 
#' \code{rshape} is not 0.
#'
#' @parameter \code{tmax} the maximum simulation time, used by 
#' \code{rexp_var}.
#'
#' @return a list of occurrences for that species, expected to be around
#' \code{\int_TS[S]^TE[S] rr(t)dt} occurrences.
#'
#' @author written by Bruno do Rosario Petrucci.
#' 
#' @examples
#' 
#' to start we need a simulation, let us generate a simple one
sim <- BDSim(1, 0.1, 0.09, 50)
#'
#' let's first sample with constant rate
r <- 1 # this means we expect t occurrences, where t is the species duration
nrep <- 1000 # sample 1000 times
sample_1 <- lapply(1:nrep, function(x) {if (x %% 100 == 0) print(x);
  return(Sample(1, sim$TE, sim$TS, r, tmax=50))})
#' after sampling the first species, we can take a look at the 
#' mean number of occurrences; it should be (TE[1]-TS[1])*1
nOcc_1 <- c()
for (i in 1:nrep) {
  nOcc_1 <- c(nOcc_1, length(sample_1[[i]]))
}
boxplot(nOcc_1)
abline(h=(r * (sim$TE[1] - sim$TS[1])))
#' as you can see, the median of the boxplot feels right where we want
#' 
#' let us try a more complicated r,
r <- function(t) 1+0.1*t
sample_2 <- lapply(1:nrep, function(x) {if (x %% 100 == 0) print(x);
  return(Sample(1, sim$TE, sim$TS, r, tmax=50))})
nOcc_2 <- c()
for (i in 1:nrep) {
  nOcc_2 <- c(nOcc_2, length(sample_2[[i]]))
}
boxplot(nOcc_2)
abline(h=integrate(r, tmax-sim$TS[1], tmax-sim$TE[1])$value)
#' r could be any function of time, but to avoid redundancy regarding tests of
#' BDSim() and rexp_var() we will stop here

Sample<-function(S,TE,TS,rr,tmax){
  TE <- tmax - TE
  TS <- tmax - TS

  # start when the species was born, end when it died
  Now<-ifelse(TS[S]<0,0,TS[S])
  End<-ifelse(TE[S]>tmax,tmax,TE[S])
  
  # make rr a function if it isn't
  r <- rr
  rr <- ifelse(is.numeric(r), Vectorize(function(t) r), r)
  
  # initialize vector
  sampled<-c()
  
  # while we have not reached the time of extinction
  while (Now<End){
    # take th waiting time for sampling, using rexp_var()
    WaitTimeR<-ifelse(rr(Now)>0,rexp_var(1,rr,Now,tmax), Inf)
    
    # advance in time
    Now<-Now+WaitTimeR
    
    # if sampling comes after extinction, we don't include this occurrence
    if (Now>End) break
    
    # add to the vector
    sampled<-c(sampled,Now)
  }
  
  # finally, we invert time so that it goes from tmax to 0
  sampled <- tmax - sampled
  
  return(sampled)
}
  