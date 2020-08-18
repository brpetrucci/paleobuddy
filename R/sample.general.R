#' Time-dependent and age-dependent rate species sampling
#' 
#' Generates a vector of occurrence times for a species in a simulation using a
#' a Poisson process. Allows for the Poisson rate to be (1) a constant, (2) a 
#' function of time, (3) age-dependent or (4) a mix of age-dependent and
#' time-dependent preservation. For more flexibility options, see \code{make.rate}
#' and \code{sample.clade}. Note that while the Poisson process occurs in forward 
#' time, we return (both in birth-death functions and here) results in backwards 
#' time, so that time is inverted using \code{tMax} both at the beginning and end
#' of \code{sample.species}, which is used in this function.
#'
#' @param sim A \code{sim} object, usually the output of \code{bd.sim}.
#' 
#' @param rr Sampling rate (per species per million years) over time. It can be
#' a \code{numeric} describing a constant rate or a \code{function(t)} describing
#' the variation in sampling over time. For more flexibility on sampling, see
#' \code{make.rate} to create more complex rates. Note that \code{rr} should
#' always be greater than or equal to zero.
#' 
#' @param tMax The maximum simulation time, used by \code{rexp.var}. A sampling
#' time greater than \code{tMax} would mean the occurrence is sampled after the
#' present, so for consistency we require this argument. This is also required
#' to ensure time follows the correct direction both in the Poisson process and
#' in the return.
#'
#' @param S A vector of species numbers to be sampled. Could be only a subset of 
#' the species if the user wishes. The default is all species in \code{sim}. 
#' Species not included in \code{S} will not be sampled by the function.
#'
#' @param adFun A density function representing the age-dependent preservation
#' model. It must be a density function, and consequently integrate to 1 (though
#' this condition is not verified by the function). If not provided, a uniform 
#' distribution will be used. The function must also have the following 
#' properties:
#'
#' \itemize{
#'
#' \item Returns a vector of preservation densities for each time in a given
#' vector \code{t} in geological time. 
#'
#' \item Be parametrized in absolute geological time (i.e. should be relative to
#' absolute geological time, in Mya, \emph{not} the lineage's age).
#'
#' \item Should be limited between \code{s} (i.e. the lineage's
#' speciation/origination in geological time) and \code{e} (i.e. the lineage's 
#' extinction in geological time), with \code{s} > \code{e}.
#'
#' \item Include the arguments \code{t}, \code{s}, \code{e} and \code{sp}. 
#' The argument sp is used to pass species-specific parameters (see examples),
#' allowing for \code{dFun} to be species-inhomogeneous.
#' }
#'
#' @param ... Additional parameters related to \code{adFun}
#'
#' @return A list of occurrences for that species, with their distribution in 
#' species relative time given by the \code{adFun} function provided by the user.
#'
#' @author Matheus Januario
#'
#' @examples
#' 
#' # sampling rate
#' rr <- function(t) {
#'   return(0.5 + 1*t) 
#' }
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # we can start with a hat-shaped increase through the duration of a species
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = tMax)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = tMax)
#' }
#' 
#' # fixing the extinction times
#' sim$TE[sim$EXTANT] <- 0
#' 
#' # here we will use the PERT function. It is described in:
#' # Silvestro et al 2014
#' 
#' # preservation function
#' dPERT <- function(t, s, e, sp, a = 3, b = 3, log = FALSE) {
#'   
#'   # check if it is a valid PERT
#'   if (e >= s) {
#'     message("There is no PERT with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # find the valid and invalid times
#'   id1 <- which(t <= e | t >= s)
#'   id2 <- which(!(t <= e | t >= s))
#'   t <- t[id2]
#'   
#'   # initialize result vector
#'   res <- vector()
#'   
#'   # if user wants a log function
#'   if (log) {
#'     # invalid times get -Inf
#'     res[id1] <- -Inf
#'     
#'     # valid times calculated with log
#'     res[id2] <- log(((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b)))
#'   }
#'   # otherwise
#'   else{
#'     res[id1] <- 0
#'     
#'     res[id2] <- ((s - t) ^ 2)*((-e + t) ^ 2)/((s - e) ^ 5*beta(a,b))
#'   }
#'   
#'   return(res)
#' }
#' 
#' 
#' # find occurrences
#' occs <- sample.general(sim, rr, tMax = tMax, S = 1:length(sim$TE), adFun = dPERT)
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE, main = "Black = no time dependence
#'  Red = with time dependence", ylab="Density of fossil occurrences",
#'      xlab = "Time (Mya)", xlim = c(max(occs[[1]]), min(occs[[1]])))
#' 
#' # expected curve - probably will not fit great because of low sample size
#' # getting expected values form age-dependence alone
#' tt <- seq(from = sim$TE[1], to = sim$TS[1], by = 0.01)
#' lines(x = tt, y = dPERT(tt, s = sim$TS[1], e = sim$TE[1], sp = 1))
#' 
#' #getting expected values:
#' Pres_time_adpp=function(t, s, e, sp) {
#'   rrmod=function(t){#correction of scale
#'     return(rr(tMax-t))
#'   } 
#'   return(rrmod(t)*dPERT(t = t, s = s, e = e, sp = sp))
#' }
#' 
#' Pres_time_adppNorm=function(t, s, e, sp) {
#'   return(Pres_time_adpp(t = t, s = s, e = e, sp = sp)/integrate(
#'     Pres_time_adpp, lower = e, upper = s, e=e, s=s, sp=sp)$value)
#' }
#' 
#' lines(x = tt, y = Pres_time_adppNorm(tt, s = sim$TS[1], e = sim$TE[1], sp = 1),
#'       col="red")
#' 
#' # now we can test the simpler scenario of uniform age-independent but
#' # time-dependent sampling (this is the same as not declaring an input to
#' # the "adFun" argument, but lets do this for the sake of understanding how
#' # the function works)
#' 
#' # sampling rate
#' rr <- function(t) {
#'   return(0.5 + 1*t) 
#' }
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # fixing the extinction times
#' sim$TE[sim$EXTANT] <- 0
#' 
#' # preservation function in respect to age
#' # occurrences are uniformly distributed
#' custom.uniform <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid uniform
#'   if (e >= s) {
#'     message("There is no uniform function with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   res <- dunif(x = t, min = e, max = s)
#'   
#'   return(res)
#' }
#' 
#' # find occurrences
#' occs <- sample.general(sim, rr, tMax = tMax, S = 1:length(sim$TE),
#'                        adFun = custom.uniform)
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE, main = "Black = no time dependence
#'  Red = with time dependence", ylab="Density of fossil occurrences",
#'      xlab = "Time (Mya)", xlim = c(max(occs[[1]]), min(occs[[1]])))
#' 
#' # expected curve - probably will not fit great because of low sample size
#' # getting expected values form age-dependence alone
#' tt <- seq(from = sim$TE[1], to = sim$TS[1], by = 0.01)
#' lines(x = tt, y = custom.uniform(tt, s = sim$TS[1], e = sim$TE[1], sp = 1))
#' 
#' #getting expected values:
#' Pres_time_adpp=function(t, s, e, sp) {
#'   rrmod=function(t){#correction of scale
#'     return(rr(tMax-t))
#'   } 
#'   return(rrmod(t)*custom.uniform(t = t, s = s, e = e, sp = sp))
#' }
#' 
#' Pres_time_adppNorm=function(t, s, e, sp) {
#'   return(Pres_time_adpp(t = t, s = s, e = e, sp = sp)/integrate(
#'     Pres_time_adpp, lower = e, upper = s, e=e, s=s, sp=sp)$value)
#' }
#' 
#' lines(x = tt, y = Pres_time_adppNorm(tt, s = sim$TS[1], e = sim$TE[1], sp = 1),
#'       col="red")
#' 
#' # now, lets use a hat-shaped increase through the duration of a species with more
#' # parameters than TS and TE
#' 
#' # sampling rate
#' rr <- function(t) {
#'   return(0.5 + 1*t) 
#' }
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.2, tMax = 10, nFinal = c(0,1))
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10, nFinal = c(0,1))
#' }
#' 
#' # fixing the extinction times
#' sim$TE[sim$EXTANT] <- 0
#' 
#' # here we will use the triangular distribution. We have some empirical evidence
#' # that taxa occurrences might present a triangular shape (in a very broad sense)
#' # see Zliobaite et al 2017
#' 
#' # preservation function
#' dTRI <- function(t, s, e, sp, md) {
#'   
#'   # please note ths function is inverted. The correspondence would be:
#'   # s = b maximum
#'   # e = a minimum
#'   # md = c distribution's mode
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # another condition we must check
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # needed to vectorize the function:
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t > md & t <= s)
#'   id4 <- which( !(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   # actually vetorizing function
#'   res <- vector()
#'   
#'   # (t >= e & t < md)
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e)*(md - e))
#'   
#'   #(t == md)
#'   res[id2] <- 2 / (s - e)
#'   
#'   #(md < t & t <= s)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e)*(s - md))
#'   
#'   #outside function's limits
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' 
#' # note we are providing the mode for the triangular sampling as an argument
#' occs <- sample.general(sim = sim, rr = rr, adFun = dTRI, 
#'                        tMax = tMax, S = 1, md = 8)
#' 
#' # please note in the original parametrization, the "md" parameter (mode) is in
#' # "absolute time", i.e. a specific number in the absolute timescale.
#' 
#' # This is not an "age-dependent" model in a strict sense. We show it here as
#' # we will construct models related to age in the following examples.
#' # This is also the reason we only plot the first lineage (note the "S"
#' # parameter) in this specific example (many lineages might not be alive in
#' # the "md" moment in time).
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE, main = "Black = no time dependence
#'  Red = with time dependence", ylab="Density of fossil occurrences",
#'      xlab = "Time (Mya)", xlim = c(max(occs[[1]]), min(occs[[1]])))
#' 
#' # expected curve - probably will not fit great because of low sample size
#' # getting expected values form age-dependence alone
#' tt <- seq(from = sim$TE[1], to = sim$TS[1], by = 0.01)
#' lines(x = tt, y = dTRI(tt, s = sim$TS[1], e = sim$TE[1], sp = 1, md = 8))
#' 
#' # getting expected values from age + time dependences:
#' #getting expected values:
#' Pres_time_adpp=function(t, s, e, sp, ...) {
#'   rrmod=function(t){#correction of scale
#'     return(rr(tMax-t))
#'   } 
#'   return(rrmod(t)*dTRI(t = t, s = s, e = e, sp = sp, ...))
#' }
#' 
#' Pres_time_adppNorm=function(t, s, e, sp, ...) {
#'   return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...)/integrate(
#'     Pres_time_adpp, lower = e, upper = s, e=e, s=s, sp=sp, ...)$value)
#' }
#' 
#' lines(x = tt, y = Pres_time_adppNorm(tt, s = sim$TS[1], e = sim$TE[1], sp = 1, md=8),
#'       col="red")
#' 
#' 
#' # we can also have a hat-shaped increase through the duration of a species
#' # with more parameters than TS and TE, but with the parameters relate to
#' # the relative age of each lineage
#' 
#' # sampling rate
#' rr <- function(t) {
#'   return(0.5 + 1*t) 
#' }
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # fixing the extinction times
#' sim$TE[sim$EXTANT] <- 0
#' 
#' # preservation function, with the "mde" of the triangle being exactly at the
#' # last quarter of the duration of EACH lineage
#' dTRImod1 <- function(t, s, e, sp) {
#'   # note that now we don't have the "md" parameter here,
#'   # but it is calculated inside the function
#'   
#'   # check if it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # calculate md
#'   md <- ((s - e) / 4) + e
#'   # md is at the last quarter of the duration of the lineage
#'   
#'   # please note that the same logic can be used to sample parameters
#'   # internally in the function, running for instance:
#'   # md<-runif (n = 1, min = e, max = s)
#'   
#'   # check it is a valid md
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times=length(t)))
#'   }
#'   
#'   # needed to vectorize function
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t>md & t <= s)
#'   id4 <- which( !(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   # vectorize the function
#'   res<-vector()
#'   
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e)*(md - e))
#'   res[id2] <- 2 / (s - e)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e)*(s - md))
#'   res[id4] <- 0
#'   
#'   return(res)
#' }
#' 
#' 
#' # note we are providing the mode for the triangular sampling as an argument
#' occs <- sample.general(sim = sim, rr = rr, adFun = dTRImod1, tMax = tMax,
#'                        S = 1:length(sim$TE))
#' 
#' # please note in this case, function dTRImod1 fixes "md" at the last quarter
#' #of the duration of the lineage
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE, main="Black = no time dependence
#'  Red = with time dependence", ylab="Density of fossil occurrences",
#'      xlab = "Time (Mya)", xlim = c(max(occs[[1]]), min(occs[[1]])))
#' 
#' # expected curve - probably will not fit great because of low sample size
#' # getting expected values form age-dependence alone
#' tt <- seq(from = sim$TE[1], to = sim$TS[1], by = 0.01)
#' lines(x = tt, y = dTRImod1(tt, s = sim$TS[1], e = sim$TE[1], sp = 1))
#' 
#' # getting expected values from age + time dependences:
#' #getting expected values:
#' Pres_time_adpp=function(t, s, e, sp, ...) {
#'   rrmod=function(t){#correction of scale
#'     return(rr(tMax-t))
#'   } 
#'   return(rrmod(t)*dTRImod1(t = t, s = s, e = e, sp = sp, ...))
#' }
#' 
#' Pres_time_adppNorm=function(t, s, e, sp, ...) {
#'   return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...)/integrate(
#'     Pres_time_adpp, lower = e, upper = s, e=e, s=s, 
#'     sp=sp, ...)$value)
#' }
#' 
#' lines(x = tt, y = Pres_time_adppNorm(tt, s = sim$TS[1], e = sim$TE[1], sp = 1),
#'       col="red")
#' 
#' # finally, a hat-shaped increase through the duration of a species with more
#' # parameters than TS and TE, but the parameters relate to each specific lineage.
#' # This is useful when the user wants to use variable parameters for each species
#' # but wants to keep track of those parameters after the sampling is over
#' 
#' # sampling rate
#' rr <- function(t) {
#'   return(0.5 + 1*t) 
#' }
#' 
#' # maximum simulation time
#' tMax <- 10
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # fixing the extinction times
#' sim$TE[sim$EXTANT] <- 0
#' 
#' # get the par and par1 vectors
#' 
#' # a random quantity
#' par <- runif(n = length(sim$TE), min = sim$TE, max = sim$TS) 
#' # its complement to the middle of the lineage's age. 
#' # Note that the interaction between these two parameters creates a 
#' # deterministic parameter, but inside the function one of them ("par") 
#' # is a random parameter
#' par1 <- (((sim$TS - sim$TE)/2) + sim$TE) - par 
#' 
#' # preservation function in respect to age, with the "mode" of the triangle
#' # being exactly at the last quarter of the duration of EACH lineage.
#' dTRImod2 <- function(t, s, e, sp) {
#'   
#'   # make sure it is a valid TRI
#'   if (e >= s) {
#'     message("There is no TRI with e >= s")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   # here is the difference from the function in example 3 and 4
#'   md <- par[sp] + par1[sp]
#'   
#'   # check that md is valid
#'   if (md < e | md > s) {
#'     message("There is no TRI with md outside [s, e] interval")
#'     return(rep(NaN, times = length(t)))
#'   }
#'   
#'   id1 <- which(t >= e & t < md)
#'   id2 <- which(t == md)
#'   id3 <- which(t > md & t <= s)
#'   id4 <- which(!(1:length(t) %in% c(id1,id2,id3)))
#'   
#'   res <- vector()
#'   
#'   res[id1] <- (2*(t[id1] - e)) / ((s - e)*(md - e))
#'   res[id2] <- 2 / (s - e)
#'   res[id3] <- (2*(s - t[id3])) / ((s - e)*(s - md))
#'   res[id4] <- 0
#'   
#'   return(res)
#'   #for more details in this function, see example 3 and 4
#' }
#' 
#' # note we are providing the mode for the triangular sampling as an argument
#' occs <- sample.general(sim = sim, rr = rr, adFun = dTRImod2, tMax = tMax,
#'                        S = 1:length(sim$TE))
#' 
#' # please note in this case, function dTRImod1 fixes "md" at the last quarter
#' #of the duration of the lineage
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE, main="Black = no time dependence
#'  Red = with time dependence", ylab="Density of fossil occurrences",
#'      xlab = "Time (Mya)", xlim = c(max(occs[[1]]), min(occs[[1]])))
#' 
#' # expected curve - probably will not fit great because of low sample size
#' # getting expected values form age-dependence alone
#' tt <- seq(from = sim$TE[1], to = sim$TS[1], by = 0.01)
#' lines(x = tt, y = dTRImod2(tt, s = sim$TS[1], e = sim$TE[1], sp = 1))
#' 
#' # getting expected values from age + time dependences:
#' #getting expected values:
#' Pres_time_adpp=function(t, s, e, sp, ...) {
#'   rrmod=function(t){#correction of scale
#'     return(rr(tMax-t))
#'   } 
#'   return(rrmod(t)*dTRImod2(t = t, s = s, e = e, sp = sp, ...))
#' }
#' 
#' Pres_time_adppNorm=function(t, s, e, sp, ...) {
#'   return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...)/integrate(
#'     Pres_time_adpp, lower = e, upper = s, e=e, s=s, 
#'     sp=sp, ...)$value)
#' }
#' 
#' 
#' lines(x = tt, y = Pres_time_adppNorm(tt, s = sim$TS[1], e = sim$TE[1], sp = 1),
#'       col="red")
#' 
#' @name sample.general
#' @rdname sample.general
#' @export

sample.general <- function(sim, rr, tMax, S = NULL, adFun = NULL, ...){
  # checking input
  # if rr is constant, make it a function
  if (is.numeric(rr)) {
    r <- rr
    rr <- Vectorize(function(t) r)
  }
  
  # if S is not given
  if (is.null(S)) {
    S = 1:length(sim$TE)
  }
  
  # get extinction and speciation times
  TE <- sim$TE
  TS <- sim$TS
  
  # fix extant extinction times
  TE[sim$EXTANT] <- 0
  
  # check if adFun is a function with the right parameters
  if(!(is.null(adFun))){
    if (sum(c("t", "s", "e", "sp") %in% names(formals(adFun))) < 3) {
      stop("adFun must have \"t\", \"s\", \"e\", and  \"sp\" parameters")
    }
  }
  
  # if adFun is null, it will be uniform
  if(is.null(adFun)){
    
    message("Preservation will be Age-independent \n")
    
    adFun <- function(t, s, e, sp) {
      
      # make sure it is a valid uniform
      if (e >= s) {
        message("There is no uniform function with e >= s")
        return(rep(NaN, times = length(t)))
      }
      
      res <- dunif(x = t, min = e, max = s)
      
      return(res)
    }
  }
  
  # creating function with the intersection of the age-dependent and the 
  # time-varying poisson process
  
  #"prior step" function (interaction between time-dependency and age-dependency):
  Pres_time_adpp=function(t, s, e, sp, ...) {
    # this correction is needed because "t" in rr means
    # "time since clade origin" and NOT absolute geologic time in Mya
    rrmod <- function(t){ 
      return(rr(tMax - t))
    }
    return(rrmod(t)*adFun(t = t, s = s, e = e, sp = sp, ...))
  }
  
  # now normalizing to get a real density function
  Pres_time_adppNorm <- function(t, s, e, sp,...) {
    return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...) / 
             integrate(Pres_time_adpp, lower = e, upper = s, s = s, 
                       e = e, sp = sp, ...)$value)
  }
  
  # getting the number of occurrences following a time-varying poisson process:
  
  Noccs <- vector()
  Noccs[S] <- unlist(x = lapply(
    lapply(S, sample.species, sim = sim, rr = rr, tMax = tMax), FUN = length))
  
  occs <- list()
  for (sp in S) {
    # getting the maximum of function for that specific species
    adFunMaxApprox <- optimize(Pres_time_adppNorm, interval = c(TE[sp],TS[sp]),
                               e = TE[sp], s = TS[sp], sp = sp, ..., 
                               maximum = T)$maximum
    
    nOccs <- Noccs[sp]
    
    # placing occurrences along lineage's age using the accept-reject method
    res <- vector()
    while (length(res) < nOccs) {
      t <- runif(n = 1, max = TS[sp], min = TE[sp])
      
      test <- runif(n = 1, min = 0, max = adFunMaxApprox*1.05) <= 
        Pres_time_adppNorm(t = t, s = TS[sp], e = TE[sp], sp = sp, ...)
      
      test <- runif(n = 1, min = 0, max = Pres_time_adppNorm(t = adFunMaxApprox, 
          s = TS[sp], e = TE[sp], sp = sp, ...)*1.05)<= Pres_time_adppNorm(t = t,
            s = TS[sp], e = TE[sp], sp = sp, ...)
      
      if (test) {
        res <- c(res, t)
      }
    }
    occs[[sp]] <- res
  }
  
  zeroOccs <- which(lapply(occs, length) == 0)
  message(paste0(length(zeroOccs), " species left no fossil"))
  return(occs)
}