#' Age-dependent Poisson Process species sampling
#'
#' Generates a list of occurrence times for each of the desired species using a
#' Poisson process with constant average rate and occurrences distributed based on
#' a model of species age. Allows for any distribution as the model of occurrences
#' during a species life, given certain requirements (see below). Allows for an
#' optional argument - the maximum of the distribution - that can make the
#' simulation faster. Also allows for extra arguments the age-dependent
#' preservation function may take. For time-varying sampling rates without 
#' age-dependency, see \code{sample.species}. For a function that unites both
#' cases and returns an organize data frame, see \code{sample.clade}.
#'
#' @param S A vector of species numbers to be sampled. Could be only a subset of 
#' the species if the user wishes. The default is all species in \code{sim}. 
#' Species not included in \code{S} will not be sampled by the function.
#'
#' @param sim A \code{sim} object, usually an output of \code{bd.sim}.
#'
#' @param rr A mean sampling rate (equivalent to the \code{lambda} of a Poisson)
#' in the Poisson process. Must be a number greater than or equal to zero.
#'
#' @param dFun A density function representing the age-dependent
#' preservation model. It must be a density function, and consequently integrate
#' to 1 (though this condition is not verified by the function). Additionally, the
#' function must have the following properties:
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
#' @param dFunMax A function that calculates the maximum (density) value of
#' \code{dFun} using its arguments, excluding \code{t}. It can also be a number 
#' representing the maximum density.
#'
#' Note: if it not provided, it will be approximated numerically, leading to
#' longer running times.
#'
#' @param ... Additional parameters related to \code{dFun} and \code{dFunMax}.
#'
#' @return A list of occurrences for that species, expected to be around 
#' \code{(s - e)*rr} occurrences, with their distribution in species relative 
#' time given by the \code{dFun} function provided by the user.
#'
#' @author Matheus Januario.
#'
#' @examples
#'
#' ###
#' # we can start with a hat-shaped increase through the duration of a species
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
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
#' # function to calculate max of the PERT
#' dPERTmax <- function(s, e, sp) {
#'   return(((s - e) / 2) + e)
#' }
#' 
#' # find occurrences
#' occs <- sample.adpp(sim = sim, rr = 1, dFun = dPERT, S = 1, dFunMax = dPERTmax)
#' 
#' # check histogram
#' hist(unlist(occs), probability = TRUE)
#' 
#' # expected curve - probably will not fit great because of low sample size, see
#' # vignettes for more rigorous tests
#' curve(dPERT(x, s = sim$TS[1], e = sim$TE[1]), 20, 0, add = TRUE, col = "red")
#' 
#' ###
#' # now we can test the simpler scenario of uniform sampling probablity
#' # through the duration of a species (= homogeneous poisson process)
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
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
#' # we will not give a dFunMax function this time. sample.adpp() will try to find
#' # the maximum density with a very simple numerical simulation
#' occs <- sample.adpp(sim = sim, rr = 2, dFun = custom.uniform, S = 1)
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE)
#' 
#' # expected curve
#' curve(dunif (x, min = sim$TE[1], max = sim$TS[1]), 10, 0,
#'       add = TRUE, col = "red")
#' 
#' ###
#' # now, a hat-shaped increase through the duration of a species with more
#' # parameters than TS and TE
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
#' 
#' # here we will use the triangular distribution. We have some empirical evidence
#' # that taxa occurrences might present triangular shape
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
#' # the dFunMax function must have the same parameters then the dFun function,
#' # even if they do not use them
#' dTRImax <- function(s, e, sp, md) {
#'   
#'   return(2 / (s - e))
#' }
#' 
#' # note we are providing the mode for the triangular sampling as an ... argument
#' occs <- sample.adpp(sim = sim, rr = 2.5, dFun = dTRI,
#'                    dFunMax = dTRImax, S = 1, md = 8)
#' 
#' # please note in the original parametrization, the "md" parameter (mode) is
#' # in "absolute time", i.e. a specific number in the absolute scale. This is not,
#' # directly, an "age-dependent" model. We show it here as we will construct
#' # models related to age in the following examples (4 and 5). This is also
#' # the reason we only plot the first lineage in this specific example
#' # (many lineages might not be alive in the "md" moment in time).
#' 
#' # check histrogram
#' hist(unlist(occs[[1]]), probability = TRUE)
#' 
#' # expected curve
#' curve(dTRI(x, e = sim$TE[1], s = sim$TS[1], md = 8), 10, 0, add = TRUE,
#'       col = "red")
#' 
#' ###
#' # we can also have a hat-shaped increase through the duration of a species
#' # with more parameters than TS and TE, but with the parameters relate to
#' # the relative age of each lineage
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
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
#' # function to calculate max
#' dTRImaxmod1<-function(s, e, sp) {
#'   return(2 / (s - e))
#' }
#' 
#' # find occurrences
#' occs <- sample.adpp(sim = sim, rr = 5,
#'                     dFun = dTRImod1, S = 1, dFunMax = dTRImaxmod1)
#' 
#' # we do not have the "md" parameter (see example 3) as it corresponds to the
#' # last quarter of the duration of each lineage
#' 
#' # check histogram
#' hist(unlist(occs[[1]]), probability = TRUE)
#' 
#' # md of dTRI
#' mid <- ((sim$TS[1] - sim$TE[1]) / 4) + sim$TE[1]
#' 
#' # expected curve
#' curve(dTRImod1(x, e = sim$TE[1], s = sim$TS[1]), 10, 0, add = TRUE, 
#'       col = "red")
#' 
#' ###
#' # finally, a hat-shaped increase through the duration of a species with more
#' # parameters than TS and TE, but the parameters relate to each specific lineage.
#' # This is useful when the user wants to use variable parameters for each species
#' # but wants to keep track of those parameters after the sampling is over
#' 
#' # simulate a group
#' sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' 
#' # in case first simulation is short-lived
#' while ((sim$TS[1] - ifelse(is.na(sim$TE[1]), 0, sim$TE[1])) < 10) {
#'   sim <- bd.sim(n0 = 1, pp = 0.1, qq = 0.1, tMax = 10)
#' }
#' 
#' # we will need to get exact durations for some examples, so
#' sim$TE[sim$EXTANT] <- 0
#' # this is necessary since the default is to have NA for extant species
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
#' # maximum function
#' dTRImaxmod2 <- function(s, e, sp) {
#'   return(2 / (s - e))
#' }
#' 
#' # get the par and par1 vectors
#' par <- runif (n = length(sim$TE), min = sim$TE, max = sim$TS)
#' par1 <- (((sim$TS - sim$TE)/2) + sim$TE) - par
#' 
#' # find occurrence list
#' occs <- sample.adpp(sim = sim, rr = 10,
#'                     dFun = dTRImod2, dFunMax = dTRImaxmod2)
#' 
#' # we do not have the "md" parameter (see example 3) as it corresponds to
#' # the last quarter of the duration of each lineage
#' 
#' # checking each species
#' for (sp in 1:length(sim$TE)) {
#'   # check histogram
#'   hist(unlist(occs[[sp]]), probability = TRUE,
#'        main = paste0("spp ", sp, " ; duration ~ ",
#'                      round(sim$TS[sp] - sim$TE[sp], digits = 2)))
#'   
#'   # calculate mid
#'   mid <- par[sp] + par1[sp]
#'   
#'   # expected curve
#'   curve(dTRImod2(x, e = sim$TE[sp], s = sim$TS[sp], sp), 10, 0, add = TRUE,
#'         col = "red", n = 100)
#'   abline(v = mid, col = "red")
#' }
#'
#' @name sample.adpp
#' @rdname sample.adpp
#' @export

sample.adpp <- function(sim, rr, dFun, S = NULL, dFunMax = NULL, ...) {
  # some error checking
  if (sum(S %in% 1:length(sim$TS)) != length(S)) {
    stop("One or more species numbers provided not in simulation")
  }
  
  # check rr is a constant greater than or equal to 0
  if (!is.numeric(rr) | (is.numeric(rr) & length(rr) != 1) |
      (is.numeric(rr) & length(rr) == 1 & rr < 0)) {
    stop("rr must be a number greater than or equal to 0")
  }
      
  
  # make S all species if it is NULL
  if (is.null(S)) {
    S = 1:length(sim$TE)
  }
  
  # get the speciation and extinction times vectors
  TE <- sim$TE
  TS <- sim$TS
  
  # make TE sensible
  TE[sim$EXTANT] <- 0
    
  # setting things and checking inputs
  printMessage <- TRUE
  
  # checking dFun has all arguments it needs
  if (sum(c("t", "s", "e", "sp") %in% names(formals(dFun))) < 3) {
    stop("dFun must have \"t\", \"s\", \"e\", and  \"sp\" parameters")
  }

  # checking whether dFunMax exists and is a function
  if (!(is.null(dFunMax))) { 
    if (!(is.numeric(dFunMax))) {
      # it must also have those parameters
      if (sum(c("t", "s", "e", "sp") %in% names(formals(dFunMax))) < 3) {
        stop("dFunMax must have \"t\", \"s\", \"e\", and  \"sp\" parameters")
      }

      # and it must also have the same parameters as dFun
      if (sum(names(formals(dFunMax)) %in% names(formals(dFun))
             [-which(names(formals(dFun)) == "t")])
         < length(names(formals(dFunMax)))) {
        stop("dFun and dFunMax must have the same arguments
             (with the exception of \"t\" argument for dFun\")")
      }

    }
  }

  # initialize occurrence list
  occs <- list()
  
  # for each lineage in the dataset
  for (sp in S) { 
    # if dFunMax is not provided by the user
    if (is.null(dFunMax)) { 
      # it is slower
      if (printMessage) {
        message("Please wait. The function will use approximate maximum point 
                for the function and that might take a while")
        printMessage <- FALSE
      }

      # a simple/lazy but useful approximation: lower threshold until it is 
      # smaller than the maximum dFun calculating at each 0.1my, then 
      # do one step back.
      
      # get the value of dFun
      aux <- dFun(seq(TE[sp], TS[sp], by=.01), e = TE[sp], s = TS[sp], ...)
      
      # start a treshold - slightly below 1
      threshold <- 0.99
      
      # whether the value anywhere is higher than threshold
      test <- aux > threshold
      
      # if the value is higher than the threshold somewhere
      if (sum(test) >= 1) {
        # get new max
        threshold <- max(aux)
      } 
      else {
        # if not
        while (sum(test) < 1) {
          # get new treshold
          threshold <- threshold - 0.01
          
          # redo test
          test <- aux > threshold
        }
      }
      
      # max approximation is slightly above treshold so we capture that point
      dFunMaxApprox <- threshold + 0.01
    }

    # number of occurrences following a poisson process:
    nOccs <- rpois(1, lambda = rr*(TS[sp] - TE[sp]))

    # time of each occurrence:
    res <- vector()

    # accept-reject method for Monte Carlo generation of random numbers from 
    # a density distribution
    
    # while we still have occurrences to check
    while (length(res) < nOccs) {
      # sample a number within species duration:
      t <- runif (n = 1, max = TS[sp], min = TE[sp])

      # calculate the density of t and check if it is smaller than a random 
      # sampled number between 0 and dFunMax
      
      # if we need to use dFunMaxApprox
      if (is.null(dFunMax)) {
        # check whether uniform variable is smaller
        test <- runif (n = 1, min = 0, max = dFunMaxApprox) <= 
          dFun(t = t, s = TS[sp], e = TE[sp], sp=sp, ...)
      } 
      
      # if we have exact max
      else {
        # if it is a number
        if (is.numeric(dFunMax)) { 
          message(paste0("dFunMax = ", dFunMax, 
                         " will be assumed as the maximum value of dFun"))
          # check whether uniform variable is smaller
          test <- runif (n = 1, min = 0, max = dFunMax) <= 
            dFun(t = t, s = TS[sp], e = TE[sp], sp = sp, ...)
        } 
        # if it is a function
        else {
          test <- runif(n = 1, min = 0, 
                        max = dFunMax(s = TS[sp], e = TE[sp], sp = sp, ...)) <= 
            dFun(t = t, s = TS[sp], e = TE[sp], sp = sp, ...)
        }
      }
      # if the sampled number is smaller or equal to the density of t, 
      # append t to the result:
      if (test) { 
        res <- c(res, t)
      }
    }
    # append in the occurence vector
    occs[[sp]] <- res 

  }

  # warning about non-sampled species:
  zeroOccs <- which(lapply(occs, length) == 0)
  message(paste0(length(zeroOccs), " species left no fossil"))

  return(occs)
}
