#' Time-dependent and age-dependent rate species sampling
#' 
#' Generates a vector of occurrence times for a species in a simulation using 
#' a Poisson process coupled with age-dependent fossil sampling. Allows for the 
#' Poisson rate to be (1) constant or (2) time-dependent, and additionally 
#' allows the distribution of occurrences throughout a species age to be any 
#' function between speciation and extinction time. For more flexibility, see 
#' \code{make.rate} and \code{sample.clade}. Note that while the simulation 
#' occurs in forward time, we return (both in birth-death functions and here)
#' results in backwards time, so that time is inverted using \code{tMax} both 
#' at the beginning and end of the function.
#' 
#' @param rho Average fossil sampling rate (per species per million years) over 
#' time. It can be a \code{numeric} describing a constant rate or a 
#' \code{function(t)} describing the variation in sampling over time. For more 
#' flexibility on sampling, see \code{make.rate} to create more complex rates.
#' If \code{adFun} is supplied, it will be used to find the number of 
#' occurrences during the species duration, and a normalized \code{rho*adFun} 
#' will determine their distribution along the species duration. Note that 
#' \code{rho} should always be greater than or equal to zero.
#'
#' @param adFun A density function representing the age-dependent preservation
#' model. It must be a density function, and consequently integrate to 1 
#' (though this condition is not verified by the function). If not provided, a 
#' uniform distribution will be used by default. The function must also 
#' have the following properties:
#'
#' \itemize{
#'
#' \item Return a vector of preservation densities for each time in a given
#' vector \code{t} in geological time. 
#'
#' \item Be parametrized in the absolute geological time associated to 
#' each moment in age (i.e. age works relative to absolute geological 
#' time, in Mya - in other words, the convention is TS > 0). The function 
#' \emph{does not} directly use the lineage's age (which would mean that
#' TS = 0 for all species whenever they are born). Because of this, it is
#' assumed to go from \code{tMax} to \code{0}, as opposed to most functions 
#' in the package.
#'
#' \item Should be limited between \code{s} (i.e. the lineage's 
#' speciation/birth) and \code{e} (i.e. the lineage's extinction/death), 
#' with \code{s} > \code{e}. It is possible to assign parameters in absolute 
#' geological time (see third example) and relative to age as long as this 
#' follows the convention of age expressed in absolute geological time (see 
#' fourth example).
#'
#' \item Include the arguments \code{t}, \code{s}, \code{e} and \code{sp}. 
#' The argument sp is used to pass species-specific parameters (see examples),
#' allowing for \code{dFun} to be species-inhomogeneous.
#' }
#'
#' @param ... Additional parameters used by \code{adFun}. See examples.
#' 
#' @inheritParams sample.time
#'
#' @return A list of vectors of occurrence times for each species in \code{S}.
#'
#' @author Matheus Januario
#' 
#' @name sample.age.time
#' @rdname sample.age.time
#' @keywords Internal

sample.age.time <- function(sim, rho, tMax, S = NULL, adFun = NULL, ...){
  # checking input
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # convert NA elements of TE to 0
  sim$TE[sim$EXTANT] <- 0
  
  # if rho is constant, make it a function
  if (is.numeric(rho)) {
    r <- rho
    rho <- Vectorize(function(t) r)
  }
  
  # if S is not given
  if (is.null(S)) {
    S <- 1:length(sim$TE)
  }
  
  # get extinction and speciation times
  TE <- sim$TE
  TS <- sim$TS
  
  # check if adFun is a function with the right parameters
  if (!(is.null(adFun))) {
    if (sum(c("t", "s", "e", "sp") %in% names(formals(adFun))) < 3) {
      stop("adFun must have \"t\", \"s\", \"e\", and  \"sp\" parameters")
    }
  }
  
  # if adFun is null, it is the same as sample.time
  if (is.null(adFun)) {
    message("Preservation will be age-independent \n")
    
    # run sample.time
    occs <- sample.time(sim, rho, tMax, S)
    
    # check how many species left no occurrences
    zeroOccs <- which(lapply(occs, length) == 0)
    message(paste0(length(zeroOccs), " species left no fossil"))
    return(occs)
  }
  
  # create function combining adFun and rho
  Pres_time_adpp <- function(t, s, e, sp, ...) {
    # this correction is needed because "t" in rho means
    # "time since clade origin" and NOT absolute geologic time in Mya
    rhoMod <- function(t){ 
      return(rho(tMax - t))
    }
    return(rhoMod(t)*adFun(t = t, s = s, e = e, sp = sp, ...))
  }
  
  # normalize it
  Pres_time_adppNorm <- function(t, s, e, sp,...) {
    return(Pres_time_adpp(t = t, s = s, e = e, sp = sp, ...) / 
             integrate(Pres_time_adpp, lower = e, upper = s, s = s, 
                       e = e, sp = sp, ...)$value)
  }
  
  # get the number of occurrences following a time-varying poisson process
  Noccs <- vector(length = length(sim$TE))
  Noccs <- unlist(lapply(sample.time(sim = sim, rho = rho, tMax = tMax,
                                     S = S), FUN = length))
  
  occs <- list()
  
  # for each species to sample
  for (sp in S) {
    # get the maximum of function for that specific species
    adFunMaxApprox <- optimize(Pres_time_adppNorm, interval = c(TE[sp],TS[sp]),
                               e = TE[sp], s = TS[sp], sp = sp, ..., 
                               maximum = T)$maximum
    
    nOccs <- Noccs[which(sp == S)]
    
    # place occurrences along lineage's age using the accept-reject method
    res <- vector()
    while (length(res) < nOccs) {
      t <- runif(n = 1, max = TS[sp], min = TE[sp])
      
      test <- runif(n = 1, min = 0, max = adFunMaxApprox*1.05) <= 
        Pres_time_adppNorm(t = t, s = TS[sp], e = TE[sp], sp = sp, ...)
      
      test <- runif(n = 1, min = 0, max = Pres_time_adppNorm(t = adFunMaxApprox, 
          s = TS[sp], e = TE[sp], sp = sp, ...)*1.05) <= 
        Pres_time_adppNorm(t = t,
            s = TS[sp], e = TE[sp], sp = sp, ...)
      
      if (test) {
        res <- c(res, t)
      }
    }
    occs[[paste0("t", sp)]] <- res
  }
  
  # check which species left no fossils
  zeroOccs <- which(lapply(occs, length) == 0)
  message(paste0(length(zeroOccs), " species left no fossil"))
  return(occs)
}
