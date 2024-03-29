#' Constant rate Birth-Death simulation
#' 
#' Simulates a species birth-death process with constant rates for any number 
#' of starting species. Allows for constraining results on the number of 
#' species at the end of the simulation, either total or extant. Returns a 
#' \code{sim} object (see \code{?sim}). It may return true extinction times or 
#' simply information on whether species lived after the maximum simulation 
#' time, depending on input. For time-varying and age-dependent simulations, 
#' see \code{bd.sim}. Please note while time runs from \code{0} to \code{tMax} 
#' in the simulation, it returns speciation/extinction times as \code{tMax} 
#' (origin of the group) to \code{0} (the "present" and end of simulation), so 
#' as to conform to other packages in the literature.
#' 
#' @inheritParams bd.sim 
#' 
#' @param lambda Speciation rate. Must be constant, and greater than or equal
#' to zero.
#'
#' @param mu Extinction rate, similar to above.
#' 
#' @param trueExt A \code{logical} indicating whether the function should return
#' true or truncated extinction times. When \code{TRUE}, time of extinction of 
#' extant species will be the true time, otherwise it will be \code{NA} if a 
#' species is alive at the end of the simulation.
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation. See 
#' \code{?sim}.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @noRd


bd.sim.constant <- function(n0, lambda, mu, tMax, 
                            nFinal = c(0, Inf), nExtant = c(0, Inf),
                            trueExt = FALSE) {
  # check that the rates are constant
  if (!(is.numeric(lambda) & length(lambda) == 1 &
        is.numeric(mu) & length(mu) == 1)) {
    stop("bd.sim.constant requires constant rates")
  }
  
  # check that the rates are non-negative
  if (lambda < 0 || mu < 0) {
    stop("rates cannot be negative")
  }
  
  # check that n0 is positive
  if (n0 <= 0) {
    stop("initial number of species must be positive")
  }
  
  # check nFinal is sensible - two numbers, maximum >=1
  if ((length(nFinal) != 2) || (typeof(nFinal) != "double")) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
  } else if (max(nFinal) < 1) {
    stop("nFinal must have a maximum number of species greater than 0")
  } else {
    # if everything is good, make sure it's sorted
    nFinal <- sort(nFinal)
  }
  
  # similarly for nExtant
  if ((length(nExtant) != 2) || (typeof(nExtant) != "double")) {
    stop("nExtant must be a vector with a minimum and maximum number 
         of species")
  } else if (max(nExtant) < 0) {
    stop("nExtant must have a maximum number of species greater 
         than or equal to 0")
  } else {
    # if everything is good, make sure it's sorted
    nExtant <- sort(nExtant)
  }
  
  # initialize test making sure while loop runs
  inBounds <- FALSE
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  while (!inBounds) {
    # initialize the vectors to hold times of speciation and extinction, parents
    # and status (extant or not)
    TS <- rep(0, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
  
    # initialize the counting variable
    sCount <- 1
  
    # while we have more species in a vector than we have analyzed,
    while (length(TE) >= sCount) {
      # start at the time of speciation of sCount
      tNow <- TS[sCount]
  
      # draw waiting times with rexp()
      waitTimeS <- ifelse(lambda > 0, rexp(1, lambda), Inf)
      waitTimeE <- ifelse(mu > 0, rexp(1, mu), Inf)
  
      # calculate when extinction will happen
      tExp <- tNow + waitTimeE
  
      # while there are speciations before the species goes extinct
      while ((tNow + waitTimeS) < min(tExp, tMax)) {
        # update time
        tNow <- tNow + waitTimeS
  
        # create a new species with corresponding TE, TS and parent
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sCount)
        isExtant <- c(isExtant, TRUE) # it is alive
  
        # take a new waiting time - if now + waitTimeS is still less than 
        # when the species goes extinct, repeat
        waitTimeS <- ifelse(lambda > 0, rexp(1, lambda), Inf)
      }
  
      # reached the time of the species extinction
      tNow <- tExp
      
      # if trueExt is true or the species went extinct before tMax,
      # record it. If both are false record it as NA
      TE[sCount] <- ifelse(tNow < tMax | trueExt, tNow, NA)
      
      # record the extinction -
      isExtant[sCount] <- ifelse(is.na(TE[sCount]) | TE[sCount] > tMax,
                                 TRUE, FALSE)
  
      # next species
      sCount <- sCount + 1
    }
  
    # finally, we invert both TE and TS to attain to the convention that time
    # runs from 0 to tMax
    TE <- tMax - TE
    TS <- tMax - TS
  
    # check whether we are in bounds
    inBounds <- (length(TE) >= nFinal[1]) &&
      (length(TE) <= nFinal[2]) &&
      (sum(isExtant) >= nExtant[1]) &&
      (sum(isExtant) <= nExtant[2])
    
    # if we have ran for too long, stop
    counter <- counter + 1
    if (counter > 100000) {
      warning("This value of nFinal took more than 100000 simulations 
              to achieve")
      return(NA)
    }
  }
  
  # create the return
  sim <- list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant)
  class(sim) <- "sim"
  
  # for more on the seem class, see ?sim

  return(sim)
}
