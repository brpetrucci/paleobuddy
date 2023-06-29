#' Non-constant rate Birth-Death simulation
#'
#' Simulates a species birth-death process with general rates for any number of
#' starting species. Allows for the speciation/extinction rate to be (1) a 
#' constant, or (2) a function of time. Allows for constraining results on the 
#' number of species at the end of the simulation, either total or extant. The 
#' function can also take an optional shape argument to generate age-dependence
#' on speciation and/or extinction, assuming a Weibull distribution as a model 
#' of age-dependence. Returns a \code{sim} object (see \code{?sim}). It may 
#' return calculated extinction times or simply information on whether species 
#' lived after the maximum simulation time, depending on input. For constant 
#' rate simulations, see \code{bd.sim.constant}. For a function that unites all
#' scenarios, see \code{bd.sim}. \code{bd.sim} also allows for extra inputs, 
#' creating a time-dependent only rate internally through \code{make.rate}. For
#' similar flexibility, use \code{make.rate} to generate the desired rate.
#' Please note while time runs from \code{0} to \code{tMax} in the simulation,
#' it returns speciation/extinction times as \code{tMax} (origin of the group) 
#' to \code{0} (the "present" and end of simulation), so as to conform to other
#' packages in the literature.
#'
#' @inheritParams bd.sim
#' 
#' @param lambda Function that defines the speciation rate over time. It will 
#' either be interpreted the rate of an exponential distribution, or the scale 
#' of a Weibull distribution if \code{lShape != NULL}. Can be constant, to 
#' allow for mixing of constant and non-constant rates. One can use constructs 
#' such as \code{ifelse()} to create rates whose underlying model change over 
#' time (see the last examples). Note that \code{lambda} should always be 
#' greater than or equal to zero.
#'
#' @param mu Similar to above, but for the extinction rate.
#' 
#' Note: rates should be considered as running from \code{0} to \code{tMax}, as
#' the simulation runs in that direction even though the function inverts 
#' speciation and extinction times before returning.
#' 
#' Note: this function is meant to be called by \code{bd.sim}, so it neither
#' allows for as much flexibility, nor calls \code{make.rate}. If the user 
#' wishes to use \code{bd.sim.general} with environmental or step-function 
#' rates, they can generate the rate with \code{make.rate} and supply it to the
#' function.
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation. See 
#' \code{?sim}.
#'
#' @author Bruno do Rosario Petrucci.
#'
#' @noRd
#' 

bd.sim.general <- function(n0, lambda, mu, tMax = Inf, N = Inf,
                         lShape = NULL, mShape = NULL,
                         nFinal = c(0, Inf), nExtant = c(0, Inf),
                         trueExt = FALSE) {
  # error check - rate cannot be negative
  if ((is.numeric(lambda))) {
    if (lambda < 0) {
      stop("speciation rate cannot be negative")
    }
  }
  
  else {
    if (optimize(lambda, interval = c(0, 1e10))$objective < 0) {
      stop("speciation rate cannot be negative at any point in time")
    }
  }
  
  if ((is.numeric(mu))) {
    if (mu < 0) {
      stop("extinction rate cannot be negative")
    }
  }
  
  else {
    if (optimize(mu, interval = c(0, 1e10))$objective < 0) {
      stop("extinction rate cannot be negative at any point in time")
    }
  }
  
  # check that n0 is not negative
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
  
  # if shape is not null, make scale a function to facilitate checking
  if (!is.null(lShape)) {
    message("since lShape is not null, lambda will be a Weibull scale and
            therefore correspond to 1/rate. See ?bd.sim or ?bd.sim.general")
    
    if (is.numeric(lambda)) {
      l <- lambda
      lambda <- Vectorize(function(t) l)
    }
    
    # check that it is never <= 0
    if (is.numeric(lShape)) {
      if (lShape <= 0) {
        stop("lShape may be nonpositive. It must always be >0")
      }
    }
    
    else {
      if (optimize(lShape, interval = c(0, 1e10))$objective < 0.01) {
        stop("lShape may be nonpositive. It must always be >0")
      }
    }
  }  
  
  if (!is.null(mShape)) {
    message("since mShape is not null, mu will be a Weibull scale and therefore
            correspond to 1/rate. See ?bd.sim or ?bd.sim.general")
    
    if (is.numeric(mu)) {
      m <- mu
      mu <- Vectorize(function(t) m)
    }
    
    # check that it is never <= 0
    if (is.numeric(mShape)) {
      if (mShape <= 0) {
        stop("mShape may be nonpositive. It must always be >0")
      }
    }
    
    else {
      if (optimize(mShape, interval = c(0, 1e10))$objective < 0.01) {
        stop("mShape may be nonpositive. It must always be >0")
      }
    }
  }
  
  # check which of N or tMax is not Inf, and condition as needed
  if ((N == Inf) && (tMax == Inf)) {
    stop("Either tMax or N must not be Inf.")
  } else if ((N != Inf) && (tMax != Inf)) {
    stop("Only one condition can be set, 
         so only one of N or tMax can be non-Inf")
  } else if (N != Inf) {
    condition = "number"
  } else {
    condition = "time"
  }
  
  # do we condition on the number of species?
  condN <- condition == "number"
  
  # whether our condition is met - rejection sampling for
  # time, or exactly number of species at the end for number
  condMet <- FALSE
  
  # counter to make sure the nFinal is achievable
  counter <- 1

  while (!condMet) {
    # create vectors to hold times of speciation, extinction, 
    # parents and status
    TS <- rep(0, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
  
    # initialize species count
    sCount <- 1
  
    # while we have species to be analyzed still
    while (length(TE) >= sCount) {
      # if sCount > nFinal[2], no reason to continue
      if (sCount > nFinal[2]) {
        # so we fail the condMet test
        sCount <- Inf
        
        # leave while
        break
      }
      
      # start at the time of speciation of sCount
      tNow <- TS[sCount]

      # find the waiting time using rexp.var if lambda is not constant
      waitTimeS <- ifelse(
        is.numeric(lambda), ifelse(lambda > 0, rexp(1, lambda), Inf),
        ifelse(lambda(tNow) > 0, 
               rexp.var(1, lambda, tNow, tMax, lShape, TS[sCount],
                        fast = !condN), Inf))
      waitTimeE <- ifelse(
        is.numeric(mu), ifelse(mu > 0, rexp(1, mu), Inf),
        ifelse(mu(tNow) > 0,
               rexp.var(1, mu, tNow, tMax, mShape, TS[sCount],
                        fast = !(trueExt || condN), Inf)))
      # fast for extinction depends on whether we want to record
      # true values when they are higher than tMax
  
      tExp <- tNow + waitTimeE
  
      # while there are fast enough speciations before the species 
      # goes extinct,
      while ((tNow + waitTimeS) < min(tExp, tMax)) {
  
        # advance to the time of speciation
        tNow <- tNow + waitTimeS
  
        # add new times to the vectors
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sCount)
        isExtant <- c(isExtant, TRUE)
  
        # get a new speciation waiting time, and include it in the vector
        waitTimeS <- ifelse(
          is.numeric(lambda), ifelse(lambda > 0, rexp(1, lambda), Inf),
          ifelse(lambda(tNow) > 0, 
                 rexp.var(1, lambda, tNow, tMax, lShape, TS[sCount],
                          fast = !condN), Inf))
        # fast is true since we do not record speciations after
        # the extinction anyway
      }
  
      # reached the time of extinction
      tNow <- tExp
  
      # if trueExt is true or the species went extinct before tMax,
      # record it. If both are false record it as NA
      TE[sCount] <- ifelse(tNow < tMax | trueExt, tNow, NA)
      
      # record the extinction -
      isExtant[sCount] <- ifelse(is.na(TE[sCount]) | TE[sCount] > tMax,
                                 TRUE, FALSE)
  
      # next species
      sCount <- sCount + 1
      
      # if we passed the limit
      if (condN && (sum(TS < tNow & (is.na(TE) | TE > tNow)) > 10*N)) {
        # function to find the excess at t
        nAlive <- Vectorize(function(t) {
          sum(TS <= t & (is.na(TE) | TE > t)) - N
        })
        
        # vector of all events
        events <- sort(c(TS, TE))
        
        # find times when we were at N
        nAliveT <- events[which(nAlive(events) == 0)]
        
        # find the times where the next event happened for each
        nextEvent <- unlist(lapply(nAliveT, function(t) events[events > t][1]))
        
        # get chosen event index
        eventChosen <- sample(1:length(nextEvent), 1, 
                              prob = (nextEvent - nAliveT))
        
        # draw uniform as above
        tMax <- runif(1, nAliveT[eventChosen], nextEvent[eventChosen])
        
        # adjust isExtant to TRUE for those alive at tMax
        isExtant[TE >= tMax] <- TRUE
        
        # adjust TE to NA for those alive at tMax
        TE[isExtant] <- NA
        
        # set condMet to true
        condMet <- TRUE
        
        # leave while
        break
      }
    }
  
    # now we invert TE and TS so time goes from tMax to 0
    TE <- tMax - TE
    TS <- tMax - TS
    
    # check which species are born after tMax
    nPrune <- which(TS <= 0)
    
    # if any, prune them
    if (length(nPrune) > 0) {
      TE <- TE[-nPrune]
      TS <- TS[-nPrune]
      isExtant <- isExtant[-nPrune]
      
      # need to be careful with parent
      parent <- c(NA, unlist(lapply(parent[-nPrune][-1], function(x)
        x - sum(nPrune < x))))
      # each species pruned lowers the parent numbering before them by 1
    }

    # check whether we are in bounds
    if (!condN) {
      condMet <- (length(TE) >= nFinal[1]) &&
        (length(TE) <= nFinal[2]) &&
        (sum(isExtant) >= nExtant[1]) &&
        (sum(isExtant) <= nExtant[2])
    }
    
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
  
  return(sim)
}
