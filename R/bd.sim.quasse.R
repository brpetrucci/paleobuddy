#' QuaSSE simulation
#'
#' Simulates a species birth-death process following the Quantitative State 
#' Speciation and Extinction (QuaSSE) model for any number of starting species. 
#' Allows for the speciation/extinction rate to be (1) a constant, or (2) a 
#' function of trait values. Traits are simulated to evolve under a Brownian
#' motion model (see references). Allows for constraining results on the 
#' number of species at the end of the simulation, either total or extant, 
#' using rejection sampling. Returns a \code{sim} object (see \code{?sim}), and
#' a list of data frames describing trait values for each interval. It may
#' return true extinction times or simply information on whether species lived
#' after the maximum simulation time, depending on input.
#' 
#' Please note while time runs from \code{0} to \code{tMax} in the simulation, 
#' it returns speciation/extinction times as \code{tMax} (origin of the group) 
#' to \code{0} (the "present" and end of simulation), so as to conform to other
#' packages in the literature.
#'
#' @inheritParams bd.sim.musse
#' 
#' @param lambda Function to hold the speciation rate over time. It should
#' either be a constant, or a function of one argument. For each species a 
#' trait evolution simulation will be run, and then used to calculate the final 
#' speciation rate. Note that \code{lambda} should always be greater than or 
#' equal to zero.
#'
#' @param mu Similar to above, but for the extinction rate.
#' 
#' @param X0 Initial trait value for original species. Can be a constant or a 
#' vector of length \code{nTraits}.
#' 
#' @param sigma Brownian motion parameter. The variance at time \code{t} after
#' the start of the simulation will be \code{sigma^2*t}
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation, and a 
#' list object with the trait functions describing the trait value for each
#' species at each time.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @references 
#' 
#' Maddison W.P., Midford P.E., Otto S.P. 2007. Estimating a binary character’s 
#' effect on speciation and extinction. Systematic Biology. 56(5):701.
#' 
#' FitzJohn R.G. 2012. Diversitree: Comparative Phylogenetic Analyses of 
#' Diversification in R. Methods in Ecology and Evolution. 3:1084–1092.
#'
#' @examples
#'
#' ###
#' 
#' @name bd.sim.quasse
#' @rdname bd.sim.quasse
#' @export

bd.sim.quasse <- function(n0, lambda, mu, condition = "time",
                         tMax = Inf, N = Inf,
                         nTraits = 1, nFocus = 1, 
                         X0 = 0, sigma = list(1),
                         drift = 0, bounds = NULL,
                         nFinal = c(0, Inf), nExtant = c(0, Inf)) {
  # check that n0 is not negative
  if (n0 <= 0) {
    stop("initial number of species must be positive")
  }
  
  # check nFinal's length
  if (length(nFinal) != 2) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
  }
  
  ### HAVE TO write error tests, including for type + tMax and N
  
  # check whether rates are trait-dependent
  tdLambda <- !is.numeric(lambda)
  tdMu <- !is.numeric(mu)
  
  # do we condition on the number of species?
  condN <- condition == "number"
  
  # if conditioning on number, need to set some 
  # tMax for the trait evolution
  if (condN) {
    # the average time it will take to get to 10*N species
    # under the average diversification rate parameters
    trTMax <- log(10*N) / (ifelse(is.numeric(lambda), lambda, lambda(X0)) -
                             ifelse(is.numeric(mu), mu, mu(X0)))
  } else trTMax <- Inf
  
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
    now <- c()
    
    # initialize species count
    sCount <- 1
    
    # initialize list of traits for first species' traits
    traits <- list(traits.quasse(tMax = min(tMax, trTMax), tStart = 0, 
                                 nTraits = nTraits, X0 = X0, sigma = sigma,
                                 bounds = bounds, drift = drift))
    
    # while we have species to be analyzed still
    while (length(TE) >= sCount) {
      # if sCount > nFinal[2], no reason to continue
      if (sCount > nFinal[2]) {
        # so we fail the condMet test
        sCount <- Inf
        
        # leave while
        break
      }
      
      # get argument for the next species
      # it will be the oldest one that hasn't lived yet
      sp <- which(TS == sort(TS)[sCount])
      
      # start at the time of speciation of sp
      tNow <- TS[sp]
      
      # get traits data set
      traitsSp <- traits[[sp]][[nFocus]]
      
      # find the waiting time using rexp.var if lambda is not constant
      waitTimeS <- ifelse(tdLambda, 
                          ifelse(sum(lambda(traitsSp$value)) > 0,
                                 rexp.traits(1, lambda(traitsSp$value), 
                                             traitsSp,
                                             tNow, tMax), Inf),
                          ifelse(lambda > 0, 
                                 rexp(1, lambda), Inf))
      
      waitTimeE <- ifelse(tdMu,
                          ifelse(sum(mu(traitsSp$value)) > 0,
                                 rexp.traits(1, mu(traitsSp$value), 
                                             traitsSp,
                                             tNow, tMax), Inf),
                          ifelse(mu > 0, 
                                 rexp(1, mu), Inf))
      
      tExp <- tNow + waitTimeE
      
      # while there are fast enough speciations before the species 
      # goes extinct,
      while ((tNow + waitTimeS) < min(tExp, tMax)) {
        # advance to the time of speciation
        tNow <- tNow + waitTimeS
        now <- c(now, tNow)
        
        # add new times to the vectors
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sp)
        isExtant <- c(isExtant, TRUE)
        
        # get trait value at tNow
        XPar <- tail(traitsSp$value[traitsSp$min < tNow], 1)
        
        # run trait evolution and append it to list
        traits[[length(traits) + 1]] <- 
          traits.quasse(tMax = min(tMax, max(trTMax, tNow + 100)), 
                        tStart = tNow, 
                        nTraits = nTraits, X0 = XPar, 
                        sigma = sigma, bounds = bounds, drift = drift)
        
        # get a new speciation waiting time
        waitTimeS <- ifelse(tdLambda, 
                            ifelse(sum(lambda(traitsSp$value)) > 0,
                                   rexp.traits(1, lambda(traitsSp$value), 
                                               traitsSp, 
                                               tNow, tMax), Inf),
                            ifelse(lambda > 0, 
                                   rexp(1, lambda), Inf))
      }
      
      # reached the time of extinction
      tNow <- tExp
      
      # if the species went extinct before tMax,
      # record it, otherwise, record NA
      TE[sp] <- ifelse(tNow < tMax, tNow, NA)
      
      # record the extinction
      isExtant[sp] <- is.na(TE[sp]) | TE[sp] > tMax
      
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
    
    # truncate traits so we only go up to tMax
    for (i in 1:length(traits)) {
      for (j in 1:nTraits) {
        # df in question
        traitsSp <- traits[[i]][[j]]
        
        # eliminate rows with min greater than tMax
        traitsSp <- traitsSp[traitsSp$min < tMax, ]
        
        # set max of last row to tMax
        traitsSp$max[nrow(traitsSp)] <- tMax
        
        # invert time for max and min
        traitsSp$max <- tMax - traitsSp$max
        traitsSp$min <- tMax - traitsSp$min
        
        # invert columns
        colnames(traitsSp) <- c("value", "max", "min")
        
        # set traits back to it
        traits[[i]][[j]] <- traitsSp
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
      traits <- traits[-nPrune]
      
      # need to be careful with parent
      parent <- c(NA, unlist(lapply(parent[-nPrune][-1], function(x)
        x - sum(nPrune < x))))
      # each species pruned lowers the parent numbering before them by 1
    }
    
    # check whether we are in bounds if rejection sampling is the thing
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
  
  res <- list("TRAITS" = traits, "SIM" = sim)
  
  return(res)
}

###
# function to run trait evolution for a given species
traits.quasse <- function(tMax, tStart = 0, nTraits = 1, 
                          X0 = 0, sigma = list(1),
                          drift = 0, bounds = NULL) {
  # create a return list
  traits <- vector(mode = "list", length = nTraits)
  
  # for each trait
  for (i in 1:nTraits) {
    # get initial state
    X0I <- ifelse(length(X0) > 1, X0[i], X0)
    
    # get drift
    driftI <- ifelse(length(drift) > 1, drift[i], drift)
    
    # get sigma
    if (length(sigma) > 1) sigmaI <- sigma[[i]] else sigmaI <- sigma[[1]]
    
    # get bounds
    if (length(bounds) > 1) boundsI <- bounds[[i]] else boundsI <- bounds[[1]]
    
    # append traits data frame to traits
    traits[[i]] <- trait.quasse(tMax, tStart, X0I, sigmaI, drift, bounds)
  }
  
  return(traits)
}

trait.quasse <- function(tMax, tStart = 0, 
                         X0 = 0, sigma = 1,
                         drift = 0, bounds = NULL,
                         nPoints = 100) {
  # make sure tMax and tStart are numbers
  if (!is.numeric(c(tMax, tStart))) {
    stop("tMax and tStart must be numeric")
  } else if (length(tMax) > 1 || length(tStart) > 1) {
    stop("tMax and tStart must be one number")
  }
  
  # make sure tMax > tStart
  if (tStart >= tMax) {
    stop("tMax must be greater than tStart")
  }
  
  # update nPoints based on the time range
  nPoints <- nPoints * (tMax - tStart)

  # calculate vector of times
  times <- seq(tStart, tMax, (tMax - tStart) / (nPoints - 1))
  
  # calculate nPoints - 1 randomly distributed normal
  # variates with variance sigma2*t
  jumps <- rnorm(nPoints - 1, mean = 0, 
                 sd = sigma*sqrt((tMax - tStart) / (nPoints - 1)))
  
  # calculate the bm vector being X0 and the cumulative sum of jumps
  bm <- c(X0, X0 + cumsum(jumps))
  
  # bound it, if bounds exists
  if (!is.null(bounds)) {
    bm <- ifelse(bm < bounds[1], bounds[1], bm)
    bm <- ifelse(bm > bounds[2], bounds[2], bm)
  }
  
  # insert drift, if necessary
  bm <- bm + drift * times
  
  # create traits data frame
  trait <- data.frame(value = bm, min = times, 
                      max = c(times[2:length(times)], Inf))
  
  return(trait)
}

values <- function(bmList, t) {
  vals <- c()
  
  for (i in 1:length(bmList)) {
    vals <- c(vals, bmList[[i]]$value[t])
  }
  
  return(vals)
}