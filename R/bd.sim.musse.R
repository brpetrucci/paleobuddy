#' MuSSE simulation
#'
#' Simulates a species birth-death process following the Multiple State 
#' Speciation and Extinction (MuSSE) model for any number of starting species. 
#' Allows for the speciation/extinction rate to be (1) a constant, or (2) a 
#' list of values for each trait state. Traits are simulated to evolve under a
#' simple Mk model (see references). Allows for constraining results on the 
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
#' @inheritParams bd.sim
#' 
#' @param lambda Function to hold the speciation rate over time. It should
#' either be a constant, or a list of size \code{nStates}. For each species a 
#' trait evolution simulation will be run, and then used to calculate the final 
#' speciation rate. Note that \code{lambda} should always be greater than or 
#' equal to zero.
#'
#' @param mu Similar to above, but for the extinction rate.
#' 
#' @param condition Whether to condition tree by total simulation running time,
#' or the number of species alive at the end of the simulation. If set to 
#' \code{"time"}, simulation will run for \code{tMax} million years. If set
#' to \code{"number"}, simulation will run until there are \code{N} species
#' alive at a given time. 
#' 
#' Note: this manner of conditioning for a number of tips conditions the 
#' phylogeny to have shorter terminal branches. While this might not be an 
#' issue for low values of extinction fraction \code{mu/lambda}, future 
#' implementation of the correct conditioning is planned.
#' 
#' @param N Number of species at the end of the simulation, if \code{condition}
#' equals \code{"number"}. End of the simulation will be set for the first time
#' a species alive at a period where \code{N} species are alive would go 
#' extinct.
#' 
#' @param nTraits The number of traits to be considered. \code{lambda} and 
#' \code{mu} need not reference every trait simulated.
#' 
#' @param nFocus Trait of focus, i.e. the one that rates depend on. If it is 
#' one number, that will be the trait of focus for both speciation and 
#' extinction rates. If it is of length 2, the first will be the focus for
#' the former, the second for the latter.
#' 
#' @param nStates Number of possible states for categorical trait. The range
#' of values will be assumed to be \code{(0, nStates - 1)}. Can be a constant
#' or a vector of length \code{nTraits}, if traits are intended to have 
#' different numbers of states.
#' 
#' @param X0 Initial trait value for original species. Must be within 
#' \code{(0, nStates - 1)}. Can be a constant or a vector of length 
#' \code{nTraits}.
#' 
#' @param Q Transition rate matrix for continuous-time trait evolution. For
#' different states \code{i} and \code{j}, the rate at which a species at
#' \code{i} transitions to \code{j} is \code{Q[i + 1, j + 1]}. Must be within
#' a list, so as to allow for different \code{Q} matrices when
#' \code{nTraits > 1}.
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation, and a 
#' list object with the trait data frames describing the trait value for each
#' species at each specified interval.
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
#' @name bd.sim.musse
#' @rdname bd.sim.musse
#' @export

bd.sim.musse <- function(n0, lambda, mu, condition = "time",
                         tMax = Inf, N = Inf,
                         nTraits = 1, nFocus = 1, nStates = 2, X0 = 0,
                         Q = list(matrix(c(0, 0.1, 0.1, 0), 
                                         ncol = 2, nrow = 2)),
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
  
  # check that rates are numeric
  if (!is.numeric(c(lambda, mu))) {
    stop("lambda and mu must be numeric")
  } else if (any(!(c(length(lambda), length(mu)) %in% c(1, nStates)))) {
    stop("lambda and mu must be of length either one or equal to nStates")
  } else {
    # if everything is good, we set flags on whether each are TD
    tdLambda <- length(lambda) == nStates
    tdMu <- length(mu) == nStates
  }
  # 
  # # if nFocus is just one number, make it a vector
  # if (length(nFocus) == 1) {
  #   nFocus <- rep(nFocus, 2)
  # }
  # 
  ### HAVE TO write error tests, including for type + tMax and N
  
  # # error checks for each trait
  # for (i in 1:nTraits) {
  #   # take the Q for this trait
  #   Qn <- Q[[i]]
  #   
  #   # check that Qn is square
  #   if (ncol(Qn) != nrow(Qn)) {
  #     stop("Q must be a square matrix")
  #   }
  #   
  #   # check that Qn has row number (and therefore column number) 
  #   # equal to the number of traits
  #   if (nrow(Qn) != nStates) {
  #     stop("Q must have transition rates for all trait combinations")
  #   }
  #   
  #   # set the diagonal of Q to 0
  #   diag(Qn) <- rep(0, nrow(Qn))
  # }
  
  # do we condition on the number of species?
  condN <- condition == "number"
  
  # if conditioning on number, need to set some 
  # tMax for the trait evolution
  if (condN) {
    # thrice the average time it will take to get to 10*N species
    # under the slowest diversification rate parameters
    trTMax <- max(log(10*N) / (lambda - mu))
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
    traits <- list(traits.musse(tMax = min(tMax, trTMax), tStart = 0, 
                                nTraits = nTraits, nStates = nStates,
                                X0 = X0, Q = Q))
    
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
                          ifelse(sum(lambda) > 0,
                                 rexp.traits(1, lambda, 
                                             traitsSp,
                                             tNow, tMax), Inf),
                          ifelse(lambda > 0, 
                                 rexp(1, lambda), Inf))
      
      waitTimeE <- ifelse(tdMu,
                          ifelse(sum(mu) > 0,
                                 rexp.traits(1, mu, 
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
          traits.musse(tMax = min(tMax, max(trTMax, tNow + 100)), 
                       tStart = tNow, 
                       nTraits = nTraits, nStates = nStates,
                       X0 = XPar, Q = Q)
        
        # get a new speciation waiting time
        waitTimeS <- ifelse(tdLambda, 
                            ifelse(sum(lambda) > 0,
                                   rexp.traits(1, lambda, 
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
traits.musse <- function(tMax, tStart = 0, nTraits = 1, nStates = 2, X0 = 0,
                         Q = list(matrix(c(0, 0.1, 0.1, 0), 2, 2))) {
  # create a return list
  traits <- vector(mode = "list", length = nTraits)
  
  # for each trait
  for (i in 1:nTraits) {
    # get number of states and initial states
    nStatesI <- ifelse(length(nStates) > 1, nStates[i], nStates)
    X0I <- ifelse(length(X0) > 1, X0[i], X0)
    
    # get Q matrix
    if (length(Q) > 1) QI <- Q[[i]] else QI <- Q[[1]]
    
    # append traits data frame to traits
    traits[[i]] <- trait.musse(tMax, tStart, nStatesI, X0I, QI)
  }
  
  return(traits)
}

trait.musse <- function(tMax, tStart = 0, nStates = 2, X0 = 0,
                         Q = matrix(c(0, 0.1, 0.1, 0), 2, 2)) {
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
  
  # make sure X0 is a possible state
  if (X0 > nStates - 1) {
    stop("X0 must be an achievable state (i.e. in c(0, nStates - 1))")
  }
  
  # create states vector from number
  states <- 0:(nStates - 1)
  
  if (length(X0) != 1 || length(tStart) != 1) {
    print(X0)
    print(tStart)
  }

  # create traits data frame
  traits <- data.frame(value = X0, min = tStart, max = NA)
  
  # start a time counter
  tNow <- tStart
  
  # and a shifts counter
  shifts <- 0
  
  # make diagonals of Q 0 
  diag(Q) <- 0
  
  # while we have not reached the end
  while (tNow < tMax) {
    # current state
    curState <- traits$value[shifts + 1]
    
    # get the total rate of transition from the current state
    rTotal <- sum(Q[curState + 1, ])
    
    # get the time until the next transition
    waitTime <- ifelse(rTotal > 0, rexp(1, rTotal), Inf)
    
    # increase time
    tNow <- tNow + waitTime
    
    # add max to traits
    traits$max[shifts + 1] <- min(tNow, tMax)
    
    # break if needed
    if (tNow >= tMax) break
    
    # increase shifts counter
    shifts <- shifts + 1
    
    # sample to find target state
    newState <- sample(states, 1, prob = Q[curState + 1, ])
    
    # add it to traits data frame
    traits[shifts + 1, ] <- c(newState, tNow, NA)
  }
  
  return(traits)
}
