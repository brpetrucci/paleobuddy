#' Constant and time-dependent rate species sampling
#' 
#' Generates a vector of occurrence times for species in a simulation using a
#' Poisson process. Allows for the Poisson rate to be (1) a constant or (2) a 
#' function of time. For fossil sampling dependent on species age in addition 
#' to absolute time, see \code{sample.clade}.
#' 
#' Note that while the Poisson process occurs in forward time, we return (both 
#' in birth-death functions and here) results in backwards time, so that time 
#' is inverted using \code{tMax} both at the beginning and end of 
#' \code{sample.time}.
#'
#' @param sim A \code{sim} object, containing extinction times, speciation 
#' times, parent, and status information (extant or extinct) for each species
#' in the simulation. See \code{?sim}.
#' 
#' @param rho Sampling rate (events per species per million years) over time. 
#' It can be a \code{numeric} describing a constant rate or a \code{function(t)} 
#' describing the variation in sampling over time. For more flexibility on 
#' sampling, see \code{make.rate} to create more complex rates. Note that 
#' \code{rho} should always be greater than or equal to zero.
#' 
#' @param tMax The maximum simulation time, used by \code{rexp.var}. A sampling
#' time greater than \code{tMax} would mean the occurrence is sampled after the
#' present, so for consistency we require this argument. This is also required
#' to ensure time follows the correct direction both in the Poisson process and
#' in the output.
#'
#' @param S A vector of species numbers to be sampled. The default is all 
#' species in \code{sim}. Species not included in \code{S} will not be sampled 
#' by the function.
#'
#' @return A list of vectors of occurrence times for each species in \code{S}.
#'
#' @author Bruno do Rosario Petrucci and Matheus Januario.
#' 
#' @name sample.time
#' @rdname sample.time
#' @keywords Internal
#' 

sample.time <- function(sim, rho, tMax, S = NULL) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # if S is NULL, make it all species in the simulation
  if (is.null(S)) {
    S <- 1:length(sim$TE)
  }
  
  # if there are species in S not in the simulation, error
  if (!(all(S %in% 1:length(sim$TE)))) {
    stop("S must contain only integers between 1 and the total 
         number of species in sim")
  }
  
  # get speciation and extinction times
  TE <- sim$TE
  TS <- sim$TS
  
  # make TE sensible
  TE[sim$EXTANT] <- 0
  
  # invert them since simulation goes from 0 to tMax
  TE <- tMax - TE
  TS <- tMax - TS
  
  # create return
  res <- list()
  
  # for each species
  for (s in S) {
    # start when the species was born, end when it died
    now <- TS[s]
    End <- TE[s]
    
    # initialize vector
    sampled <- c()
    
    # while we have not reached the time of extinction
    while (now < End) {
      # take the waiting time for sampling, using rexp.var()
      WaitTimeR <- ifelse(is.numeric(rho), 
                          ifelse(rho > 0, rexp(1, rho), Inf),
                          ifelse(rho(now) > 0, rexp.var(1, rho, now, tMax,
                                                        fast = TRUE), Inf))
      
      # advance in time
      now <- now + WaitTimeR
      
      # if sampling comes after extinction, we don't include this occurrence
      if (now > End) break
      
      # append to the vector
      sampled <- c(sampled, now)
    }
    
    # finally, we invert time so that it goes from tMax to 0
    sampled <- tMax - sampled
    
    # append to result
    res[[paste0("t", s)]] <- sampled
  }
  
  return(res)
}
