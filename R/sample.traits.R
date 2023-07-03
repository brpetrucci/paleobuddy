#' Trait-dependent fossil sampling
#' 
#' @inheritParams sample.clade
#'
#' @return A list of vectors of occurrence times for each species in \code{S}.
#'
#' @author Bruno do Rosario Petrucci and Matheus Januario.
#' 
#' @noRd
#' 

sample.traits <- function(sim, rho, tMax, traits, nFocus = 1, S = NULL) {
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
  
  # check that rho is numeric
  if (!is.numeric(rho)) {
    stop("rho must be a numeric vector for trait-dependent sampling")
  }
  
  # get number of states
  nStates <-length(unique(unlist(lapply(1:length(traits), function(x) 
    traits[[x]][[nFocus]]$value))))
  
  # make sure rho has the same length as number of states
  if (length(rho) != nStates) {
    stop("rho must have the same length as there are state values for traits.
         If rho is a constant, do not supply trait information as it would be
         useless")
  }
  
  # check that rho is nonzero
  if (sum(rho) == 0) {
    stop("rho must be nonzero for some trait value")
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
    
    # get trait data frame
    traitsSp <- traits[[s]][[nFocus]]
    
    # switch columns around to make it sensible with tMax
    minT <- traitsSp$min
    maxT <- traitsSp$max
    traitsSp$max <- tMax - maxT
    traitsSp$min <- tMax - minT
    names(traitsSp) <- c("values", "min", "max")
    
    # initialize vector
    sampled <- c()
    
    # while we have not reached the time of extinction
    while (now < End) {
      # take the waiting time for sampling, using rexp.musse
      WaitTimeR <- rexp.musse(1, rho, traitsSp, now, tMax, fast = TRUE)
      
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
