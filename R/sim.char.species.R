#' Character simulation for one species
#' 
#' Description
#' 
#' #at param...
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' #at references...
#' 
#' #at examples...
#' 
#' @name sim.char.species
#' @rdname sim.char.species
#' @export

sim.char.species <- function(tMax, tStart = 0, nStates = 2, X0 = 0,
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
  if (!is.numeric(X0) || X0 < 0 || X0 > nStates - 1) {
    stop("X0 must be an achievable state (i.e. in 0:(nStates * nHidden - 1)")
  }
  
  # create states vector from number
  states <- 0:(nStates - 1)
  
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
