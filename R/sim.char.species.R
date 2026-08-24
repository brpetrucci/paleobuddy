#' Discrete character simulation for one species
#' 
#' Simulates a discrete character following a continuous-time Markov Chain 
#' along the life of one species. In practice, we start with a given state,
#' then draw a waiting time to a state transition based on the total rate of
#' transitioning out of this state, and draw the new state based on the
#' relative transition rates from this state to each other state. We continue
#' this process until we have reached the total simulation time. This function
#' is intended to be called by other functions in the package (such as
#' \code{sim.chars} and \code{bd.sim.traits}), but can be used for simple proof
#' of concepts for one species simulations.
#' 
#' @param tMax Maximum simulation time. When called by other functions in the
#' package, this will likely be either the extinction time of this species, or
#' the maximum simulation time of the entire birth-death process. This must be
#' strictly greater than \code{tStart}.
#' 
#' @param tStart Starting time for the simulation. When called by other
#' functions in the package, this will likely be either \code{0} (the start of
#' the birth-death process) or the speciation time of this species. This must
#' be strictly less than \code{tMax}. Default is \code{0}.
#' 
#' @param nStates Number of possible states for categorical trait. The range
#' of values will be assumed to be \code{(0, nStates - 1)}. Default is \code{2}.
#' 
#' @param X0 Initial trait value for original species. Must be within 
#' \code{(0, nStates - 1)}. Default is \code{0}.
#' 
#' @param Q Transition rate matrix for continuous-time trait evolution. For
#' different states \code{i} and \code{j}, the rate at which a species at
#' \code{i} transitions to \code{j} is \code{Q[i + 1, j + 1]}. Default is a
#' 2x2 matrix with symmetrical rates equal to \code{0.1}. Note that while the
#' standard is to set the diagonal values such that rows sum to \code{0}, the
#' diagonal values are never used in our algorithm, and we generally set them
#' to \code{0}.
#' 
#' @return A data frame with columns
#' 
#' \describe{
#' \item{\code{value}}{The state value.}
#' \item{\code{min}}{The time at which the species transitioned to this state
#' value, or \code{0} if it is the state at the start of the simulation.}
#' \item{\code{max}}{The time at which the species transitioned to this state
#' value, or \code{tMax} if it is the state at the end of the simulation.}}
#' 
#' Note that the time values in the data frame are considering the time
#' since the beginning of the simulation. These values are inverted to represent
#' time before the present in the \code{bd.sim.traits} function.
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' @references 
#' 
#' Lewis, P.O.. 2001. A likelihood approach to estimating phylogeny from 
#' discrete morphological character data. Systematic Biology. 50(6):913-925.
#' 
#' @examples
#' 
#' # set seed
#' set.seed(1)
#' 
#' ###
#' # simple example using all default parameters
#' char <- sim.char.species(10)
#' # the species starts at state 0, stays there for most of its lifetime,
#' # and quickly transitions to state 1 around 7.5my in, before transitioning
#' # back to 0 shortly after, and staying there until the end
#' 
#' ###
#' # we can of course have a different number of states
#' 
#' # let's simulate a 4-state character
#' nStates <- 4
#' 
#' # we need to also set the transition matrix
#' Q <- matrix(0.1, 4, 4)
#' diag(Q) <- 0
#' # keeping it simple with all equal symmetrical transition rates,
#' # equivalent to the default for 2 states
#' 
#' # simulate character
#' char <- sim.char.species(10, nStates = nStates, Q = Q)
#' # the character just jumps to 2 early in the simulation and stays there
#' 
#' # what if we want a more labile character?
#' # adjust Q
#' Q <- matrix(1, 4, 4)
#' diag(Q) <- 0
#' 
#' # simulate character
#' char <- sim.char.species(10, nStates = nStates, Q = Q)
#' # now we see many more transitions
#' 
#' ###
#' # we can test the different effects of X0 and Q on the final state
#' 
#' \dontrun{
#' # create a data frame to hold results
#' res <- data.frame(matrix(nrow = 1000, ncol = 3))
#' 
#' # get values of asymmetrical q and X0
#' qs <- rep(seq(0.1, 0.5, 0.1), 2)
#' X0s <- c(rep(0, 5), rep(1, 5))
#' 
#' # iterate through values
#' for (i in seq_along(qs)) {
#'   # get q and X0
#'   q <- qs[i]
#'   X0 <- X0s[i]
#'   
#'   # make Q matrix
#'   Q <- matrix(c(0, q, 0.1, 0), 2, 2)
#'   
#'   # run 100 reps
#'   for (j in 1:100) {
#'     # get character
#'     char <- sim.char.species(50, X0 = X0, Q = Q)
#'   
#'     # add final character to results
#'     res[(i - 1) * 100 + j, ] <- c(q, X0, tail(char$value))
#'   }
#' }
#' 
#' # name results
#' colnames(res) <- c("q", "X0", "final_state")
#' 
#' # get the mean per q and X0 value
#' aggregate(res[[3]] ~ res[[1]] + res[[2]], data = res, FUN = mean)
#' # as expected, the proportion of simulations ending in 1 increases with
#' # q, but the effect is more pronounced when we start in that state
#' }
#' 
#' @name sim.char.species
#' @rdname sim.char.species
#' @export
#' 

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
  
  # make sure nStates is at least 2
  if (nStates < 2) {
    stop("nStates must be at least 2.")
  }
  
  # make sure X0 is a possible state
  if (!is.numeric(X0) || X0 < 0 || X0 > nStates - 1) {
    stop("X0 must be an achievable state (i.e. in 0:(nStates - 1)")
  }
  
  # make sure Q is a matrix
  if (!is.matrix(Q) && !is.data.frame(Q)) {
    stop("Q must be a matrix or data frame")
  }
  
  # make sure Q has the right number of rows and columns
  if (!(nrow(Q) == nStates && ncol(Q) == nStates)) {
    stop("Q must have transition rates for all trait combinations")
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
