#' Simulate trait evolution with a discrete number of states
#'
#' Generates trait evolution simulations for a number of traits under a continuous
#' Markov Chain process with discrete states.
#' 
#' @inheritParams traits.bm
#' 
#' @param states A numerical vector specifying possible trait values for the 
#' trait. Default is \code{c(0, 1)} for a binary trait. Note that it must not
#' contain repeated values.
#' 
#' @param X0 The initial trait value. Must be in \code{states}. Default is 
#' \code{0}, usually interpreted as the absence of a binary trait.
#' 
#' @param tStart The starting time of the simulation. Default is \code{0},
#' but this argument is required to accurately return a function based on the 
#' trait evolution of a trait for a species born after the beginning of a 
#' birth-death simulation.
#' 
#' @param Q The transition rate matrix. \code{Q[i, j]} refers to the transition
#' rate between \code{states[i]} and \code{states[j]}. These rates correspond to
#' the relative probabilities of going from the original state to the next. So
#' given that a state transition occurs from \code{states[i]}, the probability
#' that the transition result will be \code{states[j]} is
#' \code{Q[i, j] / sum(Q[i, ])}. Note that to conform to most interpretations
#' of rate matrices in the context of continuous Markov Processes, 
#' \code{Q[i, i] = 0} for all i (set inside the function).
#' 
#' @author Bruno do Rosario Petrucci.
#'
#' @return A named list of functions, where \code{traits.states$traitN(t)}
#' corresponds to the trait value of trait \code{N} at time \code{t}.
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
#' ###
#' # we can start with a simple trait with all default values
#' 
#' # set a maximum simulation time
#' tMax <- 10
#' 
#' # a starting time
#' tStart <- 0
#' 
#' # create a vector of times for plotting
#' times <- seq(tStart, tMax, 0.01)
#' 
#' # calculate trait function
#' traitFunc <- traits.states(tMax)$trait1
#' 
#' # plot it
#' plot(times, traitFunc(times), type = 'l', main = "Trait values for trait1",
#'      xlab = "Time (my)", ylab = "Trait value")
#'
#' ###
#' # we can have any number of states
#' 
#' # set a maximum simulation time
#' tMax <- 10
#' 
#' # a starting time
#' tStart <- 0
#' 
#' # create a vector of times for plotting
#' times <- seq(tStart, tMax, 0.01)
#' 
#' # set a states vector for a trait with 3 values
#' states <- c(0, 1, 2)
#' 
#' # and a rate matrix
#' Q <- matrix(c(0, 0.5, 0.25, 0.1, 0, 0.2, 0.2, 0.3, 0), nrow = 3, ncol = 3)
#' # notice the 0 diagonal
#' 
#' # calculate trait function
#' traitFunc <- traits.states(tMax, tStart = tStart, nTraits = 1, states = states,
#'                            X0 = 0, Q = Q)$trait1
#' 
#' # plot it
#' plot(times, traitFunc(times), type = 'l', main = "Trait values for trait1",
#'      xlab = "Time (my)", ylab = "Trait value")
#' 
#' ###     
#' # we can also have multiple traits
#' 
#' # set a maximum simulation time
#' tMax <- 10
#' 
#' # a starting time
#' tStart <- 0
#' 
#' # create a vector of times for plotting
#' times <- seq(tStart, tMax, 0.01)
#' 
#' # set a states vector
#' # notice states can be any number
#' states <- c(2.3, 7, pi)
#' 
#' # set a starting value for each trait
#' X0 <- c(2.3, pi, pi, 7)
#' 
#' # finally, a rate matrix
#' Q <- matrix(c(0, 0.2, 0.15, 0.05, 0, 0.1, 0.25, 0.1, 0), nrow = 3, ncol = 3)
#' # notice the 0 diagonal
#' 
#' # set seed to facilitate plotting
#' set.seed(1)
#' 
#' # calculate trait function
#' traitFuncs <- traits.states(tMax, tStart = tStart, nTraits = 4, states = states,
#'                             X0 = X0, Q = Q)
#' 
#' # plot them
#' plot(times, traitFuncs[[1]](times), type = 'l',
#'      main = "Traits with differing starting values",
#'      xlab = "Time (my)", ylab = "Trait value")
#' lines(times, traitFuncs[[2]](times), col = "RED")
#' lines(times, traitFuncs[[3]](times), col = "BLUE")
#' lines(times, traitFuncs[[4]](times), col = "GREEN")
#' legend(x = 1, y = 5, legend = c("2.3", "pi", "pi", "7"),
#'        col = c("BLACK", "RED", "BLUE", "GREEN"),
#'        lty = c(1, 1, 1, 1))
#' 
#' @name traits.states
#' @rdname traits.states
#' @export
#'

traits.states <- function(tMax, tStart = 0, nTraits = 1, 
                          states = c(0, 1), X0 = 0,
                          Q = matrix(c(0, 1/2, 1/2, 0), 2, 2)) {
  # make sure tMax > tStart
  if (tStart >= tMax) {
    stop("tMax must be greater than tStart")
  }
  
  # nTraits must be positive
  if (nTraits < 1) {
    stop("nTraits must be a positive number")
  }
  
  # X0 is either a constant, or a vector with size nTraits
  if (length(X0) != 1) {
    if (length(X0) != nTraits) {
      stop("X0 have length equal to nTraits")
    }
  }
  
  # if it is of length 1, make it a vector
  # if nTraits is not 1
  else {
    X0 <- rep(X0, nTraits)
  }
  
  # if any of X0 is not in states, problem
  if (any(!(X0 %in% states))) {
    stop("X0 must contain only values in states")
  }
  
  # states should not be repeated
  if (length(unique(states)) != length(states)) {
    stop("states must not contain repeated values")
  }
  
  # check that Q has row and column numbers
  # equal to the number of states
  if (nrow(Q) != ncol(Q)) {
    stop("Q must be a square matrix")
  } else if (nrow(Q) != length(states)) {
    stop("Q must have a number  of rows and columns equal to the 
         length of states")
  }
  
  # set the diagonal of Q to 0
  diag(Q) <- rep(0, nrow(Q))
  
  # create a results list
  res <- list()
  
  # for each trait,
  for (i in 1:nTraits) {
    # create trait vector
    traits <- c(X0[i])
    
    # start a time counter
    tNow <- tStart
    
    # create a shift times vector
    shifts <- c()
    
    # and a shifts counter
    nShifts <- 0
    
    # while we have not reached the end
    while (tNow < tMax) {
      # current state
      curState <- traits[nShifts + 1]
      
      # current state index
      stateInd <- which(states == curState)
      
      # get the total rate of transition from the current state
      rTotal <- sum(Q[stateInd, ])
      
      # get the time until the next transition
      waitTime <- rexp(1, rTotal)
      
      # increase time and break if needed
      tNow <- tNow + waitTime
      if (tNow >= tMax) break
      
      # append shift time to shifts vector
      shifts <- c(shifts, tNow)
      
      # sample to find target state
      newState <- sample(states, 1, prob = Q[stateInd, ])
      
      # add it to trait vector
      traits <- c(traits, newState)
      
      # increase nShifts
      nShifts <- nShifts + 1
    }
    
    # if there were no shifts, create a constant function
    if (nShifts == 0) {
      stateFunc <- Vectorize(function(t) {
        return(traits[1])
      })
    }
    
    # otherwise, make it a function using stepfun
    else {    
      stateFunc <- stepfun(shifts, traits)
    }
    
    # append it to the results list
    res[[paste0("trait", i)]] <- stateFunc
  }
  
  return(res)
}
