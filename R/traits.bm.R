#' Simulate trait evolution through Brownian Motion
#'
#' Generates trait evolution simulations for a number of traits under Brownian
#' Motion. Each trait starts with a value of \code{Z0} and, for each time point,
#' changes by an amount calculated stochastically by a normally distributed
#' variable with given variance.
#'
#' @param tMax The only required parameter. Ending time of the simulation, usually
#' inherited from a birth-death or sampling function.
#' 
#' @param nTraits The number of traits to be simulated. Traits are assumed to be
#' independent from one another. By default, the function simulates one trait.
#' 
#' @param sigma2 The variance of the normal distribution. This defines the 
#' expected variance of the trait values, which will be \code{sigma2 * t}, for
#' each time interval \code{t} to be considered. The default is \code{1}, so that
#' variance in trait values by time \code{nStart + t} will be \code{t}. This can 
#' optionally be a vector with size equal to \code{nTraits}, in which case each 
#' trait will have a different variance.
#' 
#' @param X0 The initial trait value. For a BM process, one can simply run the
#' process starting from \code{0} and sum \code{X0} at the end, but we add this
#' for simplicity.
#' 
#' @param tStart The starting time of the simulation. Standard Brownian Motion 
#' starts at \code{0}, so we set the default value as such, but this argument is
#' required to accurately return a function based on the trait evolution of a
#' trait for a species born after the beginning of a birth-death simulation. Note
#' that the trait value for \code{t = 0} is always \code{0}.
#' 
#' @param bounds A numeric vector with length two containing minimum and maximum
#' trait values. If trait values fall under \code{bounds[1]} or over 
#' \code{bounds[2]}, they default to the corresponding bound value. By default not
#' supplied, so that any trait value is accepted as a return. Note that setting 
#' bounds will alter the expectations of the process.
#' 
#' @param nPoints The total number of points in the discretization of the process.
#' While Brownian Motion is a continuous process, we must simulate it discretely.
#' As such, the user can provide a number of time points, by default set to 
#' \code{100}, for which the exact trait values will be calculated. Remaining
#' values of time will have values predicted by linear interpolation (see 
#' \code{?approxfun}).
#'
#' @author Bruno do Rosario Petrucci.
#'
#' @return A named list of functions, where \code{traits.bm$traitN(t)}
#' corresponds to the trait value of trait \code{N} at time \code{t}.
#' 
#' @references 
#' 
#' Felsenstein, J. (1973). Maximum-Likelihood estimation of evolutionary trees
#' from continuous chracters. American Journal of Human Genetics. 25(5):471-492.
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
#' # and a number of time points
#' nPoints <- 100
#' 
#' # which lets us find a vector of time points
#' times <- seq(tStart, tMax, (tMax - tStart) / nPoints)
#' 
#' # calculate BM function
#' bmFunc <- traits.bm(tMax)$trait1
#' 
#' # plot it
#' plot(times, bmFunc(times), type = 'l', main = "Trait values for trait1",
#'      xlab = "Time (my)", ylab = "Trait value")
#' 
#' ###
#' # a lot of customization is possible - let us try a different sigma2
#' 
#' # set a maximum simulation time
#' tMax <- 10
#' 
#' # a starting time
#' tStart <- 0
#' 
#' # and a number of time points
#' nPoints <- 100
#' 
#' # which lets us find a vector of time points
#' times <- seq(tStart, tMax, (tMax - tStart) / nPoints)
#' 
#' # try a lower value for sigma2
#' sigma2 <- 0.5
#' # note this would mean the variance for trait values at 
#' # t = tMax is 5
#' 
#' # calculate BM function
#' bmFunc <- traits.bm(tMax, sigma2 = sigma2)$trait1
#' 
#' # plot it
#' plot(times, bmFunc(times), type = 'l', main = "Trait values for trait1",
#'      xlab = "Time (my)", ylab = "Trait value")
#' 
#' ###
#' # we can calculate any number of traits at once
#' 
#' # set a maximum simulation time
#' tMax <- 10
#' 
#' # a starting time
#' tStart <- 0
#' 
#' # and a number of time points
#' nPoints <- 100
#' 
#' # which lets us find a vector of time points
#' times <- seq(tStart, tMax, (tMax - tStart) / nPoints)
#' 
#' # if sigma2 is a vector, each trait will have 
#' # a different variance, in the corresponding order
#' sigma2 <- c(0.5, 0.1, 0.05)
#' 
#' # we can also have different starting points
#' X0 <- c(0, 0.2, -0.1)
#' 
#' # set a seed so we can control ylim and legend position
#' set.seed(1) 
#' 
#' # calculate BM functions
#' bmFuncs <- traits.bm(tMax, nTraits = 3, sigma2 = sigma2, X0 = X0)
#' 
#' # plot them
#' plot(times, bmFuncs[[1]](times), type = 'l', 
#'      main = "Trait values for traits with sigma2 = 0.5, 0.1 and 0.05",
#'      xlab = "Time (my)", ylab = "Trait value", ylim = c(-2.5, 3))
#' lines(times, bmFuncs[[2]](times), col = 'RED')
#' lines(times, bmFuncs[[3]](times), col = 'BLUE')
#' legend(x = 1, y = -1, legend = c(0.5, 0.1, 0.05),
#'        col = c('BLACK', 'RED', 'BLUE'), lty = c(1, 1, 1))
#'
#' @name traits.bm
#' @rdname traits.bm
#' @export

traits.bm <- function(tMax, tStart = 0, nTraits = 1, 
                      sigma2 = 1, X0 = 0, bounds = NULL, nPoints = 100) {
  # make sure tMax > tStart
  if (tStart >= tMax) {
    stop("tMax must be greater than tStart")
  }
  
  # nTraits must be positive
  if (nTraits < 1) {
    stop("nTraits must be a positive number")
  }
  
  # sigma2 is either a constant, or a vector with size nTraits
  if (length(sigma2) == 1) {
    if (sigma2 <= 0) {
      stop("sigma2 must be a positive number")
    }
    
    # if nTraits > 1, make it a vector
    sigma2 <- rep(sigma2, nTraits)
  } else {
    if (length(sigma2) != nTraits) {
      stop("sigma2 must have length equal to nTraits")
    } else if (any(sigma2 <= 0)) {
      stop("All variance values must be positive")
    }
  }
  
  # same for X0
  if (length(X0) != 1) {
    if (length(X0) != nTraits) {
      stop("X0 must have length equal to nTraits")
    }
  }
  
  # if it is of length 1, make it a vector
  # if nTraits is not 1
  else {
    X0 <- rep(X0, nTraits)
  }
  
  # create a results list
  res <- list()
  
  # create the list for auxiliary functions
  aux <- list()
  # necessary for extending them to infinity
  
  # for each trait,
  for (i in 1:nTraits) {
    # calculate vector of times
    times <- seq(tStart, tMax, (tMax - tStart) / (nPoints - 1))
    
    # calculate nPoints - 1 randomly distributed normal
    # variates with variance sigma2
    jumps <- rnorm(nPoints - 1, mean = 0, 
                   sd = sqrt(sigma2[i] * (tMax - tStart) / (nPoints - 1)))
    
    # calculate the bm vector being X0 and the cumulative sum of jumps
    bm <- c(X0[i], cumsum(jumps))
    
    # bound it, if bounds exists
    if (!is.null(bounds)) {
      bm <- ifelse(bm < bounds[1], bounds[1], bm)
      bm <- ifelse(bm > bounds[2], bounds[2], bm)
    }

    # make it a function using linear interpolation
    bmFunc <- approxfun(times, bm, rule = 2)
  
    # append it to the results list
    res[[paste0("trait", i)]] <- bmFunc
  }

  return(res)
}
