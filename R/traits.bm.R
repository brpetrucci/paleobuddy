#' Simulate trait evolution through Brownian Motion
#'
#' Generates trait evolution simulations for a number of traits under Brownian
#' Motion. Each trait starts with a value of \code{0} and, for each time point,
#' changes by an amount calculated stochastically by a normally distributed
#' variable with given variance. Note that if one wishes to start at a nonzero
#' value, one can simply sum the desired starting value to the result.
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
#' @param tStart The starting time of the simulation. Standard Brownian Motion 
#' starts at \code{0}, so we set the default value as such, but this argument is
#' required to accurately return a function based on the trait evolution of a
#' trait for a species born after the beginning of a birth-death simulation. Note
#' that the trait value for \code{t = 0} is always \code{0}.
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
#' # set a seed so we can control ylim and legend position
#' set.seed(1) 
#' 
#' # calculate BM functions
#' bmFuncs <- traits.bm(tMax, nTraits = 3, sigma2 = sigma2)
#' 
#' # plot them
#' plot(times, bmFuncs[[1]](times), type = 'l', 
#'      main = "Trait values for traits with sigma2 = 0.5, 0.1 and 0.05",
#'      xlab = "Time (my)", ylab = "Trait value", ylim = c(-10, 10))
#' lines(times, bmFuncs[[2]](times), col = 'RED')
#' lines(times, bmFuncs[[3]](times), col = 'BLUE')
#' legend(x = 1, y = -5, legend = c(0.5, 0.1, 0.05),
#'        col = c('BLACK', 'RED', 'BLUE'), lty = c(1, 1, 1))
#'
#' @name traits.bm
#' @rdname traits.bm
#' @export

traits.bm <- function(tMax, nTraits = 1, sigma2 = 1, tStart = 0, nPoints = 100) {
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
  } else {
    if (length(sigma2) != nTraits) {
      stop("sigma2 must have length equal to nTraits")
    } else if (any(sigma2 <= 0)) {
      stop("All variance values must be positive")
    }
  }
  
  # create a results list
  res <- list()

  # for each trait,
  for (i in 1:nTraits) {
    # calculate vector of times
    times <- seq(tStart, tMax, (tMax - tStart) / (nPoints - 1))
    
    # calculate nPoints - 1 randomly distributed normal
    # variates with variance sigma2
    jumps <- rnorm(nPoints - 1, mean = 0, 
                   sd = sqrt(sigma2[i] * (tMax - tStart) / (nPoints - 1)))
    
    # calculate the bm vector being 0 and the cumulative sum of jumps
    bm <- c(0, cumsum(jumps))

    # make it a function using linear interpolation
    bmFunc <- approxfun(times, bm)
    
    # append it to the results list
    res[[paste0("trait", i)]] <- bmFunc
  }

  return(res)
}