#' Simulate trait evolution through the Ornstein-Uhlenbeck process
#'
#' Generates trait evolution simulations for a number of traits under the
#' Ornstein-Uhlenbeck process. This is a generalization of Brownian Motion that
#' has a target value. Each trait starts with a value \code{X0}, and for each time
#' interval receives an increase defined by its difference from the target value
#' and a normally distributed variate. See references for the specific mathematics
#' of the process.
#'
#' @inheritParams traits.bm
#' 
#' @param sigma2 The variance of the normal distribution associated with the
#' increase in trait value for each time step. This defines the variance of 
#' trait values, which will be 
#' \code{sigma2 / (2*theta) * (1 - exp(-2 * theta * t))} for an interval
#' \code{t}.  The default is \code{1}. This can optionally be a vector with size 
#' equal to \code{nTraits}, in which case each trait will have a different 
#' variance.
#' 
#' @param theta The intensity of drift towards the target value \code{mu}. 
#' Increasing \code{theta} leads to a higher probability of trait values being 
#' close to the target value, which decreases variance (see \code{sigma2}).
#' The default is \code{0.5}. This can optionally be a vector with size 
#' equal to \code{nTraits}, in which case each trait will have a different 
#' drift intensity.
#' 
#' @param mu The target value of the process. The higher the difference between
#' a trait value and its target, the more likely the next time step is to bringing
#' the value closer to the target. Default is \code{0}. This can optionally be a 
#' vector with size equal to \code{nTraits}, in which case each trait will have a 
#' different target value.
#' 
#' @param X0 Given the complexity of the OU process, we allow for an argument
#' defining the starting value for the trait. The default is \code{0}, so that
#' choosing all default options leads to a process starting at the same value as
#' its target. This can optionally be a vector with size equal to \code{nTraits},
#' in which case each trait will have a different starting value.
#'
#' @author Bruno do Rosario Petrucci.
#'
#' @return A named list of functions, where \code{traits.ou$traitN(t)}
#' corresponds to the trait value of trait \code{N} at time \code{t}.
#' 
#' @references 
#' 
#' Uhlenbeck, G. E., Ornstein, L. S. (1930). On the theory of Brownian Motion. 
#' Phys. Rev. 36 (5): 823–841.
#' 
#' Felsenstein, J. (1988). Phylogenies and Quantitative Characters. Annual Review
#' of Ecology and Systematics
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
#' # calculate OU function
#' ouFunc <- traits.ou(tMax)$trait1
#' 
#' # plot it
#' plot(times, ouFunc(times), type = 'l', main = "Trait values for trait1",
#'      xlab = "Time (my)", ylab = "Trait value")
#' 
#' ###
#' # a lot of customization is possible - let us try a 
#' # different sigma2 and theta
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
#' 
#' # and a higher theta
#' theta <- 1
#' 
#' # calculate OU function
#' ouFunc <- traits.ou(tMax, sigma2 = sigma2, theta = theta)$trait1
#' 
#' # plot it
#' plot(times, ouFunc(times), type = 'l', main = "Trait values for trait1",
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
#' sigma2 <- c(5, 0.5, 1)
#' 
#' # same for theta
#' theta <- c(1, 0.5, 0.25)
#' 
#' # and mu, but say we want to keep it a constant
#' mu <- 0.5
#' 
#' # finally, define the starting value for each trait
#' X0 <- c(0, 0.3, 1)
#' 
#' # set a seed so we can control ylim and legend position
#' set.seed(1)
#' 
#' # calculate BM functions
#' ouFuncs <- traits.ou(tMax, nTraits = 3, sigma2 = sigma2, theta = theta,
#'                      mu = mu, X0 = X0)
#' 
#' # plot them
#' plot(times, ouFuncs[[1]](times), type = 'l',
#'      main = "Trait values for traits with mu = 0.5",
#'      xlab = "Time (my)", ylab = "Trait value", ylim = c(-3, 3))
#' lines(times, ouFuncs[[2]](times), col = 'RED')
#' lines(times, ouFuncs[[3]](times), col = 'BLUE')
#' legend(x = 1, y = -1.5, legend = c("5, 1, 0", "0.5, 0.5, 0.3", "1, 0.25, 1"), 
#'       col = c('BLACK', 'RED', 'BLUE'), lty = c(1,1,1))
#'
#' @name traits.ou
#' @rdname traits.ou
#' @export

traits.ou <- function(tMax, nTraits = 1, sigma2 = 1, tStart = 0,
                      theta = 0.5, mu = 0, X0 = 0, nPoints = 100) {
  # make sure tMax > tStart
  if (tStart >= tMax) {
    stop("tMax must be greater than tStart")
  }
  
  # nTraits must be positive
  if (nTraits < 1) {
    stop("nTraits must be a positive number")
  }
  
  # make a list of parameters
  pars <- list(sigma2, theta, mu, X0)
  
  # sigma2 is either a constant, or a vector with size nTraits
  if (length(sigma2) == 1) {
    if (sigma2 <= 0) {
      stop("sigma2 must be a positive number")
    }
    
    # make it a vector if needed
    if (max(lengths(pars)) > 1) {
      sigma2 <- rep(sigma2, nTraits)
    }
  } else {
    if (length(sigma2) != nTraits) {
      stop("sigma2 must have length equal to nTraits")
    } else if (any(sigma2 <= 0)) {
      stop("All variance values must be positive")
    }
  }
  
  # theta is also a constant, or a vector with size nTraits
  if (length(theta) == 1) {
    if (theta <= 0) {
      stop("theta must be a positive number")
    }
    
    # make it a vector if needed
    if (max(lengths(pars)) > 1) {
      theta <- rep(theta, nTraits)
    }
  } else {
    if (length(theta) != nTraits) {
      stop("theta must have length equal to nTraits")
    } else if (any(theta <= 0)) {
      stop("All theta values must be positive")
    }
  }
  
  # finally, X0 and mu can also be constants or vectors with size nTraits
  if ((length(X0) != 1) && (length(X0) != nTraits)) {
    stop("X0 must be either a constant or a vector with size nTraits")
  }
  
  # make it a vector if needed
  if (length(X0) < max(lengths(pars))) {
    X0 <- rep(X0, nTraits)
  }
  
  if ((length(mu) != 1) && (length(mu) != nTraits)) {
    stop("mu must be either a constant or a vector with size nTraits")
  }
  
  # make it a vector if needed
  if (length(mu) < max(lengths(pars))) {
    mu <- rep(mu, nTraits)
  }
  
  # create a results list
  res <- list()

  # for each trait,
  for (i in 1:nTraits) {
    # calculate time interval
    dt <- (tMax - tStart) / (nPoints - 1)
    
    # calculate vector of times
    times <- seq(tStart, tMax, dt)
    
    # calculate nPoints - 1 randomly distributed normal
    # variates with variance (tMax - tStart) / nPoints
    dw <- rnorm(nPoints - 1, mean = 0, 
                   sd = sqrt(dt))
    
    # create the OU vector with start at X0[i]
    ou <- c(X0[i])
    
    # calculate the OU process for the remaining time points
    for (j in 2:length(times)) {
      ou[j] <- ou[j - 1] + theta[i]*(mu[i] - ou[j - 1])*dt + 
        sqrt(sigma2[i])*dw[j - 1]
    }

    # make it a function using linear interpolation
    ouFunc <- approxfun(times, ou)
    
    # append it to the results list
    res[[paste0("trait", i)]] <- ouFunc
  }

  return(res)
}
