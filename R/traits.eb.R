#' Simulate trait evolution through the Early Burst model
#' 
#' Generates trait evolution simulations following the Early Burst model,
#' an extension of Brownian Motion (see \code{traits.bm}) where the corresponding
#' variance follows an exponential decay (or increase). It is ideal for simulating
#' adaptive radiations (see References).
#' 
#' @inheritParams traits.ou
#' 
#' @param sigma2 The baseline variance of the normal distribution. The total
#' variance in trait values will be influenced by \code{b}, leading to a
#' variance of \code{sigma2 * (exp(b * t) - 1) / b}. The default is \code{1}. Can
#' also be a vector with size equal to \code{nTraits}.
#' 
#' @param b Exponential rate determining the speed of exponential decay of
#' trait value variance. Default is \code{-0.05}. Note it can have positive
#' values if one wishes to model trait variance increasing with time. Can also
#' be a vector with size equal to \code{nTraits}.
#' 
#' @param X0 Initial trait value. Default is \code{0}. Can also be a vector with
#' size equal to \code{nTraits}.
#' 
#' @author Bruno do Rosario Petrucci.
#' 
#' @return A named list of functions, where \code{traits.eb$traitN(t)}
#' corresponds to the trait value of trait \code{N} at time \code{t}.
#' 
#' @references
#' 
#' Harmon, L. J., Losos, J. B., Davies, T. D., Gillespie, R. G., Gittleman, J. L.,
#' Jennings, W. B., Hozak, K. H., McPeek, M. A., Moreno-Roark, F., Near, T. J.,
#' Purvis, A., Ricklefs, R. E., Schluter, D., Schulte II, J. A., Seehausen, O.,
#' Sidlauskas, B. L., Torres-Carvajal, O., Weir, J. T., Mooers, A. Ø. (2010).
#' Early Bursts of Body Size and Shape Evolution are Rare in Comparative Data.
#' Evolution 64-8: 2385-2396.
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
#' ebFunc <- traits.eb(tMax)$trait1
#' 
#' # plot it
#' plot(times, ebFunc(times), type = 'l', main = "Trait values for trait1",
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
#' 
#' # calculate BM function
#' ebFunc <- traits.eb(tMax, sigma2 = sigma2)$trait1
#' 
#' # plot it
#' plot(times, ebFunc(times), type = 'l', main = "Trait values for trait1",
#'      xlab = "Time (my)", ylab = "Trait value")
#' 
#' ###
#' # we can also try a different b value
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
#' # try a higher sigma2 this time
#' sigma2 <- 1.5
#' 
#' # and a faster exponential decay
#' b <- -0.15
#' 
#' # calculate BM function
#' ebFunc <- traits.eb(tMax, sigma2 = sigma2, b = b)$trait1
#' 
#' # plot it
#' plot(times, ebFunc(times), type = 'l', main = "Trait values for trait1",
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
#' # and different b, but let us just have all the same
#' b <- -0.025
#' 
#' # set a seed so we can control ylim and legend position
#' set.seed(1) 
#' 
#' # calculate BM functions
#' ebFuncs <- traits.eb(tMax, nTraits = 3, sigma2 = sigma2, b = b, X0 = X0)
#' 
#' # plot them
#' plot(times, ebFuncs[[1]](times), type = 'l', 
#'      main = "Trait values for traits with sigma2 = 0.5, 0.1 and 0.05",
#'      xlab = "Time (my)", ylab = "Trait value", ylim = c(-2.5, 3))
#' lines(times, ebFuncs[[2]](times), col = 'RED')
#' lines(times, ebFuncs[[3]](times), col = 'BLUE')
#' legend(x = 1, y = -1, legend = c(0.5, 0.1, 0.05),
#'        col = c('BLACK', 'RED', 'BLUE'), lty = c(1, 1, 1))
#'      
#' @name traits.eb
#' @rdname traits.eb
#' @export  


traits.eb <- function(tMax, tStart = 0, nTraits = 1, sigma2 = 1, 
                      b = -0.05, X0 = 0, nPoints = 100) {
  # make sure tMax > tStart
  if (tStart >= tMax) {
    stop("tMax must be greater than tStart")
  }
  
  # nTraits must be positive
  if (nTraits < 1) {
    stop("nTraits must be a positive number")
  }
  
  # make a list of parameters
  pars <- list(sigma2, b, X0)
  
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
  
  # finally, X0 and b can also be constants or vectors with size nTraits
  if ((length(X0) != 1) && (length(X0) != nTraits)) {
    stop("X0 must be either a constant or a vector with size nTraits")
  }
  
  # make it a vector if needed
  if (length(X0) < max(lengths(pars))) {
    X0 <- rep(X0, nTraits)
  }
  
  if ((length(b) != 1) && (length(b) != nTraits)) {
    stop("b must be either a constant or a vector with size nTraits")
  }
  
  # make it a vector if needed
  if (length(b) < max(lengths(pars))) {
    b <- rep(b, nTraits)
  }
  
  # create a results list
  res <- list()
  s <- c()

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
    
    # create the EB vector with start at X0[i]
    eb <- c(X0[i])
    
    # calculate the EB process for the remaining time points
    for (j in 2:length(times)) {
      eb[j] <- eb[j - 1] + sqrt(sigma2[i] * exp(b[i] * times[j]))*dw[j - 1]
    }

    s <- c(s, eb[length(eb)])
    
    # make it a function using linear interpolation
    ebFunc <- approxfun(times, eb)
    
    # append it to the results list
    res[[paste0("trait", i)]] <- ebFunc
  }
  mean(s)
  var(s)

  return(res)
}
