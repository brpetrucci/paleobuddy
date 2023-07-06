#' Trait dependent Exponential waiting times
#' 
#' Generates a waiting time following an exponential rates varying based on a
#' trait data frame. Output can be used as the waiting time to an extinction,
#' speciation, or fossil sampling event. Allows for further customization by 
#' restricting possible waiting time outputs with arguments for (1) current 
#' time, to consider only the rates starting at that time, and (2) maximum time
#' to bound the output and therefore allow for faster calculations if one only
#' cares about waiting times lower than a given amount. 
#'
#' @param n The number of waiting times to return. Usually 1, but we allow for
#' a higher \code{n} to be consistent with the \code{rexp} function.
#' 
#' @param rate The rate parameter for the exponential distribution. Can be a
#' constant, or a vector corresponding to the rate for each state of the trait.
#' 
#' @param traits The trait data frame detailing trait value for each interval 
#' of time. Its rows should be
#' 
#' \itemize{
#' 
#' \item \code{value} Trait values. These are assumed to be categorical, and 
#' therefore always in the range \code{(0, nStates - 1)}. So e.g. for two 
#' states, we have values 0 and 1. 
#' 
#' \item \code{min} The starting time for each trait value. The first row will
#' always have either 0 or the species' speciation time here. 
#' 
#' \item \code{max} The ending time for each trait value. Note that the 
#' \code{max} of a row is the \code{min} of the next row. After the max time of
#' the last row, the species is assumed to maintain the same trait value it had
#' last.}
#'
#' @param now The current time. Needed if one wants to consider only the 
#' interval starting at the current time time-varying nate. Note this does 
#' means the waiting time is \code{>= now}, so we also subtract \code{now} from
#' the result before returning. The default is \code{0}.
#'
#' @param tMax The simulation ending time. If the waiting time would be too
#' high, we return \code{tMax + 0.01} to signify the event never happens, if
#' \code{fast == TRUE}. Otherwise we return the true waiting time. By default,
#' \code{tMax} will be \code{Inf}, but if \code{FAST == TRUE} one must supply
#' a finite value.
#' 
#' @param fast If set to \code{FALSE}, waiting times larger than \code{tMax} will
#' not be thrown away. This argument is needed so one can test the function
#' without bias.
#'
#' @return A vector of waiting times following the exponential distribution 
#' with the given parameters.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @examples
#' 
#' ###
#' # this function, similarly to rexp.var, can also work as rexp
#' 
#' # rate
#' rate <- 0.1
#' 
#' # traits
#' traits <- data.frame(value = 0, min = 0, max = 100)
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the waiting time
#' t <- rexp.musse(1, rate, traits)
#' t
#' # this is the same as t <- rexp(1, rate)
#' 
#' ###
#' # we can also use it as intended, though, and use trait-dependent rates
#' 
#' # rates
#' rate <- c(0.1, 0.2)
#' 
#' # traits
#' traits <- data.frame(value = c(0, 1, 0), min = c(0, 10, 20), 
#'                      max = c(10, 20, 30))
#'                      
#' # set seed
#' set.seed(1)
#' 
#' # find waiting time
#' t <- rexp.musse(1, rate, traits)
#' t
#' 
#' ###
#' # traits can have more states of course
#' 
#' # rates
#' rate <- c(0.1, 0.2, 0.15)
#' 
#' # traits
#' traits <- data.frame(value = c(0, 1, 2, 0),
#'                      min = c(0, 5, 10, 15),
#'                      max = c(5, 10, 15, 20))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find waiting time
#' t <- rexp.musse(1, rate, traits)
#' t
#' 
#' ###
#' # sometimes the traits data frame doesn't include all traits
#' 
#' # rates
#' rate <- c(0.1, 0.2, 0.15)
#' 
#' # traits
#' traits <- data.frame(value = c(0, 1),
#'                      min = c(0, 5),
#'                      max = c(5, 10))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find waiting time
#' t <- rexp.musse(1, rate, traits)
#' t
#'
#' @import stats
#'
#' @noRd

rexp.musse <- function(n, rate, traits,
                        now = 0, tMax = Inf, fast = FALSE) {
  # some error checking
  if (tMax == Inf & fast) {
    stop("Need a valid tMax for fast computation")
  }
  
  # make a vector to hold the results
  vars <- rep(0, n)
  
  # if tMax is not supplied, need another upper for uniroot
  upper <- ifelse(tMax == Inf, now + 100, tMax)
  # time will still be truncated without tMax, but only if fast = TRUE
  
  # manipulate traits to include only relevant portion
  if (now < max(traits$max)) {
    # delete the rows in traits that are lower than now
    traits <- traits[traits$max > now, ]
    
    # substitute the first min by now
    traits$min[1] <- now
  } else {
    # make last max now
    traits$max[length(traits$max)] <- now
    
    # delete all other rows
    traits <- traits[nrow(traits), ]
  }
  
  
  # if rate has length 1, make it a vector
  if (length(rate) == 1) {
    rate <- rep(rate, max(traits$value) + 1)
  }
  
  # if only one rate is to be considered, return rexp
  if (length(unique(rate)) == 1 || nrow(traits) == 1) {
    if (rate[traits$value[1] + 1] == 0) {
      return(Inf)
    } else {
      return(rexp(n, rate[traits$value[1] + 1]))
    }
  }
  # this happens when rates are all the same, or 
  # when there are no rate shifts after now
  
  # make a rates data frame
  rates <- data.frame(rate[traits$value + 1], traits$min, traits$max)
  colnames(rates) <- colnames(traits)
  
  for (i in 1:n) {
    # draw an uniform random variable from 0 to 1
    p <- runif (1)
    
    # if f(t) = 0, t is exponentially distributed
    f <- Vectorize(function(t) {
      rates <- rates[rates$min < t, ]
      rates$max[length(rates$max)] <- t
      1 - p - exp(-sum(rates$value*(rates$max - rates$min)))
    })
    
    if ((rates$value[length(rates$value)] == 0) && 
        (f(rates$max[length(rates$max)]) < 0)) {
      vars[i] <- Inf
      next
    }
    
    # change tolerance if p is really small
    tol <- ifelse(p < 1e-5, .Machine$double.eps^0.5, .Machine$double.eps^0.25)
    
    # find the root and subtract now
    vars[i] <- uniroot(f, c(now, now + 100), extendInt = "yes", 
                       tol = tol)$root - now
    
    # when p is really small, uniroot might find now as the root since
    # the real root is really small, so we increase the tolerance
    while (vars[i] == 0) {
      tol <- tol / 2
      vars[i] <- uniroot(f, c(now, now + 100), extendInt = "yes", 
                         tol = tol)$root - now
    }
  }
  return(vars)
}