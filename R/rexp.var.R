#' General rate exponential and Weibull waiting times
#' 
#' Generates a waiting time following an exponential or Weibull distribution 
#' with constant or varying rates. Output can be used as the waiting time to an
#' extinction, speciation, or fossil sampling event. Allows for an optional 
#' shape parameter, in which case \code{rate} will be taken as a Weibull scale.
#' Allows for further customization by restricting possible waiting time 
#' outputs with arguments for (1) current time, to consider only the rates
#' starting at that time, (2) maximum time, to bound the output and therefore
#' allow for faster calculations if one only cares about waiting times lower 
#' than a given amount, and (3) speciation time, necessary to scale rates in
#' the case where the output is to follow a Weibull distribution, i.e. for 
#' age-dependent processes. This function is used in birth-death and sampling
#' functions, but can also be convenient for any user looking to write their
#' own code requiring exponential or Weibull distributions with varying rates.
#'
#' @param n The number of waiting times to return. Usually 1, but we allow for
#' a higher \code{n} to be consistent with the \code{rexp} function.
#' 
#' @param rate The rate parameter for the exponential distribution. If
#' \code{shape} is not \code{NULL}, \code{rate} is a scale for a Weibull
#' distribution. In both cases we allow for any time-varying function. Note
#' \code{rate} can be constant.
#'
#' @param now The current time. Needed if one wants to consider only the 
#' interval between the current time and the maximum time for the time-varying
#' rate. Note this does means the waiting time is \code{>= now}, so we also 
#' subtract \code{now} from the result before returning. The default is 
#' \code{0}.
#'
#' @param tMax The simulation ending time. If the waiting time would be too
#' high, we return \code{tMax + 0.01} to signify the event never happens, if
#' \code{fast == TRUE}. Otherwise we return the true waiting time. By default,
#' \code{tMax} will be \code{Inf}, but if \code{FAST == TRUE} one must supply
#' a finite value.
#'
#' @param shape Shape of the Weibull distribution. Can be a \code{numeric} for 
#' constant shape or a \code{function(t)} for time-varying. When smaller than
#' one, rate will decrease along species' age (negative age-dependency). When 
#' larger than one, rate will increase along each species' age (positive 
#' age-dependency). Default is \code{NULL}, so the function acts as an 
#' exponential distribution. For \code{shape != NULL} (including when equal to 
#' one), \code{rate} will be considered a scale (= 1/rate), and \code{rexp.var} 
#' will draw a Weibull distribution instead of an exponential. This means 
#' Weibull(rate, 1) = Exponential(1/rate). Notice even when 
#' \code{Shape != NULL}, \code{rate} may still be time-dependent. 
#'
#' Note: Time-varying shape is within expectations for most cases, but if it is
#' lower than 1 and varies too much (e.g. \code{0.5 + 0.5*t}), it can be 
#' slightly biased for higher waiting times due to computational error. Slopes 
#' (or equivalent, since it can be any function of time) of the order of 0.01 
#' are advisable. It rarely also displays small biases for abrupt variations. 
#' In both cases, error is still quite low for the purposes of the package.
#' 
#' Note: We do not test for shape < 0 here since as we allow shape to be a 
#' function this would severely slow the rest of the package. It is tested on
#' the birth-death functions, and the user should make sure not to use any
#' functions that become negative eventually.
#'
#' @param TS Speciation time, used to account for the scaling between simulation 
#' and species time. The default is \code{0}. Supplying a \code{TS > now} will
#' return an error.
#'
#' @param fast If set to \code{FALSE}, waiting times larger than \code{tMax} will
#' not be thrown away. This argument is needed so one can test the function
#' without bias.
#'
#' @return A vector of waiting times following the exponential or Weibull 
#' distribution with the given parameters.
#'
#' @author Bruno do Rosario Petrucci.
#'
#' @import stats
#'
#' @examples
#'
#' ###
#' # let us start by checking a simple exponential variable
#' 
#' # rate
#' rate <- 0.1
#'
#' # set seed
#' set.seed(1)
#' 
#' # find the waiting time
#' t <- rexp.var(n = 1, rate)
#' t
#' 
#' # this is the same as t <- rexp(1, rate)
#' 
#' ###
#' # now let us try a linear function for the rate
#' 
#' # rate
#' rate <- function(t) {
#'   return(0.01*t + 0.1)
#' }
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the waiting time
#' t <- rexp.var(n = 1, rate)
#' t
#' 
#' ###
#' # what if rate is exponential?
#' 
#' # rate
#' rate <- function(t) {
#'   return(0.01 * exp(0.1*t) + 0.02)
#' }
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the waiting time
#' t <- rexp.var(n = 1, rate)
#' t
#' 
#' ###
#' # if we give a shape argument, we have a Weibull distribution
#' 
#' # scale - note that this is equivalent to 1/rate if shape were NULL
#' rate <- 2
#' 
#' # shape
#' shape <- 1
#' 
#' # speciation time
#' TS <- 0
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the vector of waiting time
#' t <- rexp.var(n = 1, rate, shape = shape, TS = TS)
#' t
#' 
#' ###
#' # when shape = 1, the Weibull is an exponential, we could do better
#' 
#' # scale
#' rate <- 10
#' 
#' # shape
#' shape <- 2
#' 
#' # speciation time
#' TS <- 0
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the vector of waiting times - it doesn't need to be just one
#' t <- rexp.var(n = 5, rate, shape = shape, TS = TS)
#' t
#' 
#' ###
#' # shape can be less than one, of course
#' 
#' # scale
#' rate <- 10
#' 
#' # shape
#' shape <- 0.5
#' 
#' # note we can supply now (default 0) and tMax (default Inf)
#' 
#' # now matters when we wish to consider only waiting times
#' # after that time, important when using time-varying rates
#' now <- 3
#' 
#' # tMax matters when fast = TRUE, so that higher times will be
#' # thrown out and returned as tMax + 0.01
#' tMax <- 40
#' 
#' # speciation time - it will be greater than 0 frequently during a 
#' # simulation, as it is used to represent where in the species life we 
#' # currently are and rescale accordingly
#' TS <- 2.5
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the vector of waiting times
#' t <- rexp.var(n = 10, rate, now, tMax,
#'               shape = shape, TS = TS, fast = TRUE)
#' t
#' 
#' # note how some results are tMax + 0.01, since fast = TRUE
#' 
#' ###
#' # both rate and shape can be time varying for a Weibull
#' 
#' # scale
#' rate <- function(t) {
#'   return(0.25*t + 5)
#' }
#' 
#' # shape
#' shape <- 3
#' 
#' # current time
#' now <- 0
#' 
#' # maximum time to check
#' tMax <- 40
#' 
#' # speciation time
#' TS <- 0
#' 
#' # set seed
#' set.seed(1)
#' 
#' # find the vector of waiting times
#' t <- rexp.var(n = 5, rate, now, tMax,
#'               shape = shape, TS = TS, fast = TRUE)
#' t
#'
#' @name rexp.var
#' @rdname rexp.var
#' @export

rexp.var <- function(n, rate, 
                     now = 0, tMax = Inf, 
                     shape = NULL, TS = 0, fast = FALSE) {
  # some error checking
  if (tMax == Inf & fast) {
    stop("Need a valid tMax for fast computation")
  }
  
  # TS must be greater than current time
  if (now < TS) {
    stop("TS must be greater than the current time")
  }
  
  # make a vector to hold the results
  vars <- rep(0, n)

  # default is not age dependent, will change this later
  AD <- FALSE

  # if rate is numeric, check error and call rexp if shape is null
  if (is.numeric(rate)) {
    # check that it is one number
    if (length(rate) != 1) {
      stop("If rate is numeric, it must be one number")
    }
    
    # if shape is null, we could call rexp
    if (is.null(shape)) {
      return(rexp(n, rate))
    }
  }
  
  # if shape is a constant, make it a function
  if (!is.null(shape)) {
    # it is AD
    AD <- TRUE
    
    # test whether scale and shape are numeric
    scaleNumeric <- is.numeric(rate)
    shapeNumeric <- is.numeric(shape)
    
    # if both are numeric, we can call Weibull
    if (scaleNumeric && shapeNumeric) {
      return(rweibull(n, shape, rate))
    }
    
    # if one of them isn't, make the other a function as well
    if (shapeNumeric) {
      # make sure it is of length one
      if (length(shape) != 1) {
        stop("If shape is numeric, it must be one number")
      }
      s <- shape
      shape <- Vectorize(function(t) s)
    }
    
    # and the rate
    if (scaleNumeric) {
      r <- rate
      rate <- Vectorize(function(t) r)
    }
  }

  # if tMax is not supplied, need another upper for uniroot
  upper <- ifelse(tMax == Inf, now + 100, tMax)
  # time will still be truncated without tMax, but only if fast = TRUE
  
  for (i in 1:n) {
    # draw an uniform random variable from 0 to 1
    p <- runif (1)

    if (AD) {
      # if it is age dependent, find the current species time
      spnow <- now - TS
      
      # and adjust upper as well
      upper <- upper - TS

      # calculate the probability that the event will happen at all
      total <- 1 - exp(-(integrate(
        Vectorize(function(x) 1/rate(x + TS)), lower = spnow, 
        upper = upper, subdivisions = 2000, 
        stop.on.error = FALSE)$value) ^ shape(upper))

      # if the probability is lower than p, the event will not happen
      if (total < p & fast) {
        vars[i] <- tMax + 0.01
      }

      # if fast is false or the probability is higher than p, we find the time
      else {
        # give upper a buffer in case it is too low for uniroot to find a root
        if (total < p) {
          upper <- 1.5 * upper
        }

        # create a function to hold the CDF of the distribution minus the
        # uniform variable - if we find t where this is 0, this t is
        # distributed as a weibull
        f <- Vectorize(function(t) 
          {
          1 - p - exp(-(integrate(Vectorize(function(x) 1/rate(x + TS)), 
                                  lower = spnow, upper = t, 
                                  subdivisions = 2000, 
                                  stop.on.error = FALSE)$value) ^ shape(t))})

        # if f(t) = 0, t is distributed as a Weibull
        vars[i] <- suppressWarnings(uniroot(f, c(spnow, upper), 
                                            extendInt = "yes"))$root - spnow
        
        # if rate is really high and the integral goes to +-infinity 
        # (computationally speaking), uniroot substitutes it for a really
        # high/low value instead. Since this does not change our results, we 
        # accept it and simply suppress the warning
      }
    }
    else {
      # calculate the probability that the event will happen at all
      total <- 1 - exp(-integrate(Vectorize(function(x) rate(x)), 
                                  lower = now, upper = upper, 
                                  subdivisions = 2000, 
                                  stop.on.error = FALSE)$value)
      
      # if the probability is lower than p, the event will not happen
      if (total < p & fast) {
        vars[i] <- tMax + 0.01
      }
      
      # if fast is false or the probability is higher than p, we find the time
      else {
        # give upper a buffer in case it is too low for uniroot to find a root
        if (total < p) {
          upper <- 1.5 * upper
        }
        
        # if f(t) = 0, t is exponentially distributed
        f <- Vectorize(function(t) {
          1 - p - exp(-integrate(Vectorize(function(x) rate(x)), lower = now, 
                                 upper = t, subdivisions = 2000, 
                                 stop.on.error = FALSE)$value)})

        # if rate is really high and the integral goes to +-infinity 
        # (computationally speaking), uniroot substitutes it for a really
        # high/low value instead. Since this does not change our results, we 
        # accept it and simply suppress the warning
        vars[i] <- max(suppressWarnings(uniroot(f, c(now, upper), 
                                            extendInt = "yes"))$root - now,
                       1e-6)
        # max so that we never have a 0 due to numerical issues when p is small
      }
    }
  }
  return(vars)
}
