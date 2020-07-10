#' General rate exponential and Weibull waiting times
#'
#' \code{rexp_var} uses a rate (that could be constant or time-varying), a range
#' of time and optionally a shape for age dependent rates. It also takes a
#' number of waiting times to return, and whether the user wishes to throw away
#' waiting times that pass \code{tMax}. It returns a time from an exponential
#' distribution, or a weibull distribution if \code{shape != NULL}. 
#'
#' @param n the number of waiting times to return. The default is 1, but we
#'  allow for a higher \code{n} to be consistent with the \code{rexp} function.
#' 
#' @param lambda the rate parameter for the exponential
#' distribution. If shape is not \code{NULL}, \code{lambda} is a scale for a
#' Weibull distribution. In both cases we allow for any time-varying function. 
#' If one wants a constant \code{c}, please use 
#' \code{lambda <- Vectorize(function(t) c)}.
#'
#' @param now the current time. Needed so that we consider only the interval
#' between the current time and the maximum time for the time-varying rate.
#' Notice this does means the waiting time is \code{>= now}, so we also
#' subtract \code{now} from the result before returning.
#'
#' @param tMax the simulation ending time. If the waiting time would be too
#' high, we return \code{2*tMax} to signify the event never happens, if 
#' \code{FAST == TRUE}. The function only considers the rate between \code{now}
#' and \code{tMax}.
#'
#' @param shape the shape of a weibull distribution. If not \code{NULL}, the 
#' distribution dis taken to be a weibull. Otherwise, it is considered an
#' exponential.
#'
#' Notes: if \code{shape} is really low it may be impossible
#' for \code{uniroot} to find a root; Time-varying shape is implemented, but
#' not yet thoroughly tested.
#'
#' @param TS if shape is given, there must be a \code{TS} parameter to account
#' for the scaling between simulation and species time. Supplying one without
#' the other leads to an error.
#'
#' @param fast if set to \code{FALSE}, waiting times larger than \code{tMax} 
#' will not be thrown away. This argument is needed so one can testt he function
#' without bias.
#'
#' @return a vector of waiting times for the exponential or weibull
#' distribution with the given rates.
#'
#' @author written by Bruno do Rosario Petrucci.
#'
#' @import stats
#'
#' @examples
#'
#' # Note: to effectively test, you must pass \code{fast=FALSE} to the function.
#' 
#' # use the \code{fitdistrplus} package to see if our distributions fit what
#' # we expect
#' if (requireNamespace("fitdistrplus", quietly = TRUE)) {
#'   # let us start by checking a simple exponential variable
#'   
#'   # rate
#'   lambda <- 0.1
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 10000, lambda, now, tMax, fast = FALSE)
#'   
#'   # find rate
#'   rate <- unname(fitdistrplus::fitdist(dis, distr = "exp")$estimate)
#'   
#'   # check
#'   print(paste("rate: ", paste(rate, paste(", expected: ", lambda))))
#' }
#'   
#'   # another simple exponential
#'   
#'   # rate
#'   lambda <- 0.5
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 10000, lambda, now, tMax, fast = FALSE)
#'   
#'   # find rate
#'   rate <- unname(fitdistrplus::fitdist(dis, distr = "exp")$estimate)
#'   
#'   # check
#'   print(paste("rate: ", paste(rate, paste(", expected: ", lambda))))
#' }
#'   
#'   # now let us try a linear function
#'   
#'   # rate
#'   lambda <- function(t) {
#'     return(0.01*t + 0.1)
#'   }
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 10000, lambda, now, tMax, fast = FALSE)
#'   
#'   # find the rate
#'   rate <- unname(fitdistrplus::fitdist(dis, distr = "exp")$estimate)
#'   
#'   # we can integrate to find the mean rate
#'   
#'   # the integrate of 1/this is the mean rate
#'   f <- Vectorize(function(t) {
#'     ifelse(lambda(t) == Inf, 0,
#'            exp(-integrate(function(x) lambda(x),
#'                           lower = 0, upper = t, subdivisions = 2000)$value))})
#'   
#'   # calculate the mean
#'   mean <- 1/integrate(f, 0, Inf)$value
#'   
#'   # check
#'   print(paste("rate: ", paste(rate, paste(", expected: ", mean))))
#'   
#'   # and we can also integrate to find the mean directly
#'   
#'   # get distribution function
#'   p <- Vectorize(function(t) {
#'     return(
#'       ifelse(lambda(t) == Inf, 0,
#'              lambda(t)*exp(-integrate(Vectorize(
#'                function(x) lambda(x)), lower = 0, upper = t)$value))
#'     )
#'   })
#'   
#'   # integrate it on the whole space
#'   mean<-integrate(function(t) t*p(t), 0, Inf)$value
#'   
#'   # check
#'   print(paste("mean: ", paste(mean(dis), paste(", expected: ", mean))))
#' }
#'   
#'   # what if lambda is exponential?
#'   
#'   # rate
#'   lambda <- function(t) {
#'     return(0.01 * exp(0.1*t) + 0.02)
#'   }
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 10000, lambda, now, tMax, fast = FALSE)
#'   
#'   # we can integrate to find the mean rate
#'   
#'   # the integrate of 1/this is the mean rate
#'   f <- Vectorize(function(t) {
#'     ifelse(lambda(t) == Inf, 0,
#'            exp(-integrate(function(x) lambda(x),
#'                           lower = 0, upper = t, subdivisions = 2000)$value))})
#'   
#'   # calculate the mean
#'   mean <- 1/integrate(f, 0, Inf)$value
#'   
#'   # check
#'   print(paste("rate: ", paste(rate, paste(", expected: ", mean))))
#'   
#'   # and we can also integrate to find the mean directly
#'   
#'   # get distribution function
#'   p <- Vectorize(function(t) {
#'     return(
#'       ifelse(lambda(t) == Inf, 0,
#'              lambda(t)*exp(-integrate(Vectorize(
#'                function(x) lambda(x)), lower = 0, upper = t)$value))
#'     )
#'   })
#'   
#'   # integrate it on the whole space
#'   mean<-integrate(function(t) t*p(t), 0, Inf)$value
#'   
#'   # check
#'   print(paste("mean: ", paste(mean(dis), paste(", expected: ", mean))))
#' }
#'   
#'   # now we can also test the age dependency
#'   
#'   # scale
#'   lambda <- 2
#'   
#'   # shape
#'   shape <- 1
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 5000, lambda, now, tMax, 
#'                   shape = shape, TS = 0, fast = FALSE)
#'   
#'   # find the estimates
#'   
#'   # shape
#'   eShape <- unname(fitdistrplus::fitdist(dis, distr = "weibull")$estimate)[1]
#'   
#'   # scale
#'   eScale <- unname(fitdistrplus::fitdist(dis, distr = "weibull")$estimate)[2]
#'   
#'   # check
#'   print(paste("shape: ", paste(eShape, paste(", expected: ", shape))))
#'   print(paste("scale: ", paste(eScale, paste(", expected: ", lambda))))
#'   
#' }
#'   
#'   # shape = 1 is an exponential, we could do better
#'   
#'   # scale
#'   lambda <- 10
#'   
#'   # shape
#'   shape <- 2
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{   
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 5000, lambda, now, tMax, 
#'                   shape = shape, TS = 0, fast = FALSE)
#'   
#'   # find the estimates
#'   
#'   # shape
#'   eShape <- unname(fitdistrplus::fitdist(dis, distr = "weibull")$estimate)[1]
#'   
#'   # scale
#'   eScale <- unname(fitdistrplus::fitdist(dis, distr = "weibull")$estimate)[2]
#'   
#'   # check
#'   print(paste("shape: ", paste(eShape, paste(", expected: ", shape))))
#'   print(paste("scale: ", paste(eScale, paste(", expected: ", lambda))))
#' }
#'   
#'   # fitdist gets a bit weird with shape less than 1,
#'   # so we need some extra arguments to get a fit
#'   
#'   # scale
#'   lambda <- 10
#'   
#'   # shape
#'   shape <- 0.5
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n=10000, lambda, now, tMax, 
#'                   shape = shape, TS = 0, fast = FALSE)
#'   
#'   # find the estimates
#'   
#'   # shape
#'   eShape <- unname(fitdistrplus::fitdist(dis, distr = "weibull",
#'                                          start = list(shape = 1, scale = 1),
#'                                          method = "mge", gof = "CvM")$estimate)[1]
#'   
#'   # scale
#'   eScale <- unname(fitdistrplus::fitdist(dis, distr = "weibull",
#'                                          start = list(shape = 1, scale = 1),
#'                                          method = "mge", gof = "CvM")$estimate)[2]
#'   
#'   # check
#'   print(paste("shape: ", paste(eShape, paste(", expected: ", shape))))
#'   print(paste("scale: ", paste(eScale, paste(", expected: ", lambda))))
#' }
#'   
#'   # when lambda varies, we might need to do some more to test
#'   
#'   # scale
#'   lambda <- function(t) {
#'     return(0.25*t + 5)
#'   }
#'   
#'   # shape
#'   shape <- 3
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 5000, lambda, now, tMax, 
#'                   shape = shape, TS = 0, fast = FALSE)
#'   
#'   # we will integrate the distribution to find the mean directly
#'   
#'   # distribution function
#'   p <- Vectorize(function(t) {
#'     res<-shape/lambda(t)*(integrate(
#'       Vectorize(function(x) 1/lambda(x)), lower = 0, upper = t,
#'       subdivisions = 2000)$value) ^ (shape - 1)*
#'       exp(-(integrate(Vectorize(function(x) 1/lambda(x)), 0, t,
#'                       subdivisions = 2000)$value)^shape)
#'     # na simply means R doesn't know how to multiply 0 by infinity, but we just
#'     # need to make it 0 since t*exp(-t) goes to 0 when t goes to infinity
#'     return(res)
#'   })
#'   
#'   # find the mean
#'   mean <- integrate(function(t) t*p(t), 0, Inf)$value
#'   
#'   # check
#'   print(paste("rate: ", paste(mean(dis), paste(", expected: ", mean))))
#' }
#'   
#'   # lambda can be any function of time, remember
#'   
#'   # scale
#'   lambda <- function(t) {
#'     return(0.2*exp(0.1*t) + 5)
#'   }
#'   
#'   # shape
#'   shape <- 3
#'   
#'   # current time
#'   now <- 0
#'   
#'   # maximum time to check 
#'   tMax <- 40
#'   
#' \dontrun{
#'   # find the list of waiting times
#'   dis <- rexp_var(n = 5000, lambda, now, tMax, 
#'                   shape = shape, TS = 0, fast = FALSE)
#'   
#'   # we will integrate the distribution to find the mean directly
#'   
#'   # distribution function
#'   p <- Vectorize(function(t) {
#'     res<-shape/lambda(t)*(integrate(
#'       Vectorize(function(x) 1/lambda(x)), lower = 0, upper = t,
#'       subdivisions = 2000)$value) ^ (shape-1)*
#'       exp(-(integrate(Vectorize(function(x) 1/lambda(x)), 0, t,
#'                       subdivisions = 2000)$value) ^ shape)
#'     # na simply means R doesn't know how to multiply 0 by infinity, but we just
#'     # need to make it 0 since t*exp(-t) goes to 0 when t goes to infinity
#'     return(ifelse(is.na(res), 0, res))
#'   })
#'   
#'   # find the mean
#'   mean <- integrate(function(t) t*p(t), 0, Inf)$value
#'   
#'   # check
#'   print(paste("rate: ", paste(mean(dis), paste(", expected: ", mean))))
#' }
#' }
#'
#' @name rexp_var
#' @rdname rexp_var
#' @export

rexp_var<-function(n = 1, lambda, now, tMax, shape = NULL, TS = NULL, fast = TRUE) {
  # make a vector to hold the results
  vars <- rep(0, n)

  # default is not age dependent, will change this later
  AD <- FALSE

  # make lambda a function if it is a constant
  l <- lambda
  lambda <- ifelse(is.numeric(l), Vectorize(function(t) l), l)

  # same for shape
  if (!is.null(shape)) {
    AD <- TRUE
    s <- shape
    shape <- ifelse(is.numeric(s), Vectorize(function(t) s), s)
  }

  for (i in 1:n) {
    # draw an uniform random variable from 0 to 1
    p <- runif (1)

    if (AD) {
      # if it is age dependent, find the current species time
      spnow <- now - TS

      # calculate the probability that the event will happen at all
      total <- 1 - exp(-(integrate(
        Vectorize(function(x) 1/lambda(x + TS)), lower = spnow, 
        upper = tMax - TS, subdivisions = 2000)$value) ^ shape(tMax))

      # if the probability is lower than p, the event will not happen
      if (total < p & fast) {
        vars[i] <- 2*tMax + 0.01
      }

      else {

        # create a function to hold the CDF of the distribution minus the
        # uniform variable - if we find t where this is 0, this t is
        # distributed as a weibull
        f <- Vectorize(function(t) 
          {
          1 - p - exp(-(integrate(Vectorize(function(x) 1/lambda(x + TS)), 
                                  lower = spnow, upper = t, 
                                  subdivisions = 2000)$value) ^ shape(t))})

        # if f(t) = 0, t is distributed as a Weibull
        vars[i] <- suppressWarnings(uniroot(f, c(spnow, tMax), 
                                            extendInt="yes"))$root - spnow
        
        # if lambda is really high and the integral goes to +-infinity (computationally
        # speaking), uniroot substitutes it for a really high/low value instead. Since
        # this does not change our results, we accept it and simply suppress the warning
      }
    }
    else {
      # calculate the probability that the event will happen at all
      total <- 1 - exp(-integrate(Vectorize(function(x) lambda(x)), 
                                  lower = now, upper = tMax, 
                                  subdivisions = 2000)$value)
      
      # if the probability is lower than p, the event will not happen
      if (total < p & fast) {
        vars[i] <- 2*tMax + 0.01
      }
      
      else {
        # if f(t) = 0, t is exponentially distributed
        f <- Vectorize(function(t) {
          1 - p - exp(-integrate(Vectorize(function(x) lambda(x)), lower = now, 
                                 upper = t, subdivisions = 2000)$value)})

        # if lambda is really high and the integral goes to +-infinity (computationally
        # speaking), uniroot substitutes it for a really high/low value instead. Since
        # this does not change our results, we accept it and simply suppress the warning
        vars[i] <- suppressWarnings(uniroot(f, c(0, tMax), 
                                            extendInt="yes"))$root - now
      }
    }
  }
  return(vars)
}
