#' Expected diversity for general exponential rates
#'
#' Calculates the expected species diversity on an interval given a (possibly time
#' dependent) exponential rate.  Takes as the base rate (1) a constant, (2) a 
#' function of time, (3) a function of time interacting with an environmental 
#' variable, or (4) a vector of numbers describing rates as a step function. 
#' Requires information regarding the maximum simulation time, and allows for 
#' optional extra parameters to tweak the baseline rate. For more information on
#' the creation of the final rate, see \code{make.rate}.
#'
#' @inheritParams make.rate
#' 
#' @param t A time vector over which to consider the distribution.
#'
#' @param n0 The initial number of species is by default 1, but one can change to
#' any nonnegative number.
#'
#' Note: \code{var.rate.div} will find the expected number of daughters given a
#' rate \code{rate} and an initial number of parents \code{n0}, so in a
#' biological context \code{rate} is diversification rate, not speciation (unless
#' extinction is \code{0}).
#'
#' @return A vector of the expected number of species per time point supplied
#' in \code{t}, which can then be used to plot vs. \code{t}.
#'
#' @examples
#'
#' # let us first create a vector of times to use in these examples.
#' time <- seq(0, 50, 0.1)
#' 
#' ###
#' # we can start simple: create a constant rate
#' rate <- 0.1
#' 
#' # set this up so we see rates next to diversity
#' par(mfrow = c(1,2))
#' 
#' # see how the rate looks
#' r <- make.rate(0.5)
#' plot(time, rep(r, length(time)), type = 'l')
#' 
#' # get the diversity and plot it
#' div <- var.rate.div(rate, time)
#' plot(time, div, type = 'l')
#' 
#' ###
#' # something a bit more complex: a linear rate
#' rate <- function(t) {
#'   return(0.1 - 0.005*t)
#' }
#' 
#' # visualize the rate
#' r <- make.rate(rate)
#' plot(time, r(time), type = 'l')
#' 
#' # get the diversity and plot it
#' div <- var.rate.div(rate, t = time)
#' plot(time, div, type = 'l')
#' 
#' ###
#' # remember: rate is diversity!
#' 
#' # we can create speciation...
#' pp <- function(t) {
#'   return(-0.01*t + 0.2)
#' }
#' 
#' # ...and extinction...
#' qq <- function(t) {
#'   return(0.01*t)
#' }
#' 
#' # ...and code rate as diversification
#' rate <- function(t) {
#'   return(pp(t) - qq(t))
#' }
#' 
#' # visualize the rate
#' r <- make.rate(rate)
#' plot(time, r(time), type = 'l')
#' 
#' # get diversity and plot it
#' div <- var.rate.div(rate, time, n0 = 2)
#' plot(time, div, type = 'l')
#' 
#' ###
#' # remember: rate can be any time-varying function!
#' 
#' # such as a sine
#' rate <- function(t) {
#'   return(sin(t)*0.5)
#' }
#' 
#' # visualize the rate
#' r <- make.rate(rate)
#' plot(time, r(time), type = 'l')
#' 
#' # we can have any number of starting species
#' div <- var.rate.div(rate, time, n0 = 2)
#' plot(time, div, type = 'l')
#' 
#' ###
#' # we can use ifelse() to make a step function like this
#' rate <- function(t) {
#'   return(ifelse(t < 2, 0.1,
#'           ifelse(t < 3, 0.3,
#'            ifelse(t < 5, -0.2, 0.05))))
#' }
#' 
#' # change t so things are faster
#' time <- seq(0, 10, 0.1)
#' 
#' # visualize the rate
#' r <- make.rate(rate)
#' plot(time, r(time), type = 'l')
#' 
#' # get the diversity and plot it
#' div <- var.rate.div(rate, time)
#' plot(time, div, type = 'l')
#' 
#' # important note: this method of creating a step function might be annoying,
#' # but when running thousands of simulations it will provide a much faster
#' # integration than when using our method of transforming a rates and shifts
#' # vector into a function of time
#' 
#' ###
#' # ...which we can do as follows
#' 
#' # rates vector
#' rate <- c(0.1, 0.3, -0.2, 0.05)
#' 
#' # rate shifts vector
#' rateShifts <- c(0, 2, 3, 5)
#' 
#' # visualize the rate
#' r <- make.rate(rate, tMax = 10, rateShifts = rateShifts)
#' plot(time, r(time),type = 'l')
#'  
#' \dontrun{
#' # get the diversity and plot it
#' div <- var.rate.div(rate, time, tMax = 10, rateShifts = rateShifts)
#' plot(time, div, type = 'l')
#' }
#' 
#' # we set it to not run due to the time this method takes. integrating a step
#' # function created using the methods in make.rate() is slow, as explained in
#' # the make.rate documentation)
#' 
#' # it is also impractical to supply a rate and a shifts vector and
#' # have an environmental dependency, so in cases where one looks to run
#' # more than a couple dozen simulations, and when one is looking to have a
#' # step function modified by an environmental variable, consider using ifelse()
#' 
#' # finally let us see what we can do with environmental variables
#' 
#' # get the temperature data
#' data(temp)
#' 
#' # diversification
#' rate <- function(t, env) {
#'   return(0.002*env)
#' }
#' 
#' # visualize the rate
#' r <- make.rate(rate, tMax = tMax, envRate = temp)
#' plot(time, r(time), type = 'l')
#' 
#' 
#' # get diversity and plot it
#' div <- var.rate.div(rate, time, tMax = tMax, envRate = temp)
#' plot(time, div, type = 'l')
#' 
#' ###
#' # we can also have a function that depends on both time AND temperature
#' 
#' # diversification
#' rate <- function(t, env) {
#'   return(0.02 * env - 0.001 * t)
#' }
#' 
#' # visualize the rate
#' r <- make.rate(rate, tMax = tMax, envRate = temp)
#' plot(time, r(time), type = 'l')
#' 
#' # get diversity and plot it
#' div <- var.rate.div(rate, time, tMax = tMax, envRate = temp)
#' plot(time, div, type = 'l')
#'   
#' ###
#' # as mentioned above, we could also use ifelse() to construct a step function
#' # that is modulated by temperature
#' 
#' # diversification
#' rate <- function(t, env) {
#'   return(ifelse(t < 2, 0.1 + 0.01*env,
#'           ifelse(t < 5, 0.2 - 0.005*env,
#'            ifelse(t < 8, 0.1 + 0.005*env, 0.08))))
#' }
#' 
#' # visualize the rate
#' r <- make.rate(rate, tMax = tMax, envRate = temp)
#' plot(time, r(time), type = 'l')
#' 
#' \dontrun{
#' # get diversity and plot it
#' div <- var.rate.div(rate, time, tMax = tMax, envRate = temp)
#' plot(time, div, type = 'l')
#' }
#' 
#' # takes a bit long so we set it to not run, but the user
#' # should feel free to explore this and other scenarios
#' 
#' @name var.rate.div
#' @rdname var.rate.div
#' @export

var.rate.div <- function(rate, t, n0 = 1, tMax = NULL, 
                         envRate = NULL, rateShifts = NULL) {
  # check that n0 is nonnegative
  if (n0 < 0) {
    stop("initial number of species n0 must be nonnegative")
  }
  
  # check that t is a numeric vector
  if (!is.numeric(t)) {
    stop("t must be a numeric vector")
  }
  
  # get the corresponding rate
  r <- make.rate(rate, tMax = tMax, envRate = envRate, rateShifts = rateShifts)

  if (!is.numeric(r)) {
    # integrate the rate for each t
    integral <- lapply(t, function(x) {
      integrate(Vectorize(r), 0, x, subdivisions = 2000)$value
      })

    # make the integral numerical so we can plot
    for (i in 1:length(integral)) {
      integral[i] <- as.numeric(integral[[i]])
      }
    integral <- as.numeric(integral)
  }

  else {
    integral <- r*t
  }

  return(n0*exp(integral))
}
