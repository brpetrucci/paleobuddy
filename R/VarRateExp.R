#' Calculates the expected number of species for a given time varying rate
#' using the exponential distribution for variable rates.
#'
#' \code{VarRateExp} takes a function, an initial number of species and a time
#' vector and calculates the predicted exponential with that function as rate
#' on that interval. This allows for efficient testing of the diversity curves
#' produced by \code{PaleoBuddy} simulations.
#'
#' @param ff a rate for the exponential distribution that can be
#' any function of time. One can also supply data for an environmental
#' variable (see below for the \code{env_f} param) and get the expected
#' number of species for a hybrid function of time and said variable. Finally,
#' one can instead supply a vector of rates to \code{ff} and a vector of shifts
#' to \code{fshifts} and get a step function. It is more efficient to create a
#' stepfunction using \code{ifelse} however (see examples below).
#'
#' @param N0 the initial number of species is by default 1, but one
#' can change to any positive number. We allow for negative initial values as
#' well, but of course that will not help in testing the package.
#'
#' NOTE: \code{VarRateExp} will find the expected number of daughters given a
#' rate \code{ff} and an initial number of parents \code{N0}, so in a
#' biological context \code{ff} is diversification rate, not speciation (unless
#' extinction is 0, of course).
#'
#' @param t a time vector over which to consider the distribution.
#'
#' @param env_f a two dimensional dataframe with time as a first
#' column and the desired environmental variable as a second. Note that
#' supplying a function with one argument and a non-NULL \code{env_f}, and vice
#' versa, will return an error.
#'
#' @param fshifts a vector of rate shifts. Then used with the rates
#' vector to create a step function for the rates. If supplied without a rates
#' vector, and vice versa, will return an error.
#'
#' @return a vector of the expected number of species per time point supplied
#' in \code{t}, which can then be used to plot vs. \code{t}.
#'
#' @examples
#'
#' # let us create a vector of times to use in these examples.
#' T <- seq(0, 50, 0.1)
#'
#' # we will proceed with the same examples as the \code{MakeRate} function.
#' # note that here \code{ff} means the diversification rate, if we want to
#' # interpret this in a species birth-death context.
#'
#' # let us start simple: create a constant rate
#' ff <- 0.1
#'
#' # set this up so we see rates next to diversity
#' par(mfrow = c(1,2))
#'
#' # see how the rate looks
#' r <- MakeRate(0.5)
#' plot(T, r(T), type='l')
#'
#' # get the diversity and plot it
#' div <- VarRateExp(ff, t=T)
#' plot(T, div, type='l')
#'
#' # something a bit more complex: a linear rate
#' ff <- function(t) {
#'   return(0.01*t)
#' }
#' r <- MakeRate(ff)
#' plot(T, r(T), type='l')
#'
#' div <- VarRateExp(ff, t=T)
#' plot(T, div, type='l')
#'
#' # remember: ff is diversity!
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
#' # ...and code ff as diversification
#' ff <- function(t) {
#'   return(pp(t) - qq(t))
#' }
#' r <- MakeRate(ff)
#' plot(T, r(T), type='l')
#'
#' # not a good time for this poor clade
#' div <- VarRateExp(ff, N0 = 2, T)
#' plot(T, div, type='l')
#'
#' # remember: ff can be any time-varying function!
#' ff <- function(t) {
#'   return(sin(t)*0.5)
#' }
#' r <- MakeRate(ff)
#' plot(T, r(T), type='l')
#'
#' div <- VarRateExp(ff, N0=2, T)
#' plot(T, div, type='l')
#'
#' # we can use ifelse() to make a step function like this
#' ff <- function(t) {
#'   return(ifelse(t < 10, 0.1,
#'                 ifelse(t < 20, 0.3,
#'                        ifelse(t < 30, -0.2,
#'                               ifelse(t < 40, 0.05, 0)))))
#' }
#' r <- MakeRate(ff)
#' plot(T, r(T), type='l')
#'
#' div <- VarRateExp(ff, t=T)
#' plot(T, div, type='l')
#' # important note: this method of creating a step function might be annoying,
#' # but when running thousands of simulations it will provide a much faster
#' # integration than when using our method of transforming a rates and shifts
#' # vector into a function of time...
#'
#' # ...which we can do as follows
#' ff <- c(0.2, -0.1, 0.1, -0.05)
#' fshifts <- c(0, 10, 20, 35)
#' r <- MakeRate(ff, fshifts = fshifts)
#' plot(T, r(T),type='l')
#'
#' div <- VarRateExp(ff, t=T, fshifts=fshifts)
#' plot(T, div, type='l')
#'
#' # note the delay in running VarRateExp using this method. Integrating a step
#' # function created using the methods in MakeRate() is slow, as explained in
#' # the above and below notes (reproduced from the MakeRate documentation)
#'
#' # it is also impractical to supply a rate and a shifts vector and
#' # have an environmental dependency, so in cases where one looks to run
#' # more than a couple dozen simulations, and when one is looking to have a
#' # step function modified by an environmental variable, consider using ifelse()
#'
#' # finally let us see what we can do with environmental variables
#'
#' # RPANDA supplies us with some really useful environmental dataframes
#' library(RPANDA)
#'
#' # to use as an example, let us try temperature
#' data(InfTemp)
#'
#' ff <- function(t, env) {
#'   return(0.002*env)
#' }
#' r <- MakeRate(ff, env_f = InfTemp)
#' plot(T, r(T), type='l')
#'
#' div <- VarRateExp(ff, t=T, env_f=InfTemp)
#' plot(T, div, type='l')
#'
#' # we can also have a function that depends on both time AND temperature
#' ff <- function(t, env) {
#'   return(0.03 * env - 0.01 * t)
#' }
#' r <- MakeRate(ff, env_f = InfTemp)
#' plot(T, r(T), type='l')
#'
#' div <- VarRateExp(ff, t=T, env_f=InfTemp)
#' plot(T, div, type='l')
#'
#' # as mentioned above, we could also use ifelse() to construct a step function
#' # that is modulated by temperature
#' ff <- function(t, env) {
#'   return(ifelse(t < 10, 0.1 + 0.01*env,
#'                 ifelse(t < 30, 0.2 - 0.005*env,
#'                        ifelse(t <= 50, 0.1 + 0.005*env, 0))))
#' }
#' r <- MakeRate(ff, env_f = InfTemp)
#' plot(T, r(T), type='l')
#'
#' div <- VarRateExp(ff, t=T, env_f = InfTemp)
#' plot(T, div, type='l')
#'
#' @name VarRateExp
#' @rdname VarRateExp
#' @export

VarRateExp<-function(ff,N0=1,t,env_f=NULL,fshifts=NULL){
  # get the corresponding rate
  f <- MakeRate(ff, env_f = env_f, fshifts = fshifts)

  if (!is.numeric(f)) {
    # integrate the rate for each t
    Integral<-lapply(t,function(x){integrate(Vectorize(f),0,x,subdivisions=2000)[1]})

    # make the integral numerical so we can plot
    for (i in 1:length(Integral)){Integral[i]<-as.numeric(Integral[[i]])}
    Integral<-as.numeric(Integral)
  }

  else {
    Integral <- f*t
  }

  return(N0*exp(Integral))
}
