#' Create a flexible rate for birth-death or sampling simulations
#' 
#' Generates a function determining the variation of a rate (speciation,
#' extinction, sampling) with respect to time. To be used on birth-death or 
#' sampling functions, it takes as the base rate (1) a constant, (2) a function
#' of time, (3) a function of time and a time-series (usually an environmental
#' variable), or (4) a vector of numbers describing rates as a step function. 
#' Requires information regarding the maximum simulation time, and allows for 
#' optional extra parameters to tweak the baseline rate.
#'
#' @param rate The baseline function with which to make the rate.
#' It can be a
#'
#' \describe{
#' \item{A number}{For constant birth-death rates.}
#'
#' \item{A function of time}{For rates that vary with time. Note that this can 
#' be any function of time.} 
#' 
#' \item{A function of time and an environmental variable}{For rates varying 
#' with time and an environmental variable, such as temperature. Note that 
#' supplying a function on more than one variable without an accompanying 
#' \code{envRate} will result in an error.}
#'
#' \item{A numeric vector}{To create step function rates. Note this must be
#'  accompanied by a corresponding vector of rate shift times, 
#'  \code{rateShifts}.}}
#'
#' @param tMax Ending time of simulation, in million years after the clade's 
#' origin. Needed to ensure \code{rateShifts} runs the correct way.
#'
#' @param envRate A \code{data.frame} representing a time-series, usually an 
#' environmental variable (e.g. CO2, temperature, etc) varying with time. The 
#' first column of this \code{data.frame} must be time, and the second column 
#' must be the values of the variable. The function will return an error if 
#' the user supplies \code{envRate} without \code{rate} being a function of two
#' variables. \code{paleobuddy} has two environmental data frames, \code{temp}
#' and \code{co2}. One can check \code{RPANDA} for more examples.
#' 
#' Note that, since simulation functions are run in forward-time (i.e. with 
#' \code{0} being the origin time of the simulation), the time column of
#' \code{envRate} is assumed to do so as well, so that the row corresponding to
#' \code{t = 0} is assumed to be the value of the time-series when the 
#' simulation starts, and \code{t = tMax} is assumed to be its value when the
#' simulation ends (the present).
#' 
#' Acknowledgements: The strategy to transform a function of \code{t} and 
#' \code{envRate} into a function of \code{t} only using \code{envRate} was 
#' adapted from \code{RPANDA}.
#'
#' @param rateShifts A vector indicating the time of rate shifts in a step 
#' function. The first element must be the first or last time point for the 
#' simulation, i.e. \code{0} or \code{tMax}. Since functions in paleobuddy run
#' from \code{0} to \code{tMax}, if \code{rateShifts} runs from past to present
#' (meaning \code{rateShifts[2] < rateShifts[1]}), we take 
#' \code{tMax - rateShifts} as the shifts vector. Note that supplying 
#' \code{rateShifts} when \code{rate} is not a numeric vector of the same 
#' length will result in an error.
#'
#' @return A constant or time-varying function (depending on input) that can
#' then be used as a rate in the other \code{paleobuddy} functions. 
#'
#' @author Bruno do Rosario Petrucci
#' 
#' @references 
#' 
#' Morlon H. et al (2016) RPANDA: an R package for macroevolutionary analyses on 
#' phylogenetic trees. \emph{Methods in Ecology and Evolution} 7: 589-597.
#'
#' @examples
#' 
#' # first we need a time vector to use on plots
#' time <- seq(0, 50, 0.1)
#' 
#' ###
#' # we can have a step function rate
#' 
#' # vector of rates
#' rate <- c(0.1, 0.2, 0.3, 0.2)
#' 
#' # vector of rate shifts
#' rateShifts <- c(0, 10, 20, 35)
#' # this could be c(50, 40, 30, 15) for equivalent results
#' 
#' # make the rate
#' r <- make.rate(rate, tMax = 50, rateShifts = rateShifts)
#' 
#' # plot it
#' plot(time, rev(r(time)),type = 'l', xlim = c(max(time), min(time)))
#' 
#' # note that this method of generating a step function rate is slower to
#' # numerically integrate
#' 
#' # it is also not possible a rate and a shifts vector and a time-series 
#' # dependency, so in cases where one looks to run many simulations, or have a
#' # step function modified by an environmental variable, consider 
#' # using ifelse() (see below)
#' 
#' ###
#' # we can have an environmental variable (or any time-series)
#' 
#' # temperature data
#' data(temp)
#' 
#' # function
#' rate <- function(t, env) {
#'   return(0.05*env)
#' }
#' 
#' # make the rate
#' r <- make.rate(rate, envRate = temp)
#' 
#' # plot it
#' plot(time, rev(r(time)), type = 'l', xlim = c(max(time), min(time)))
#' 
#' ###
#' # we can have a rate that depends on time AND temperature
#' 
#' # temperature data
#' data(temp)
#' 
#' # function
#' rate <- function(t, env) {
#'   return(0.001*exp(0.1*t) + 0.05*env)
#' }
#' 
#' # make a rate
#' r <- make.rate(rate, envRate = temp)
#' 
#' # plot it
#' plot(time, rev(r(time)), type = 'l', xlim = c(max(time), min(time)))
#' 
#' ###
#' # as mentioned above, we could also use ifelse() to 
#' # construct a step function that is modulated by temperature
#' 
#' # temperature data
#' data(temp)
#' 
#' # function
#' rate <- function(t, env) {
#'   return(ifelse(t < 10, 0.1 + 0.01*env,
#'                 ifelse(t < 30, 0.2 - 0.005*env,
#'                        ifelse(t <= 50, 0.1 + 0.005*env, 0))))
#' }
#' 
#' # rate
#' r <- make.rate(rate, envRate = temp)
#' 
#' # plot it
#' plot(time, rev(r(time)), type = 'l', xlim = c(max(time), min(time)))
#' 
#' # while using ifelse() to construct a step function is more
#' # cumbersome, it leads to much faster numerical integration,
#' # so in cases where the method above is proving too slow,
#' # consider using ifelse() even if there is no time-series dependence
#' 
#' ###
#' # make.rate will leave some types of functions unaltered
#' 
#' # constant rates
#' r <- make.rate(0.5)
#' 
#' # plot it
#' plot(time, rep(r, length(time)), type = 'l', 
#'      xlim = c(max(time), min(time)))
#' 
#' ###
#' # linear rates
#' 
#' # function
#' rate <- function(t) {
#'   return(0.01*t)
#' }
#' 
#' # create rate
#' r <- make.rate(rate)
#' 
#' # plot it
#' plot(time, rev(r(time)), type = 'l', xlim = c(max(time), min(time)))
#' 
#' ###
#' # any time-varying function, really
#' 
#' # function
#' rate <- function(t) {
#'   return(abs(sin(t))*0.1 + 0.05)
#' }
#' 
#' # create rate
#' r <- make.rate(rate)
#' 
#' # plot it
#' plot(time, r(time), type = 'l')
#'
#' @name make.rate
#' @rdname make.rate
#' @export

make.rate <- function(rate, tMax = NULL, envRate = NULL, rateShifts = NULL) {
  # use this to check for length and how many arguments
  nargs <- ifelse(is.numeric(rate), length(rate), length(formals(rate)))

  # let us first check for some errors
  
  # if the given rate is a number or vector of numbers
  if (is.numeric(rate)) {
    # if rate is constant, we should not see any envRate or rateShifts
    if (nargs == 1) {
      if (!is.null(envRate) || !is.null(rateShifts)) {
        stop("Constant rate with environmental variable or shifts")
      } 
      
      # if it is just the number, return it
      else {
        return(rate)
      }
    }

    # if length(rate) > 1 we have a rates vector, so we 
    # must have a shifts vector
    else if (is.null(rateShifts)) {
      stop("Rate vector supplied without a shift vector")
    }

    # if they are not the same size, we have a problem
    else if (length(rate) != length(rateShifts)) {
      stop("Rate vector and shifts vector must have the same length")
    }

    # if we have a rates vector and shifts vector, should not have envRate
    else if (!is.null(envRate)) {
      stop("Rates and shifts supplied with environmental variable;
           use ifelse()")
    }
  }

  # if we have a rates vector, we need a shifts vector
  else if (!is.null(rateShifts)) {
    stop("Shifts vector supplied without a rates vector")
  }

  # some sanity checks on the number of arguments
  else if (nargs > 1 || (nargs == 0 && !is.numeric(rate))) {
    # rate should not have more than two arguments or less than 1
    if (nargs > 2 || nargs == 0) {
      stop("Function can only depend on time and an environmental variable")
    }

    # if it has two, we must also have a non-null envRate
    else if (nargs == 2) {
      if (is.null(envRate)) {
        stop("Environmental function supplied with no environmental variable")
      }
      
      # and envRate needs to have two columns
      if (ncol(envRate) != 2) {
        stop("Environmental variable must be a dataframe with two columns")
      }
    }
  }

  # need a function of environmental var if we have said var
  else if (nargs == 1 && !is.null(envRate)) {
    stop("Environmental variable supplied with one argument function")
  }

  # check if there are shifts - i.e. if the rate is a step function
  if (!is.null(rateShifts)) {
    # check that we have tMax in case we want to invert time
    if (is.null(tMax)) {
      stop("Need tMax for creating step functions")
    }

    fList <- rate

    # if user gave a vector from past to present, invert it
    if (rateShifts[2] < rateShifts[1]) {
      rateShifts <- tMax - rateShifts
    }
    
    if (min(rateShifts) != 0) {
      rateShifts[which(rateShifts == min(rateShifts))] <- 0
    }

    # create the step function
    r <- stepfun(rateShifts[2:length(rateShifts)], fList)

    # vectorize the function so we can integrate it
    r <- Vectorize(r)
  }

  # if we want it to be dependent on environmental variables
  else if (!is.null(envRate)) {
    # find degrees of freedom
    df <- smooth.spline(x = envRate[, 1], envRate[, 2])$df

    # now that we have the degrees of freedom, perform the spline
    spline_result <- smooth.spline(envRate[, 1], envRate[, 2], df = df)

    # use predict to find the rate at all times
    envFunc <- function(t) {
      ifelse(t < max(envRate[, 1]), predict(spline_result, t)$y,
             predict(spline_result, max(envRate[, 1]))$y)
    }

    # make it a function of time only
    r <- function(t) {
      return(rate(t, envFunc(t)))
    }
  }

  # otherwise it is a normal function of time, so return itself
  else {
    r <- rate
  }

  return(r)
}
