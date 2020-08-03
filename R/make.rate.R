#' Create a flexible rate for birth-death or sampling simulations
#' 
#' Generates a rate to be used on birth-death or sampling functions. Takes as the 
#' base rate (1) a constant, (2) a function of time, (3) a function of time 
#' interacting with an environmental variable, or (4) a vector of numbers 
#' describing rates as a step function. Requires information regarding the maximum
#' simulation time, and allows for optional extra parameters to tweak the baseline
#' rate.
#'
#' @param ff The baseline function with which to make the rate.
#' It can be a
#'
#' \describe{
#' \item{A number}{For constant birth-death rates.}
#'
#' \item{A function of time}{For rates that vary with time. Note that this can be
#' any function of time.} 
#' 
#' \item{A function of time and an environmental variable}{For rates varying with
#' time and an environmental variable, such as temperature. Note that supplying a
#' function on more than one variable without an accompanying \code{envF} will
#' result in an error.}
#'
#' \item{A numeric vector}{To create step function rates. Note this must be
#'  accompanied by a corresponding vector of shifts \code{fShifts}.}}
#'
#' @param tMax Ending time of simulation, in million years after the clade's 
#' origin. Needed to ensure \code{fShifts} runs the correct way.
#'
#' @param envF A \code{data.frame} representing the variation of an environmental 
#' variable (e.g. CO2, temperature, available niches, etc) with time. The first 
#' column of this \code{data.frame} must be time, and the second column must be 
#' the values of the variable. The function will return an error if supplying 
#' \code{envF} without \code{ff} being a function of two variables. Note 
#' \code{paleobuddy} has two environmental data frames, \code{temp} and 
#' \code{co2}. One can check \code{RPANDA} for more examples.
#' 
#' Acknowledgements: The strategy to transform a function of \code{t} and 
#' \code{env} into a function of \code{t} only using \code{envF} was adapted from
#' RPANDA (see below).
#'
#' @param fShifts A vector indicating the time placement of rate shifts in a step 
#' function. The first element must be the first time point for the simulation. 
#' This may be \code{0} or \code{tMax}. Since functions in paleobuddy run from 
#' \code{0} to \code{tMax}, if \code{fShifts} runs from past to present 
#' (\code{fShifts[2] < fShifts[1]}), we take \code{tMax - fShifts} as the shifts 
#' vector. Note that supplying \code{fShifts} when \code{ff} is not a numeric 
#' vector of the same length will result in an error.
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
#' t <- seq(0, 50, 0.1)
#' 
#' # make.rate will leave some types of functions unaltered, like the following
#' 
#' ###
#' # let us start simple: create a constant rate
#' r <- make.rate(0.5)
#' 
#' # plot it
#' plot(t, rep(r, length(t)), type = 'l')
#' 
#' ###
#' # something a bit more complex: a linear rate
#' 
#' # function
#' ff <- function(t) {
#'   return(0.01*t)
#' }
#' 
#' # create rate
#' r <- make.rate(ff)
#' 
#' # plot it
#' plot(t, r(t), type = 'l')
#' 
#' ###
#' # remember: this can be any time-varying function!
#' 
#' # function
#' ff <- function(t) {
#'   return(sin(t)*0.01)
#' }
#' 
#' # create rate
#' r <- make.rate(ff)
#' 
#' # plot it
#' plot(t, r(t), type = 'l')
#' 
#' ###
#' # we can use ifelse() to make a step function like this
#' ff <- function(t) {
#'   return(ifelse(t < 10, 0.1,
#'                 ifelse(t < 20, 0.3,
#'                        ifelse(t < 30, 0.2,
#'                               ifelse(t < 40, 0.05, 0)))))
#' }
#' 
#' # and make it into a rate - in this case, as the previous, it does not alter
#' # ff. We put it here as a contrast to the other way to make a step function
#' r <- make.rate(ff)
#' 
#' # plot it
#' plot(t, r(t), type = 'l')
#' 
#' # important note: this method of creating a step function might be annoying,
#' # but when running thousands of simulations it will provide a much faster
#' # integration than when using our method of transforming a rates and shifts
#' # vector into a function of time
#' 
#' # this is a good segway into the cases where make.rate actually makes a rate!
#' # note that while the previous ones seemed useless, we need that implementation
#' # so that the birth-death functions work
#' 
#' ###
#' # now we can demonstrate the other way of making a step function
#' 
#' # vector of rates
#' ff <- c(0.1, 0.2, 0.3, 0.2)
#' 
#' # vector of rate shifts
#' fShifts <- c(0, 10, 20, 35)
#' # this could be c(50, 40, 30, 15) for equivalent results
#' 
#' # make the rate
#' r <- make.rate(ff, tMax = 50, fShifts = fShifts)
#' 
#' # plot it
#' plot(t, r(t),type = 'l')
#' 
#' # as mentioned above, while this works well it will be a pain to integrate.
#' # Furthermore, it is impractical to supply a rate and a shifts vector and
#' # have an environmental dependency, so in cases where one looks to run
#' # more than a couple dozen simulations, or when one is looking to have a
#' # step function modified by an environmental variable, consider using ifelse()
#' 
#' ###
#' # finally let us see what we can do with environmental variables
#' 
#' # temperature data
#' data(temp)
#' 
#' # function
#' ff <- function(t, env) {
#'   return(0.05*env)
#' }
#' 
#' # make the rate
#' r <- make.rate(ff, envF = temp)
#' 
#' # plot it
#' plot(t, r(t), type = 'l')
#' 
#' ###
#' # we can also have a function that depends on both time AND temperature
#' 
#' # function
#' ff <- function(t, env) {
#'   return(0.001*exp(0.1*t) + 0.05*env)
#' }
#' 
#' # make a rate
#' r <- make.rate(ff, envF = temp)
#' 
#' # plot it
#' plot(t, r(t), type = 'l')
#' 
#' ###
#' # as mentioned above, we could also use ifelse() to construct a step function
#' # that is modulated by temperature
#' 
#' # function
#' ff <- function(t, env) {
#'   return(ifelse(t < 10, 0.1 + 0.01*env,
#'                 ifelse(t < 30, 0.2 - 0.005*env,
#'                        ifelse(t <= 50, 0.1 + 0.005*env, 0))))
#' }
#' 
#' # rate
#' r <- make.rate(ff, envF = temp)
#' 
#' # plot it
#' plot(t, r(t), type = 'l')
#'
#' @name make.rate
#' @rdname make.rate
#' @export

make.rate <- function(ff, tMax = NULL, envF = NULL, fShifts = NULL) {
  # use this to check for length and how many arguments
  nargs = ifelse(is.numeric(ff), length(ff), length(formals(ff)))

  # let us first check for some errors
  
  # if the given rate is a number or vector of numbers
  if (is.numeric(ff)) {
    # if ff is constant, we should not see any envF or fShifts
    if (nargs == 1) {
      if (!is.null(envF) || !is.null(fShifts)) {
        stop("constant rate with environmental variable or shifts")
      } 
      
      # if it is just the number, return it
      else {
        return(ff)
      }
    }

    # if length(ff) > 1 we have a rates vector, so we must have a shifts vector
    else if (is.null(fShifts)) {
      stop("rate vector supplied without a shift vector")
    }

    # if they are not the same size, we have a problem
    else if (length(ff) != length(fShifts)) {
      stop("rate vector and shifts vector must have the same length")
    }

    # if we have a rates vector and shifts vector, should not have envF
    else if (!is.null(envF)) {
      stop("rates and shifts supplied with environmental variable;
           use ifelse()")
    }
  }

  # if we have a rates vector, we need a shifts vector
  else if (!is.null(fShifts)) {
    stop("shifts vector supplied without a rates vector")
  }

  # some sanity checks on the number of arguments
  else if (nargs > 1 || (nargs == 0 && !is.numeric(ff))) {
    # ff should not have more than two arguments or less than 1
    if (nargs > 2 || nargs == 0) {
      stop("function can only depend on time and environmental var")
    }

    # if it has two, we must also have a non-null envF
    else if (nargs == 2) {
      if (is.null(envF)) {
        stop("environmental function supplied with no environmental variable")
      }
      
      # and envF needs to have two columns
      if (ncol(envF) != 2) {
        stop("environmental variable must be a dataframe with two columns")
      }
    }
  }

  # need a function of environmental var if we have said var
  else if (nargs == 1 && !is.null(envF)) {
    stop("environmental variable supplied with one argument function")
  }

  # check if there are shifts - i.e. if the rate is a step function
  if (!is.null(fShifts)) {
    # check that we have tMax in case we want to invert time
    if (is.null(tMax)) {
      stop("need tMax for creating step functions")
    }

    fList <- ff

    # if user gave a vector from past to present, make it from present to past
    if (fShifts[2] < fShifts[1]) {
      fShifts <- tMax - fShifts
    }

    # create the step function
    f <- function(t) {

      if (t < 0) {
        return(0)
      }

      else {
        # get the rate for this time by subtracting the shifts
        # and finding where the subtraction is positive
        return(fList[utils::tail(which(t - fShifts >= 0), n = 1)])
      }
    }

    # vectorize the function so we can integrate it
    f <- Vectorize(f)
  }

  # if we want it to be dependent on environmental variables
  else if (!is.null(envF)) {
    # find degrees of freedom
    df <- smooth.spline(x = envF[, 1], envF[, 2])$df

    # now that we have the degrees of freedom, perform the spline
    spline_result <- smooth.spline(envF[, 1], envF[, 2], df = df)

    # use predict to find the rate at all times
    envFunc <- function(t) {
      predict(spline_result, t)$y
    }

    # make it a function of time only
    f <- function(t) {
      return(ff(t, envFunc(t)))
    }
  }

  # otherwise it is either constant or a function of time, so return itself
  else {
    f <- ff
  }

  return(f)
}
