#' Returns a rate based on a time-varying function, environmental variable
#' and/or vectors of rates and shifts
#'
#' \code{MakeRate} takes a function \code{ff}, which could be a constant, a
#' function of time or a vector of rates. If it is a constant or a time-varying
#' function, nothing else need be supplied. Otherwise, if \code{ff} is a vector
#' of rates, the user must supply the accompanying vector of rate shifts
#' \code{fShifts} to create a step function rate. Finally, if \code{ff} takes
#' an environmental variable as well, the user must supply a dataframe of time
#' vs. the environmental param, \code{envF}. Note that an error is
#' returned if the user suplies \code{fShifts} AND \code{envF}. If one wants a
#' step function modified by an environmental variable, use \code{ifelse} to
#' give \code{ff(t, env)} (see examples below).
#'
#' @param ff the baseline function with which to make the rate.
#' It can be a
#'
#' \describe{
#' \item{\code{Constant}}{For constant birth-death rates}
#'
#' \item{\code{Function of time}}{For rates that vary with time. Note that this
#' can be any function of time, but one should not supply a function that
#' depends on more than one variable without an accompanying \code{envF} -
#' that will result in an error}
#'
#' \item{\code{Vector of rates}}{To create step function rates. Note this must
#' be accompanied by a corresponding vector of shifts \code{fShifts}}}
#'
#' @param tMax a number corresponding to the maximum simulation time.
#' Needed to ensure \code{fShifts} runs the correct way.
#'
#' @param envF a dataframe representing an environmental variable
#' (time, CO2 etc) with time. The first column must be time, second column the
#' values of the variable. See below; one good resource for these dataframes is
#' \href{https://cran.r-project.org/web/packages/RPANDA/}{RPANDA}.
#'
#' @param fShifts a vector of rate shifts. The first element must
#' be the first time point for the simulation. This may be 0 or tMax. Since
#' functions in paleobuddy run from 0 to tMax, if \code{fShifts} runs from past
#' to present, in other words \code{fShifts[2] < fShifts[1]}, we take
#' \code{tMax-fShifts} as the shifts vector. Note that supplying \code{fShifts}
#' when \code{ff} is not a rates vector will return an error.
#'
#' @return returns a constant or time-varying function (depending on input)
#' that can then be used as a rate in the other \code{paleobuddy} functions.
#' The returned function will invariably be either a number or a function of
#' one variable only, usually set as time.
#'
#' @author written by Bruno do Rosario Petrucci
#'
#' @examples
#'
#' # let us start simple: create a constant rate
#' r <- MakeRate(0.5)
#' plot(1:50, rep(r, 50), type='l')
#'
#' # something a bit more complex: a linear rate
#' ff <- function(t) {
#'   return(0.01*t)
#' }
#' r <- MakeRate(ff)
#' plot(1:50, r(1:50), type='l')
#'
#' # remember: this can be any time-varying function!
#' ff <- function(t) {
#'   return(sin(t)*0.01)
#' }
#' r <- MakeRate(ff)
#' plot(1:50, r(1:50), type='l')
#'
#' # we can use ifelse() to make a step function like this
#' ff <- function(t) {
#'   return(ifelse(t < 10, 0.1,
#'                 ifelse(t < 20, 0.3,
#'                        ifelse(t < 30, 0.2,
#'                               ifelse(t < 40, 0.05, 0)))))
#' }
#' r <- MakeRate(ff)
#' plot(1:50, r(1:50), type='l')
#'
#' # important note: this method of creating a step function might be annoying,
#' # but when running thousands of simulations it will provide a much faster
#' # integration than when using our method of transforming a rates and shifts
#' # vector into a function of time...
#'
#' # ...which we can do as follows
#' ff <- c(0.1, 0.2, 0.3, 0.2)
#' fShifts <- c(0, 10, 20, 35)
#' r <- MakeRate(ff, fShifts = fShifts)
#' plot(1:50, r(1:50),type='l')
#'
#' # as mentioned above, while this works well it will be a pain to integrate.
#' # Furthermore, it is impractical to supply a rate and a shifts vector and
#' # have an environmental dependency, so in cases where one looks to run
#' # more than a couple dozen simulations, and when one is looking to have a
#' # step function modified by an environmental variable, consider using ifelse()
#'
#' # finally let us see what we can do with environmental variables
#'
#' # RPANDA supplies us with some really useful environmental dataframes
#' # to use as an example, let us try temperature
#' if (requireNamespace("RPANDA", quietly = TRUE)) {
#'   data(InfTemp, package="RPANDA")
#'
#'   ff <- function(t, env) {
#'     return(0.05*env)
#'   }
#'   r <- MakeRate(ff, envF = InfTemp)
#'   plot(1:50, r(1:50), type='l')
#' }
#' # we can also have a function that depends on both time AND temperature
#' ff <- function(t, env) {
#'   return(0.001*exp(0.1*t) + 0.05*env)
#' }
#' r <- MakeRate(ff, envF = InfTemp)
#' plot(1:50, r(1:50), type='l')
#'
#' # as mentioned above, we could also use ifelse() to construct a step function
#' # that is modulated by temperature
#' ff <- function(t, env) {
#'   return(ifelse(t < 10, 0.1 + 0.01*env,
#'                 ifelse(t < 30, 0.2 - 0.005*env,
#'                        ifelse(t <= 50, 0.1 + 0.005*env, 0))))
#' }
#' r <- MakeRate(ff, envF = InfTemp)
#' plot(1:50, r(1:50), type='l')
#'
#' @name MakeRate
#' @rdname MakeRate
#' @export

MakeRate<-function(ff,tMax, env_f=NULL,fShifts=NULL) {
  # may use this soon
  nargs = ifelse(is.numeric(ff), length(ff), length(formals(ff)))

  # let us first check for some errors
  if (is.numeric(ff)) {
    # if ff is constant, we should not see any env_f or fShifts
    if (nargs == 1) {
      if (!is.null(env_f) || !is.null(fShifts)) {
        stop("constant rate with environmental variable or shifts")
      } else {
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

    # if we have a rates vector and shifts vector, should not have env_f
    else if (!is.null(env_f)) {
      stop("rates and shifts supplied with environmental variable;
           use ifelse()")
    }
  }


  else if (!is.null(fShifts)) {
    stop("shifts vector supplied without a rates vector")
  }

  else if (nargs > 1 || (nargs == 0 && !is.numeric(ff))) {
    # ff should not have more than two arguments or less than 1
    if (nargs > 2 || nargs == 0) {
      stop("function can only depend on time and environmental var")
    }

    # if it has two, we must also have a non-null env_f
    else if (nargs == 2) {
      if (is.null(env_f)) {
        stop("environmental function supplied with no environmental variable")
      }
    }
  }

  else if (nargs == 1 && !is.null(env_f)) {
    stop("environmental variable supplied with one argument function")
  }

  # check if there are shifts - i.e. if the rate is a step function
  if (!is.null(fShifts)) {

    fList<-ff

    # if user gave a list from past to present, make it from present to past
    if (fShifts[2] < fShifts[1]) {
      fShifts<-tMax - fShifts
    }

    # create the step function
    f<-function(t) {

      if (t<0) {
        return(0)
      }

      else {
        # get the rate for this time by subtracting the shifts
        # and finding where the subtraction is positive
        return(fList[utils::tail(which(t-fShifts>=0), n=1)])
      }
    }

    # vectorize the function so we can integrate it
    f <- Vectorize(f)
  }

  # if we want it to be dependent on environmental variables
  else if (!is.null(env_f)) {
    # find degrees of freedom
    df <- smooth.spline(x=env_f[,1], env_f[,2])$df

    # now that we have the degrees of freedom, perform the spline
    spline_result <- smooth.spline(env_f[,1],env_f[,2], df=df)

    # use predict to find the rate at all times
    env_func <- function(t) {
      predict(spline_result,t)$y
    }

    # make it a function of time only
    f<-function(t) {
      return(ff(t,env_func(t)))
    }
  }

  # otherwise it is either constant or a function of time, so return itself
  else {
    f <- ff
  }

  return(f)
}
