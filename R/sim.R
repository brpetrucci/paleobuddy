#' Details, generics, and methods for the \code{sim} class
#' 
#' @description The \code{sim} class is a frequent return and input argument for functions in
#' paleobuddy. It contains the following four elements.
#' 
#' \describe{
#' \item{\code{TE}}{List of extinction times, with \code{NA} as the time of
#' extinction for extant species.}
#'
#' \item{\code{TS}}{List of speciation times, with \code{tMax} as the time of
#' speciation for species that started the simulation.}
#'
#' \item{\code{PAR}}{List of parents. Species that started the simulation have
#' \code{NA}, while species that were generated during the simulation have their
#' parent's number. Species are numbered as they are born.}
#'
#' \item{\code{EXTANT}}{List of booleans representing whether each species is
#' extant.}}
#' 
#' Here we declare useful generics and methods for \code{sim} objects.
#' 
#' @param sim,x,object Object of class "sim"
#' 
#' @param t Time t (in Mya). Used for counting and/or plotting births, deaths
#' and species number.
#' 
#' @param ... Further arguments inherited from generics.
#'
#' @name sim
#' 
#' @importFrom graphics plot par
#' @importFrom utils head
#' 
NULL

#' @rdname sim
#' 
#' @details \code{is.sim} A \code{sim} object must contain 4 members (usually 
#' vectors for extinction times, speciation times, species' parents and status), 
#' and all of these must have the correct length (i.e. same as all the others) and
#' types. We do not utilize the members' order inside \code{sim} for our tests, 
#' since they are accessed with the $ operator and therefore the order is 
#' irrelevant.
#' 
#' @export
#' 

is.sim <- function(sim) {
  # checks that the class is sim
  cla <- class(sim) == "sim"
  
  # if it is not, should be
  if (!cla) {
    return(FALSE)
  }
  
  # checks that there are four members
  len <- (length(sim) == 4)
  
  # if there are not, the other tests are redundant
  if (!len) {
    return(FALSE)
  }
  
  # checks that the members have same size
  siz <- (length(unique(c(length(sim[[1]]), length(sim[[2]]), length(sim[[3]]),
                    length(sim[[4]])))) == 1)
  
  # checks that there are three double vectors and one logical, or two of each
  # when there is only one species
  types <- unlist(lapply(1:4, function(x) typeof(sim[[x]])))
  typ <- (sum(types == "double") == 3 && sum(types == "logical") == 1) ||
    (sum(types == "double") == 2 && sum(types == "logical") == 2 &&
       length(sim[[1]]) == 1)
  

  # check that, if typ is false, it is because there are NA-only 
  return(siz && typ)
}

#' @rdname sim
#' 
#' @details \code{print.sim} The printing of a sim object is formatted into a more
#' straightforward and informative sequence manner. We provide details only for 
#' the first few species, since otherwise this print could be overwhelming for 
#' simulations with 10+ species.
#' 
#' @export
#' 

print.sim <- function(x, ...) {
  # change name just for clarity of the object
  sim <- x
  
  # first check that it is a valid sim
  if (!is.sim(sim)) {
    stop("Invalid sim object, see ?sim")
  }
  
  # make the extant info more straightforward
  status <- rep("extinct", length(sim$EXTANT))
  status[sim$EXTANT] <- "extant"
  
  # then print some basic information about the sim
  cat(paste0("\nBirth-death simulation object with ", length(sim$TE), 
             " species and ", sum(sim$EXTANT), " extant species\n"))
  cat("\nDetails for some species:\n")
  
  # and then some details for first five
  cat("\nExtinction times (NA means extant)\n")
  print(utils::head(sim$TE))
  
  cat("\n\nSpeciation times \n")
  print(utils::head(sim$TS))
  
  cat("\n\nSpecies parents (NA for initial)\n")
  print(utils::head(sim$PAR))
  
  cat("\n\nSpecies status (extinct or extant)\n")
  print(utils::head(status))
  
  # to see the whole vector
  cat("\n\nFor more details on vector y, try sim$y, with y one of\n")
  cat(names(sim))
}

#' @rdname sim
#' 
#' @details \code{summary.sim} Quantitative details on the \code{sim} object. 
#' Prints the number of species, number of extant species, summary of durations
#' and speciation waiting times, in case there are more than one species.
#' 
#' @export
#'

summary.sim <- function(object, ...) {
  # change name just for clarity of the object
  sim <- object
  
  # check that it is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid sim object, see ?sim")
  }
  
  # state the name and nature of the object
  cat("\nBirth-death simulation object: ", deparse(substitute(object)))
  
  # some quick numbers
  cat("\n\n  Total number of species: ", length(sim$TE))
  cat("\n  Number of extant species: ", sum(sim$EXTANT))
  
  # get durations (for extinct species)
  dur <- sim$TS[!sim$EXTANT] - sim$TE[!sim$EXTANT]
  
  # summarize durations
  if (length(dur) > 0) {
    cat("\n  Durations (for extinct species):")
    cat("\n    mean: ", mean(dur))
    cat("\n    standard deviation: ", sd(dur))
    cat("\n    summary:\n")
    print(summary(dur)[-4])
  }
  
  else {
    cat("\n  No extinct species, so no duration information.")
  }
  
  # get speciation waiting times
  spec <- sim$TS[sim$PAR][2:length(sim$TS)] - sim$TS[2:length(sim$TS)]
  
  # summarize time to speciation
  if (length(spec) > 0) {
    cat("\n  Times to speciation:")
    cat("\n    mean: ", mean(spec))
    cat("\n    standard deviation: ", sd(spec))
    cat("\n    summary:\n")
    print(summary(spec)[-4])
  }
  
  else {
    cat("\n  Only one species, so no time to speciation information.")
  }
}

#' @rdname sim
#' 
#' @details \code{plot.sim.ltt} Plots births, deaths, and diversity through time for
#' the sim object.
#' 
#' @export
#' 

plot.sim.ltt <- function(x, ...) {
  # change name just for clarity of the object
  sim <- x
  
  # check that it is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid sim object, see ?sim")
  }
  
  # set up three plots
  par(mfrow = c(3, 1))
  
  # make TE sensible
  sim$TE[sim$EXTANT] <- 0
  
  # create a time vector
  t <- seq(max(sim$TS) + 0.1, min(sim$TE) - 0.1, -0.001)
  
  # count births, deaths, and diversity
  counts <- sim.counts(sim, t)
  
  # plot births
  plot(t, counts$births, xlim = c(max(t), 0), ylim = c(0, length(sim$TE) + 0.2), 
       type = 'l', main = "Speciation number through time",
       ylab = "Births", xlab = "Time (mya)")
  
  # deaths
  plot(t, counts$deaths, xlim = c(max(t), 0), ylim = c(0, sum(!sim$EXTANT) + 0.2),
       type = 'l', main = "Extinction number through time", 
       ylab = "Deaths", xlab = "Time (mya)")
  
  # and diversity
  plot(t, counts$div, xlim = c(max(t), 0), ylim = c(0, max(counts$div) + 0.2), 
       type = 'l', main = "Species number through time", 
       ylab = "Species", xlab = "Time (mya)")
}

#' @rdname sim
#' 
#' @details \code{sim.counts} Calculates the births, deaths, and diversity for a 
#' sim at time t.
#' 
#' @export
#' 

sim.counts <- function(sim, t) {
  # check that it is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid sim object, see ?sim")
  }
  
  # make TE -Inf if extant (so that the number of extinctions
  # does not artificially jump at the ending time of sim)
  sim$TE[sim$EXTANT] <- -Inf
  
  # calculates births at t
  births <- unlist(lapply(t, function(t) sum(sim$TS >= t)))
  
  # calculates deaths at t
  deaths <- unlist(lapply(t, function(t) sum(sim$TE >= t)))
  
  # calculates diversity at t
  div <- unlist(lapply(t, function(t) sum(sim$TS >= t & sim$TE <= t)))
  
  return(list(births = births, deaths = deaths, div = div))
}
