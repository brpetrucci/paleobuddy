#' MuSSE simulation
#'
#' Simulates a species birth-death process following the Multiple 
#' State-dependent Speciation and Extinction (MuSSE) or the Hidden 
#' State-dependent Speciation and Extinction (HiSSE) model for any number of 
#' starting species. Allows for the speciation/extinction rate to be (1) a 
#' constant, or (2) a list of values for each trait state. Traits are simulated 
#' to evolve under a simple Mk model (see references). Results can be 
#' conditioned on either total simulation time, or total number of extant 
#' species at the end of the simulation. Also allows for constraining results 
#' on a range of number of species at the end of the simulation, either total 
#' or extant, using rejection sampling. Returns a \code{sim} object 
#' (see \code{?sim}), and a list of data frames describing trait values for 
#' each interval. It may  return true extinction times or simply information 
#' on whether species lived after the maximum simulation time, depending on 
#' input. Can simulate any number of traits, but rates need to depend on only
#' one (each, so speciation and extinction can depend on different traits). 
#' 
#' Please note while time runs from \code{0} to \code{tMax} in the simulation, 
#' it returns speciation/extinction times as \code{tMax} (origin of the group) 
#' to \code{0} (the "present" and end of simulation), so as to conform to other
#' packages in the literature.
#'
#' @inheritParams bd.sim
#' 
#' @param lambda Vector to hold the speciation rate over time. It should either
#' be a constant, or a list of size \code{nStates}. For each species a trait
#' evolution simulation will be run, and then used to calculate the final 
#' speciation rate. Note that \code{lambda} should always be greater than or 
#' equal to zero.
#'
#' @param mu Similar to above, but for the extinction rate.
#' 
#' @param nTraits The number of traits to be considered. \code{lambda} and 
#' \code{mu} need not reference every trait simulated.
#' 
#' @param nFocus Trait of focus, i.e. the one that rates depend on. If it is 
#' one number, that will be the trait of focus for both speciation and 
#' extinction rates. If it is of length 2, the first will be the focus for
#' the former, the second for the latter.
#' 
#' @param nStates Number of possible states for categorical trait. The range
#' of values will be assumed to be \code{(0, nStates - 1)}. Can be a constant
#' or a vector of length \code{nTraits}, if traits are intended to have 
#' different numbers of states.
#' 
#' @param nHidden Number of hidden states for categorical trait. Default is 
#' \code{1}, in which case there are no added hidden traits. Total number of
#' states is then \code{nStates * nHidden}. States will then be set to a value
#' in the range of \code{(0, nStates - 1)} to simulate that hidden states are
#' hidden. This is done by setting the value of a state to the remainder of
#' \code{state / nStates}. E.g. if \code{nStates = 2} and \code{nHidden = 3},
#' possible states during simulation will be in the range \code{(0, 5)}, but
#' states \code{(2, 4)} (corresponding to \code{(0B, 0C)} in the nomenclature
#' of the original HiSSE reference) will be set to \code{0}, and states 
#' \code{(3, 5)} (corresponding to \code{(1B, 1C)}) to \code{1}.
#' 
#' @param X0 Initial trait value for original species. Must be within 
#' \code{(0, nStates - 1)}. Can be a constant or a vector of length 
#' \code{nTraits}. 
#' 
#' @param Q Transition rate matrix for continuous-time trait evolution. For
#' different states \code{i} and \code{j}, the rate at which a species at
#' \code{i} transitions to \code{j} is \code{Q[i + 1, j + 1]}. Must be within
#' a list, so as to allow for different \code{Q} matrices when
#' \code{nTraits > 1}.
#' 
#' Note that for all of \code{nStates}, \code{nHidden}, \code{X0} and \code{Q},
#' if \code{nTraits > 1} and any of those is of length \code{1}, they will be
#' considered to apply to all traits equally. This might lead to problems if,
#' e.g., two traits have different states but the same \code{Q}, so double
#' check that you are providing all parameters for the required traits.
#'
#' @return A \code{sim} object, containing extinction times, speciation times,
#' parent, and status information for each species in the simulation, and a 
#' list object with the trait data frames describing the trait value for each
#' species at each specified interval.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @references 
#' 
#' Maddison W.P., Midford P.E., Otto S.P. 2007. Estimating a binary character’s 
#' effect on speciation and extinction. Systematic Biology. 56(5):701.
#' 
#' FitzJohn R.G. 2012. Diversitree: Comparative Phylogenetic Analyses of 
#' Diversification in R. Methods in Ecology and Evolution. 3:1084–1092.
#' 
#' Beaulieu J.M., O'Meara, B.C. 2016. Detecting Hidden Diversification Shifts 
#' in Models of Trait-Dependent Speciation and Extinction. Systematic Biology.
#' 65(4):583-601.
#'
#' @examples
#'
#' ###
#' # first, it's good to check that it can work with constant rates
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation
#' lambda <- 0.1
#' 
#' # extinction
#' mu <- 0.03
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation, making sure we have more than one species in the end
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nFinal = c(2, Inf))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   ape::plot.phylo(phy)
#' }
#' 
#' ###
#' # now let's actually make it trait-dependent, a simple BiSSE model
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 40
#' 
#' # speciation, higher for state 1
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, trait-independent
#' mu <- 0.03
#' 
#' # number of traits and states (1 binary trait)
#' nTraits <- 1
#' nStates <- 2
#' 
#' # initial value of the trait
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical transition rates
#' Q <- list(matrix(c(0, 0.1,
#'                    0.1, 0), ncol = 2, nrow = 2))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits, 
#'                     nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = c("red", "blue")[traits + 1])
#' }
#' 
#' ###
#' # extinction can be trait-dependent too, of course
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # number of species at the end of the simulation
#' N <- 20
#' 
#' # speciation, higher for state 1
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, highe for state 0
#' mu <- c(0.06, 0.03)
#' 
#' # number of traits and states (1 binary trait)
#' nTraits <- 1
#' nStates <- 2
#' 
#' # initial value of the trait
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical transition rates
#' Q <- list(matrix(c(0, 0.1,
#'                    0.1, 0), ncol = 2, nrow = 2))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, N = N, nTraits = nTraits,
#'                     nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = c("red", "blue")[traits + 1])
#' }
#' 
#' ###
#' # we can complicate the model further by making transition rates asymmetric
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 20
#' 
#' # speciation, higher for state 1
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, lower for state 1
#' mu <- c(0.03, 0.01)
#' 
#' # number of traits and states (1 binary trait)
#' nTraits <- 1
#' nStates <- 2
#' 
#' # initial value of the trait
#' X0 <- 0
#' 
#' # transition matrix, with q01 higher than q10
#' Q <- list(matrix(c(0, 0.1,
#'                    0.25, 0), ncol = 2, nrow = 2))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits, 
#'                     nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = c("red", "blue")[traits + 1])
#' }
#' 
#' ###
#' # MuSSE is BiSSE but with higher numbers of states
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # number of species at the end of the simulation
#' N <- 20
#' 
#' # speciation, higher for state 1, highest for state 2
#' lambda <- c(0.1, 0.2, 0.3)
#' 
#' # extinction, higher for state 2
#' mu <- c(0.03, 0.03, 0.06)
#' 
#' # number of traits and states (1 trinary trait)
#' nTraits <- 1
#' nStates <- 3
#' 
#' # initial value of the trait
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical, fully reversible transition rates
#' Q <- list(matrix(c(0, 0.1, 0.1,
#'                    0.1, 0, 0.1,
#'                    0.1, 0.1, 0), ncol = 3, nrow = 3))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, N = N, nTraits = nTraits, 
#'                     nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # 0 tips = red, 1 tips = blue, 2 tips = green
#'   ape::plot.phylo(phy, tip.color = c("red", "blue", "green")[traits + 1])
#' }
#' 
#' ###
#' # HiSSE is like BiSSE, but with the possibility of hidden traits
#' # here we have 4 states, representing two states for the observed trait
#' # (0 and 1) and two for the hidden trait (A and B), i.e. 0A, 1A, 0B, 1B
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 20
#' 
#' # speciation, higher for state 1A, highest for 1B
#' lambda <- c(0.1, 0.2, 0.1, 0.3)
#' 
#' # extinction, lowest for 0B
#' mu <- c(0.03, 0.03, 0.01, 0.03)
#' 
#' # number of traits and states (1 binary observed trait, 
#' # 1 binary hidden trait)
#' nTraits <- 1
#' nStates <- 2
#' nHidden <- 2
#' 
#' # initial value of the trait
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical transition rates. Only one transition
#' # is allowed at a time, i.e. 0A can go to 0B and 1A,
#' # but not to 1B, and similarly for others
#' Q <- list(matrix(c(0, 0.1, 0.1, 0,
#'                    0.1, 0, 0, 0.1,
#'                    0.1, 0, 0, 0.1,
#'                    0, 0.1, 0.1, 0), ncol = 4, nrow = 4))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits,
#'                     nStates = nStates, nHidden = nHidden,
#'                     X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#' 
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = c("red", "blue")[traits + 1])
#' }
#' 
#' ###
#' # we can also increase the number of traits, e.g. to have a neutral trait 
#' # evolving with the real one to compare the estimates of the model for each
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 20
#' 
#' # speciation, higher for state 1
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, lowest for state 0
#' mu <- c(0.01, 0.03)
#' 
#' # number of traits and states (2 binary traits)
#' nTraits <- 2
#' nStates <- 2
#' 
#' # initial value of both traits
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical transition rates for trait 1,
#' # and asymmetrical (and higher) for trait 2
#' Q <- list(matrix(c(0, 0.1,
#'                    0.1, 0), ncol = 2, nrow = 2),
#'           matrix(c(0, 1,
#'                    0.5, 0), ncol = 2, nrow = 2))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits, 
#'                     nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits1 <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' traits2 <- unlist(lapply(sim$TRAITS, function(x) tail(x[[2]]$value, 1)))
#' 
#' # make index for coloring tips
#' index <- ifelse(!(traits1 | traits2), "red", 
#'                 ifelse(traits1 & !traits2, "purple",
#'                        ifelse(!traits1 & traits2, "magenta", "blue")))
#' # 00 = red, 10 = purple, 01 = magenta, 11 = blue
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = index)
#' }
#' 
#' ###
#' # we can then do the same thing, but with the
#' # second trait controlling extinction
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 20
#' 
#' # speciation, higher for state 10 and 11
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, lowest for state 00 and 01
#' mu <- c(0.01, 0.03)
#' 
#' # number of traits and states (2 binary traits)
#' nTraits <- 2
#' nStates <- 2
#' nFocus <- c(1, 2)
#' 
#' # initial value of both traits
#' X0 <- 0
#' 
#' # transition matrix, with symmetrical transition rates for trait 1,
#' # and asymmetrical (and higher) for trait 2
#' Q <- list(matrix(c(0, 0.1,
#'                    0.1, 0), ncol = 2, nrow = 2),
#'           matrix(c(0, 1,
#'                    0.5, 0), ncol = 2, nrow = 2))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits, 
#'                     nStates = nStates, nFocus = nFocus,
#'                     X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits1 <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' traits2 <- unlist(lapply(sim$TRAITS, function(x) tail(x[[2]]$value, 1)))
#' 
#' # make index for coloring tips
#' index <- ifelse(!(traits1 | traits2), "red", 
#'                 ifelse(traits1 & !traits2, "purple",
#'                        ifelse(!traits1 & traits2, "magenta", "blue")))
#' # 00 = red, 10 = purple, 01 = magenta, 11 = blue
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = index)
#' }
#' 
#' ###
#' # as a final level of complexity, let us change the X0
#' # and number of states of the trait controlling extinction
#' 
#' # initial number of species
#' n0 <- 1
#' 
#' # maximum simulation time
#' tMax <- 20
#' 
#' # speciation, higher for state 10 and 11
#' lambda <- c(0.1, 0.2)
#' 
#' # extinction, lowest for state 00, 01, and 02
#' mu <- c(0.01, 0.03, 0.03)
#' 
#' # number of traits and states (2 binary traits)
#' nTraits <- 2
#' nStates <- c(2, 3)
#' nFocus <- c(1, 2)
#' 
#' # initial value of both traits
#' X0 <- c(0, 2)
#' 
#' # transition matrix, with symmetrical transition rates for trait 1,
#' # and asymmetrical, directed, and higher rates for trait 2
#' Q <- list(matrix(c(0, 0.1,
#'                    0.1, 0), ncol = 2, nrow = 2),
#'           matrix(c(0, 1, 0,
#'                    0.5, 0, 0.75,
#'                    0, 1, 0), ncol = 3, nrow = 3))
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run the simulation
#' sim <- bd.sim.traits(n0, lambda, mu, tMax, nTraits = nTraits, 
#'                     nStates = nStates, nFocus = nFocus,
#'                     X0 = X0, Q = Q, nFinal = c(2, Inf))
#' 
#' # get trait values for all tips
#' traits1 <- unlist(lapply(sim$TRAITS, function(x) tail(x[[1]]$value, 1)))
#' traits2 <- unlist(lapply(sim$TRAITS, function(x) tail(x[[2]]$value, 1)))
#' 
#' # make index for coloring tips
#' index <- ifelse(!(traits1 | (traits2 != 0)), "red", 
#'                 ifelse(traits1 & (traits2 == 0), "purple",
#'                        ifelse(!traits1 & (traits2 == 1), "magenta", 
#'                               ifelse(traits1 & (traits2 == 1), "blue",
#'                                      ifelse(!traits1 & (traits2 == 2), 
#'                                             "orange", "green")))))
#'                                      
#' # 00 = red, 10 = purple, 01 = magenta, 11 = blue, 02 = orange, 12 = green
#' 
#' # we can plot the phylogeny to take a look
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   phy <- make.phylo(sim$SIM)
#'   
#'   # color 0 valued tips red and 1 valued tips blue
#'   ape::plot.phylo(phy, tip.color = index)
#' }
#' # one could further complicate the model by adding hidden states
#' # to each trait, each with its own number etc, but these examples
#' # include all the tools necessary to make these or further extensions
#' 
#' @name bd.sim.traits
#' @rdname bd.sim.traits
#' @export

bd.sim.traits <- function(n0, lambda, mu,
                         tMax = Inf, N = Inf,
                         nTraits = 1, nFocus = 1, nStates = 2, nHidden = 1,
                         X0 = 0, Q = list(matrix(c(0, 0.1, 0.1, 0), 
                                                 ncol = 2, nrow = 2)),
                         nFinal = c(0, Inf), nExtant = c(0, Inf)) {
  # check that n0 is not negative
  if (n0 <= 0) {
    stop("initial number of species must be positive")
  }
  
  # check nFinal's length
  if (length(nFinal) != 2) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
  }
  
  # if nFocus is just one number, make it a vector
  if (length(nFocus) == 1) {
    nFocus <- rep(nFocus, 2)
  }
  
  # if nStates and nHidden is just one number, make them vectors
  if (length(nStates) == 1) {
    nStates <- rep(nStates, nTraits)
  }
  if (length(nHidden) == 1) {
    nHidden <- rep(nHidden, nTraits)
  }
  
  # finally, same for Q
  if (length(Q) == 1) {
    Q <- rep(list(Q[[1]]), nTraits)
  }

  # check that rates are numeric
  if (!is.numeric(c(lambda, mu))) {
    stop("lambda and mu must be numeric")
  } else if (any(!(c(length(lambda), length(mu)) %in% 
                   c(1, nStates * nHidden)))) {
    stop("lambda and mu must be of length either one or equal to nStates")
  } else {
    # if everything is good, we set flags on whether each are TD
    tdLambda <- length(lambda) == nStates[nFocus[1]] * nHidden[nFocus[1]]
    tdMu <- length(mu) == nStates[nFocus[2]] * nHidden[nFocus[2]]
  }
  
  # check that rates are non-negative
  if (any(lambda < 0) || any(mu < 0)) {
    stop("rates cannot be negative")
  }
  
  # check nFinal is sensible - two numbers, maximum >=1
  if ((length(nFinal) != 2) || (typeof(nFinal) != "double")) {
    stop("nFinal must be a vector with a minimum and maximum number 
         of species")
  } else if (max(nFinal) < 1) {
    stop("nFinal must have a maximum number of species greater than 0")
  } else {
    # if everything is good, make sure it's sorted
    nFinal <- sort(nFinal)
  }
  
  # similarly for nExtant
  if ((length(nExtant) != 2) || (typeof(nExtant) != "double")) {
    stop("nExtant must be a vector with a minimum and maximum number 
         of species")
  } else if (max(nExtant) < 0) {
    stop("nExtant must have a maximum number of species greater 
         than or equal to 0")
  } else {
    # if everything is good, make sure it's sorted
    nExtant <- sort(nExtant)
  }

  # error checks for each trait
  for (i in 1:nTraits) {
    # take the Q for this trait
    Qn <- Q[[i]]

    # check that Qn is square
    if (ncol(Qn) != nrow(Qn)) {
      stop("Q must be a square matrix")
    }

    # check that Qn has row number (and therefore column number)
    # equal to the number of traits
    if (nrow(Qn) != nStates[i] * nHidden[i]) {
      stop("Q must have transition rates for all trait combinations")
    }
  }
  
  # check which of N or tMax is not Inf, and condition as needed
  if ((N == Inf) && (tMax == Inf)) {
    stop("Either tMax or N must not be Inf.")
  } else if ((N != Inf) && (tMax != Inf)) {
    stop("Only one condition can be set, 
         so only one of N or tMax can be non-Inf")
  } else if (N != Inf) {
    condition = "number"
  } else {
    condition = "time"
  }
  
  # do we condition on the number of species?
  condN <- condition == "number"
  
  # if conditioning on number, need to set some 
  # tMax for the trait evolution
  if (condN) {
    # thrice the average time it will take to get to 10*N species
    # under the slowest diversification rate parameters
    trTMax <- max(log(10*N) / (lambda - mu))
  } else trTMax <- Inf
  
  # whether our condition is met - rejection sampling for
  # time, or exactly number of species at the end for number
  condMet <- FALSE
  
  # counter to make sure the nFinal is achievable
  counter <- 1
  
  while (!condMet) {
    # create vectors to hold times of speciation, extinction, 
    # parents and status
    TS <- rep(0, n0)
    TE <- rep(NA, n0)
    parent <- rep(NA, n0)
    isExtant <- rep(TRUE, n0)
    now <- c()

    # initialize species count
    sCount <- 1
    
    # initialize list of traits for first species' traits
    traits <- list(traits.musse(tMax = min(tMax, trTMax), tStart = 0, 
                                nTraits = nTraits, nStates = nStates,
                                nHidden = nHidden,
                                X0 = X0, Q = Q))
    
    # while we have species to be analyzed still
    while (length(TE) >= sCount) {
      # if sCount > nFinal[2], no reason to continue
      if (sCount > nFinal[2]) {
        # so we fail the condMet test
        sCount <- Inf
        
        # leave while
        break
      }

      # get argument for the next species
      # it will be the oldest one that hasn't lived yet
      sp <- which(TS == sort(TS)[sCount])
      
      # start at the time of speciation of sp
      tNow <- TS[sp]

      # get traits data set for each rate
      traitsSpLambda <- traits[[sp]][[nFocus[1]]]
      traitsSpMu <- traits[[sp]][[nFocus[2]]]
      
      # find the waiting time using rexp.var if lambda is not constant
      waitTimeS <- ifelse(tdLambda, 
                          ifelse(sum(lambda) > 0,
                                 rexp.musse(1, lambda, 
                                             traitsSpLambda,
                                             tNow, tMax), Inf),
                          ifelse(lambda > 0, 
                                 rexp(1, lambda), Inf))
      
      waitTimeE <- ifelse(tdMu,
                          ifelse(sum(mu) > 0,
                                 rexp.musse(1, mu, 
                                             traitsSpMu,
                                             tNow, tMax), Inf),
                          ifelse(mu > 0, 
                                 rexp(1, mu), Inf))
      
      tExp <- tNow + waitTimeE
      
      # while there are fast enough speciations before the species 
      # goes extinct,
      while ((tNow + waitTimeS) < min(tExp, tMax)) {
        # advance to the time of speciation
        tNow <- tNow + waitTimeS
        now <- c(now, tNow)
        
        # add new times to the vectors
        TS <- c(TS, tNow)
        TE <- c(TE, NA)
        parent <- c(parent, sp)
        isExtant <- c(isExtant, TRUE)
        
        # get trait value at tNow
        XPar <- unlist(lapply(traits[[sp]], function(x) 
          tail(x$value[x$min < tNow], 1)))

        # run trait evolution and append it to list
        traits[[length(traits) + 1]] <- 
          traits.musse(tMax = min(tMax, max(trTMax, tNow + 100)), 
                       tStart = tNow, 
                       nTraits = nTraits, nStates = nStates, nHidden = nHidden,
                       X0 = XPar, Q = Q)
        
        # get a new speciation waiting time
        waitTimeS <- ifelse(tdLambda, 
                            ifelse(sum(lambda) > 0,
                                   rexp.musse(1, lambda, 
                                               traitsSpLambda, 
                                               tNow, tMax), Inf),
                            ifelse(lambda > 0, 
                                   rexp(1, lambda), Inf))
      }
      
      # reached the time of extinction
      tNow <- tExp

      # if the species went extinct before tMax,
      # record it, otherwise, record NA
      TE[sp] <- ifelse(tNow < tMax, tNow, NA)
      
      # record the extinction
      isExtant[sp] <- is.na(TE[sp]) | TE[sp] > tMax
      
      # next species
      sCount <- sCount + 1
      
      # if we passed the limit
      if (condN && (sum(TS < tNow & (is.na(TE) | TE > tNow)) > 10*N)) {
        # function to find the excess at t
        nAlive <- Vectorize(function(t) {
          sum(TS <= t & (is.na(TE) | TE > t)) - N
        })
        
        # vector of all events
        events <- sort(c(TS, TE))
        
        # find times when we were at N
        nAliveT <- events[which(nAlive(events) == 0)]
        
        # find the times where the next event happened for each
        nextEvent <- unlist(lapply(nAliveT, function(t) events[events > t][1]))
        
        # get chosen event index
        eventChosen <- sample(1:length(nextEvent), 1, 
                              prob = (nextEvent - nAliveT))
        
        # draw uniform as above
        tMax <- runif(1, nAliveT[eventChosen], nextEvent[eventChosen])

        # adjust isExtant to TRUE for those alive at tMax
        isExtant[TE >= tMax] <- TRUE
        
        # adjust TE to NA for those alive at tMax
        TE[isExtant] <- NA
        
        # set condMet to true
        condMet <- TRUE
        
        # leave while
        break
      }
    }
    
    # truncate traits so we only go up to tMax
    for (i in 1:length(traits)) {
      for (j in 1:nTraits) {
        # df in question
        traitsSp <- traits[[i]][[j]]

        # eliminate rows with min greater than tMax
        traitsSp <- traitsSp[traitsSp$min < tMax, ]
        
        # set max of last row to tMax
        traitsSp$max[nrow(traitsSp)] <- tMax
        
        # invert time for max and min
        traitsSp$max <- tMax - traitsSp$max
        traitsSp$min <- tMax - traitsSp$min
        
        # invert columns
        colnames(traitsSp) <- c("value", "max", "min")
        
        # set traits back to it
        traits[[i]][[j]] <- traitsSp
      }
    }
    
    # now we invert TE and TS so time goes from tMax to 0
    TE <- tMax - TE
    TS <- tMax - TS
    
    # check which species are born after tMax
    nPrune <- which(TS <= 0)
    
    # if any, prune them
    if (length(nPrune) > 0) {
      TE <- TE[-nPrune]
      TS <- TS[-nPrune]
      isExtant <- isExtant[-nPrune]
      traits <- traits[-nPrune]
      
      # need to be careful with parent
      parent <- c(NA, unlist(lapply(parent[-nPrune][-1], function(x)
        x - sum(nPrune < x))))
      # each species pruned lowers the parent numbering before them by 1
    }
    
    # check whether we are in bounds if rejection sampling is the thing
    if (!condN) {
      condMet <- (length(TE) >= nFinal[1]) &&
      (length(TE) <= nFinal[2]) &&
      (sum(isExtant) >= nExtant[1]) &&
      (sum(isExtant) <= nExtant[2])
    }
    
    # if we have ran for too long, stop
    counter <- counter + 1
    if (counter > 100000) {
      warning("This value of nFinal took more than 100000 simulations 
              to achieve")
      return(NA)
    }
  }
  
  # truncate traits so that last min time is extinction time
  for (sp in 1:length(TE)) {
    # iterate through number of traits
    for (tr in 1:length(traits[[sp]])) {
      # get traits data frame
      traitsSp <- traits[[sp]][[tr]]

      # if there are hidden states
      if (nHidden[tr] > 1) {
        # set them to normal states
        traitsSp$value <- traitsSp$value %% nStates[tr]
        
        if (nrow(traitsSp) > 1) {
          # duplicate rows
          dup <- c()
          
          # counting how many duplicates in a row
          count <- 0
          
          #iterate through rows to make sure there are no duplicates
          for (i in 2:nrow(traitsSp)) {
            if (traitsSp$value[i] == traitsSp$value[i - 1]) {
              # add to dup
              dup <- c(dup, i)
              
              # increase count of dups
              count <- count + 1
            } else {
              # if count > 0, change the max of count rows ago to max of last row
              if (count > 0) {
                traitsSp$min[i - 1 - count] <- traitsSp$min[i - 1]
                
                # return count to 0
                count <- 0
              }
            }
          }
          
          # need to do a last check in case the last row is a duplicate
          if (count > 0) {
            traitsSp$min[i - count] <- traitsSp$min[i]
          }
          
          # delete duplicates
          if (!is.null(dup)) traitsSp <- traitsSp[-dup, ]
        }
      }

      # check if it is extinct
      if (!is.na(TE[sp])) {
        # if so, find last row with maximum higher than TE
        lastRow <- max(which(traitsSp$max > TE[sp]))
        
        # set min of last row to TE
        traitsSp$min[lastRow] <- TE[sp]
        
        # delete other rows
        if (lastRow < nrow(traitsSp)) {
          # row(s) to delete
          if (lastRow + 1 < nrow(traitsSp)) {
            delRow <- (lastRow + 1):nrow(traitsSp)
          } else {
            delRow <- lastRow + 1
          }
          traitsSp <- traitsSp[-delRow, ]
        }
      }
      
      # reassign traits data frame
      traits[[sp]][[tr]] <- traitsSp
    }
  }
  
  # create the return
  sim <- list(TE = TE, TS = TS, PAR = parent, EXTANT = isExtant)
  class(sim) <- "sim"
  
  res <- list("TRAITS" = traits, "SIM" = sim)

  return(res)
}

###
# function to run trait evolution for a given species
traits.musse <- function(tMax, tStart = 0, nTraits = 1, nStates = 2,
                         nHidden = 1, X0 = 0,
                         Q = list(matrix(c(0, 0.1, 0.1, 0), 2, 2))) {
  # create a return list
  traits <- vector(mode = "list", length = nTraits)
  
  # for each trait
  for (i in 1:nTraits) {
    # get number of states and initial states
    nStatesI <- nStates[i]
    nHiddenI <- nHidden[i]
    X0I <- ifelse(length(X0) > 1, X0[i], X0)
    
    # get Q matrix
    if (length(Q) > 1) QI <- Q[[i]] else QI <- Q[[1]]
    
    # append traits data frame to traits
    traits[[i]] <- trait.musse(tMax, tStart, nStatesI, nHiddenI, X0I, QI)
  }
  
  return(traits)
}

trait.musse <- function(tMax, tStart = 0, nStates = 2, nHidden = 1, X0 = 0,
                         Q = matrix(c(0, 0.1, 0.1, 0), 2, 2)) {
  # make sure tMax and tStart are numbers
  if (!is.numeric(c(tMax, tStart))) {
    stop("tMax and tStart must be numeric")
  } else if (length(tMax) > 1 || length(tStart) > 1) {
    stop("tMax and tStart must be one number")
  }

  # make sure tMax > tStart
  if (tStart >= tMax) {
    stop("tMax must be greater than tStart")
  }
  
  # make sure X0 is a possible state
  if (X0 > nStates * nHidden - 1) {
    stop("X0 must be an achievable state (i.e. in 0:(nStates * nHidden - 1)")
  }
  
  # create states vector from number
  states <- 0:(nStates * nHidden - 1)

  # create traits data frame
  traits <- data.frame(value = X0, min = tStart, max = NA)
  
  # start a time counter
  tNow <- tStart
  
  # and a shifts counter
  shifts <- 0
  
  # make diagonals of Q 0 
  diag(Q) <- 0
  
  # while we have not reached the end
  while (tNow < tMax) {
    # current state
    curState <- traits$value[shifts + 1]
    
    # get the total rate of transition from the current state
    rTotal <- sum(Q[curState + 1, ])
    
    # get the time until the next transition
    waitTime <- ifelse(rTotal > 0, rexp(1, rTotal), Inf)
    
    # increase time
    tNow <- tNow + waitTime
    
    # add max to traits
    traits$max[shifts + 1] <- min(tNow, tMax)
    
    # break if needed
    if (tNow >= tMax) break
    
    # increase shifts counter
    shifts <- shifts + 1
    
    # sample to find target state
    newState <- sample(states, 1, prob = Q[curState + 1, ])
    
    # add it to traits data frame
    traits[shifts + 1, ] <- c(newState, tNow, NA)
  }
  
  return(traits)
}
