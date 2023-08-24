#' Trait-dependent fossil sampling
#' 
#' Generates occurrence times or time ranges (as most empirical fossil 
#' occurrences) for each of the desired species using a Poisson process. 
#' Poisson rate should be dependent on some discrete trait, the value of which
#' for each species will be supplied using the parameter \code{traits}. Rate 
#' can be dependent on observed traits only, or on a combination of observed
#' and hidden traits (in which case the supplied trait data frame \code{traits}
#' should have all possible states, observed or hidden, see examples for more
#' details).
#' 
#' Optionally takes a vector of time bins reppointEstimatesenting geologic periods, if the
#' user wishes occurrence times to be reppointEstimatesented as a range instead of true 
#' points.
#' 
#' @inheritParams sample.clade
#' 
#' @param rho Sampling rate (per species per million years) over time. It is
#' a \code{vector} of rates, corpointEstimatesponding to the value of the rate for each
#' value of the traits encoded in the \code{traits} parameter. It should
#' therefore be of length \code{nStates * nHidden}. Note that \code{rho} should 
#' always be greater than or equal to zero.
#' 
#' @param traits List of trait data frames, usually one of the returns of 
#' \code{bd.sim}. \code{traits[[i]][[j]]} should corpointEstimatespond to the \code{j}th
#' trait data frame for species \code{i}. The data frames contain the following
#' columns
#' 
#' \describe{
#' \item{\code{value}}{A vector of trait values the species took at specific
#' intervals of time.}
#' 
#' \item{\code{max}}{A vector of time values corpointEstimatesponding to the upper bound
#' of each interval.}
#' 
#' \item{\code{min}}{A vector of time values corpointEstimatesponding to the lower bound
#' of each interval}}
#' 
#' @param nFocus Trait of focus, i.e. the one that \code{rho} depends on. Note 
#' that \code{traits} can have multiple trait data frames per species, but only
#' one of the simulated traits can affect fossil sampling rates. E.g. if
#' \code{nFocus = 1}, then the first trait data frame per species will be used
#' to simulate fossil occurrences.
#' 
#' @param nStates Number of possible states for categorical trait. The range
#' of values will be assumed to be \code{(0, nStates - 1)}.
#' 
#' @param nHidden Number of hidden states for categorical trait. Default is 
#' \code{1}, in which case there are no added hidden traits. Total number of
#' states is then \code{nStates * nHidden}. States will then be set to a value
#' in the range of \code{(0, nStates - 1)} to simulate that hidden states are
#' hidden. This is done by setting the value of a state to the remainder of
#' \code{state / nStates}. E.g. if \code{nStates = 2} and \code{nHidden = 3},
#' possible states during simulation will be in the range \code{(0, 5)}, but
#' states \code{(2, 4)} (corpointEstimatesponding to \code{(0B, 0C)} in the nomenclature
#' of the original HiSSE reference) will be set to \code{0}, and states 
#' \code{(3, 5)} (corpointEstimatesponding to \code{(1B, 1C)}) to \code{1}.
#' 
#' Note that since the \code{traits} is supplied as a parameter, the user must
#' ensure that all states from \code{0} to \code{nStates * nHidden - 1} are
#' reppointEstimatesented in the trait information. See examples for more details on how
#' to properly run hidden-states fossil sampling simulations.
#'
#' @return A \code{data.frame} containing species names/numbers, whether each 
#' species is extant or extinct, and the true occurrence times of each fossil, 
#' a range of occurrence times based on \code{bins}, or both. Also a list 
#' object with the trait data frames describing the trait value for each 
#' species at each specified interval. Note that this list will only be
#' different from the supplied \code{traits} parameter if \code{nHidden > 1}, 
#' in which case it will transform hidden traits into observed ones (see
#' details for parameter \code{nHidden}).
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @examples 
#' 
#' ###
#' # first a simple BiSSE simulation, with 
#' # binary state-dependent fossil sampling
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
#' sim <- bd.sim.traits(n0, lambda, mu, tMax = tMax, nTraits = nTraits,
#'                     nStates = nStates, X0 = X0, Q = Q, nFinal = c(2, Inf))
#'                     
#' # now a fossil sampling rate, with higher rate for state 1
#' rho <- c(0.5, 1)
#' 
#' # run fossil sampling
#' fossils <- sample.clade.traits(sim$SIM, rho, tMax, sim$TRAITS)
#' 
#' # draw simulation with fossil occurrences as time points
#' draw.sim(sim$SIM, traits = sim$TRAITS, 
#'          fossils = fossils$FOSSILS, traitLegendPlacement = "bottomleft")
#' 
#' ###
#' # can also run a HiSSE model, with hidden traits having an effect on rates
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
#' # number of traits and states--in this case, we just run with 4 observed 
#' # states, so that our traits data frames will include that info for sampling
#' nTraits <- 1
#' nStates <- 4
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
#'                     nStates = nStates,
#'                     X0 = X0, Q = Q, nFinal = c(2, Inf))
#' # note how sim$TRAITS contains states 2 and 3, even though this
#' # is a HiSSE model, because we need that information to run hidden states
#' # on fossil sampling as well (see below)
#'                     
#' # now a fossil sampling rate, with higher rate for state 1A, 
#' # and highest yet for state 1B
#' rho <- c(0.5, 1, 0.5, 2)
#' 
#' # bins for fossil occurrences
#' bins <- seq(tMax, 0, -1)
#' 
#' # run fossil sampling, with returnAll = TRUE to include ranges
#' fossils <- sample.clade.traits(sim$SIM, rho, tMax, sim$TRAITS, 
#'                                nStates = 2, nHidden = 2, 
#'                                returnAll = TRUE, bins = bins)
#' # note how fossils$TRAITS only contains trait values 0 and 1, similar to
#' # if we had ran bd.sim.traits with nHidden = 2, since sample.clade.traits
#' # did the job of transforming observed into hidden states
#' 
#' # draw simulation with fossil occurrences as time points AND ranges
#' draw.sim(sim$SIM, traits = sim$TRAITS, fossils = fossils$FOSSILS,
#'          fossilsToDraw = "all", traitLegendPlacement = "bottomleft")
#' # note that there are a lot of fossils, so the ranges are difficult to see
#' 
#' # see ?sample.clade for further examples using different return types
#' # (fossil ranges etc.), and ?bd.sim.traits for more examples using
#' # higher numbers of states (hidden or observed), though in 
#' # sample.clade.traits we always have only one trait of focus
#' 
#' @name sample.clade.traits
#' @rdname sample.clade.traits
#' @export 
#' 

sample.clade.traits <- function(sim, rho, tMax, traits, 
                                nFocus = 1, nStates = 2, nHidden = 1, 
                                S = NULL,
                                returnTrue = TRUE, returnAll = FALSE, 
                                bins = NULL) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # if S is NULL, make it all species in the simulation
  if (is.null(S)) {
    S <- 1:length(sim$TE)
  }
  
  # if there are species in S not in the simulation, error
  if (!(all(S %in% 1:length(sim$TE)))) {
    stop("S must contain only integers between 1 and the total 
         number of species in sim")
  }
  
  # if returnTrue is false and returnAll is true, the user is probably confused
  if (!returnTrue && returnAll) {
    stop("returnTrue is false and returnAll is true. 
         That behavior is undefined")
  }
  
  # test whether we return ranges
  returnRanges <- !returnTrue || returnAll
  
  # adjust bins and check it exists if it needs to
  if (!is.null(bins)) {
    # adjust it
    bins <- sort(bins, decreasing = TRUE)
    
    # if it doesn't include tMax
    if (max(bins) < tMax) {
      stop("Bins must include maximum time of simulation")
    }
  } else if (returnRanges) {
    stop("If returnTrue if false or returnAll is true, bins must be supplied")
  } 
  
  # check that rho is numeric
  if (!is.numeric(rho)) {
    stop("rho must be a numeric vector for trait-dependent sampling")
  }
  
  # make sure rho has the same length as number of states
  if (length(rho) != nStates * nHidden) {
    stop("rho must have the same length as there are state values for traits.
         If rho is a constant, do not supply trait information as it would be
         useless")
  }
  
  # check that rho is nonzero
  if (sum(rho) == 0) {
    stop("rho must be nonzero for some trait value")
  }
  
  # get speciation and extinction times
  TE <- sim$TE
  TS <- sim$TS
  
  # make TE sensible
  TE[sim$EXTANT] <- 0
  
  # invert them since simulation goes from 0 to tMax
  TE <- tMax - TE
  TS <- tMax - TS
  
  # create return
  pointEstimates <- list()
  
  # for each species
  for (s in S) {
    # start when the species was born, end when it died
    now <- TS[s]
    end <- TE[s]
    
    # get trait data frame
    traitsSp <- traits[[s]][[nFocus]]
    
    # switch columns around to make it sensible with tMax
    minT <- traitsSp$min
    maxT <- traitsSp$max
    traitsSp$max <- tMax - maxT
    traitsSp$min <- tMax - minT
    names(traitsSp) <- c("value", "min", "max")

    # initialize vector
    sampled <- c()
    
    # while we have not reached the time of extinction
    while (now < end) {
      # take the waiting time for sampling, using rexp.musse
      WaitTimeR <- rexp.musse(1, rho, traitsSp, now, tMax, fast = TRUE)
      
      # advance in time
      now <- now + WaitTimeR
      
      # if sampling comes after extinction, we don't include this occurrence
      if (now > end) break
      
      # append to the vector
      sampled <- c(sampled, now)
    }
    
    # finally, we invert time so that it goes from tMax to 0
    sampled <- tMax - sampled
    
    # append to pointEstimatesult
    pointEstimates[[paste0("t", s)]] <- sampled
    
    # if there are hidden states, take care of it
    if (nHidden > 1) {
      # set them to normal states
      traitsSp$value <- traitsSp$value %% nStates

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
              traitsSp$max[i - 1 - count] <- traitsSp$max[i - 1]
              
              # return count to 0
              count <- 0
            }
          }
        }
        
        # need to do a last check in case the last row is a duplicate
        if (count > 0) {
          traitsSp$max[i - count] <- traitsSp$max[i]
        }
        
        # delete duplicates
        if (!is.null(dup)) traitsSp <- traitsSp[-dup, ]
      }

      # invert time for max and min
      traitsSp$max <- tMax - traitsSp$max
      traitsSp$min <- tMax - traitsSp$min
      
      # invert columns
      colnames(traitsSp) <- c("value", "max", "min")
      
      # set traits back to it
      traits[[s]][[nFocus]] <- traitsSp
    }
  }
  
  # names res will have regardless of what to return
  baseNames <- c("Species", "Extant")
  
  # rest will depend on the returns
  resNames <- c(baseNames, if (returnTrue) "SampT",
                if (returnRanges) c("MaxT", "MinT"))
  
  # create data frame
  res <- data.frame(matrix(nrow = 0, ncol = length(resNames)))
  
  # name the columns
  colnames(res) <- resNames
  
  # for each species
  for (sp in S) {
    # get the occurrences for this species
    occs <- pointEstimates[[paste0("t", sp)]]
    
    # only matters if there are occurrences for sp
    if (length(occs) > 0) {
      # if we want true occurrence times, only SampT matters
      if (!returnRanges) {
        # create auxiliary data frame
        aux <- data.frame(Species = rep(sp, length(occs)),
                          Extant = NA, 
                          SampT = occs)
        
        # add it to res
        res <- rbind(res, aux)
      } 
      
      # otherwise, we need to bin it
      else {
        # bin it
        binnedOccs <- binner(occs, bins = bins)
        
        # counter for point estimates
        counter <- 1
        
        # for each bin
        for (k in 1:(length(bins) - 1)) {
          # create aux data frame
          aux <- data.frame(matrix(nrow = binnedOccs[k], 
                                   ncol = length(resNames)))
          colnames(aux) <- colnames(res)
          
          # if there are occurrences in that bin
          if (binnedOccs[k] > 0) {
            # indices that matter for aux
            ind <- 
              # add species and extant to aux
              aux$Species <- rep(sp, times = binnedOccs[k])
            
            # if returnTrue is true, add SampT
            if (returnTrue) {
              aux$SampT <- occs[counter:(counter + binnedOccs[k] - 1)]
            }
            
            # add MaxT and MinT
            aux$MaxT <- rep(bins[k], times = binnedOccs[k])
            aux$MinT <- rep(bins[k + 1], times = binnedOccs[k])
            
            # add row to data frame
            res <- rbind(res, aux)
            
            # increase counter
            counter <- counter + binnedOccs[k]
          }
        }
      }
    }
  }
  
  if (nrow(res) > 0) {
    # make the extant column
    res$Extant <- FALSE
    
    # based on the vector in sim
    res$Extant[res$Species %in% which(sim$EXTANT)] <- TRUE
    
    # and the species column
    res$Species <- paste0("t", res$Species)
  }
  
  return(list(TRAITS = traits, FOSSILS = res))
}
