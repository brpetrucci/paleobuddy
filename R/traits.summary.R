#' Summarizing trait data
#'
#' Summarizes trait data from a \code{sim} object, usually the ouput of
#' \code{bd.sim} in the case where diversification rates are trait-dependent.
#' Returns a list of trait values at the present or the time of extinction
#' (depending on whether the species is alive at present), and optionally
#' returns values at the time of fossil sampling if provided with a fossil
#' record object \code{fossils}, usually the output of \code{sample.clade}. 
#' Does not make assumptions on the number of traits described in the
#' \code{traits} parameter, so that if that list has more than one trait per
#' species, multiple vectors will be returned by the function.
#' 
#' @inheritParams make.phylo
#' 
#' @param traits List of trait data frames, usually one of the returns of 
#' \code{bd.sim}. \code{traits[[i]][[j]]} should correspond to the \code{j}th
#' trait data frame for species \code{i}. The data frames contain the following
#' columns
#' 
#' \describe{
#' \item{\code{value}}{A vector of trait values the species took at specific
#' intervals of time.}
#' 
#' \item{\code{max}}{A vector of time values corresponding to the upper bound
#' of each interval.}
#' 
#' \item{\code{min}}{A vector of time values corresponding to the lower bound
#' of each interval}}
#' 
#' @param selection Which subset of species to collect trait data for. If set
#' to \code{"all"}, it will return every trait value it has access to, i.e.
#' either all species, living or dead, or all species plus fossils if
#' \code{fossils} is supplied. If set to \code{"extant"}, it will return only
#' trait values for living species. If set to \code{"extinct"}, it will return
#' only trait values for extinct species, and fossils if \code{fossils} is 
#' supplied. If set to \code{"fossil"}, it will return values for only the
#' fossil species (and therefore requires a \code{fossils} parameter). If set 
#' to \code{"sampled"}, it will function the same as in the case for 
#' \code{"extant"}, except it will also return values for the fossils if 
#' \code{fossils} is supplied.
#' 
#' @return A named list of named vectors of trait values. List element names
#' refer to each trait, so i.e. \code{res$traitN} will correspond to the vector
#' of trait values for trait \code{N}. Vector element names refer to the
#' species, using the default naming convention of the package (\code{tN} is
#' the \code{N}th species in the simulation, and \code{tN.M} is the \code{M}th
#' sampled fossil of that species). 
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' @examples
#' 
#' ###
#' # need a simple simulation to use as an example
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
#' # get all trait values
#' traitSummary <- traits.summary(sim$SIM, sim$TRAITS)
#' traitSummary
#' 
#' # could get only the extant values, instead
#' traitSummary <- traits.summary(sim$SIM, sim$TRAITS, selection = "extant")
#' traitSummary
#' 
#' # or all the extinct values
#' traitSummary <- traits.summary(sim$SIM, sim$TRAITS, selection = "extinct")
#' traitSummary
#' 
#' # set seed
#' set.seed(1)
#' 
#' # maybe we want to take a look at the traits of fossil records too
#' fossils <- sample.clade(sim$SIM, rho = 0.5, tMax = max(sim$SIM$TS))
#' 
#' # get the trait values for all extinct species, including fossil samples
#' traitSummary <- traits.summary(sim$SIM, sim$TRAITS, 
#'                                fossils = fossils, selection = "extinct")
#' traitSummary
#' 
#' # can also get the values for all sampled species, i.e. extant or fossils
#' traitSummary <- traits.summary(sim$SIM, sim$TRAITS, 
#'                                fossils = fossils, selection = "sampled")
#' traitSummary
#' 
#' # or just the fossil species
#' traitSummary <- traits.summary(sim$SIM, sim$TRAITS, 
#'                                fossils = fossils, selection = "fossil")
#' traitSummary
#' 
#' @name traits.summary
#' @rdname traits.summary
#' @export
#' 

traits.summary <- function(sim, traits, fossils = NULL, selection = "all") {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # check that traits has the same length as the number of species in sim
  if (length(traits) != length(sim$TE)) {
    stop("traits must have trait values for all species")
  }
  
  # check that selection is within our bounds
  if (!(selection %in% c("all", "extant", "extinct", "fossil", "sampled"))) {
    stop("selection parameter must be 'all', 'extant', 
         'extinct', 'fossil', or 'sampled'")
  }
  
  # create final list return
  res <- vector(mode = "list", length = length(traits[[1]]))

  # iterate through traits
  for (t in 1:length(traits[[1]])) {
    # create return
    traitList <- c()
    traitNames <- c()
    traitStatus <- c()
    
    # iterate through all species
    for (sp in 1:length(sim$TE)) {
      # get traits data frame for sp
      traitsSp <- traits[[sp]][[t]]
      
      # get trait value at present or TE
      traitList <- c(traitList, tail(traitsSp$value, 1))
      
      # add name and status
      traitNames <- c(traitNames, paste0("t", sp))
      traitStatus <- c(traitStatus, c("extinct", "extant")[sim$EXTANT[sp] + 1])

      # check whether we need to go through fossils or not
      if (!is.null(fossils) && sum(fossils$Species == paste0("t", sp)) > 0) {
        # select fossil rows that are for this species
        fossilsSp <- fossils[fossils$Species == paste0("t", sp), ]
        
        # iterate through each fossil
        for (i in 1:nrow(fossilsSp)) {
          # find time to sample trait value
          traitTime <- fossilsSp$SampT[i]
          
          # find trait value at that time
          traitList <- c(traitList, 
                         traitsSp$value[which(traitsSp$max >= traitTime & 
                                                traitsSp$min <= traitTime)])
          
          # add name and status
          traitNames <- c(traitNames, 
                          paste0("t", sp +
                                   i/10^ceiling(log(nrow(fossilsSp) + 1, 10)),
                                 ifelse(i %% 10 == 0, "0", "")))
          traitStatus <- c(traitStatus, "fossil")
        }
      }
    }

    # name traitList
    names(traitList) <- traitNames
    
    # get final result based on selection
    if (selection == "all") {
      res[[t]] <- traitList
    } else if (selection == "extant") {
      res[[t]] <- traitList[traitStatus == "extant"]
    } else if (selection == "extinct") {
      res[[t]] <- traitList[traitStatus %in% c("extinct", "fossil")]
    } else if (selection == "fossil") {
      res[[t]] <- traitList[traitStatus == "fossil"]
    } else if (selection == "sampled") {
      res[[t]] <- traitList[traitStatus %in% c("extant", "fossil")]
    }
  }
  
  # name list
  names(res) <- paste0("trait", 1:length(traits[[1]]))
  
  return(res)
}
