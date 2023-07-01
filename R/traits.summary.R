#' Summarizing trait data
#'
#' ...
#'
#' @param sim ...
#' 
#' @param traits ...
#' 
#' @param nFocus ...
#' 
#' @param fossils ...
#' 
#' @param selection ...
#' 
#' @return 
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' @references
#' 
#' @examples
#' 
#' @name traits.summary
#' @rdname traits.summary
#' @export
#' 

traits.summary <- function(sim, traits, nFocus = 1,
                           fossils = NULL, selection = "all") {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid argument, must be a sim object. See ?sim")
  }
  
  # check that traits has the same length as the number of species in sim
  if (length(traits) != length(sim$TE)) {
    stop("traits must have trait values for all species")
  }
  
  # check that selection is within our bounds
  if (!(selection %in% c("all", "extant", "extinct", "sampled"))) {
    stop("selection parameter must be 'all', 'extant', 
         'extinct', or 'sampled'")
  }
  
  # create return
  traitList <- c()
  traitNames <- c()
  traitStatus <- c()
  
  # iterate through all species
  for (sp in 1:length(sim$TE)) {
    # get traits data frame for sp
    traitsSp <- traits[[sp]][[nFocus]]

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
    res <- traitList
  } else if (selection == "extant") {
    res <- traitList[traitStatus == "extant"]
  } else if (selection == "extinct") {
    res <- traitList[traitStatus %in% c("extinct", "fossil")]
  } else if (selection == "sampled") {
    res <- traitList[traitStatus %in% c("extant", "fossil")]
  }
  
  return(res)
}
