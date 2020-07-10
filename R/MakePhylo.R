#' Converts a paleobuddy simulation into a phylogeny
#'
#' \code{MakePhylo} generates a \code{phylo} object using a simulation from the
#' \code{BDSim} function. The phylogeny follows a "Hennigian" (sensu Ezard et
#' al 2011) format. If the simulation has only one lineage, the function
#' returns \code{NA} as there is no phylogeny for a simulation with only one
#' lineage.
#'
#' @param sim a simulation from the \code{BDSim} function.
#'
#' @author Function written by Matheus Januario. 
#' 
#' Reference: Ezard, T. H.,
#' Pearson, P. N., Aze, T., & Purvis, A. (2012). The meaning of birth and death
#' (in macroevolutionary birth-death models). Biology letters, 8(1), 139-142.
#'
#' @return A \code{phylo} object.
#'
#' @examples
#'
#' # generating a phylogeny using constant rates
#' sim <- BDSim(n0 = 1, pp = 0.2, qq = 0.05, tMax = 10)
#' 
#' # in case first simulation has only one species
#' while (length(sim$TE) < 2) {
#'   sim <- BDSim(n0 = 1, pp = 0.2, qq = 0.05, tMax = 10)
#' }
#' 
#' phy <- MakePhylo(sim)
#' 
#' # we need ape to plot it
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   ape::plot.phylo(phy)
#'   
#'   # we can also plot only the molecular phylogeny
#'   ape::plot.phylo(ape::drop.fossil(phy))
#' }
#' 
#' # this works for sim generated with any of the scenarios in BDSim, of course
#' sim <- BDSim(n0 = 1, pp = function(t) 0.12 + 0.01*t, qq = 10, 
#'              tMax = 10, qShape = 1.3)
#' 
#' # in case first simulation has only one species
#' while (length(sim$TE) < 2) { 
#'   sim <- BDSim(n0 = 1, pp = function(t) 0.12 + 0.01*t, qq = 10, 
#'                tMax = 10, qShape = 1.3)
#' }
#' phy <- MakePhylo(sim)
#' 
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   ape::plot.phylo(phy)
#'   ape::plot.phylo(ape::drop.fossil(phy))
#' }
#' @name MakePhylo
#' @rdname MakePhylo
#' @export

MakePhylo <- function(sim) {

  # simulations with just one species do not have a phylogeny
  if (length(sim$TE) < 2) {
    message("There is no phylogeny for a simulation with only one lineage")
    return(NA)
  }

  if (sum(is.na(sim$PAR)) > 1) {
    stop("Multiple starting species. Use function findClades()")
  }

  all.dir.daughter <- function(lin, x) {
    # all.dir.daughters returns the name of each direct daughter species
    # x = a simulation from paleobuddy
    # lin = a numeric specyfing the name of a lineage
    return(which(x$PAR == lin))
  }

  # current node
  curNode <- length(sim$TE) + 1 
  
  # create the edge matrix
  edge <- matrix(nrow = 1, ncol = 2, data = c(curNode, NA)) 
  
  # lineages which the function already put in the phylogeny
  passed <- vector()
  
  # current lineage
  i <- 2 
  
  # lineages which the function still has to solve (at least)
  lins <- c(1, 2)
  
  # internal variable to help control the node function
  jump <- 0 
  
  # number of nodes in the phylogeny
  nNode <- length(sim$TE) - 1
  
  # vector storing the node corresponding to each birth
  birthsNode <- rep(NA, times = length(sim$TE)) 
  birthsNode[2] <- curNode
  
  # needed for debugging
  counter <- 0

  # while some tip does not have a place in the phylogeny
  while (length(lins) > 0) {
    # find daughters
    dau <- all.dir.daughter(lin = i, x = sim)
    dau <- dau[!(dau %in% passed)]

    # if lineage has daughters
    if (is.numeric(dau) & length(dau) > 0) {

      # if a whole clade has very recently been put in the phylogeny
      if (jump == 1) {
        curNode <- max(edge) + 1

        # append it to the edge matrix
        if (is.na(edge[nrow(edge), 2])) {
          # if there is no edge there currently
          edge[nrow(edge), 2] <- curNode
        } else {
          # if there is
          edge <- rbind(edge,
                        matrix(nrow = 1, ncol = 2, data = c(prevNode, curNode)))
        }
        
        # update birthsNode
        birthsNode[dau[1]] <- curNode
        
        # update jump
        jump<-0

      # if the current lineage is a non-monophyletic branch
      } else { 
        # update curNode
        curNode <- curNode + 1
        
        # append to edge matrix, as above
        if (is.na(edge[nrow(edge), 2])) {
          edge[nrow(edge), 2] <- curNode
        } else {
          edge <- rbind(edge,
                        matrix(nrow = 1, ncol = 2, 
                               data = c(curNode - 1, curNode)))
        }
        
        # update birthsNode
        birthsNode[dau[1]] <- curNode
      }

      # update edge
      edge <- rbind(edge,
                    matrix(nrow = 1, ncol = 2, data = c(curNode, NA)))
      
      # update lineage list and current lineage
      lins <- c(lins, dau[1])
      i <- lins[length(lins)]
    }

    # if lineage has no daughters
    if (is.numeric(dau) & length(dau) == 0) {

      # append lineage to the edge matrix
      if (is.na(edge[nrow(edge), 2])) {
        # if there is no edge there currently
        edge[nrow(edge), 2] <- i
      } else {
        # if there is
        edge <- rbind(edge, 
                    matrix(nrow = 1, ncol = 2, 
                           data = c(max(
                             edge[!(duplicated(edge[,1]) | 
                                duplicated(edge[,1], fromLast=TRUE)), 1]), i)))
      }
      
      # we put the lineage on the phylogeny
      passed <- c(passed, i)
      
      # update lineage list and current lineage
      lins <- lins[-length(lins)]
      i <- lins[length(lins)]
    }
    
    # this means that the function reached the end of the lineage of the curNode
    if (sum(edge[, 1] %in% curNode) > 1) {
      # the warning here only "affects" a a condition which is never satisfied
      # (jump when there is previous opened edge).
      suppressWarnings(
        {prevNode <- 
          max(edge[!(duplicated(edge[, 1]) | 
                       duplicated(edge[, 1], fromLast = TRUE)), 1])})
      
      # update jump
      jump <- 1
    }

    # registering bugs (if any)
    counter <- counter + 1
    
    # if the function ran for too long
    if (counter > 10*dim(edge)[1]) {
      return("The function is lost and seems that it will not find a phylogeny.
             Please report this error and provide this simulation for debugging")}

  }

  # calculating edge length
  edgeLength <- vector()
  for (i in 1:nrow(edge)) {
    # make auxiliary variables
    aux1 <- edge[i, 1]
    aux2 <- edge[i, 2]

    # if the branch is a tip
    if (aux2 <= length(sim$TE)) {
      # calculate length
      edgeLength[i] <- sim$TS[which(birthsNode == aux1)] - sim$TE[aux2]
    } else {
      # calculate length
      edgeLength[i] <- sim$TS[which(birthsNode == aux1)] -
                        sim$TS[which(birthsNode == aux2)]
    }

  }

  # Tyding all together to create the phylo object
  phy <- list(
    tip.label = paste0("t", 1:length(sim$TE)),
    edge = edge,
    edge.length = edgeLength,
    Nnode = nNode)
  class(phy) <- "phylo"

  return(phy)
}
