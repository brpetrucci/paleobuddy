#' Converting a phylogeny in a paleobuddy object
#'
#' Generates a \code{sim} object using a \code{phylo} object and some additional 
#' information (depending on other inputs). It is the inverse of the 
#' \code{make.phylo} function. Input is (1) a phylogeny, following a "Hennigian" 
#' (sensu Ezard et al 2011) format, (2) information on the "mother lineage" of 
#' each tip in the phylogeny (see "mothers" section in "details", below), (3) the 
#' status ("extant" or "extinct") of each lineage, (4) the stem age (or age of
#' origination of the clade), and (5) the "stem length" (or time interval 
#' between the stem age and the appearance of the first "daughter"). The user can
#' also choose if the event dating should be done from root to tips or from tips 
#' to root (this choice is important - see "dating" section in "details" below. 
#' The function returns a \code{sim} object (which contains speciation and
#' extinction times, parent, and status information). The function does not accept
#' more than one species having \code{NA} as parent (which is interpreted as if 
#' there were no single common ancestor in the phylogeny). 
#'
#' @param phy A \code{phylo} object, which may contain only extant or extant and
#' extinct lineages.
#'
#' @param mothers Vector containing the mother of each tip in the phylogeny. 
#' First species' mother should be \code{NA}.
#' 
#' @param extant Logical vetor indicating which lineages are extant and extinct.
#' 
#' @param dateFromPresent Logical vector indicating if TS/TE events should be 
#' dated from present to root (\code{TRUE}, default value) of from root to 
#' present. Please see "dating" section in "details", below. It is internally set 
#' to \code{FALSE} and prints a message in the prompt if there are no extant
#' species in the \code{extant} vector.
#' 
#' @param stemAge Numeric vetor indicating the age, in absolute geological time
#' (Mya), when the first lineage of the clade originated. It is not needed when
#' \code{dateFromPresent} is \code{TRUE} and \code{stemLength} is provided, or 
#' when \code{phy} has a \code{root.edge}. This argument is required if 
#' \code{dateFromPresent} is \code{FALSE}.
#' 
#' @param stemLength Numeric vector indicating the time difference between the 
#' \code{stemAge} and the appearance of its first "daughter" lineage (that is, the
#' second lineage to originate in the phylogeny). This argument is required if 
#' \code{dateFromPresent} is \code{FALSE}, but users have no need to assign values
#' in this parameter if \code{phy} have a \code{$root.edge}, which is taken by the
#' function as the \code{stemLength} value.
#' 
#' @return A \code{sim} object organized in the same format as the output of 
#' \code{bd.sim}. Items in the object follow their tip assignment in the 
#' phylogeny.
#' 
#' @details 
#' 
#' Mothers:
#' 
#' The function needs the indication of a mother lineage for every tip in the 
#' phylogeny but one (which is interpreted as the first known lineage in the 
#' clade, which should have \code{NA} as the mother). This assignment might be 
#' straightforward for simulations (as in the examples section below), but is 
#' evidently a non-trivial task for empirical phylogenies. As there are many 
#' ways to assign impossible combinations of motherthood, the function does not
#' return any specific error message if the provided motherhood does not map to
#' possible lineages given the phylogeny. Instead, simulations conducted by the 
#' author showed the function tends to crash when an "impossible" motherhood is 
#' assigned, but is not guaranteed that this will happen because of the enormous
#' universe of "impossible" ways to assign motherhood. However, if the function 
#' crashes when all lineages have reasonable motherhood, users are invited to 
#' contact the author.
#' 
#' Dating:
#' 
#' Phylogenies store relative distances between speciation (and possibly 
#' extinction) times of each lineage. However, to get absolute times for those 
#' events (which are required to construct the output of this function), users 
#' should provide a moment in absolute geological time to position the phylogeny. 
#' This could be (1) the present, in the case at least one lineage in the 
#' phylogeny is extant, or (2) some time in the past, which is the \code{stemAge} 
#' parameter. Those two possible dating methods are used by setting 
#' \code{dateFromPresent} to \code{TRUE} or \code{FALSE}, respectively (see 
#' \code{dateFromPresent} above). If users do not have a reasonable value for 
#' \code{stemAge}, they are encouraged to use present to root dating 
#' (\code{dateFromPresent = TRUE}), as deviations in the value of 
#' \code{stemLength} will only affect the speciation time of the first lineage of
#' the clade. When \code{dateFromPresent} is set to \code{FALSE}, eventual errors 
#' in \code{stemAge} or \code{stemLength} will bias the dating of all nodes in the
#' phylogeny.
#' 
#' @author Matheus Januario. 
#' 
#' @references
#' 
#' Ezard, T. H., Pearson, P. N., Aze, T., & Purvis, A. (2012). The meaning of 
#' birth and death (in macroevolutionary birth-death models). Biology letters, 
#' 8(1), 139-142.
#'
#' @examples
#'
#' # to check the usage of the function, let us make sure it transforms a 
#' # phylogeny generated with make.phylo back into the original simulation
#' 
#' ### 
#' # birth-death process
#' 
#' # simulate the clade
#' tmax <- 10
#' sim <- bd.sim(1, lambda = 0.3, mu = 0.1, tMax = tmax, nFinal = c(10, Inf))
#' 
#' # convert birth-death into phylo
#' phy <- make.phylo(sim)
#' 
#' # convert phylo into birth-death again
#' res <- phylo.to.sim(phy = phy, extant = sim$EXTANT, mothers = sim$PAR)
#' 
#' # test if simulation and converted object are the same
#' all.equal(sim, res)
#' 
#' 
#' ###
#' # pure birth process
#' 
#' # simulate the clade
#' tmax <- 10
#' sim <- bd.sim(1, lambda = 0.2, mu = 0, tMax = tmax, nFinal = c(10, Inf))
#' 
#' # convert birth-death into phylo
#' phy <- make.phylo(sim)
#' 
#' # convert phylo into birth-death again
#' # note we can supply optional arguments, see description above
#' res <- phylo.to.sim(phy = phy, extant = sim$EXTANT, mothers = sim$PAR, 
#'                 stemAge = tmax, stemLength = (tmax-sim$TS[2]), 
#'                 dateFromPresent = TRUE)
#' 
#' # testing if simulation and converted object are the same
#' all.equal(sim, res)
#' 
#' @name phylo.to.sim
#' @rdname phylo.to.sim
#' @export
#' 

phylo.to.sim <- function(phy, mothers, extant, dateFromPresent = TRUE,
                         stemAge = NULL, stemLength = NULL){
  
  # checking inputs
  if (is.null(stemLength)) {
    if("root.edge" %in% names(phy)){
      stemLength <- phy$root.edge
    }
  }
  
  # if more than one species has NA as a mother, we have a problem
  if (sum(is.na(mothers) > 1)) {
    stop("\n Function assumes all lineages (except one) are monophyletic, so 
         only one lineage is allowed to have \"no mother\".")
  }
  
  # if the user wants dateFromPresent
  if (dateFromPresent) {
    # it can only work if there are extant lineages
    if (sum(extant) < 1) {
      message("\n No extant lineages, \"dateFromPresent\" will be set 
              to FALSE")
      
      # changing date mode
      dateFromPresent <- FALSE 
    }
    
    else {
      # we need either stemAge or stemLength to date from the present
      if (sum(c(is.null(stemAge), is.null(stemLength))) > 1) {
        stop("\n Please provide \"stemAge\" OR \"stemLength\" if
             \"dateFromPresent\" is TRUE")
      }  
    }
  }
  
  # if the user wants to date from the root, need stemAge and stemLength
  if (!(dateFromPresent)) { 
    if (is.null(stemAge) & is.null(stemAge)) {
      stop("\n Please provide \"stemAge\" AND \"stemLength\" if
           \"dateFromPresent\" is FALSE")
    }
  }
  
  # declaring functions
  
  # dating functions
  #dates from furthest to closest to the present
  date.nodes.forward <- function(phy, stemAge, stemLength) {
    dating <- vector()
    
    #dates first node (first birth of sim):
    dating[which(phy$edge[, 1] == min(phy$edge[, 1]))] <- stemAge - stemLength
    
    # choose a focal node (the closest to StemAge without dating)
    fnode <- min(phy$edge[, 1]) + 1 
    
    # while dating is not finished
    while (fnode <= max(phy$edge)) { 
      # find where the nodes connects
      ids <- which(phy$edge[, 1] == fnode) 
      
      # dates it
      dating[ids] <- dating[which(phy$edge[, 2] == fnode)] - 
        phy$edge.length[which(phy$edge[, 2] == fnode)]
      
      # changes fnode to posterior node
      fnode <- fnode + 1 
    }
    
    # wrapping it all together
    res <- unique(cbind(phy$edge[, 1], dating))
    colnames(res) <- c("node", "dating")
    res <- as.data.frame(res)
    return(res)
  }
  
  
  # dates from to present to the past
  date.nodes.rewind <- function(phy, extant) {
    # error if function not apply
    if (sum(extant) == 0) { 
      stop("\n No extant lineages. \"dateFromPresent\" is impossible")
    }
    
    dating <- vector()
    
    # dating extant (0 - [each extant edge length])
    dating[phy$edge[, 2] %in% which(extant)] <- 
      phy$edge.length[phy$edge[, 2] %in% which(extant)] 
    
    # gets which is missing  dating
    tab <- table(phy$edge[!(is.na(dating)), 1])
    
    ids <- as.numeric(names(tab[tab == 1]))
    
    # dating the extinct sisters of an extant
    for (i in 1:length(ids)) {
      dating[which(phy$edge[, 1] == ids[i] & is.na(dating))] <- 
        dating[which(phy$edge[, 1] == ids[i] & !(is.na(dating)))]
    }
    
    aux <- min(phy$edge[, 1])
    
    # dating nodes in the insides of the phylo, from superficial to deep nodes
    while ((sum(is.na(dating[phy$edge[, 1] == min(phy$edge[, 1])]))) > 0) { 
      # gets a nono-dated node
      aux <- min(phy$edge[!(is.na(dating)), 1])
      
      # dates node
      new_date <- unique(dating[phy$edge[,1] == aux]) + 
        phy$edge.length[which(phy$edge[, 2] == aux)]
      aux <- phy$edge[which(phy$edge[, 2] == aux), 1]
      
      # assign dating
      dating[which(phy$edge[, 1] == aux)] <- new_date 
    }
    
    # dating from root to tips
    aux <- min(phy$edge[, 1]) 
    
    # dating in direction of root
    while (sum(is.na(dating)) > 0){ 
      aux <- min(phy$edge[((is.na(dating))), 1])
      new_date <- dating[phy$edge[, 2] == aux] - 
        phy$edge.length[phy$edge[, 2] == aux]
      dating[which(phy$edge[, 1] == aux)] <- new_date
    }
    
    # wrapping it all together
    res <- unique(cbind(phy$edge[, 1], dating))
    colnames(res) <- c("node", "dating")
    res <- as.data.frame(res)
    return(res)
  }
  
  # coalescence function
  
  # "coalesces" lineage until first node in phy
  # returns list of nodes until that event
  coal.lin <- function(lin, phy) { 
    # if lienage has no mother
    if (is.na(lin)) { 
      return(NA) 
    }
    
    lin.coal <- lin
    leng.aft <- 1
    stop <- FALSE
    
    # coalescing branches until ancestor
    while (!(stop)) {
      
      # append lins between mother and phy
      lin.coal <- c(lin.coal, 
                    phy$edge[(phy$edge[, 2] == lin.coal[length(lin.coal)]), 1])
      len.bef <- length(lin.coal)
      
      # if not coalescing anymore
      if (len.bef == leng.aft) { 
        # stop while
        stop <- TRUE
      }
      
      # updating leng for testing next time
      else { 
        leng.aft <- len.bef 
      }
      
    }
    return(lin.coal)
  }
  
  # choose dating method given inputs
  if (dateFromPresent) {
    dated_nodes <- date.nodes.rewind(phy, extant)
  }
  
  else {
    dated_nodes <- date.nodes.forward(phy, stemAge, stemLength)
  }
  
  
  # dating births and deaths
  res <- list(TE = vector(), TS = vector(), PAR = mothers, EXTANT = extant)
  
  for (i in 1:length(phy$tip.label)) {
    
    mot.coal <- coal.lin(mothers[i], phy)
    dau.coal <- coal.lin(i, phy)
    
    # coalescence between mother and child
    ids <- which(!(dau.coal %in% mot.coal)) 
    
    # correcting
    ids <- c(ids, max(ids) + 1) 
    dau.coal <- dau.coal[ids]
    
    # first lineage
    if (is.na(mothers[i])) {
      if (is.null(stemAge)) {
        stemAge <- stemLength+max(dated_nodes$dating)
      }
      res$TS[i] <- stemAge 
    }
    
    # other lineages
    else {  
      res$TS[i] <- dated_nodes$dating[dated_nodes$node == 
                                        dau.coal[length(dau.coal)]]
    }
    
    # if it is extant, it has no TE
    if (res$EXTANT[i]) {
      res$TE[i] <- NA
    }
    
    else {
      # if extinct, calculates TE summing edgelength w/ dated node
      res$TE[i] <-  dated_nodes$dating[dated_nodes$node == dau.coal[2]] -
        phy$edge.length[which(apply(phy$edge, 1, function(x) all.equal(
          x, sort(dau.coal[1:2], decreasing = TRUE))) == "TRUE")]
    }
  }
  
  # one simple transformation if all lineages are extant:
  if (class(res$TE) == "logical") {
    res$TE <- as.numeric(res$TE)
  }
  
  return(res)
}