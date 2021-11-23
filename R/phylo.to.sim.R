#' Converting a phylogeny in a paleobuddy object
#'
#' Generates a \code{sim} object using a \code{phylo} object and some 
#' additional information (depending on other arguments). It is the inverse of 
#' the \code{make.phylo} function. Input is (1) a phylogeny, following an 
#' evolutionary Hennigian (sensu Ezard et al 2011) format (i.e., a fully 
#' bifurcating phylogeny), (2) information on the "mother lineage" of each tip 
#' in the phylogeny, (3) the status ("extant" or "extinct") of each lineage, 
#' (4) the stem age (or age of origination of the clade), and (5) the stem 
#' length (or time interval between the stem age and the first speciation 
#' event). The user can also choose if the event dating should be done from 
#' root to tips or from tips-to-root. The function returns a \code{sim} object 
#' (see \code{?sim}). The function does not accept more than one species having
#' \code{NA} as parent (which is interpreted as if there were no single common
#' ancestor in the phylogeny). In that case, use \code{find.lineages} first. 
#' 
#' See Details below for more information on each argument.
#'
#' @param phy A \code{phylo} object, which may contain only extant or extant 
#' and extinct lineages.
#'
#' @param mothers Vector containing the mother of each tip in the phylogeny. 
#' First species' mother should be \code{NA}. See details below.
#' 
#' @param extant Logical vetor indicating which lineages are extant and 
#' extinct.
#' 
#' @param dateFromPresent Logical vector indicating if speciation/extinction 
#' events should be dated from present-to-root (\code{TRUE}, default value) or
#' from root-to-present. As it is impossible to date "from present" without a
#' living lineage, it is internally set to \code{FALSE} and prints a message in
#' the prompt if there are no extant species.
#' 
#' @param stemAge Numeric vetor indicating the age, in absolute geological time
#' (million years ago), when the first lineage of the clade originated. It is 
#' not needed when \code{dateFromPresent} is \code{TRUE} and \code{stemLength} 
#' is provided, or when \code{phy} has a \code{root.edge}. This argument is 
#' required if \code{dateFromPresent} is \code{FALSE}.
#' 
#' @param stemLength Numeric vector indicating the time difference between the 
#' \code{stemAge} and the first speciation event of the group. This argument is
#' required if \code{dateFromPresent} is \code{FALSE}, but users have no need 
#' to assign values in this parameter if \code{phy} has a \code{$root.edge},
#' which is taken by the function as the \code{stemLength} value.
#' 
#' @return A \code{sim} object. For details, see \code{?sim}. Items in the 
#' object follow their tip assignment in the phylogeny.
#' 
#' @details 
#' 
#' Mothers:
#' 
#' The function needs the indication of a mother lineage for every tip in the 
#' phylogeny but one (which is interpreted as the first known lineage in the 
#' clade, and should have \code{NA} as the mother). This assignment might be 
#' straightforward for simulations (as in the examples section below), but is a
#' non-trivial task for empirical phylogenies. As there are many ways to assign
#' impossible combinations of motherthood, the function does not return any 
#' specific error message if the provided motherhood does not map to possible 
#' lineages given the phylogeny. Instead, the function tends to crash when an 
#' "impossible" motherhood is assigned, but this is not guaranteed to happen 
#' because the set of "impossible" ways to assign motherhood is vast, and 
#' therefore has not allowed for a test of every possibility. If the function 
#' crashes when all lineages have reasonable motherhood, users should submit an
#' issue report at 
#' \url{https://https://github.com/brpetrucci/paleobuddy/issues}.
#' 
#' Dating:
#' 
#' Phylogenies store the relative distances between speciation (and possibly 
#' extinction) times of each lineage. However, to get absolute times for those 
#' events (which are required to construct the output of this function), users 
#' should provide a moment in absolute geological time to position the 
#' phylogeny. This could be (1) the present, which is used as reference 
#' in the case at least one lineage in the phylogeny is extant (i.e., default 
#' behavior of the function), or (2) some time in the past, which is the 
#' \code{stemAge} parameter. Those two possible dating methods are used by 
#' setting \code{dateFromPresent} to \code{TRUE} or \code{FALSE}. If users have
#' extant lineages in their phylogeny but do not have a reasonable value for 
#' \code{stemAge}, they are encouraged to use present-to-root dating 
#' (\code{dateFromPresent = TRUE}), as in that case deviations in the value of 
#' \code{stemLength} will only affect the speciation time of the first lineage 
#' of the clade. In other words, when \code{dateFromPresent} is set to 
#' \code{FALSE}, user error in \code{stemAge} or \code{stemLength} will bias
#' the absolute (but not the relative) dating of all nodes in the phylogeny.
#' 
#' @author Matheus Januario.
#' 
#' @references
#' 
#' Ezard, T. H., Pearson, P. N., Aze, T., & Purvis, A. (2012). The meaning of 
#' birth and death (in macroevolutionary birth-death models). 
#' \emph{Biology letters}, 8(1), 139-142.
#'
#' @examples
#'
#' # to check the usage of the function, let us make sure it transforms a 
#' # phylogeny generated with make.phylo back into the original simulation
#' 
#' ### 
#' # birth-death process
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run simulation
#' sim <- bd.sim(1, lambda = 0.3, mu = 0.1, tMax = 10, nFinal = c(10, Inf))
#' 
#' # convert birth-death into phylo
#' phy <- make.phylo(sim)
#' 
#' # convert phylo into a sim object again
#' res <- phylo.to.sim(phy = phy, extant = sim$EXTANT, mothers = sim$PAR)
#' 
#' # test if simulation and converted object are the same
#' all.equal(sim, res)
#' 
#' ### 
#' # birth-death process with extinct lineages:
#' # set seed
#' set.seed(1)
#' 
#' # run simulation
#' sim <- bd.sim(1, lambda = 0.1, mu = 0.3, tMax = 10, nFinal = c(2, 4))
#' 
#' # convert birth-death into phylo
#' phy <- make.phylo(sim)
#' 
#' # convert phylo into a sim object again
#' res <- phylo.to.sim(phy = phy, extant = sim$EXTANT, mothers = sim$PAR, stemAge = max(sim$TS))
#' 
#' # test if simulation and converted object are the same
#' all.equal(sim, res)
#' 
#' ###
#' # pure birth process
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run simulation
#' sim <- bd.sim(1, lambda = 0.2, mu = 0, tMax = 10, nFinal = c(10, Inf))
#' 
#' # convert birth-death into phylo
#' phy <- make.phylo(sim)
#' 
#' # convert phylo into birth-death again
#' # note we can supply optional arguments, see description above
#' res <- phylo.to.sim(phy = phy, extant = sim$EXTANT, mothers = sim$PAR, 
#'                 stemAge = 10, stemLength = (10 - sim$TS[2]))
#' 
#' # testing if simulation and converted object are the same
#' all.equal(sim, res)
#' 
#' @name phylo.to.sim
#' @rdname phylo.to.sim
#' @export
#' 

phylo.to.sim <- function(phy, mothers, extant, dateFromPresent = TRUE,
                         stemAge = NULL, stemLength = NULL) {
  
  # if stemLength is not supplied, get it from root.edge if possible
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
  
  # auxiliary functions
  
  # dates from furthest to closest to the present
  date.nodes.forward <- function(phy, stemAge, stemLength) {
    # create return
    dating <- vector()
    
    # date first node (first birth of sim)
    dating[which(phy$edge[, 1] == min(phy$edge[, 1]))] <- stemAge - stemLength
    
    # choose a focal node (the closest to stemAge without dating)
    fnode <- min(phy$edge[, 1]) + 1 
    
    # while dating is not finished
    while (fnode <= max(phy$edge)) { 
      # find where the nodes connects
      ids <- which(phy$edge[, 1] == fnode) 
      
      # date it
      dating[ids] <- dating[which(phy$edge[, 2] == fnode)] - 
        phy$edge.length[which(phy$edge[, 2] == fnode)]
      
      # change fnode to posterior node
      fnode <- fnode + 1 
    }
    
    # wrap it all together
    res <- unique(cbind(phy$edge[, 1], dating))
    colnames(res) <- c("node", "dating")
    res <- as.data.frame(res)
    return(res)
  }
  
  
  # date from to present to the past
  date.nodes.rewind <- function(phy, extant) {
    # error if function not apply
    if (sum(extant) == 0) { 
      stop("\n No extant lineages. \"dateFromPresent\" is impossible")
    }
    
    # create return
    dating <- vector()
    
    # date extant (0 - [each extant edge length])
    dating[phy$edge[, 2] %in% which(extant)] <- 
      phy$edge.length[phy$edge[, 2] %in% which(extant)] 
    
    # get which is missing  dating
    tab <- table(phy$edge[!(is.na(dating)), 1])
    
    ids <- as.numeric(names(tab[tab == 1]))
    
    # date the extinct sisters of an extant
    for (i in 1:length(ids)) {
      dating[which(phy$edge[, 1] == ids[i] & is.na(dating))] <- 
        dating[which(phy$edge[, 1] == ids[i] & !(is.na(dating)))]
    }
    
    # auxiliary vector for dating
    aux <- min(phy$edge[, 1])
    
    # date nodes in the insides of the phylo, from superficial to deep nodes
    while ((sum(is.na(dating[phy$edge[, 1] == min(phy$edge[, 1])]))) > 0) { 
      # get a nono-dated node
      aux <- min(phy$edge[!(is.na(dating)), 1])
      
      # date node
      new_date <- unique(dating[phy$edge[,1] == aux]) + 
        phy$edge.length[which(phy$edge[, 2] == aux)]
      aux <- phy$edge[which(phy$edge[, 2] == aux), 1]
      
      # assign dating
      dating[which(phy$edge[, 1] == aux)] <- new_date 
    }
    
    # date from root to tips
    aux <- min(phy$edge[, 1]) 
    
    # date in direction of root
    while (sum(is.na(dating)) > 0){ 
      aux <- min(phy$edge[((is.na(dating))), 1])
      new_date <- dating[phy$edge[, 2] == aux] - 
        phy$edge.length[phy$edge[, 2] == aux]
      dating[which(phy$edge[, 1] == aux)] <- new_date
    }
    
    # wrap it all together
    res <- unique(cbind(phy$edge[, 1], dating))
    colnames(res) <- c("node", "dating")
    res <- as.data.frame(res)
    return(res)
  }
  
  # "coalesce" lineage until first node in phy
  # returns list of nodes until that event
  coal.lin <- function(lin, phy) { 
    # if lienage has no mother
    if (is.na(lin)) { 
      return(NA) 
    }
    
    # initialize aux variables
    lin.coal <- lin
    leng.aft <- 1
    stop <- FALSE
    
    # coalesce branches until ancestor
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
      
      # update leng for testing next time
      else { 
        leng.aft <- len.bef 
      }
      
    }
    return(lin.coal)
  }
  
  # choose dating method given inputs
  if (dateFromPresent) {
    dated_nodes <- date.nodes.rewind(phy, extant)
  } else {
    dated_nodes <- date.nodes.forward(phy, stemAge, stemLength)
  }
  
  
  # date births and deaths
  res <- list(TE = vector(), TS = vector(), PAR = mothers, EXTANT = extant)
  
  for (i in 1:length(phy$tip.label)) {
    # aux variables
    mot.coal <- coal.lin(mothers[i], phy)
    dau.coal <- coal.lin(i, phy)
    
    # coalescence between mother and child
    ids <- which(!(dau.coal %in% mot.coal)) 
    
    # correct it
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
  
  # set res as a sim object
  class(res) <- "sim"
  
  # make sure it is valid
  if (!is.sim(res)) {
    stop("Result is not a valid sim object")
  }
  
  return(res)
}
