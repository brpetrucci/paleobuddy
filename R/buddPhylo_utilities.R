#' @rdname buddPhylo
#'
#' @inheritParams buddPhylo
#'
#' @details \code{getParents.buddPhylo} Find all the parental lineages of a
#' \code{focalLineage} from a \code{buddPhylo$name} object, from first to
#' last ancestor. Note this does not use the information in the
#' \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario
#'
#' @export

getParents.buddPhylo <- function(buddPhylo, focalLineage) {
  res <- vector()
  res <- c(res, buddPhylo$parent[buddPhylo$name == focalLineage])
  
  while (!is.na(res[length(res)])) {
    focalLineage <- res[length(res)]
    
    res <- c(res, buddPhylo$parent[buddPhylo$name == focalLineage])
  }
  
  # removing the NA at the end:
  res <- res[-length(res)]
  return(res)
}

#' @rdname buddPhylo
#'
#' @inheritParams buddPhylo
#' 
#' @details \code{getChildren.buddPhylo} returns the immediate descendant
#' (daughter) lineages of a given \code{focalLineage} in a \code{buddPhylo}
#' object (using the \code{name} column). The output is a character vector of
#' length two containing the \code{name} of each child, equivalent to
#' \code{getDescendants.buddPhylo(..., onlyImmediates = TRUE)}.
#'
#' By default (\code{returnInfo = TRUE}), the function also annotates each
#' descendant based on its status. If \code{focalLineage} leads to a sampled
#' ancestor, one child is labeled accordingly and the other as a "lineage". If
#' it instead \code{focalLineage} leads to a speciation event, descendants are
#' labeled according to their \code{orientation}. Note that this function does
#' not use information from the \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario
#'
#' @export

getChildren.buddPhylo <- function(buddPhylo, focalLineage, returnInfo = TRUE) {
  desc <- buddPhylo$name[buddPhylo$parent == focalLineage]
  desc <- desc[!is.na(desc)]
  res <- desc
  
  if (returnInfo) {
    if (0 %in% buddPhylo$length[buddPhylo$name %in% desc]) {
      aux <- buddPhylo$length[match(desc, buddPhylo$name)]
      aux <- ifelse(aux == 0, yes = "sampAnc", no = "lineage")
      names(res) <- aux
      res <- res[order(names(res), decreasing = TRUE)]
    } else { # if no children is sampled ancesotr:
      aux <- buddPhylo$orientation[match(desc, buddPhylo$name)]
      names(res) <- aux
      # make sure ancestor is always first
      res <- res[order(names(res))]
    }
  }
  
  return(res)
}

#' @rdname buddPhylo
#'
#' @inheritParams buddPhylo
#'
#' @details \code{getDescendants.buddPhylo} Find all the descendant lineages of a
#' \code{focalLineage} from a \code{buddPhylo$name} object, from most immediate
#' to most distal descendant. \code{getChildren.buddPhylo} returns only the
#' most immediate descendants (same effect as \code{getDescendants.buddPhylo(..., onlyImmediates = T)}). Note this does not use the information in the
#' \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario
#'
#' @export

getDescendants.buddPhylo <- function(buddPhylo, focalLineage, 
                                     onlyImmediates = FALSE, all = FALSE) {
  checked <- vector()
  toCheck <- focalLineage
  res <- vector()
  
  while (length(toCheck) > 0) {
    focalLineage <- toCheck[1]
    desc <- buddPhylo$name[buddPhylo$parent == focalLineage]
    desc <- desc[!is.na(desc)]
    
    if (length(desc) > 0) {
      toCheck <- c(toCheck, desc)
      res <- c(res, desc)
    }
    
    if (onlyImmediates) {
      break
    } else {
      toCheck <- toCheck[!toCheck %in% focalLineage]
    }
  }
  
  if (!all) {
    internal_nodes <- which(res %in% buddPhylo$name[buddPhylo$type == "node"])
    if (length(internal_nodes) > 0) res <- res[-internal_nodes]
  }
  
  return(res)
}


#' @rdname buddPhylo
#'
#' @inheritParams buddPhylo
#'
#' @details \code{getDescendants.buddPhylo} Find all the descendant lineages of a
#' \code{focalLineage} from a \code{buddPhylo$name} object, from most immediate
#' to most distal descendant. Note this does not use the information in the
#' \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario
#'
#' @export

getMRCA.buddPhylo <- function(buddPhylo, nameList) {
  # First,get the LUCA of the list:
  ParentList <- list()
  for (i in 1:length(nameList)) {
    ParentList[[i]] <- getParents.buddPhylo(buddPhylo, nameList[i])
  }
  # find longest par-desc length:
  lengs <- unlist(lapply(ParentList, length))
  ref <- which.max(lengs)
  
  listChecker <- function(pars, ca) {
    return(ca %in% pars)
  }
  
  UCA <- toCheck <- rev(unlist(ParentList[ref])) # check tipwards from LUCA
  test <- all(unlist(lapply(ParentList, listChecker, toCheck[1])))
  while (test) {
    toCheck <- toCheck[-1]
    test <- all(unlist(lapply(ParentList, listChecker, toCheck[1])))
  }
  
  UCA <- UCA[!UCA %in% toCheck]
  MRCA <- rev(UCA)[1]
  return(MRCA)
}


#' @rdname buddPhylo
#'
#' @inheritParams buddPhylo
#'
#' @details \code{getDescendants.buddPhylo} Find all the descendant lineages of a
#' \code{focalLineage} from a \code{buddPhylo$name} object, from most immediate
#' to most distal descendant. Note this does not use the information in the
#' \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario
#'
#' @export

is.monophyletic.buddPhylo <- function(buddPhylo, nameList,
                                      excludeSampAnc = TRUE) {
  MRCA <- getMRCA.buddPhylo(buddPhylo, nameList)
  
  alldesc <- getDescendants.buddPhylo(buddPhylo, MRCA)
  
  if (excludeSampAnc) {
    alldesc <- alldesc[alldesc %in% buddPhylo$name[buddPhylo$type == "sampAnc"]]
    res <- all(alldesc %in% nameList)
  } else {
    res <- all(alldesc %in% nameList)
  }
  
  return(res)
}

#' Genertaing a phlyo object from a buddPhylo object
#'
#' Generates a phylogeny ina \code{phylo} format based on the ifnromation within 
#' a \code{buddPhylo} object.
#' 
#' @param buddPhylo Object of class "buddPhylo"
#' 
#' @return A \code{phylo} object compatible with fucntions from package \code{ape}
#' 
#' @author Matheus Januario
#' 

as.phylo.buddPhylo <- function(buddPhylo) {
  ########################################################
  # aux fucnctions:
  # To count resembalnces in node paths:
  countEqual <- function(x, ref) {
    return(sum(ref %in% x))
  }
  
  # function to re-order edges so the phlyo object does not break:
  reorder_edges <- function(edge) { #input is any edge matrix
    parents <- edge[, 1]
    children <- edge[, 2]
    
    root <- setdiff(parents, children)
    
    adj <- split(edge[, 2], edge[, 1])
    
    new_edge <- matrix(ncol = 2, nrow = 0)
    
    dfs <- function(node) {
      if (!is.null(adj[[as.character(node)]])) {
        for (child in adj[[as.character(node)]]) {
          # record edge in traversal order
          new_edge <<- rbind(new_edge, c(node, child))
          
          # recurse
          dfs(child)
        }
      }
    }
    
    dfs(root)
    
    return(new_edge)
  }
  ########################################################
  
  ########### Part 1: Labels
  # labels
  tip.label <- buddPhylo$name[buddPhylo$type %in% c("tip", "sampAnc")]
  node.label <- buddPhylo$name[!buddPhylo$type %in% c("tip", "sampAnc")]
  
  ########### Part 2: Edge matrix
  # Make a node dictionary:
  node.id <- rep(NA, times = length(node.label))
  names(node.id) <- node.label
  # First, make a list with all descendats:
  pars <- list()
  for (tt in 1:length(tip.label)) {
    pars[[tt]] <- getParents.buddPhylo(buddPhylo,
                                       focalLineage = tip.label[tt]
    )
  }
  names(pars) <- tip.label
  
  # we can start from anywhere, so choose the tip w/ most descs:
  plen <- unlist(lapply(pars, length))
  
  refDesc <- tip.label[which.max(plen)]
  refID <- which(tip.label == refDesc)
  nextNode <- length(tip.label) + 1
  while (length(pars) > 0) {
    aux <- pars[[refID]]
    tempIDS <- !aux %in% names(node.id[!is.na(node.id)])
    sumTempIDS <- sum(tempIDS)
    if (sumTempIDS > 0) {
      #print(paste0("N = ", length(pars)))
      aux2 <- aux[tempIDS]
      node.id[match(rev(aux2), names(node.id))] <- nextNode:(nextNode + sumTempIDS - 1)
      nextNode <- max(node.id, na.rm = T) + 1
    }
    
    # choose & update next ref
    pars[[refID]] <- NULL
    refID <- which.max(unlist(lapply(pars, countEqual, ref = aux)))
    refDesc <- names(pars[refID])
  }
  
  # aux objects"
  all.labels <- c(tip.label, names(node.id))
  all.ids <- c(1:length(tip.label), node.id)
  
  # edge matrix
  edge1 <- node.id[match(buddPhylo$parent, node.label)]
  edge2 <- all.ids[match(buddPhylo$name, all.labels)]
  edge <- cbind(edge1, edge2)
  edge <- edge[-which(is.na(edge[, 1])), ] # NA element is the root, ignore it
  
  # Finally, we re-order so phylo() doesnt break:
  edge <- reorder_edges(edge)
  
  ########### Part 3: Edge lengths
  edge_names <- all.labels[match(edge[, 2], all.ids)] # aux obj
  edge.length <- buddPhylo$length[match(edge_names, buddPhylo$name)]
  
  ########### Part 4: Assemble the tree:
  phy <- list(
    edge = edge, edge.length = edge.length, Nnode = length(node.label),
    node.label = node.label, tip.label = tip.label,
    root.edge = buddPhylo$length[buddPhylo$orientation == "uca"]
  )
  class(phy) <- "phylo"
  return(phy)
}

#' Adjust the timescale of a buddPhylo pbject
#'
#' Given a vector of fossil occurrences and time bins to represent geological
#' ranges, returns the occurrence counts in each bin.
#'
#' @param buddPhylo Object of class "buddPhylo"
#'
#' @param timeFromRoot A logical indicating whether time should be expressed as
#' time since the root (\code{timeFromRoot=TRUE}). If set to \code{FALSE}
#' (default), time is expressed as millions of years ago (Mya), taking the tip
#' extant tip that is further from the root as a reference for "present".
#'
#' @return A \code{buddPhylo} with the time adjusted accordingly.
#'
#' @author Matheus Januario
#'

adjust.timescale.buddPhylo <- function(buddPhylo, timeFromRoot = FALSE) {
  # first check that it is a valid buddPhylo
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  if (timeFromRoot) {
    needsCorrection1 <- (diff(c(0, max(buddPhylo$x_coord[buddPhylo$extant]))) < 10E-10)
    if (needsCorrection1) {
      buddPhylo$x_coord <- max(buddPhylo$x_coord) - buddPhylo$x_coord
    }
  } else { # if time from present:
    needsCorrection2 <- !(diff(c(0, max(buddPhylo$x_coord[buddPhylo$extant]))) < 10E-10)
    if (needsCorrection2) {
      refMax <- max(c(buddPhylo$x_coord, buddPhylo$x_par))
      buddPhylo$x_coord <- refMax - buddPhylo$x_coord
      buddPhylo$x_par <- refMax - buddPhylo$x_par
    }
  }
  return(buddPhylo)
}

#' List the monophyletic clades in a budding phylogeny
#'
#' Get a list of monophyletic clades in this budding phylogeny.
#'
#' @param buddPhylo Object of class "buddPhylo"
#' 
#' @param returnInfo Whether to return the MRCA (if the defining node of the
#' clade is a sampled-ancestor node) or orientation (if the defining node of
#' the clade is a speciation node) for each clade. Defaults to \code{TRUE}.
#'
#' @return A list of length equal to the number of nodes, with each element
#' being a vector of strings (listing the tips/SAs that are descended
#' from each node). If \code{returnInfo = TRUE}, each element contains two 
#' vectors of strings: either the SA and lineage, or the ancestor clade and
#' the descendant clade (similar to \code{getChildren.buddPhylo}).
#'
#' @author Bruno do Rosario Petrucci
#'

subtrees.buddPhylo <- function(buddPhylo, returnInfo = TRUE) {
  # extract information from buddPhylo
  nms   <- buddPhylo$name
  lins <- buddPhylo$lineage
  pars  <- buddPhylo$parent
  types <- buddPhylo$type
  oris  <- buddPhylo$orientation
  
  # if names and lineage differ, get lineage instead
  if (any(nms != lins, na.rm = TRUE)) 
    nms[!is.na(lins)] <- lins[!is.na(lins)]

  # build children lookup
  children <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_along(nms)) {
    # get parent of nms
    p <- pars[i]
    
    # check that there is a parent
    if (!is.na(p)) {
      # check if the parent is already in the environment
      if (exists(p, envir = children, inherits = FALSE))
        # if so, add this species to its descendants
        children[[p]] <- c(children[[p]], nms[i])
      else
        # otherwise, just start it as this species
        children[[p]] <- nms[i]
    }
  }
  
  # the descendants of all tips (i.e. the tips themselves)
  tip_set <- new.env(hash = TRUE, parent = emptyenv())
  
  # iterate through all names
  for (i in seq_along(nms)) {
    # if this is not a node, i.e. if it is a tip or SA, add it to tip_set
    if (types[i] != "node") tip_set[[nms[i]]] <- nms[i]
  }
  
  # get root, i.e. UCA
  root  <- nms[is.na(pars) & types == "node"]
  
  # empty character vector, we will use to keep track of the order
  order <- character(sum(types == "node"))
  
  # start at root
  stack <- root
  
  # how many nodes we've already checked
  k <- 0
  
  # while we have something on the stack
  while (length(stack)) {
    # get the first stack element
    n <- stack[1]
    
    # remove first element from stack
    stack <- stack[-1]
    
    # increase k
    k <- k + 1
    
    # make kth order equal to the current node
    order[k] <- n
    
    # if n is a node, get its children, otherwise nothing
    ch    <- if (exists(n, envir = children, inherits = FALSE)) children[[n]]
    else character(0L)
    
    # add children of this node to the stack, if any
    stack <- c(ch[types[match(ch, nms)] == "node"], stack)
  }
  # now we know the order to fill tip_sets
  
  # fill tip sets in reverse (post-order = children before parents)
  for (n in rev(order[seq_len(k)])) {
    # get the childrenm of this node, if they exist
    ch <- if (exists(n, envir = children, inherits = FALSE)) children[[n]]
    else character(0L)
    
    # add to tip_set
    tip_set[[n]] <- unlist(lapply(ch, function(x) tip_set[[x]]))
  }
  # now we have a hash table for the descendants of each node
  
  # get all the nodes
  nodes <- nms[types == "node"]
  
  # iterate through nodes
  lapply(nodes, function(n) {
    # get the children of this node
    ch <- children[[n]]
    
    # get the types of children
    type <- types[match(ch, nms)]
    
    # if one of the children is an SA, this is an SA node
    if (any(type == "sampAnc")) {
      # get which one is the SA
      sa  <- ch[type == "sampAnc"]
      
      # and which one is the lineage
      lin <- ch[type != "sampAnc"]
      
      # return the named list with the SA and the descendants of lin
      if (returnInfo) list(sampAnc = sa,
                           lineage = sort(unlist(lapply(lin, 
                                                        function(x) tip_set[[x]]))))
      else sort(c(sa, unlist(lapply(lin, 
                                         function(x) tip_set[[x]]))))
    } else {
      # otherwise, get the ancestor and descendant lineage
      ori <- oris[match(ch, nms)]
      anc <- ch[ori == "ancestor"]
      des <- ch[ori == "descendant"]
      
      # and return the named list with each and their descendants
      if (returnInfo) list(ancestor = sort(unlist(lapply(anc, function(x) tip_set[[x]]))),
           descendant = sort(unlist(lapply(des, function(x) tip_set[[x]]))))
      else sort(c(unlist(lapply(anc, function(x) tip_set[[x]])),
                       unlist(lapply(des, function(x) tip_set[[x]]))))
    }
  })
}


# subtrees.buddPhylo <- function(buddPhylo, returnInfo = TRUE) {
#   # internal nodes
#   nodes <- buddPhylo$name[buddPhylo$type == "node"]
#   
#   # number of subtrees
#   n_subtrees <- length(nodes)
#   
#   # start the results list
#   res <- vector("list", n_subtrees)
#   
#   # iterate through number of subtrees
#   for (i in 1:n_subtrees) {
#     # check if we want the info
#     if (returnInfo) {
#       # get the children of this node
#       children <- getChildren.buddPhylo(buddPhylo, nodes[i])
#       
#       res[[i]] <- lapply(children, function(x) {
#         # get the descendants
#         desc <- sort(getDescendants.buddPhylo(buddPhylo, x))
#         
#         # if it has any descendants, return those, otherwise just return it
#         if (length(desc) > 0) desc else x
#       })
#       
#       # name it
#       names(res[[i]]) <- names(children)
#     } else {
#       # get the descendants of this node
#       desc <- sort(getDescendants.buddPhylo(buddPhylo, nodes[i]))
#       
#       # add to res
#       res[[i]] <- desc
#     }
#   }
#   
#   # return result
#   return(res)
# }

