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
  nms <- buddPhylo$name
  lins <- buddPhylo$lineage
  pars <- buddPhylo$parent
  types <- buddPhylo$type
  oris <- buddPhylo$orientation
  ranges <- buddPhylo$range
  
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
    ch <- if (exists(n, envir = children, inherits = FALSE)) children[[n]]
    else character(0L)
    
    # add children of this node to the stack, if any
    stack <- c(ch[types[match(ch, nms)] == "node"], stack)
  }
  # now we know the order to fill tip_sets
  
  # fill tip sets in reverse (post-order = children before parents)
  for (n in rev(order[seq_len(k)])) {
    # get the childrenm of this node, if they exist
    ch <- if (exists(n, envir = children, inherits = FALSE)) children[[n]]
    else character(0)
    
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
    
    # get the range
    ran <- ranges[match(n, nms)]
    
    # if one of the children is an SA, this is an SA node
    if (any(type == "sampAnc")) {
      # get which one is the SA
      sa  <- ch[type == "sampAnc"]
      
      # and which one is the lineage
      lin <- ch[type != "sampAnc"]

      
      # return the named list with the SA and the descendants of lin
      if (returnInfo) list(sampAnc = sa,
                           lineage = sort(unlist(lapply(lin, 
                                                        function(x) tip_set[[x]]))),
                           range = ran)
      else sort(c(sa, unlist(lapply(lin, 
                                    function(x) tip_set[[x]]))))
    } else {
      # otherwise, get the ancestor and descendant lineage
      ori <- oris[match(ch, nms)]
      anc <- ch[ori == "ancestor"]
      des <- ch[ori == "descendant"]
      
      # and return the named list with each and their descendants
      if (returnInfo) list(ancestor = sort(unlist(lapply(anc, 
                                                         function(x) tip_set[[x]]))),
                           descendant = sort(unlist(lapply(des, 
                                                           function(x) tip_set[[x]]))),
                           range = ran)
      else sort(c(unlist(lapply(anc, function(x) tip_set[[x]])),
                  unlist(lapply(des, function(x) tip_set[[x]]))))
    }
  })
}

#' Get clades from buddPhylos and phylos
#'
#' Get a list of monophyletic clades from a phylogenetic tree. Could be budding
#' or bifurcating. 
#'
#' @param tree Object of class "buddPhylo" or "phylo".
#' 
#' @param considerSAs Logical indicating whether sampled ancestors should be
#' considered in the definition of a clade.
#' 
#' @param considerOrientation Logical indicating whether speciation node 
#' orientation should be considered in the definition of a clade.
#'
#' @return A list of length equal to the number of internal nodes, with each
#' element being a string with taxa in a given clade separated by commas.
#' If \code{considerSAs = TRUE}, clades defined by an SA rather than speciation
#' node will end with \code{"\001SA:X"}, where X is the taxon in question.
#' If \code{considerOrientation = TRUE}, clades defined by a speciation node
#' will end with \code{"\001ANC:X;Y;..."}, where X, Y, etc. are the taxa in the
#' ancestral branch.
#'
#' @author Bruno do Rosario Petrucci

extract.clades <- function(tree,
                           considerSAs = FALSE,
                           considerOrientation = FALSE) {
  # check if this is a phylo
  is_phylo <- inherits(tree, "phylo")
  
  # check type
  if (is_phylo) {
    # tip labels
    tip_labels <- tree$tip.label
    
    # number of tips
    N <- length(tip_labels)
    
    # number of nodes
    Mn <- tree$Nnode
    
    # total number of nodes
    total <- N + Mn
    
    # post-order traversal edge order
    edge_po <- ape::reorder.phylo(tree, "postorder")$edge
    
    # start SA and anc nodes (the latter will never be filled)
    node_sa <- NULL

    # check if we want to fill SA nodes
    if (considerSAs && !is.null(tree$edge.length)) {
      # tip edge length, indexed by tip id
      tip_el <- tree$edge.length[match(seq_len(N), tree$edge[, 2])]
      
      # sa tips
      sa_tips <- which(tip_el < 1e-8)
      
      # check if there are any
      if (length(sa_tips) > 0) {
        # internal-node parent of each SA tip
        sa_pars <- tree$edge[match(sa_tips, tree$edge[, 2]), 1]
        
        # initialize full-length vector so it can be indexed by internal
        node_sa <- rep(NA_character_, total)
        
        # if a node has multiple SA children, take the first
        keep_first <- !duplicated(sa_pars)
        node_sa[sa_pars[keep_first]] <- tip_labels[sa_tips[keep_first]]
      }
    }
  } else if (inherits(tree, "buddPhylo")) {
    # get each column
    nm <- tree$name
    par <- tree$parent
    typ <- tree$type
    ori <- tree$orientation
    
    # get the number of rows
    N_row <- length(nm)
    
    # find what's a tip and what's an internal node
    is_tip  <- typ %in% c("tip", "sampAnc")
    is_node <- typ == "node"
    
    # store those values
    N  <- sum(is_tip)
    Mn <- sum(is_node)
    
    # total  number of nodes
    total <- N + Mn
    
    # ids for each node
    new_id <- integer(N_row)
    new_id[is_tip]  <- seq_len(N)
    new_id[is_node] <- N + seq_len(Mn)
    
    # tip labels
    tip_labels <- nm[is_tip]
    
    # parent index per row (NA for root)
    par_idx <- match(par, nm)
    
    # ensure order is correct
    stopifnot("buddPhylo rows are not in pre-order" =
                all(par_idx < seq_len(N_row), na.rm = TRUE))
    
    # traversal order
    po_order <- rev(seq_len(N_row))
    
    # make an SA map, if we want to consider SAs
    node_sa <- NULL
    if (considerSAs) {
      # get rows representing SA tips
      sa_rows <- which(typ == "sampAnc")
      
      # check if there are any
      if (length(sa_rows)) {
        # get the parents of each SA
        sa_par_rows <- par_idx[sa_rows]
        
        # avoid duplicated SA parents
        keep <- !duplicated(sa_par_rows)
        
        # set node_sa to which nodes are defined by SAs
        node_sa <- rep(NA_character_, total)
        node_sa[new_id[sa_par_rows[keep]]] <- nm[sa_rows[keep]]
      }
    }
  } else {
    stop("`tree` must be of class 'phylo' or 'buddPhylo'.")
  }
  
  # precompute order of tip labels
  rank_of_tip <- integer(N)
  rank_of_tip[order(tip_labels)] <- seq_len(N)
  
  # get members of each clade
  members <- vector("list", total)
  if (is_phylo) {
    # make trivial clades to make things easier
    members[seq_len(N)] <- as.list(rank_of_tip)
    
    # get the clades
    for (k in seq_len(nrow(edge_po))) {
      # parent and children node for this edge
      p <- edge_po[k, 1]
      c <- edge_po[k, 2]
      
      # add members of c clade to the p clade
      members[[p]] <- c(members[[p]], members[[c]])
    }
  } else {
    # get rows of tips
    tip_rows <- which(is_tip)
    
    # make trivial clades to make things easier
    for (r in tip_rows) members[[new_id[r]]] <- rank_of_tip[new_id[r]]
    
    # get the clades
    for (r in po_order) {
      # get parent for this clade
      pr <- par_idx[r]
      
      # we don't need to check the root since it has no parent
      if (is.na(pr)) next
      
      # add members of child clade to parent clade
      members[[new_id[pr]]] <- c(members[[new_id[pr]]], members[[new_id[r]]])
    }
  }
  
  # get sorted labels
  sorted_labels <- tip_labels[order(tip_labels)]
  
  # get internal labels
  internal <- N + seq_len(Mn)
  
  # make clade strings
  clades <- vapply(internal, function(nd) {
    paste(sorted_labels[sort.int(members[[nd]])], collapse = ", ")
  }, character(1))
  
  # check if we want to consider SAs
  if (considerSAs && !is.null(node_sa)) {
    # get the SA nodes
    sa <- node_sa[internal]
    
    # nodes for which we want SAs
    has <- !is.na(sa)
    
    # add to clade
    clades[has] <- paste0(clades[has], "\x01SA:", sa[has])
  }
  
  # check if we want to consider ancestors
  if (inherits(tree, "buddPhylo") && considerOrientation) {
    # start node list
    node_anc <- vector("list", total)
    
    # children of each parent
    ch_by_par <- split(seq_len(N_row), par_idx)
    
    # iterate through internal nodes
    for (r in which(is_node)) {
      # get the children of this node
      ch <- ch_by_par[[as.character(r)]]
      
      # if there are none, skip
      if (!length(ch)) next
      
      # if any of the children are sampled ancestors, skip (SA node)
      if (any(typ[ch] == "sampAnc")) next
      
      # get the ancestor child
      anc_ch <- ch[ori[ch] == "ancestor"]
      
      # if there are none, skip
      if (!length(anc_ch)) next
      
      # get the children of the ancestral branch
      anc_ranks <- members[[new_id[anc_ch[1]]]]
      
      # add to node_anc
      node_anc[[new_id[r]]] <- sorted_labels[sort.int(anc_ranks)]
    }
    
    # reorder ancestor nodes
    anc <- node_anc[internal]
    
    # check if there are any
    has <- lengths(anc) > 0
    if (any(has)) {
      # get the string of ancestral clade
      anc_str <- vapply(anc[has],
                        function(a) paste(a, collapse = "; "),
                        character(1))
      
      # add to clade
      clades[has] <- paste0(clades[has], "\x01ANC:", anc_str)
    }
  }
  
  # return clades
  return(clades)
}
