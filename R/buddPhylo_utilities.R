#' @rdname buddPhylo
#'
#' @details \code{getParents.buddPhylo} Finds all the parental lineages of a
#' \code{focalLineage} from a \code{buddPhylo$name} object, from first to
#' last ancestor. Note this does not use the information in the
#' \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#' 

getParents.buddPhylo <- function(buddPhylo, focalLineage) {
  # make a map of parent to child index
  par_idx <- match(buddPhylo$parent, buddPhylo$name)
  
  # get the index of the focal lineage
  i <- match(focalLineage, buddPhylo$name)
  
  # start the vector of parents
  res <- character(0)
  
  # while we haven't reached the root...
  while (!is.na(par_idx[i])) {
    # make the parent of the focal lineage the focal lineage
    i <- par_idx[i]
    
    # add number to result
    res <- c(res, buddPhylo$name[i])
  }
  
  return(res)
}

#' @rdname buddPhylo
#' 
#' @details \code{getChildren.buddPhylo} Returns the immediate descendant
#' (daughter) lineages of a given \code{focalLineage} in a \code{buddPhylo}
#' object (using the \code{name} column). The output is a character vector of
#' length two containing the \code{name} of each child, equivalent to
#' \code{getDescendants.buddPhylo(..., onlyImmediates = TRUE)}.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#' 

getChildren.buddPhylo <- function(buddPhylo, focalLineage, returnInfo = TRUE) {
  # find the lineages with this lineage as parent
  desc <- buddPhylo$name[buddPhylo$parent == focalLineage]
  desc <- desc[!is.na(desc)]
  res <- desc
  
  # if res is of length 0, this is a node with no children, return NULL
  if (!length(res)) return(NULL)
  
  # check if we want to return the status of each child
  if (returnInfo) {
    # check if any of the children have 0 length branches
    if (0 %in% buddPhylo$length[buddPhylo$name %in% desc]) {
      # if so, name the one with 0 length branch the sampled ancestor
      desc_len <- buddPhylo$length[match(desc, buddPhylo$name)]
      desc_len <- ifelse(desc_len == 0, yes = "sampAnc", no = "lineage")
      
      # set names
      names(res) <- desc_len
      
      # order with SA first
      res <- res[order(names(res), decreasing = TRUE)]
    } else {
      # otherwise, just order based on orientation (ancestor first)
      desc_ori <- buddPhylo$orientation[match(desc, buddPhylo$name)]
      names(res) <- desc_ori
      res <- res[order(names(res))]
    }
  }
  
  return(res)
}

#' @rdname buddPhylo
#'
#' @details \code{getDescendants.buddPhylo} Find all the descendant lineages of
#' a \code{focalLineage} from a \code{buddPhylo$name} object, from most 
#' immediate to most distal descendant. \code{getChildren.buddPhylo} returns 
#' only the most immediate descendants (same effect as
#' \code{getDescendants.buddPhylo(..., onlyImmediates = TRUE)}). Note this does 
#' not use the information in the \code{buddPhylo$taxon} column.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#' 

getDescendants.buddPhylo <- function(buddPhylo, focalLineage, 
                                     onlyImmediates = FALSE, 
                                     internalNodes = FALSE) {
  # initialize variables for bookkeeping
  toCheck <- focalLineage
  res <- vector()
  
  # while we still have lineages to check
  while (length(toCheck) > 0) {
    # get the first one
    focalLineage <- toCheck[1]
    
    # get the children of this lineage
    desc <- buddPhylo$name[buddPhylo$parent == focalLineage]
    desc <- desc[!is.na(desc)]
    
    # check if there are any children
    if (length(desc) > 0) {
      # if so, add those to the toCheck and result
      toCheck <- c(toCheck, desc)
      res <- c(res, desc)
    }
    
    # check if we only want the immediate children
    if (onlyImmediates) {
      break
    } else {
      # otherwise, remove focalLineage from toCheck
      toCheck <- toCheck[!toCheck %in% focalLineage]
    }
  }
  
  # if we don't want internal nodes, remove them from result
  if (!internalNodes) {
    internal_nodes <- which(res %in% buddPhylo$name[buddPhylo$type == "node"])
    if (length(internal_nodes) > 0) res <- res[-internal_nodes]
  }
  
  # if res is of length 0, this is a node with no children, return NULL
  if (!length(res)) return(NULL)
  
  return(res)
}

#' @rdname buddPhylo
#'
#' @details \code{getMRCA.buddPhylo} Finds the node representing the most recent
#' common ancestor of the tips in \code{tipList}.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#' 

getMRCA.buddPhylo <- function(buddPhylo, tipList) {
  # check that tipList is in the phylogeny
  if (any(!(tipList %in% buddPhylo$name))) {
    stop("Tips tipList must be in buddPhylo$name.")
  }
  
  # if tipList is of length one, MRCA is itself
  if (length(tipList) == 1) return(tipList)
  
  # create the list of parents
  ParentList <- list()
  
  # iterate through list of tips
  for (i in 1:length(tipList)) {
    # get the parents of each of the tips in tipList
    ParentList[[i]] <- getParents.buddPhylo(buddPhylo, tipList[i])
  }
  
  # find longest parent-descendant length
  lengs <- unlist(lapply(ParentList, length))
  ref <- which.max(lengs)
  
  # function to check if an element is in a list
  listChecker <- function(pars, ca) {
    return(ca %in% pars)
  }
  
  # start at the root, check all the elements 
  # in longest parent-descendant length
  path_to_root <- toCheck <- rev(unlist(ParentList[ref]))
  
  # check that the first node (i.e. the root) is in the path
  node_in_path <- all(unlist(lapply(ParentList, listChecker, toCheck[1])))
  
  # while the node in question is in the path
  while (node_in_path) {
    # remove checked node
    toCheck <- toCheck[-1]
    
    # check again
    node_in_path <- all(unlist(lapply(ParentList, listChecker, toCheck[1])))
  }
  
  # remove the checked nodes from the path to the root
  path_to_root <- path_to_root[!(path_to_root %in% toCheck)]
  
  # the MRCA is the last node in the remaining path
  # i.e. the last node that was in both parent lists
  MRCA <- rev(path_to_root)[1]
  return(MRCA)
}


#' @rdname buddPhylo
#'
#' @details \code{is.monophyletic.buddPhylo} Checks if the tips in
#' \code{tipList} constitute a monophyletic group in \code{buddPhylo}.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#' 

is.monophyletic.buddPhylo <- function(buddPhylo, tipList,
                                      excludeSampAnc = TRUE) {
  # get the MRCA of this list of tips
  MRCA <- getMRCA.buddPhylo(buddPhylo, tipList)
  
  # get all the descendants of the MRCA
  alldesc <- getDescendants.buddPhylo(buddPhylo, MRCA)
  
  # check if we want to exclude sampled ancestors
  if (excludeSampAnc) {
    # exclude sampled ancestors
    alldesc <- alldesc[!(alldesc %in% 
                           buddPhylo$name[buddPhylo$type == "sampAnc"])]
  }
  
  # ensure that all descendants of the MRCA and tipList are the same
  res <- identical(sort(alldesc), sort(tipList))
  
  return(res)
}

#' Generating a phylo object from a buddPhylo object
#'
#' Generates a phylogeny in a \code{phylo} format based on the information
#' within a \code{buddPhylo} object.
#' 
#' @param x A \code{buddPhylo} object.
#' 
#' @param ... Further arguments inherited from generics.
#' 
#' @return A \code{phylo} object compatible with functions from the \code{ape}
#' package.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#' 
#' @export
#' 
#' @importFrom ape as.phylo
#' 

as.phylo.buddPhylo <- function(x, ...) {
  # assign to buddPhylo for clarity
  buddPhylo <- x
  
  # function to re-order edges so the phylo object does not break
  reorder_edges <- function(edge) {
    # get parents and children from edge matrix
    parents <- edge[, 1]
    children <- edge[, 2]
    
    # find the root
    root <- setdiff(parents, children)
    
    # separate children edges into parent edges they connect to
    adj <- split(children, parents)
    
    # make a matrix for new edges
    new_edge <- matrix(ncol = 2, nrow = 0)
    
    # depth-first search
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
    
    # get the dfs for the root
    dfs(root)
    
    return(new_edge)
  }

  # get tip and node labels
  tip.label <- buddPhylo$name[buddPhylo$type %in% c("tip", "sampAnc")]
  node.label <- buddPhylo$name[!buddPhylo$type %in% c("tip", "sampAnc")]
  
  # match parents and children
  par_idx  <- match(buddPhylo$parent, buddPhylo$name)
  kids <- split(seq_len(nrow(buddPhylo)), par_idx)
  
  # start a node index map
  node.id <- integer(length(node.label))
  names(node.id) <- node.label
  nextNode <- length(tip.label) + 1
  
  # start stack at root
  stack <- which(is.na(par_idx))
  
  # while the stack has something in it
  while (length(stack)) {
    # get first element of stack and remove it
    r <- stack[1]
    stack <- stack[-1]
    
    # number this row if it is an internal node
    if (buddPhylo$type[r] == "node") {
      node.id[[buddPhylo$name[r]]] <- nextNode
      nextNode <- nextNode + 1
    }
    
    # push children, preserving row order
    ch <- kids[[as.character(r)]]
    if (!is.null(ch)) stack <- c(ch, stack)
  }
  
  # get labels and ids for these labels
  all.labels <- c(tip.label, names(node.id))
  all.ids <- c(1:length(tip.label), node.id)
  
  # create edge matrix
  edge1 <- node.id[match(buddPhylo$parent, node.label)]
  edge2 <- all.ids[match(buddPhylo$name, all.labels)]
  edge <- cbind(as.integer(edge1), as.integer(edge2))
  edge <- edge[-which(is.na(edge[, 1])), ] 
  # NA element is the root, ignore it
  
  # re-order phylo edges
  edge <- reorder_edges(edge)
  
  # set edge lengths
  edge_names <- all.labels[match(edge[, 2], all.ids)]
  edge.length <- buddPhylo$length[match(edge_names, buddPhylo$name)]
  
  # assemble the tree
  phy <- list(
    edge = edge, edge.length = edge.length, Nnode = length(node.label),
    node.label = node.label, tip.label = tip.label,
    root.edge = buddPhylo$length[buddPhylo$orientation == "uca"]
  )
  
  # set class
  class(phy) <- "phylo"
  
  return(phy)
}

#' Adjust the timescale of a buddPhylo object
#'
#' @details \code{adjust.timescale.buddPhylo} Adjusts the time coordinates of a
#' buddPhylo object, either with time going from the root or from the present 
#' (based on an argument).
#'
#' @param buddPhylo Object of class \code{buddPhylo}.
#'
#' @return A \code{buddPhylo} with the time adjusted accordingly.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @noRd
#'

adjust.timescale.buddPhylo <- function(buddPhylo, timeFromRoot = FALSE) {
  # first check that it is a valid buddPhylo
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  # get root time
  root_x <- buddPhylo$x_coord[buddPhylo$orientation == "uca"]
  
  # check if time is from root
  is_time_from_root <- root_x < 10E-10
  
  # check if we need correction
  if (is_time_from_root != timeFromRoot) {
    # find maximum time
    refMax <- max(c(buddPhylo$x_coord, buddPhylo$x_par), na.rm = TRUE)
    
    # reverse time
    buddPhylo$x_coord <- refMax - buddPhylo$x_coord
    buddPhylo$x_par  <- refMax - buddPhylo$x_par
  }
  
  return(buddPhylo)
}

#' All subtrees of a budding phylogenetic tree
#'
#' This function returns a list of all the subtrees of a budding phylogenetic
#' tree, i.e. a list of \code{buddPhylo} objects.
#'
#' @param buddPhylo Object of class "buddPhylo".
#' 
#' @return A list of \code{buddPhylo} objects, one per internal node of 
#' \code{buddPhylo}, each rooted at that node. Each subtree will be rerooted,
#' i.e. the node that is the root will have \code{parent=NA},
#' \code{orientation = "uca"}, and \code{length=0}.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @export
#'

subtrees.buddPhylo <- function(buddPhylo) {
  # check it is a valid buddPhylo object
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  # get number of rows
  N <- nrow(buddPhylo)
  
  # get parents for each node
  par_idx <- match(buddPhylo$parent, buddPhylo$name)
  
  # make sure it is in pre-order (should be if built by paleobuddy)
  stopifnot("buddPhylo rows are not in pre-order" =
              all(par_idx < seq_len(N), na.rm = TRUE))
  
  # get internal node rows
  node_rows <- which(buddPhylo$type == "node")
  
  # descendant row indices per row, including self
  desc <- as.list(seq_len(N))
  for (i in rev(seq_len(N))) {
    # get parent of this node
    p <- par_idx[i]
    
    # add to descendants
    if (!is.na(p)) desc[[p]] <- c(desc[[p]], desc[[i]])
  }
  
  # build list of trees
  lapply(node_rows, function(r) {
    # get the subtree
    sub <- buddPhylo[desc[[r]], , drop = FALSE]
    
    # find root
    root_i <- match(buddPhylo$name[r], sub$name)
    
    # fix root values
    sub$parent[root_i] <- NA_character_
    sub$orientation[root_i] <- "uca"
    
    # rename rows
    rownames(sub) <- NULL
    
    # zero root length
    sub$length[root_i] <- 0
    sub$x_par[root_i] <- sub$x_coord[root_i]
    
    # set classes
    class(sub) <- c("buddPhylo", "data.frame")
    
    # return
    sub
  })
}

#' Get clades from buddPhylos and phylos
#'
#' Get a list of monophyletic clades from a phylogenetic tree. Could be budding
#' or bifurcating. 
#'
#' @param tree Object of class \code{buddPhylo} or \code{phylo}.
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
#' @author Bruno do Rosario Petrucci.
#' 
#' @export
#' 

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

#' 
#' @details \code{fix.coords} Fixes the x and y coordinates of a buddPhylo 
#' object based on its associated phylo object.
#' 
#' @author Bruno do Rosario Petrucci.
#' 
#' @noRd
#' 

fix.coords <- function(buddPhylo, phylo, fix_x = TRUE) {
  # ladderize phylo
  phylo <- ape::ladderize(phylo)
  
  # ensure length is numeric
  buddPhylo$length <- as.numeric(buddPhylo$length)
  
  # build row->name lookup to make things faster later
  parent_idx <- match(buddPhylo$parent, buddPhylo$name)
  
  # function to find the y coordinates of the tips of a phylo object
  get_tip_y_coords <- function(phylo) {
    # collect useful variables
    n_tips <- length(phylo$tip.label)
    n_nodes <- phylo$Nnode
    n_total <- n_tips + n_nodes
    edge <- phylo$edge
    edge_lengths <- phylo$edge.length
    
    # zero-length edges
    zero_len_edges <- which(!is.null(edge_lengths) & edge_lengths == 0)
    
    # true tips: tip indices NOT reachable only via zero-length edges
    true_tips <- setdiff(1:n_tips, edge[zero_len_edges, 2])
    
    # assign y positions in the order they appear in the tip.label vector,
    # respecting the ladderized traversal order
    yy <- numeric(n_total)
    
    # make map of children per parent
    children_by_parent <- split(edge[, 2], edge[, 1])
    
    # iterative pre-order traversal (left to right = first child first)
    stack <- phylo$edge[phylo$edge[, 1] == (n_tips + 1), 2]
    
    # reverse so first child is on top of stack
    stack <- rev(stack)
    
    # create vector to hold visit order
    visit_order <- integer(n_tips)
    
    # start at the beginning of the stack
    k <- 1
    stk <- stack
    
    # while we have a stack
    while (length(stk) > 0) {
      # get this last node on the stack
      node <- stk[length(stk)]
      
      # remove it from the stack
      stk <- stk[-length(stk)]
      
      # check if it is a tip
      if (node <= n_tips) {
        # if so, add it to visit order
        visit_order[k] <- node
        
        # increment order
        k <- k + 1
      } else {
        # if not, get its children
        children <- children_by_parent[[as.character(node)]]
        
        # add children to stack
        stk <- c(stk, rev(children))
      }
    }
    
    # assign tip number in traversal order
    for (i in seq_along(visit_order)) {
      yy[visit_order[i]] <- i
    }
    
    # process in reverse of edge order (post-order)
    node_order <- rev(unique(edge[, 1]))
    for (nd in node_order) {
      # get children
      children <- children_by_parent[[as.character(nd)]]
      
      # set y based on the mean of children ys
      yy[nd] <- mean(yy[children])
    }
    
    # filter to only true tips if zero-length branches exist
    result <- yy[true_tips]
    names(result) <- phylo$tip.label[true_tips]
    
    # return result
    return(result)
  }
  
  # if fix_x is false, need to save x coordinates
  if (!fix_x) {
    x_coords_bak <- buddPhylo$x_coord
    x_par_bak <- buddPhylo$x_par
  }
  
  # get the y coordinates for the tips of the object
  yy <- get_tip_y_coords(phylo)
  
  # get x coordinates for these tips
  xx <- ape::node.depth.edgelength(phylo)[which(phylo$tip.label %in% names(yy))]
  
  # start coordinates as NA
  buddPhylo$y_coord <- NA
  buddPhylo$x_coord <- NA
  
  # get the match between phylo and buddPhylo nodes
  matching_ids <- match(names(yy), buddPhylo$name)
  
  # set the coordinates for tips
  buddPhylo$y_coord[matching_ids] <- yy
  buddPhylo$x_coord[matching_ids] <- xx
  
  # start the parent coordinates as NA
  buddPhylo$x_par <- NA
  buddPhylo$y_par <- NA
  
  # check which already have a set y coordinate
  buddPhylo$y_set <- !is.na(buddPhylo$y_coord)
  
  # precompute row ids of the children for each parent
  children_by_par <- split(seq_len(nrow(buddPhylo)), parent_idx)
  
  # check which are tips
  tip_ids <- which(buddPhylo$type == "tip")
  
  # iterate through tips
  for (ttID in tip_ids) {
    # set todo
    todo <- TRUE
    
    # while we need to set these coordinates...
    while (todo) {
      # get parent and grandparent IDs
      ParID <- parent_idx[ttID]
      grandParID <- parent_idx[ParID]
      
      # check if parent ID is set
      if (buddPhylo$y_set[ParID]) {
        # if so, parent x-coordinate will just be 
        # based on this tip's x-coordinate and length
        buddPhylo$x_par[ttID] <- buddPhylo$x_coord[ttID] - buddPhylo$length[ttID]
        todo <- FALSE
      } else {
        # get siblings of this tip
        sibs <- children_by_par[[as.character(ParID)]]
        sibs <- sibs[sibs != ttID]
        
        # find non-SA siblings that are ancestors
        anc_sibs <- sibs[buddPhylo$orientation[sibs] == "ancestor" &
                           buddPhylo$type[sibs] != "sampAnc"]
        
        # check if this tip has an ancestor sibling we have not hit yet
        hasUnwalkedAncestorSibling <-
          length(anc_sibs) > 0 && any(is.na(buddPhylo$x_par[anc_sibs]))
        
        # check if this is a descendant tip with an unwalked ancestor sib
        if (buddPhylo$orientation[ttID] == "descendant" && 
            hasUnwalkedAncestorSibling) {
          # if so, set parent coordinates based on this tip
          buddPhylo$x_coord[ParID] <- buddPhylo$x_coord[ttID] - 
            buddPhylo$length[ttID]
          buddPhylo$x_par[ttID] <- buddPhylo$x_coord[ttID] - 
            buddPhylo$length[ttID]
          todo <- FALSE
        } else {
          # otherwise, set y coordinates too
          buddPhylo$y_coord[ParID] <- buddPhylo$y_coord[ttID]
          buddPhylo$x_coord[ParID] <- buddPhylo$x_coord[ttID] - 
            buddPhylo$length[ttID]
          buddPhylo$x_par[ttID] <- buddPhylo$x_coord[ttID] - buddPhylo$length[ttID]
          buddPhylo$y_set[ParID] <- TRUE
        }
        # in this way we avoid setting parent coordinates to
        # the coordinates of a descendant child
      }
      
      # check if there is a grandparent
      if (is.na(grandParID) || length(grandParID) == 0) {
        # if not, we have reached the node for which the parent is the root
        # set the parent coordinates and stop
        buddPhylo$x_par[ParID] <- buddPhylo$x_coord[ParID] - 
          buddPhylo$length[ParID]
        buddPhylo$y_par[ParID] <- buddPhylo$y_coord[ParID]
        todo <- FALSE
      } else {
        # otherwise, continue to the parent
        ttID <- ParID
      }
    }
  }
  
  # reset y_set to NULL
  buddPhylo$y_set <- NULL
  
  # now get the sampled ancestor nodes
  sa_ids <- which(buddPhylo$type == "sampAnc")
  sa_par_ids <- parent_idx[sa_ids]
  
  # set their coordinates and parent coordinates
  buddPhylo$y_coord[sa_ids] <- buddPhylo$y_coord[sa_par_ids]
  buddPhylo$x_coord[sa_ids] <- buddPhylo$x_coord[sa_par_ids]
  buddPhylo$x_par[sa_ids] <- buddPhylo$x_coord[sa_par_ids]
  
  # get minimum x (present)
  min_x <- min(c(buddPhylo$x_coord, buddPhylo$x_par), na.rm = TRUE)
  
  # set time scale based on min_x
  buddPhylo$x_coord <- buddPhylo$x_coord - min_x
  buddPhylo$x_par <- buddPhylo$x_par - min_x
  
  # match y_par to parent coordinates
  buddPhylo$y_par <- buddPhylo$y_coord[match(buddPhylo$parent, buddPhylo$name)]
  buddPhylo$y_par[buddPhylo$orientation == "uca"] <- 
    buddPhylo$y_coord[buddPhylo$orientation == "uca"]
  
  # restore x coords if we wanted to not mess with it
  if (!fix_x) {
    buddPhylo$x_coord <- x_coords_bak
    buddPhylo$x_par <- x_par_bak
  }
  
  # fix row types
  buddPhylo$length <- as.numeric(buddPhylo$length)
  buddPhylo$x_coord <- as.numeric(buddPhylo$x_coord)
  buddPhylo$y_coord <- as.numeric(buddPhylo$y_coord)
  
  # return
  return(buddPhylo)
}

#' @rdname buddPhylo
#' 
#' @details \code{node.time.buddPhylo} Finds the time of a node with a given
#' \code{name} field.
#' 
#' @author Bruno do Rosario Petrucci.
#' 
#' @export
#' 

node.time.buddPhylo <- function(buddPhylo, node) {
  # check that the node exists in buddPhylo
  if (!(node %in% buddPhylo$name)) {
    stop("node must be in buddPhylo$name")
  }
  
  # get the time
  buddPhylo$x_coord[buddPhylo$name == node]
}