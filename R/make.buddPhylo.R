#' Generate a budding phylogenetic tree from simulations
#'
#' @description \code{make.buddPhylo} generates a budding phylogenetic tree
#' from a birth-death simulation (e.g. an output from \code{bd.sim}) and 
#' (optionally) a fossil record (e.g. an output from \code{sample.clade}).
#' 
#' @param returnTrueExt A logical indicating whether to include in tree the
#' tips representing the true extinction time of extinct species. If set to
#' \code{TRUE} (default), the returned tree will include those tips. If 
#' \code{FALSE}, they will be dropped and instead the last sampled fossil of
#' a given species will be the last sampled tip of that species. Note that if
#' a species was not sampled it will then not appear in the tree. If no fossils
#' have been added to the tree with the parameter \code{fossils}, this will
#' have the equivalent effect as the ape function \code{drop.fossil}, returning
#' an ultrametric tree.
#' 
#' @inheritParams make.phylo
#' 
#' @details .... TODO
#' 
#' @return ... TODO
#' 
#' @author Bruno do Rosario Petrucci
#' 
#' @name make.buddPhylo
#' @rdname make.buddPhylo
#' @export

make.buddPhylo <- function(sim, fossils = NULL, returnTrueExt = TRUE) {
  # check that sim is a valid sim object
  if (!is.sim(sim)) {
    stop("Invalid sim object")
  } 
  
  # start buddPhylo
  buddPhylo <- data.frame(matrix(nrow = 0, ncol = 13))
  colnames(buddPhylo) <- c("lineage", "taxon", "orientation", "length", "name",
                           "range", "parent", "type", "y_coord", "x_coord", 
                           "x_par", "y_par", "extant")
  
  # start internal node number
  n_node <- 1
  
  # flag to indicate we haven't set the UCA yet
  uca <- TRUE
  
  # iterate through lineages
  for (i in 1:length(sim$TS)) {
    # get all the speciation events that came from this lineage
    # exclude first since that's always NA (first lineage has no parent)
    lin_children_speciation <- sim$TS[sim$PAR == i][-1]
    
    # get all fossil samples for this lineage
    lin_samples <- fossils$SampT[fossils$Species == paste0("t", i)]
    
    # get the extinction time for this lineages
    lin_ext <- ifelse(is.na(sim$TE[i]), 0, sim$TE[i])
    
    # arrange it all in a data frame with all this lineage's events
    lin_events <- data.frame(times = c(lin_children_speciation,
                                       lin_samples,
                                       lin_ext),
                             sp = c(which(sim$TS %in% lin_children_speciation),
                                    rep(i, length(lin_samples)),
                                    i),
                             type = c(rep("speciation", length(lin_children_speciation)),
                                      rep("sampling", length(lin_samples)),
                                      "extinction"))
    
    # and sort it by time
    lin_events <- lin_events[order(lin_events$times, decreasing = TRUE), ]
    
    # if we don't want the true extinction time, remove
    # extinction event if lineage is extinct
    if (!returnTrueExt && lin_ext != 0) {
      lin_events <- lin_events[-nrow(lin_events), ]
      
      # if we have an unsampled lineage, skip
      if (nrow(lin_events) == 0) next
    }
    
    # starting orientation
    ori <- ifelse(uca, "uca", 
                  ifelse(i == 1, "ancestor", "descendant"))
    
    # never start within a range
    range <- NA
    
    # get parent lineage time
    t_parent <- sim$TS[i]
    
    # and parent node
    parent_node <- ifelse(uca, NA, 
                          ifelse(i == 1, n_tips + 1, 
                                 buddPhylo$name[buddPhylo$x_coord == t_parent]))
    
    # next fossil sample will be the first
    sample_id <- 1
    
    # iterate through all lineage events
    for (j in 1:nrow(lin_events)) {
      # get this line
      event <- lin_events[j, ]
      
      # and the time for the event
      time <- event$time
      
      # if this will be the UCA, we have to set the parent time
      if (uca) t_parent <- time
      
      # get which taxon this is
      tax <- paste0("t", i)
      
      # check the event type
      if (event$type == "speciation") {
        # check if this speciation would lead to no samples
        if (!returnTrueExt && 
            !check.speciation.sampled(which(sim$TS == time), sim, fossils)) {
          next
        }
        
        # add speciation node to buddPhylo
        buddPhylo <- rbind(buddPhylo,
                           as.data.frame(as.list(
                             c(lineage = n_node, taxon = n_node, 
                               orientation = ori, length = t_parent - time, 
                               name = n_node, range = range, parent = parent_node, 
                               type = "node", y_coord = i, 
                               x_coord = time, 
                               x_par = t_parent, y_par = i, 
                               extant = FALSE))))
        
        # set the parent node for the next event as this node
        parent_node <- n_node
        
        # increase node number by 1
        n_node <- n_node + 1
        
        # if this lineage generated a species, it is now an ancestor
        ori = "ancestor"
        
        # set parent time for the next event as this time
        t_parent <- time
        
        # if uca flag was TRUE, now it isn't
        if (uca) uca <- FALSE
      } else if (event$type == "sampling") {
        # check if this is the last event, or if it is the last sampled event
        last_sampled_event <- j == nrow(lin_events) ||
          (all(lin_events$type[(j + 1):nrow(lin_events)] == "speciation") &&
             all(!unlist(lapply(lin_events$sp[(j + 1):nrow(lin_events)],
                                check.speciation.sampled, sim = sim, 
                                fossils = fossils))))
        
        # if this is not the last event, we need an auxiliary node for the
        # sampled ancestor branch to be attached to
        if (!last_sampled_event) {
          buddPhylo <- rbind(buddPhylo,
                             as.data.frame(as.list(
                               c(lineage = n_node, taxon = n_node, 
                                 orientation = ori,  length = t_parent - time,
                                 name = n_node, range = range, parent = parent_node, 
                                 type = "node", y_coord = i, x_coord = time, 
                                 x_par = t_parent, y_par = i, 
                                 extant = FALSE))))
          
          # set this as the SA parent node
          parent_node <- n_node
          
          # next internal node numbers
          n_node <- n_node + 1
          
          # parent time for the next lineage is this time
          t_parent <- time
          
          # if uca flag is set, need to reset it and orientation
          if (uca) {
            uca <- FALSE
            ori <- "ancestor"
          }
        }
        
        # number of fossils for this species
        nFossilsSp <- sum(fossils$Species == tax)
        
        # numerical name
        numName <- i + sample_id / 10 ^ ceiling(log(nFossilsSp + 1, 10))
        
        # fossil name
        fossil_name <- paste0("t", numName, 
                              ifelse(sample_id %% 10 == 0, "0", ""))
        
        # get y coordinate for this lineage
        sp_y_coord <- ifelse(last_sampled_event, i, i - 0.5)
        
        # add sampled ancestor or fossil tip to buddPhylo
        buddPhylo <- rbind(buddPhylo,
                           as.data.frame(as.list(
                             c(lineage = fossil_name, 
                               taxon = tax, orientation = ori, 
                               length = t_parent - time,
                               name = fossil_name, 
                               range = range, parent = parent_node,
                               type = ifelse(last_sampled_event, 
                                             "tip", "sampAnc"), 
                               y_coord = sp_y_coord, 
                               x_coord = time, x_par = t_parent, y_par = i, 
                               extant = FALSE))))
        
        # if this is one of many sampling events, we're now in a range
        range <- ifelse(sum(lin_events$type %in% c("sampling", "extinction")) == 1, NA, 
                        paste0("t", i))
        
        # increase the id for the next sample
        sample_id <- sample_id + 1
      } else if (event$type == "extinction" && (returnTrueExt || time == 0)) {
        # if returnTrueExt is true or this is an extant lineage, include tip
        buddPhylo <- rbind(buddPhylo,
                           as.data.frame(as.list(
                             c(lineage = tax, taxon = tax, orientation = ori,
                               length = t_parent - time, name = tax,
                               range = range, parent = parent_node,
                               type = "tip", y_coord = i, x_coord = time,
                               x_par = t_parent, y_par = i, extant = (time == 0)))))
        # note that if returnTrueExt is false, the tip for this species
        # will be the last sampling event (if any)
      }
    }
  }
  
  # set length to numeric
  buddPhylo$length <- as.numeric(buddPhylo$length)
  
  # if returnTrueExt if false, we need to drop some rows
  if (!returnTrueExt) {
    rows_to_drop <- c()
    
    # iterate through buddPhylo
    for (i in 1:nrow(buddPhylo)) {
      # get the children of this node
      children <- getChildren.buddPhylo(buddPhylo, buddPhylo$name[i])
      
      # if there is only one child, we need to remove this node,
      # as nodes should have either one or zero children
      if (length(children) == 1) {
        # get the indices for these children
        child_indices <- which(buddPhylo$name == children)
        
        # correct the coordinates and parent identities
        buddPhylo$y_coord[child_indices] <- buddPhylo$y_coord[i]
        buddPhylo$y_par[child_indices] <- buddPhylo$y_par[i]
        buddPhylo$parent[child_indices] <- buddPhylo$parent[i]
        
        # correct length and orientation
        buddPhylo$length[child_indices] <- buddPhylo$length[child_indices] + 
          buddPhylo$length[i]
        buddPhylo$x_par[child_indices] <- buddPhylo$x_par[i]
        buddPhylo$orientation[child_indices] <- buddPhylo$orientation[i]
        
        # get the children of the child
        next_children <- getChildren.buddPhylo(buddPhylo, children)
        
        # if there are any, we need to correct their properties too
        while (length(next_children) > 0) {
          # get indices for these children
          nxt_child_indices <- buddPhylo$name %in% next_children
          
          # correct y coordinates and orientation 
          if (all(names(next_children) == c("sampAnc", "lineage"))) {
            buddPhylo$orientation[nxt_child_indices] <- buddPhylo$orientation[i]
            buddPhylo$y_coord[nxt_child_indices] <- buddPhylo$y_coord[i]
            
            # get the next children
            next_children <- getChildren.buddPhylo(buddPhylo, next_children[2])
          } else if (all(names(next_children) == c("ancestor", "descendant"))) {
            buddPhylo$y_coord[buddPhylo$name == next_children[1]] <- buddPhylo$y_coord[i]
            
            # continue only for ancestor side
            next_children <- getChildren.buddPhylo(buddPhylo, next_children[1])
          } else {
            # means we've found another 1-child node, it will be
            # taken care of later
            break
          }
        }
        
        # add to rows to drop
        rows_to_drop <- c(rows_to_drop, i)
      }
    }
    
    # if there are any unsampled nodes, drop them
    if (length(rows_to_drop) > 0) buddPhylo <- buddPhylo[-rows_to_drop, ]
  }
  
  # correct y_par
  buddPhylo$y_par <- c(1, unlist(lapply(buddPhylo$parent[-1],
                                        function(x) buddPhylo$y_coord[buddPhylo$name == x])))
  
  # preserve y_coords and y_par
  y_coords <- buddPhylo$y_coord
  y_par <- buddPhylo$y_par
  
  # start at the uca
  cur_node <- 1
  new_cur_node <- sum(buddPhylo$type %in% c("tip", "sampAnc")) + 1
  
  # start with the root/uca node renaming
  node_rename <- data.frame(old = cur_node, new = new_cur_node)
  
  # get all the nodes we need to rename
  nodes_to_rename <- buddPhylo$name[buddPhylo$type == "node"][-1]
  
  # while that's not empty
  while (length(nodes_to_rename) > 0) {
    # get children of current node
    children <- getChildren.buddPhylo(buddPhylo, cur_node)
    
    # check how many descendants each child has
    n_desc <- unlist(lapply(children, function(x) 
      length(getDescendants.buddPhylo(buddPhylo, x))))
    
    # if both are 0, cur node will be the next node to rename
    if (sum(n_desc) == 0) {
      cur_node <- nodes_to_rename[1]
    } else {
      # otherwise, it will be the one with the most descendants
      cur_node <- children[which.max(n_desc)]
      
      # if the other child has descendants, send it to the beginning
      if (min(n_desc) > 0) {
        nodes_to_rename <- c(children[which.min(n_desc)],
                             nodes_to_rename[nodes_to_rename != children[which.min(n_desc)]])
      }
    }
    
    # new current node
    new_cur_node <- new_cur_node + 1
    
    # add to rename data frame
    node_rename <- rbind(node_rename,
                         c(old = cur_node, new = new_cur_node))
    
    # remove from list
    nodes_to_rename <- nodes_to_rename[-which(nodes_to_rename == cur_node)]
  }
  
  # get a lookup
  lookup <- setNames(node_rename$new, node_rename$old)
  
  # match
  m <- match(as.matrix(buddPhylo), node_rename$old)
  
  # substitute values
  buddPhylo[] <- ifelse(is.na(m),
                        as.matrix(buddPhylo),
                        node_rename$new[m])
  
  # redo y coords
  buddPhylo$y_coord <- y_coords
  buddPhylo$y_par <- y_par
  
  # make equivalent bifurcating phylogeny
  if (is.null(fossils)) {
    if (returnTrueExt) {
      # if we want true extinction times and don't have fossils,
      # just get full phylogeny
      phylo <- make.phylo(sim)
    } else {
      # if we want don't want true extinction times, just get the
      # extant phylogeny
      phylo <- ape::drop.fossil(make.phylo(sim))
    }
  } else {
    # otherwise, get phylogeny based on other parameters
    phylo <- make.phylo(sim, fossils, returnTrueExt = returnTrueExt)
  }
  
  # fix y coordinates
  buddPhylo <- fix.coords(buddPhylo, phylo,
                          fix_x = FALSE)
  
  # correct type of some columns
  buddPhylo$extant <- as.logical(buddPhylo$extant)
  buddPhylo$length <- as.numeric(buddPhylo$length)
  buddPhylo$x_coord <- as.numeric(buddPhylo$x_coord)
  buddPhylo$y_coord <- as.numeric(buddPhylo$y_coord)
  buddPhylo$x_par <- as.numeric(buddPhylo$x_par)
  buddPhylo$y_par <- as.numeric(buddPhylo$y_par)
  
  # set classes
  class(buddPhylo) <- c("buddPhylo", "data.frame")  
  
  # return
  return(buddPhylo)
}

# auxiliary function to check which speciation events have left samples
check.speciation.sampled <- function(lin, sim, fossils = NULL) {
  # get the collection of all species that are descended from this lineage
  descendants <- getDescendants.sim(sim, lin)
  
  # get all fossil samples from these lineages
  samples <- which(fossils$Species %in% paste0("t", c(lin, descendants)))
  
  # check which of these lineages are extant
  ext_samples <- which(sim$EXTANT[c(lin, descendants)])
  
  # if either of those are higher than 0, there are descendants of this lineage
  if (length(c(samples, ext_samples)) > 0) TRUE else FALSE
}