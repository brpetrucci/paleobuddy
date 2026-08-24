#' Details, generics, and methods for the \code{buddPhylo} class
#'
#' @description The \code{buddPhylo} class is a format that stores the
#' information from a phylogeny that contains budding speciation events. Every
#' \code{buddPhylo} is also a \code{data.frame}, and every line refers to a
#' branch (and the node at the tip of said branch) in the phylogeny. Every 
#' node (i.e., "line") may contain a myriad of types of information ("columns"),
#' but the following elements organize the core information regarding every node
#' within a \code{buddPhylo} object:
#'
#' \describe{
#' \item{\code{name}}{Name of the node. Must be a character column used to
#' trace evolutionary relationships among edges. It can be any character
#' (including numbers), as long as it is consistently named across each
#' relevant branches/edge connection.}
#' \item{\code{length}}{Numeric column containing the length of the edge that
#' lead to the node. All functions assume time is measured in million years
#' since the present.}
#' \item{\code{orientation}}{Whether this is an ancestor or descendant branch.
#' This is where these trees differ from traditional bifurcating trees. After
#' every speciation event, one of the branches is assigned an orientation, with
#' ancestor branches assigned to the progenitor species, and descendant branches
#' assigned to the generated species.}
#' \item{\code{taxon}}{Name of the species that should be associated with the 
#' node (e.g., in plotting functions). Could be the same as \code{name}, but
#' only terminal branches (tips and sampled ancestors) will have a \code{taxon}
#' field.}
#' \item{\code{parent}}{\code{name} of the parent node of the node in question.
#' Note this refers to \code{name}, not to \code{taxon}. The root node has
#' parent \code{NA}.}
#' \item{\code{type}}{Character assigiong whether the node is a \code{node},
#' a \code{tip}, a sampled ancestor (\code{sampAnc}), or the root (\code{uca}).}
#' \item{\code{y_coord}}{Y-coordinate to be used as a reference for that node
#' for the plotting functions.}
#' \item{\code{x_coord}}{X-coordinate to be used as a reference for that node
#' for the plotting functions.}
#' \item{\code{extant}}{A logical indicating whether the tip is alive at the
#' present (\code{TRUE}) or not (\code{FALSE}).}}
#'
#' Here we declare useful generics and methods for \code{buddPhylo} objects.
#'
#' @param buddPhylo,x,object Object of class "buddPhylo".
#' 
#' @param focalLineage A node/lineage to be used as reference for the tree operation.
#'
#' @param tipList A character vector with a list of \code{name}s (i.e., not
#' "taxon") to use as a reference for the node set operation. These could be
#' tips or sampled ancestors in the tree, but not internal nodes.
#'
#' @param excludeSampAnc Logical indicating whether sampled ancestors should be
#' excluded when considering a monophyletic group. Default is \code{TRUE}.
#'
#' @param onlyImmediates Logical indicating if only the most immediate
#' descendant (the two daughter) lineages of a given \code{focalLineage} should
#' be considered.
#' 
#' @param node A character representing the name of a node.
#' 
#' @param returnInfo Logical on whether to return information on each
#' child. If \code{TRUE}, the function also annotates each descendant based on 
#' its status. If \code{focalLineage} is a sampled ancestor node, one child is
#' labeled "sampAnc" and the other "lineage". If it instead is a speciation 
#' node, descendants are labeled according to their \code{orientation}
#' ("ancestor" or "descendant"). Note that this function does not use 
#' information from the \code{buddPhylo$taxon} column.
#' 
#' @param internalNodes Logical on whether to return internal node descendants.
#' If \code{FALSE} (default), only tips and sampled ancestors will be returned.
#'
#' @param ... Further arguments inherited from generics.
#' 
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#' 
#' @examples
#' 
#' # set seed
#' set.seed(1)
#' 
#' # run a simple birth-death simulation
#' sim <- bd.sim(1, 0.2, 0.05, 20)
#' 
#' # and sample fossils on it
#' sample <- suppressMessages(sample.clade(sim, 0.5, 20))
#' 
#' ###
#' # plotting, printing, summarizing
#' 
#' # build and plot the tree
#' budd <- make.buddPhylo(sim, sample, returnTrueExt = FALSE)
#' plot(budd, show.SA.label = TRUE)
#' 
#' # check out print and summary
#' print(budd)
#' summary(budd)
#' 
#' ###
#' # learning about the relationships in the tree
#' 
#' # check the MRCA of t13 and t1
#' mrca <- getMRCA.buddPhylo(budd, c("t1", "t13"))
#' # should be the first speciation node
#' 
#' # what about the parent of this node?
#' par <- getParents.buddPhylo(budd, mrca)[1]
#' 
#' # we can check what kind of node it is by looking at its children
#' children <- getChildren.buddPhylo(budd, par)
#' # children[1] is an SA, the third t1 occurrence
#' 
#' # which species are descended from the descendant lineage
#' # in that first budding speciation?
#' desc_node <- getChildren.buddPhylo(budd, mrca)["descendant"]
#' descendants <- getDescendants.buddPhylo(budd, desc_node)
#' 
#' # since this is all descendants from a node, it makes a monophyletic group
#' is.monophyletic.buddPhylo(budd, descendants)
#' # FALSE because by default we only consider tips
#' is.monophyletic.buddPhylo(budd, descendants, excludeSampAnc = FALSE)
#' 
#' # if you'd like to compare the budding and bifurcating trees,
#' # you can make use of the conversion function
#' phy <- ape::as.phylo(budd)
#' 
#' ###
#' # advanced clade operations
#' 
#' # similarly to ape::subtrees, you can extract all subtrees of a buddPhylo
#' sub_budd <- subtrees.buddPhylo(budd)
#' 
#' # of course it's hard to comb through all clades of a tree
#' # by hand, so we have a nice wrapper to easily extract the
#' # relevant info for each clade implied by your tree
#' clades <- extract.clades(budd)
#' 
#' # you can get more info by using the arguments within
#' sa_clades <- extract.clades(budd, 
#'                             considerSAs = TRUE)
#' oriented_clades <- extract.clades(budd, 
#'                                   considerSAs = TRUE, 
#'                                   considerOrientation = TRUE)
#' # the SA species or ancestor lineage is included in the output now
#'
#' @name buddPhylo
#'
#' @importFrom graphics plot par
#' @importFrom utils head tail
#'
NULL

#' @rdname buddPhylo
#'
#' @details \code{is.buddPhylo} A \code{buddPhylo} object must contain all 9
#' members: \code{name}, \code{length}, \code{orientation}, \code{taxon},
#' \code{parent}, \code{type}, \code{y_coord}, \code{x_coord}. All of these 
#' must have the correct length (i.e. \code{2 * Ntip - 1)}, and correct class
#' (\code{character}, \code{logical}, or \code{numeric} depending on the 
#' member).
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#'

is.buddPhylo <- function(x) {
  # call the object buddPhylo for clarity
  buddPhylo <- x
  
  # checks that the class is buddPhylo
  cla <- "buddPhylo" %in% class(buddPhylo)
  
  # if it is not, should be
  if (!cla) {
    return(FALSE)
  }
  
  # checks that there is minimal information
  mustHaveCols <- c("name", "length", "orientation", "taxon", 
                    "parent", "type", "y_par", "x_par",
                    "y_coord", "x_coord")
  cols <- all(mustHaveCols %in% colnames(buddPhylo))
  
  # if there are not, the other tests are redundant
  if (!cols) {
    message(paste0(
      "\"buddPhylo\" is lacking the columns: ",
      paste0(mustHaveCols[!mustHaveCols %in% colnames(buddPhylo)], 
             collapse = ", ")
    ))
    return(FALSE)
  }
  
  # check the classes of each column
  cmatch <- match(mustHaveCols, colnames(buddPhylo))
  chr_cols <- c("name", "orientation", "taxon", "parent", "type")
  num_cols <- c("length", "x_coord", "y_coord")
  
  cc1 <- vapply(buddPhylo[chr_cols], is.character, logical(1))
  cc2 <- vapply(buddPhylo[num_cols],  is.numeric,  logical(1))
  
  if (!all(cc1)) {
    message("Columns: ", paste(chr_cols[!cc1], collapse = ", "),
            " should be character")
    return(FALSE)
  }
  
  if (!all(cc2)) {
    message("Columns: ", paste(num_cols[!cc2], collapse = ", "),
            " should be numeric")
    return(FALSE)
  }
  
  if (!"extant" %in% colnames(buddPhylo)) {
    message("There is no \"extant\" in \"buddPhylo\"")
    return(FALSE)
  }
  
  if (!is.logical(buddPhylo$extant)) {
    message("Column \"extant\" is not logical in in \"buddPhylo\"")
    return(FALSE)
  }
  
  # ensure names are unique
  if (anyDuplicated(buddPhylo$name)) {
    message("\"name\" values must be unique")
    return(FALSE)
  }
  
  # ensure there is only one root
  if (sum(buddPhylo$orientation == "uca") != 1) {
    message("buddPhylo must have exactly one \"uca\" row")
    return(FALSE)
  }
  
  # endure every parent is also in the name column
  par_idx <- match(buddPhylo$parent, buddPhylo$name)
  if (any(is.na(par_idx) & buddPhylo$orientation != "uca")) {
    message("some \"parent\" values are not in \"name\"")
    return(FALSE)
  }
  
  # ensure we are in pre-order
  if (!all(par_idx < seq_len(nrow(buddPhylo)), na.rm = TRUE)) {
    message("rows must be in pre-order (each parent before its children)")
    return(FALSE)
  }
  
  # if passed in all of this, return true
  return(TRUE)
}

#' @rdname buddPhylo
#'
#' @details \code{print.buddPhylo} The printing of a buddPhylo object is 
#' formatted into a more straightforward and informative sequence manner. We 
#' provide details only for the first few species, since otherwise this print 
#' could be overwhelming for simulations with 10+ species.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#'

print.buddPhylo <- function(x, ...) {
  # change name just for clarity of the object
  buddPhylo <- x
  
  # first check that it is a valid buddPhylo
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  # important variables
  sp <- unique(buddPhylo$taxon)
  n_sp <- sum(!is.na(sp))
  n_tip <- sum(buddPhylo$type == "tip")
  n_sa <- n_sp - n_tip
  
  # then print some basic information about the buddPhylo
  cat(paste0(
    "\nBudding phylogenetic tree with ", 
    n_tip - 1,
    " budding speciation node",
    ifelse((n_tip - 1) == 1, "", "s"),
    " and ", 
    sum(buddPhylo$type %in% c("tip", "sampAnc")),
    " occurrences for ",
    n_sp,
    " species, of which ",
    n_tip,
    ifelse(n_tip == 1, " is a tip", " are tips"),
    " and ",
    n_sa,
    ifelse(n_sa == 1, " is a sampled ancestor.", " are sampled ancestors. ")
  ))
  
  # get the first occurrence of each species
  first_nodes <- buddPhylo[buddPhylo$taxon %in% sp[!is.na(sp)], ]
  first_nodes <- first_nodes[order(first_nodes$x_coord, decreasing = TRUE), ]
  first_nodes <- first_nodes[!duplicated(first_nodes$taxon), ]
  
  # get the parent ranges of these species
  parent_ranges <- buddPhylo$range[buddPhylo$name %in% first_nodes$parent]
  
  # and then some details for first five
  cat("\nSpecies names:\n")
  print(head(first_nodes$taxon))
  
  cat("\nNames of progenitor taxa:\n")
  print(head(parent_ranges))
  
  # to see the whole vector
  cat("\n\nFor more details on vector y, try buddPhylo$y, with y one of:\n")
  cat(names(buddPhylo))
}

#'
# @details \code{head.buddPhylo} Selects only a number of species from the
# #' beginning of a \code{buddPhylo} object.
# #'#' @param subType A character telling which "type" should be used to subset the
# #' output, for instance only tips \code{"tip"} or only sampled ancestors
# #' \code{"sampAnc"}. If set to \code{NULL} (default), it shows all elements.

# @author Matheus Januario.
#'

# head.buddPhylo <- function(x, subType = NULL, ...) {
#   # change name for clarity
#   buddPhylo <- x
#   
#   # first check that it is a valid buddPhylo
#   if (!is.buddPhylo(buddPhylo)) {
#     stop("Invalid buddPhylo object, see ?buddPhylo")
#   }
#   
#   # check subType
#   if (is.null(subType)) {
#     # if null, just do head
#     res <- head(as.data.frame(buddPhylo), ...)
#   } else {
#     # otherwise, check if it is a type
#     if (subType %in% buddPhylo$type) {
#       # if so, get head of just that type
#       res <- head(as.data.frame(buddPhylo[buddPhylo$type == subType, ]), ...)
#     } else {
#       # otherwise, error
#       stop("\"subType\" is not a type in \"buddPhylo\"")
#     }
#   }
#   
#   # make it a buddPhylo
#   class(res) <- c("buddPhylo", "data.frame")
#   
#   # return it
#   return(res)
# }

#'
# @details \code{tail.buddPhylo} Selects only a number of species from the end of a
# \code{buddPhylo} object.
#'
# @author Matheus Januario
#'

# tail.buddPhylo <- function(x, subType = NULL, ...) {
#   # change name for clarity
#   buddPhylo <- x
#   
#   # first check that it is a valid buddPhylo
#   if (!is.buddPhylo(buddPhylo)) {
#     stop("Invalid buddPhylo object, see ?buddPhylo")
#   }
#   
#   # check subType
#   if (is.null(subType)) {
#     # if null, just do tail
#     res <- tail(as.data.frame(buddPhylo), ...)
#   } else {
#     # otherwise, check if it is a type
#     if (subType %in% buddPhylo$type) {
#       # if so, get tail of just that type
#       res <- tail(as.data.frame(buddPhylo[buddPhylo$type == subType, ]), ...)
#     } else {
#       # otherwise, error
#       stop("\"subType\" is not a type in \"buddPhylo\"")
#     }
#   }
#   
#   # make it a buddPhylo
#   class(res) <- c("buddPhylo", "data.frame")
#   
#   # return it
#   return(res)
# }

#' @rdname buddPhylo
#'
#' @details \code{summary.buddPhylo} Quantitative details on the 
#' \code{buddPhylo} object. Prints the number of species, number of extant 
#' species, summary of durations and speciation waiting times.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#'

summary.buddPhylo <- function(object, ...) {
  # change name just for clarity of the object
  buddPhylo <- object
  
  # check that it is a valid buddPhylo object
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  # important variables
  sp <- unique(buddPhylo$taxon)
  n_sp <- sum(!is.na(sp))
  n_tip <- sum(buddPhylo$type == "tip")
  n_sa <- n_sp - n_tip
  
  # some quick numbers
  cat("\n\n  Total number of species: ", n_sp)
  cat("\n  Number of sampled ancestors: ", n_sa)
  
  # get durations (for extinct species)
  dur <- stats::aggregate(buddPhylo$length[!buddPhylo$extant & 
                                             buddPhylo$length > 0],
                          by = list(buddPhylo$taxon[!buddPhylo$extant &
                                                      buddPhylo$length > 0]), 
                          sum, na.rm = TRUE
  )
  
  # summarize durations
  if (nrow(dur) > 0) {
    cat("\n  Durations (only extinct species):")
    cat("\n    mean: ", mean(dur$x))
    cat("\n    standard deviation: ", sd(dur$x))
    cat("\n    summary:\n")
    print(summary(dur$x)[-4])
  } else {
    cat("\n  No extinct species, so no duration information.")
  }
  
  # get speciation nodes
  spec <- buddPhylo[
    sapply(buddPhylo$name, function(x) {
      children <- getChildren.buddPhylo(buddPhylo, x)
      if (length(children) > 0) {
        if (!("sampAnc" %in% names(children))) {
          return(which(buddPhylo$name == children[names(children) == 
                                                    "descendant"]))
        }
      }
      return(FALSE)
    }), 
  ]
  
  # start waiting times vector
  s_wt <- numeric(nrow(spec))
  
  # iterate through nodes to find waiting times
  for (i in seq_len(nrow(spec))) {
    # get this node
    node <- spec[i, ]
    
    # and start at its parent
    cur_node <- node
    
    # start checker for whether we found the previous speciation
    prev_sp <- FALSE
    
    # while we haven't found it
    while (!prev_sp) {
      # get parent of this node
      parent <- buddPhylo[buddPhylo$name == cur_node$parent, ]
      
      # get the children of this parent
      children <- getChildren.buddPhylo(buddPhylo, parent$name)
      
      # set new current node
      cur_node <- parent
      
      # check if it has SA children
      if (!("sampAnc" %in% names(children))) {
        # if not, we reached the next speciation, so we break
        prev_sp <- TRUE
      }
    }
    
    # add waiting time
    s_wt[i] <- cur_node$x_coord - node$x_coord
  }
  
  # summarize time to speciation
  if (length(s_wt) > 0) {
    cat("\n  Speciation waiting times:")
    cat("\n    mean: ", mean(s_wt))
    cat("\n    standard deviation: ", sd(s_wt))
    cat("\n    summary:\n")
    print(summary(s_wt)[-4])
  } else {
    cat("\n  Only one species, so no time to speciation information.")
  }
}

#' Plot a budding phylogenetic tree
#'
#' @details \code{plot.buddPhylo} Plots a \code{buddPhylo} objects as a budding
#' phylogenetic tree. The function was coded to work similarly, thought not 
#' exactly like, \code{ape::plot.phylo()}.
#'
#' @param x Object of class "buddPhylo".
#'
#' @param type a character string specifying the type of phylogeny to be drawn;
#' it must be one of "phylogram" (the default) or "fan".
#'
#' @param cex a numeric value giving the factor scaling of the tip and node
#' labels (character expansion). The default is to take the current value from
#' the graphical parameters.
#'
#' @param nodeCex a numeric value giving the factor scaling of the node
#' labels (character expansion). The default is to take the
#' current value from the graphical parameters.
#'  
#' @param tipCex a numeric value giving the factor scaling of the tip
#' labels (character expansion). The default is to take the
#' current value from the graphical parameters.
#'
#' @param tip.color the colors used for the tip labels, eventually recycled.
#'  
#' @param lineageAsNodeLabels A logical indicating whether the information on 
#' \code{buddPhylo$lineage} must be plotted as node labels, instead of 
#' \code{buddPhylo$name} (default). 
#'
#' @param sampAncJitterMean Bias (average) of the jitter associated to each 
#' sampled ancestor in the tree. The jitter is added to the y-coordinate of the
#' parent lineage as a normally-distributed random deviation.
#'
#' @param sampAncJitterSD Spreading (standard deviation) of the jittering
#' associated to each sampled ancestor in the tree. The jitter is added to the 
#' y-coordinate of the parent lineage as a normally-distributed random 
#' deviation.
#'
#' @param show.tip.label A logical indicating whether to show the tip labels on
#'  the phylogeny (defaults to \code{TRUE}, i.e. the labels are shown).
#'  
#' @param show.SA.label A logical indicating whether to show the sampled 
#' ancestor labels on the phylogeny (defaults to \code{FALSE}, i.e. the labels 
#' are not shown).
#'  
#' @param show.node.label A logical indicating whether to show the node labels
#' on the phylogeny (defaults to \code{FALSE}, i.e. the labels are not shown).
#'
#' @param edge.color A vector of mode character giving the colors used to draw
#' the branches of the plotted phylogeny. These are taken to be in the same
#' order than the component edge of the \code{buddPhylo}. If fewer colors are 
#' given than the number of edges, then the colors are recycled.
#'
#' @param edge.width A numeric vector giving the width of the branches of the
#' plotted phylogeny. These are taken to be in the same order than the
#' component edge of the \code{buddPhylo}. If fewer widths are given than the
#' number of edges, then these are recycled.
#'
#' @param edge.lty Same as the previous argument but for line types; 1: plain,
#' 2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash.
#'
#' @param node.color A vector of mode character giving the colors used to draw
#' the perpendicular lines associated with each node of the plotted phylogeny.
#' These are taken to be in the same order than the component node of the
#' \code{buddPhylo}. If fewer colors are given than the number of nodes, then
#' the colors are recycled.
#'
#' @param node.width Same as the previous argument, but for line widths.
#'
#' @param node.lty Same as the previous argument, but for line types; 1: plain,
#' 2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash.
#'
#' @param font An integer specifying the type of font for the labels: 1 (plain
#' text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#'
#' @param adj A numeric specifying the justification of the text strings of
#' the labels: 0 (left-justification), 0.5 (centering), or 1
#' (right-justification). The value is set with respect of direction.
#'
#' @param srt A numeric giving how much the labels are rotated in degrees
#' (negative values are allowed resulting in clock-like rotation); the value
#' has an effect respectively to the value of direction. This option has no 
#' effect if type = "fan".
#'
#' @param label.offset A numeric giving the space between the nodes and the tips
#'  of the phylogeny and their corresponding labels.
#'
#' @param underscore A logical specifying whether the underscores in tip
#' labels should be written as spaces (the default) or left as are (if 
#' \code{TRUE}).
#'
#' @param align.tip.label A logical value or an integer. If \code{TRUE}, the
#' tips are aligned and dotted lines are drawn between the tips of the tree and 
#' the labels. If an integer, the tips are aligned and this gives the type of 
#' the lines (\code{lty}).
#'
#' @param x.lim A 2-characters long numeric vector of length one or two giving
#' the limit(s) of the x-axis. If \code{NUL}L, this is computed with respect to
#' various parameters such as the string lengths of the labels and the branch 
#' lengths.
#'
#' @param y.lim Same as the previous argument, for the y-axis.
#'
#' @param direction A character string specifying the direction of the tree.
#' Four values are possible: "rightwards" (the default), or "leftwards".
#'
#' @param rotate.tree For "fan" trees: the rotation of the whole tree in d
#' egrees (negative values are accepted).
#'
#' @param open.angle If type = "fan", the angle in degrees left blank. Use a 
#' non-zero value if you want to add a time axis after the tree is plotted.
#'
#' @param add.fan.tip.guide A logical indicating if a dotted line connecting 
#' every tip to its label should be added. Has no effect if 
#' \code{align.tip.label = FALSE} or if \code{type = "phylogram"}.
#'
#' @param ... Further arguments inherited from graphics::plot()
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @export
#'

plot.buddPhylo <- function(x, type = "phylogram",
                           cex = par("cex"), tipCex = par("cex"), 
                           nodeCex = par("cex"), tip.color = par("col"), 
                           lineageAsNodeLabels = FALSE, 
                           sampAncJitterMean = 0, sampAncJitterSD = 0.3,
                           show.tip.label = TRUE, show.SA.label = FALSE, 
                           show.node.label = FALSE,
                           edge.color = "black", edge.width = NULL, 
                           edge.lty = NULL,
                           node.color = "grey30", node.width = NULL, 
                           node.lty = NULL,
                           font = 3, adj = NULL, srt = 0,
                           label.offset = 0, underscore = FALSE, 
                           align.tip.label = FALSE,
                           x.lim = NULL, y.lim = NULL, direction = "rightwards",
                           rotate.tree = 0, open.angle = 0, 
                           add.fan.tip.guide = TRUE, ...) {
  # change name just for clarity of the object
  buddPhylo <- x
  
  # check that it is a valid buddPhylo object
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  # make sure to reset par options after functions
  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar))
  
  # find the names of the sampled ancestor tips
  SAncNames <- unique(buddPhylo$name[buddPhylo$type == "sampAnc"])
  
  # iterate through sampled ancestors
  for (ss in SAncNames) {
    # get row for this node
    SAID <- which(buddPhylo$name == ss)
    
    # and the row of its parent
    SAParID <- which(buddPhylo$name == buddPhylo$parent[SAID])
    
    # set the y_coord as the parent y_coord with a tiny jitter
    buddPhylo$y_coord[SAID] <- buddPhylo$y_coord[SAParID] + rnorm(
      n = 1,
      mean = sampAncJitterMean, sd = sampAncJitterSD
    )
    
    # set x_coord and x_par as before
    buddPhylo$x_coord[SAID] <- buddPhylo$x_coord[SAParID]
    buddPhylo$x_par[SAID] <- buddPhylo$x_coord[SAParID]
  }
  
  # set parental coordinates
  buddPhylo$y_par <- buddPhylo$y_coord[match(buddPhylo$parent, buddPhylo$name)]
  buddPhylo$y_par[buddPhylo$orientation == "uca"] <- 
    buddPhylo$y_coord[buddPhylo$orientation == "uca"]
  buddPhylo$x_par <- buddPhylo$x_coord[match(buddPhylo$parent, buddPhylo$name)]

  # set colors to default values if not given
  if (is.null(tip.color)) {
    tip.color <- rep("black", times = sum(buddPhylo$type %in% "tip"))
  }
  if (is.null(edge.color)) {
    edge.color <- rep("black", times = nrow(buddPhylo))
  }
  if (is.null(node.color)) {
    node.color <- rep("black", times = nrow(buddPhylo))
  }
  if (is.null(edge.width)) {
    edge.width <- rep(1, times = nrow(buddPhylo))
  }
  if (is.null(node.width)) {
    node.width <- rep(0.5, times = nrow(buddPhylo))
  }
  if (is.null(edge.lty)) {
    edge.lty <- rep(1, times = nrow(buddPhylo))
  }
  if (is.null(node.lty)) {
    node.lty <- rep(1, times = nrow(buddPhylo))
  }
  
  # set x limits
  if (is.null(x.lim)) {
    xrange <- range(c(buddPhylo$x_coord, buddPhylo$x_par), na.rm = TRUE)
  } else {
    xrange <- x.lim
  }
  
  # set y limits
  if (is.null(y.lim)) {
    yrange <- range(buddPhylo$y_coord, na.rm = TRUE)
  } else {
    yrange <- y.lim
  }
  
  # check which labels to show
  if (show.SA.label) {
    row_tip_labels <- which(buddPhylo$type %in% c("tip", "sampAnc"))
  } else {
    row_tip_labels <- which(buddPhylo$type == "tip")
  }
  
  # check if we want to keep the underscore in tip names
  if (!underscore) {
    labels <- sub("_", " ", buddPhylo$taxon[row_tip_labels])
  } else {
    labels <- buddPhylo$taxon[row_tip_labels]
  }
  
  # check type
  if (type == "phylogram") {
    # if phylogram, go by direction
    if (direction == "leftwards") {
      # reverse yrange if we want leftwards
      yrange <- rev(yrange)
      
      # reverse label offset
      label.offset <- label.offset * -1
      
      # and set adjustment
      if (is.null(adj)) {
        adjust <- 1
      } else {
        adjust <- adj
      }
    } else {
      # if not leftwards, just need to set adjustment
      if (is.null(adj)) {
        adjust <- 0
      } else {
        adjust <- adj
      }
    }
    
    # now decide if time is from the present (default)
    timeFromPresent <- min(buddPhylo$x_coord[buddPhylo$type == "tip"]) <
      min(buddPhylo$x_coord[buddPhylo$type == "node"])
    if (timeFromPresent) {
      # if so, need to reverse xrange
      xrange <- rev(xrange)
      
      # start plot window
      plot(NA,
           xlim = xrange + c(0.5, diff(xrange) * .1), ylim = yrange, 
           frame.plot = FALSE, xaxt = "n", yaxt = "n", xaxs = "i",
           ylab = "", xlab = ""
      )
      
      # set ticks
      ticks <- pretty(xrange)
      
      # adjust ticks if needed
      if (any(abs(ticks - xrange[1]) < 3)) ticks <- 
        ticks[-which(abs(ticks - xrange[1]) < 3)]
      
      # get first tick
      ticks[1] <- xrange[1]
      
      # sort ticks
      ticks <- sort(unique(c(ticks, 0)), decreasing = TRUE)
      
      # create x-axis
      graphics::axis(1,
                     at = ticks,
                     labels = round(ticks, 2)
      )
      graphics::mtext("Time (Mya)", side = 1, line = 2.5)
    } else {
      # start plot window
      plot(NA,
           xlim = xrange * c(1, 1.1), ylim = yrange, 
           frame.plot = FALSE, xaxt = "n", yaxt = "n", xaxs = "i",
           ylab = "", xlab = ""
      )
      
      # set ticks
      ticks <- pretty(xrange)
      
      # end ticks at the end of the range
      ticks[length(ticks)] <- xrange[2]
      
      # create x-axis
      graphics::axis(1,
                     at = ticks, xrange[2],
                     labels = round(ticks, 2)
      )
      graphics::mtext("Time elapsed since root (My)", side = 1, line = 2.5)
    }
    
    # create segments for nodes
    segments(
      col = edge.color, lwd = edge.width, lty = edge.lty,
      x0 = buddPhylo$x_par, x1 = buddPhylo$x_par,
      y0 = buddPhylo$y_par, y1 = buddPhylo$y_coord
    )
    
    # and for edges
    segments(
      col = edge.color, lwd = edge.width, lty = edge.lty,
      x0 = buddPhylo$x_par, x1 = buddPhylo$x_coord,
      y0 = buddPhylo$y_coord, y1 = buddPhylo$y_coord
    )
    
    # check if we want to show tip labels
    if (show.tip.label) {
      # get the tip x-coordinates
      xtip <- buddPhylo$x_coord[row_tip_labels]
      
      # # check if we want to align
      # if (align.tip.label) {
      #   if (timeFromPresent) {
      #     xtip <- rep(min(xtip, na.rm = TRUE))
      #   } else {
      #     xtip <- rep(max(xtip, na.rm = TRUE))
      #   }
      # }
      
      # write labels
      text(
        x = xtip,
        y = buddPhylo$y_coord[row_tip_labels] - 0.01,
        labels = labels, adj = adjust, srt = srt,
        cex = tipCex, font = font
      )
    }
    
    # check if we want to show node label
    if (show.node.label) {
      # check if we want lineages as the node labels
      if (lineageAsNodeLabels) {
        # write labels
        text(
          x = buddPhylo$x_coord[buddPhylo$type != "tip"], col = node.color,
          y = buddPhylo$y_coord[buddPhylo$type != "tip"],
          labels = buddPhylo$lineage[buddPhylo$type != "tip"],
          cex = nodeCex, font = font, offset = label.offset
        )
      } else {
        # write labels
        text(
          x = buddPhylo$x_coord[buddPhylo$type != "tip"], col = node.color,
          y = buddPhylo$y_coord[buddPhylo$type != "tip"],
          labels = buddPhylo$name[buddPhylo$type != "tip"],
          cex = nodeCex, font = font, offset = label.offset
        )
      }
    }
  } else if (type == "fan") {
    # map ycoord in angle + ope angle
    theta_map <- function(y) {
      (y - min(buddPhylo$y_coord, na.rm = TRUE)) /
        (max(buddPhylo$y_coord, na.rm = TRUE) - 
           min(buddPhylo$y_coord, na.rm = TRUE)) *
        (2 * pi - open.angle * pi / 180)
    }
    
    # angles + rotation
    th <- theta_map(buddPhylo$y_coord) + rotate.tree * pi / 180
    
    # radius
    r0 <- buddPhylo$x_coord - min(buddPhylo$x_coord, na.rm = T)
    r1 <- (buddPhylo$x_coord - buddPhylo$length) - 
      min(buddPhylo$x_coord, na.rm = TRUE)
    
    # convert to Cartesian
    x0 <- r0 * cos(th)
    y0 <- r0 * sin(th)
    x1 <- r1 * cos(th)
    y1 <- r1 * sin(th)
    
    # set x limits
    if (is.null(x.lim)) {
      xrange <- range(c(x0, x1), na.rm = TRUE)
    } else {
      xrange <- x.lim
    }
    
    # set y limites
    if (is.null(y.lim)) {
      yrange <- range(c(y0, y1), na.rm = TRUE)
    } else {
      yrange <- y.lim
    }
    
    # stretch plot size if plotting tips
    if (show.tip.label) {
      yrange <- yrange * c(2, 1.5)
    }
    
    # draw empty plot
    plot(NA,
         xlim = xrange, ylim = yrange,
         frame.plot = F, yaxt = "n", xaxt = "n",
         asp = 1, xlab = "", ylab = ""
    )
    
    # iterate through angles
    for (i in 1:length(th)) {
      # if there is no parent, skip
      if (is.na(buddPhylo$y_par[i])) {
        next
      }
      
      # make a vector of angles
      th_seq <- seq(buddPhylo$y_par[i], buddPhylo$y_coord[i], length.out = 100)
      th_seq <- theta_map(th_seq)
      
      # r parameters
      r_par <- buddPhylo$x_coord[i] - buddPhylo$length[i] - min(
        buddPhylo$x_coord,
        na.rm = TRUE
      )
      
      # draw the angular segments
      graphics::lines(
        x = r_par * cos(th_seq), y = r_par * sin(th_seq),
        col = node.color, lwd = node.width, lty = node.lty
      )
    }
    
    # draw edges
    segments(
      col = edge.color, lwd = edge.width, lty = edge.lty,
      x0 = x0, x1 = x1, y0 = y0, y1 = y1
    )
    
    # find positions of labels
    srts <- (srt + (th * 180 / pi))
    adjs <- rep(0, times = length(srts))
    adjs[srts > 90 & srts < 270] <- 1
    srts[srts > 90 & srts < 270] <- srts[srts > 90 & srts < 270] + 180
    
    # check if we want node labels
    if (show.node.label) {
      # find the nodes that are not tips
      ids2plot <- which(buddPhylo$type != "tip")
      
      # iterate through them
      for (i in ids2plot) {
        # check if we want lineages as the labels
        if (lineageAsNodeLabels) {
          # write labels
          text(
            x = x0[i], y = y0[i], srt = srts[i], col = node.color,
            labels = buddPhylo$lineage[i], adj = adjs[i],
            cex = nodeCex, font = font, offset = label.offset
          )
        } else {
          # write labels
          text(
            x = x0[i], y = y0[i], srt = srts[i], col = node.color,
            labels = buddPhylo$name[i], adj = adjs[i],
            cex = nodeCex, font = font, offset = label.offset
          )
        }
      }
    }
    
    # check if we want tip labels
    if (show.tip.label) {
      # get angular position of tip labels
      x0 <- (r0) * cos(th)
      y0 <- (r0) * sin(th)
      
      # get tip labels
      ids2plot <- row_tip_labels
      
      # check if we want to align them
      if (align.tip.label) {
        # backup x0 and y0 
        rx0 <- x0
        ry0 <- y0
        
        # adjust them based on cosine
        x0 <- (max(r0, na.rm = TRUE)) * cos(th)
        y0 <- (max(r0, na.rm = TRUE)) * sin(th)
        
        # check if we want the fan tip guide
        if (add.fan.tip.guide) {
          # draw guides
          segments(
            x0 = x0[ids2plot], x1 = rx0[ids2plot],
            y0 = y0[ids2plot], y1 = ry0[ids2plot],
            lty = 3, col = "grey80", lwd = .5
          )
        }
      }
      
      # iterate through tip labels
      for (i in 1:length(ids2plot)) {
        # write labels
        text(
          x = x0[ids2plot[i]], y = y0[ids2plot[i]], srt = srts[ids2plot[i]],
          offset = label.offset, labels = labels[i],
          cex = tipCex, font = font, adj = adjs[ids2plot[i]]
        )
      }
    }
  } else {
    stop("wrong \"type\"")
  }
}


