#' Details, generics, and methods for the \code{buddPhylo} class
#'
#' @description The \code{buddPhylo} class is a format that stores the
#' information from a phylogeny that contains budding speciation events. Every
#' \code{buddPhylo} is also a \code{data.frame}, and every line refers to a
#' node in the phylogeny. Every node (i.e., "line") may contain a myriad of
#' types of information ("columns"), but the following elements organize the
#' core information regarding every node within a \code{buddPhylo} object:
#'
#' \describe{
#' \item{\code{name}}{Name of the node. Must be a character column used to
#' trace evolutionary relationships among edges. It can be any character
#' (including numbers), as long as it is consistently named across each
#' relevant branches/edge connection}.
#'
#' \item{\code{length}}{Numeric column containing the length of the edge that
#' lead to the node. All functions assume time is measured in millions of years}
#'
#' \item{\code{orientation}}{}
#'
#' \item{\code{taxon}}{Name that should be associated with the node (e.g., in
#' plotting functions)}
#'
#' \item{\code{parent}}{ \code{name} of the parental edge of the node in
#' question. Note this refer to \code{name}, not to \code{taxon}. The edge that
#' lead to the \code{uca} has its parent set as \code{NA}.}
#'
#' \item{\code{type}}{Character assigiong whether the node is a \code{node},
#' a \code{tip}, or a sampled ancestor (\code{sampAnc}).}
#'
#' \item{\code{y_coord}}{Y coordinate to be used as a reference for that node
#' for the plotting functions}.
#'
#' \item{\code{x_coord}}{X coordinate to be used as a reference for that node
#' for the plotting functions}.
#'
#' \item{\code{extant}}{A logical indicating whether the tip is alive at the
#' present \code{extant=TRUE} or not}.
#'
#' Here we declare useful generics and methods for \code{buddPhylo} objects.
#'
#' @param buddPhylo,x,object Object of class "buddPhylo"
#'
#' @param subType A character telling which "type" should be used to subset the
#' output, for instance only tips \code{subType="tip"} or only sampled ancestors
#' \code{subType="sampAnc"}. If set to \code(NULL) (default), it shows all elements.
#'
#' @param focalLineage A node/lineage to be used as reference for the tree operation.
#'
#' @param nameList A character vector with a list of \code{name} (i.e., not
#' "taxon") to use as a reference for the node set operation.
#'
#' @param excludeSampAnc Logical indicating whether sampled ancestors should be
#' exlcuded from the output.
#'
#' @param onlyImmediates Logical indicating if only the most immediate
#' descendant (the two daughter) lineages of a given \code{focalLineage} should
#' be considered.
#'
#' @param ... Further arguments inherited from generics.
#' 
#' @author Matheus Januario
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
#' members: \code{"name}, \code{length}, \code{orientation}, \code{taxon},
#' \code{parent}, \code{type}, \code{y_coord}, \code{x_coord"}. And all of
#' these must have the correct length (i.e. \code{= (2*Ntip)-1)}.
#'
#' @author Matheus Januario
#'
#' @export
#'

is.buddPhylo <- function(x) {
  buddPhylo <- x
  
  # checks that the class is buddPhylo
  cla <- "buddPhylo" %in% class(buddPhylo)
  
  # if it is not, should be
  if (!cla) {
    return(FALSE)
  }
  
  # checks that there is enough minimal information
  mustHaveCols <- c("name", "length", "orientation", "taxon", "parent", "type", "y_coord", "x_coord")
  cols <- all(mustHaveCols %in% colnames(buddPhylo))
  
  # if there are not, the other tests are redundant
  if (!cols) {
    message(paste0(
      "\"buddPhylo\" is lacking the columns: ",
      paste0(mustHaveCols[!mustHaveCols %in% colnames(buddPhylo)], collapse = ", ")
    ))
    return(FALSE)
  }
  
  # Checks the classes of each colmun:
  cmatch <- match(mustHaveCols, colnames(buddPhylo))
  cc1 <- (apply(buddPhylo[, cmatch[1:6]], 2, class) == "character")
  cc2 <- (apply(buddPhylo[, cmatch[7:8]], 2, class) == "numeric")
  
  if (!all(cc1)) {
    message(paste0(
      "Columns: ",
      paste0(mustHaveCols[!cc1], collapse = ", "),
      " have the wrong class"
    ))
    return(FALSE)
  }
  
  if (!all(cc2)) {
    message(paste0(
      "Columns: ",
      paste0(mustHaveCols[!cc2], collapse = ", "),
      " have the wrong class"
    ))
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
  
  # if passed in all of this:
  return(TRUE)
}

#' @rdname buddPhylo
#'
#' @details \code{print.buddPhylo} The printing of a buddPhylo object is formatted into a more
#' straightforward and informative sequence manner. We provide details only for
#' the first few species, since otherwise this print could be overwhelming for
#' simulations with 10+ species.
#'
#' @author Matheus Januario
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
  
  # then print some basic information about the buddPhylo:
  cat(paste0(
    "\nPhylogenetic tree with ", sum(buddPhylo$orientation == "descendant"),
    " budding bifurcations and ", sum(buddPhylo$type == "tip"),
    " species"
  ))
  
  # and then some details for first five
  cat("\nSpecies names:\n")
  print(head(buddPhylo$taxon[buddPhylo$type == "tip"]))
  
  cat("\nNames of Progenitor taxa:\n")
  print(
    buddPhylo$taxon[na.omit(match(buddPhylo$name, head(buddPhylo$parent[buddPhylo$type == "tip"])))]
  )
  
  cat("\nNode names:\n")
  print(head(buddPhylo$name[buddPhylo$type == "node"]))
  
  # to see the whole vector
  cat("\n\nFor more details on vector y, try buddPhylo$y, with y one of:\n")
  cat(names(buddPhylo))
}

#' @rdname buddPhylo
#'
#' @details \code{head.buddPhylo} Selects only a number of species from the
#' beginning of a \code{buddPhylo} object.
#'
#' @author Matheus Januario
#'
#' @export

head.buddPhylo <- function(x, subType = NULL, ...) {
  # change name for clarity
  buddPhylo <- x
  
  # first check that it is a valid buddPhylo
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  if (is.null(subType)) {
    res <- head(buddPhylo, ...)
  } else {
    if (subType %in% buddPhylo$type) {
      res <- head(buddPhylo[buddPhylo$type == subType, ], ...)
    } else {
      stop("\"subType\" is not a type in \"buddPhylo\"")
    }
  }
  
  # make it a buddPhylo
  class(res) <- c("buddPhylo", "data.frame")
  
  # return it
  return(res)
}

#' @rdname buddPhylo
#' @details \code{tail.buddPhylo} Selects only a number of species from the end of a
#' \code{buddPhylo} object.
#'
#' @author Matheus Januario
#'
#' @export
#'

tail.buddPhylo <- function(x, subType = NULL, ...) {
  # change name for clarity
  buddPhylo <- x
  
  # first check that it is a valid buddPhylo
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  if (is.null(subType)) {
    res <- tail(buddPhylo, ...)
  } else {
    if (subType %in% buddPhylo$type) {
      res <- tail(buddPhylo[buddPhylo$type == subType, ], ...)
    } else {
      stop("\"subType\" is not a type in \"buddPhylo\"")
    }
  }
  
  # make it a buddPhylo
  class(res) <- c("buddPhylo", "data.frame")
  
  # return it
  return(res)
}

#' @rdname buddPhylo
#'
#' @details \code{summary.buddPhylo} Quantitative details on the \code{buddPhylo} object.
#' Prints the number of species, number of extant species, summary of durations
#' and speciation waiting times, in case there are more than one species.
#'
#' @author Matheus Januario
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
  
  # some quick numbers
  cat("\n\n  Total number of species: ", sum(buddPhylo$type == "tip"))
  cat("\n  Number of sampled ancestors: ", sum(buddPhylo$type == "sampAnc"))
  
  # get durations (for extinct species)
  dur <- stats::aggregate(buddPhylo$length[!buddPhylo$extant],
                          by = list(buddPhylo$taxon[!buddPhylo$extant]), sum, na.rm = TRUE
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
  
  # get speciation waiting times
  spec <- diff(sort(buddPhylo$x_coord[buddPhylo$orientation == "descendant"] -
                      buddPhylo$length[buddPhylo$orientation == "descendant"]))
  
  # summarize time to speciation
  if (length(spec) > 0) {
    cat("\n  Speciation waiting times:")
    cat("\n    mean: ", mean(spec))
    cat("\n    standard deviation: ", sd(spec))
    cat("\n    summary:\n")
    print(summary(spec)[-4])
  } else {
    cat("\n  Only one species, so no time to speciation information.")
  }
}

#' Visualize a buddPhylo object
#'
#' @details \code{plot.buddPhylo} Plots a buddPhylo objects as phylogenetic
#' trees. The function was coded to work similarly, thought not exactly like,
#' \code{ape::plot.phylo()}.
#'
#' @inheritParams buddPhylo
#'
#' @param x Object of class "buddPhylo"
#'
#' @param type a character string specifying the type of phylogeny to be drawn;
#' it must be one of "phylogram" (the default) or "fan".
#'
#' @param cex a numeric value giving the factor scaling of the tip and node
#'  labels (Character EXpansion). The default is to take the current value from
#'  the graphical parameters.
#'
#' @param nodeCex a numeric value giving the factor scaling of the tip
#'  labels (Character EXpansion). The default is to take the
#'  current value from the graphical parameters.
#'  
#' @param tipCex a numeric value giving the factor scaling of the tip
#'  labels (Character EXpansion). The default is to take the
#'  current value from the graphical parameters.
#'
#' @param tip.color the colours used for the tip labels,
#'  eventually recycled.
#'  
#' @param lineageAsNodeLabels A logical indicating whether the inferomation on 
#'  \code{buddPhylo$lineage} must be plotted as node labels, instead of 
#'  \code{buddPhylo$name} (default). 
#'
#' @param sampAncJitterMean Bias (average) of the jitter
#'  associated to each sampled ancestor in the tree. The jitter is added with a
#'  normally-distributed random deviation from the parental lineage.
#'
#' @param sampAncJitterSD Spreading (standard deviation) of the jittering
#'  associated to each sampled ancestor in the tree. The jitter is added with a
#'  normally-distributed random deviation from the parental lineage and is
#'  "turned off" when this parameter is set to zero.
#'
#' @param show.tip.label a logical indicating whether to show the tip labels on
#'  the phylogeny (defaults to TRUE, i.e. the labels are shown).
#'
#' @param show.node.label a logical indicating whether to show the node labels
#'  on the phylogeny (defaults to FALSE, i.e. the labels are not shown).
#'
#' @param edge.color a vector of mode character giving the colours used to draw
#'  the branches of the plotted phylogeny. These are taken to be in the same
#'  order than the component edge of phy. If fewer colors are given than the
#'  length of edge, then the colours are recycled.
#'
#' @param edge.width a numeric vector giving the width of the branches of the
#'  plotted phylogeny. These are taken to be in the same order than the
#'  component edge of phy. If fewer widths are given than the length of edge,
#'  then these are recycled.
#'
#' @param edge.lty same as the previous argument but for line types; 1: plain,
#'  2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash.
#'
#' @param node.color a vector of mode character giving the colors used to draw
#'  the perpendicular lines associated with each node of the plotted phylogeny.
#'  These are taken to be in the same order than the component node of phy. If
#'  fewer colours are given than the length of node, then the colors are recycled.
#'
#' @param node.width as the previous argument, but for line widths.
#'
#' @param node.lty  as the previous argument, but for line types; 1: plain,
#'  2: dashed, 3: dotted, 4: dotdash, 5: longdash, 6: twodash.
#'
#' @param font an integer specifying the type of font for the labels: 1 (plain
#'  text), 2 (bold), 3 (italic, the default), or 4 (bold italic).
#'
#' @param adj 	a numeric specifying the justification of the text strings of
#'  the labels: 0 (left-justification), 0.5 (centering), or 1
#'  (right-justification). The value is set with respect of direction.
#'
#' @param srt a numeric giving how much the labels are rotated in degrees
#'  (negative values are allowed resulting in clock-like rotation); the value
#'  has an effect respectively to the value of direction (see Examples). This
#'  option has no effect if type = "fan".
#'
#' @param label.offset a numeric giving the space between the nodes and the tips
#'  of the phylogeny and their corresponding labels. This option has no effect if
#'  type = "unrooted".	a numeric giving the space between the nodes and the tips
#'  of the phylogeny and their corresponding labels. This option has no effect
#'  if type = "unrooted".
#'
#' @param underscore a logical specifying whether the underscores in tip
#'  labels should be written as spaces (the default) or left as are (if TRUE).
#'
#' @param align.tip.label a logical value or an integer. If TRUE, the tips are
#'  aligned and dotted lines are drawn between the tips of the tree and the
#'  labels. If an integer, the tips are aligned and this gives the type of the
#'  lines (lty).
#'
#' @param x.lim a 2-characters long numeric vector of length one or two giving
#'  the limit(s) of the x-axis. If NULL, this is computed with respect to
#'  various parameters such as the string lengths of the labels and the
#'  branch lengths.
#'
#' @param y.lim same than above for the y-axis.
#'
#' @param direction a character string specifying the direction of the tree.
#'  Four values are possible: "rightwards" (the default), "leftwards",
#'  "upwards", and "downwards".
#'
#' @param rotate.tree for "fan" trees: the rotation
#'  of the whole tree in degrees (negative values are accepted).
#'
#' @param open.angle if type = "fan" , the angle in degrees left blank.
#'  Use a non-zero value if you want to add a time axis after the tree is plotted.
#'
#' @param add.fan.tip.guide a logical indicating if a dotted line connecting every tip to its label show be added. Has no effect if \code{align.tip.label = FALSE}
#'  or if \code{type = "phylogram"}.
#'
#' @param ... Further arguments inherited from graphics::plot()
#'
#' @author Matheus Januario
#'
#' @export
#'

plot.buddPhylo <- function(x, type = "phylogram",
                           cex = par("cex"), tipCex = par("cex"), nodeCex = par("cex"),
                           tip.color = par("col"), lineageAsNodeLabels = F,
                           sampAncJitterMean = 0, sampAncJitterSD = 0.3,
                           show.tip.label = TRUE, show.node.label = FALSE,
                           edge.color = "black", edge.width = NULL, edge.lty = NULL,
                           node.color = "grey30", node.width = NULL, node.lty = NULL,
                           font = 3, adj = NULL, srt = 0,
                           label.offset = 0, underscore = FALSE, align.tip.label = FALSE,
                           x.lim = NULL, y.lim = NULL, direction = "rightwards",
                           rotate.tree = 0, open.angle = 0, add.fan.tip.guide = TRUE, ...) {
  # change name just for clarity of the object
  buddPhylo <- x
  
  # check that it is a valid buddPhylo object
  if (!is.buddPhylo(buddPhylo)) {
    stop("Invalid buddPhylo object, see ?buddPhylo")
  }
  
  # make sure to reset par options after functions
  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar))
  
  # Jitter sampling ancestor coordinates:
  SAncNames <- unique(buddPhylo$name[buddPhylo$type == "sampAnc"])
  ss <- SAncNames[1]
  for (ss in SAncNames) {
    SAID <- which(buddPhylo$name == ss)
    SAParID <- which(buddPhylo$name == buddPhylo$parent[SAID])
    
    buddPhylo$y_coord[SAID] <- buddPhylo$y_coord[SAParID] + rnorm(
      n = 1,
      mean = sampAncJitterMean, sd = sampAncJitterSD
    )
    buddPhylo$x_coord[SAID] <- buddPhylo$x_coord[SAParID]
    buddPhylo$x_par[SAID] <- buddPhylo$x_coord[SAParID]
  }
  
  # Get parental coordinates:
  buddPhylo$y_par <- buddPhylo$y_coord[match(buddPhylo$parent, buddPhylo$name)]
  buddPhylo$y_par[buddPhylo$orientation == "uca"] <- buddPhylo$y_coord[buddPhylo$orientation == "uca"]
  buddPhylo$x_par <- buddPhylo$x_coord[match(buddPhylo$parent, buddPhylo$name)]
  
  
  # Edge properties
  # If there is no color for every branch, assign all as black:
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
  
  # Plotting:
  if (is.null(x.lim)) {
    xrange <- range(c(0, buddPhylo$x_coord, buddPhylo$x_par), na.rm = T)
  } else {
    xrange <- x.lim
  }
  
  if (is.null(y.lim)) {
    yrange <- range(buddPhylo$y_coord, na.rm = T)
  } else {
    yrange <- y.lim
  }
  
  if (!underscore) {
    labels <- sub("_", " ", buddPhylo$taxon[buddPhylo$type == "tip"], )
  } else {
    labels <- buddPhylo$taxon[buddPhylo$type == "tip"]
  }
  
  if (type == "phylogram") {
    if (direction == "leftwards") {
      yrange <- rev(yrange)
      label.offset <- label.offset * -1
      if (is.null(adj)) {
        adjust <- 1
      } else {
        adjust <- adj
      }
    } else {
      if (is.null(adj)) {
        adjust <- 0
      } else {
        adjust <- adj
      }
    }
    
    # now decide if time is form root
    timeFromPresent <- (diff(c(0, max(buddPhylo$x_coord[buddPhylo$extant]))) < 10E-10)
    if (timeFromPresent) {
      xrange <- rev(xrange)
    }
    
    if (timeFromPresent) {
      plot(NA,
           xlim = xrange + c(0, diff(xrange) * .1), ylim = yrange, frame.plot = F,
           xaxt = "n", yaxt = "n", xaxs = "i",
           ylab = "", xlab = ""
      )
    } else {
      plot(NA,
           xlim = xrange * c(1, 1.1), ylim = yrange, frame.plot = F,
           xaxt = "n", yaxt = "n", xaxs = "i",
           ylab = "", xlab = ""
      )
    }
    
    
    if (timeFromPresent) {
      ticks <- pretty(xrange)
      ticks[1] <- xrange[1]
      ticks <- unique(c(ticks, 0))
      
      graphics::axis(1,
           at = ticks, xrange[1],
           labels = round(ticks, 2)
      )
      graphics::mtext("Time (Mya)", side = 1, line = 2.5)
    } else {
      ticks <- pretty(xrange)
      ticks[length(ticks)] <- xrange[2]
      
      graphics::axis(1,
           at = ticks, xrange[2],
           labels = round(ticks, 2)
      )
      graphics::mtext("Time elapsed since root (My)", side = 1, line = 2.5)
    }
    
    # nodes (parenthood connectors):
    segments(
      col = edge.color, lwd = edge.width, lty = edge.lty,
      x0 = buddPhylo$x_par, x1 = buddPhylo$x_par,
      y0 = buddPhylo$y_par, y1 = buddPhylo$y_coord
    )
    # edges (branches):
    
    segments(
      col = edge.color, lwd = edge.width, lty = edge.lty,
      x0 = buddPhylo$x_par, x1 = buddPhylo$x_coord,
      y0 = buddPhylo$y_coord, y1 = buddPhylo$y_coord
    )
    
    # Names:
    if (show.tip.label) {
      xtip <- buddPhylo$x_coord[buddPhylo$type == "tip"]
      
      if (align.tip.label) {
        if (timeFromPresent) {
          xtip <- rep(min(xtip, na.rm = T))
        } else {
          xtip <- rep(max(xtip, na.rm = T))
        }
      }
      
      text(
        x = xtip,
        y = buddPhylo$y_coord[buddPhylo$type == "tip"] - .05,
        labels = labels, adj = adjust, srt = srt,
        cex = tipCex, font = font
      )
    }
    if (show.node.label) {
      if (lineageAsNodeLabels) {
        text(
          x = buddPhylo$x_coord[buddPhylo$type != "tip"], col = node.color,
          y = buddPhylo$y_coord[buddPhylo$type != "tip"],
          labels = buddPhylo$lineage[buddPhylo$type != "tip"],
          cex = nodeCex, font = font, offset = label.offset
        )
      } else {
        text(
          x = buddPhylo$x_coord[buddPhylo$type != "tip"], col = node.color,
          y = buddPhylo$y_coord[buddPhylo$type != "tip"],
          labels = buddPhylo$name[buddPhylo$type != "tip"],
          cex = nodeCex, font = font, offset = label.offset
        )
      }
    }
  } else if (type == "fan") { # if phylo should be plotted as a fan:
    
    # map ycoord in angle + ope angle
    theta_map <- function(y) {
      (y - min(buddPhylo$y_coord, na.rm = T)) /
        (max(buddPhylo$y_coord, na.rm = T) - min(buddPhylo$y_coord, na.rm = T)) *
        (2 * pi - open.angle * pi / 180)
    }
    
    # angles + rotation
    th <- theta_map(buddPhylo$y_coord) + rotate.tree * pi / 180
    
    # radius
    r0 <- buddPhylo$x_coord - min(buddPhylo$x_coord, na.rm = T)
    r1 <- (buddPhylo$x_coord - buddPhylo$length) - min(buddPhylo$x_coord, na.rm = T)
    
    # convert to Cartesian
    x0 <- r0 * cos(th)
    y0 <- r0 * sin(th)
    x1 <- r1 * cos(th)
    y1 <- r1 * sin(th)
    
    # set axis lims:
    if (is.null(x.lim)) {
      xrange <- range(c(x0, x1), na.rm = T)
    } else {
      xrange <- x.lim
    }
    
    if (is.null(y.lim)) {
      yrange <- range(c(y0, y1), na.rm = T)
    } else {
      yrange <- y.lim
    }
    
    if (show.tip.label) { # stretch plot size if plotting tips:
      yrange <- yrange * c(2, 1.5)
    }
    
    # Empty plot
    plot(NA,
         xlim = xrange, ylim = yrange,
         frame.plot = F, yaxt = "n", xaxt = "n",
         asp = 1, xlab = "", ylab = ""
    )
    
    # Draw arcs for the node connections:
    for (i in 1:length(th)) {
      if (is.na(buddPhylo$y_par[i])) {
        next
      }
      
      th_seq <- seq(buddPhylo$y_par[i], buddPhylo$y_coord[i], length.out = 100)
      th_seq <- theta_map(th_seq)
      
      r_par <- buddPhylo$x_coord[i] - buddPhylo$length[i] - min(
        buddPhylo$x_coord,
        na.rm = T
      )
      
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
    
    # Now all labels:
    srts <- (srt + (th * 180 / pi))
    adjs <- rep(0, times = length(srts))
    adjs[srts > 90 & srts < 270] <- 1
    srts[srts > 90 & srts < 270] <- srts[srts > 90 & srts < 270] + 180
    
    if (show.node.label) {
      ids2plot <- which(buddPhylo$type != "tip")
      for (i in ids2plot) {
        if (lineageAsNodeLabels) {
          text(
            x = x0[i], y = y0[i], srt = srts[i], col = node.color,
            labels = buddPhylo$lineage[i], adj = adjs[i],
            cex = nodeCex, font = font, offset = label.offset
          )
        } else {
          text(
            x = x0[i], y = y0[i], srt = srts[i], col = node.color,
            labels = buddPhylo$name[i], adj = adjs[i],
            cex = nodeCex, font = font, offset = label.offset
          )
        }
      }
    }
    
    # Add tip labels:
    if (show.tip.label) {
      x0 <- (r0) * cos(th)
      y0 <- (r0) * sin(th)
      
      ids2plot <- which(buddPhylo$type == "tip")
      
      if (align.tip.label) {
        rx0 <- x0
        ry0 <- y0
        
        x0 <- (max(r0, na.rm = T)) * cos(th)
        y0 <- (max(r0, na.rm = T)) * sin(th)
        
        if (add.fan.tip.guide) {
          segments(
            x0 = x0[ids2plot], x1 = rx0[ids2plot],
            y0 = y0[ids2plot], y1 = ry0[ids2plot],
            lty = 3, col = "grey80", lwd = .5
          )
        }
      }
      
      # plot tip labels:
      for (i in 1:length(ids2plot)) {
        text(
          x = x0[ids2plot[i]], y = y0[ids2plot[i]], srt = srts[ids2plot[i]],
          offset = label.offset, labels = labels[i],
          cex = tipCex, font = font, adj = adjs[ids2plot[i]]
        )
      }
    }
  } else {
    stop("wrong \"type\"")
  } # end of fan phylo
}


