#' Subset and parse a nexus file
#'
#' This function reads and parses one or several trees in a NEXUS file.
#'
#' @param file 	A character containing a nexus file name.
#'
#' @param treeNumbers A numeric vector containing which trees should be read.
#' If \code{NULL} (default), the function reads all trees within \code{file}.
#'
#' @return A character with the parenthetic format within the nexus file.
#' If more than one tree is read, a list of characters is returned in the
#' order given by \code{treeNumbers}
#'
#' @author Matheus Januario
#'
#' @details The does not take note on  whichever name was assigned for any tree
#' in the nexus file.
#'
#' @export
#'

parse.nexus <- function(text = NULL, file = NULL, treeNumbers = NULL) {
  if (is.null(text)) {
    # Read file
    lines <- readLines(file)
  } else {
    lines <- text
  }
  
  # Find the lines with the tree-indicating pattern
  target_line <- grep("tree ", lines, value = TRUE)
  
  # restrict to the treeNumbers lines if it is not NULL
  if (!is.null(treeNumbers)) target_line <- target_line[treeNumbers]
  trees <- sub('^\\s*tree\\s+\\S+\\s*=\\s*"?', "", target_line, perl = TRUE)
  trees <- sub('"\\s*$', "", trees, perl = TRUE)
  
  return(trees)
}

#' Extract the "translation" block in a nexus file
#'
#' This function only the tip translations (numbers to characters) in NEXUS file.
#'
#' @param file 	A character containing a nexus file name.
#'
#' @return A data frame containing the number \code{id} and the \code{taxon}
#' associated with it.
#'
#' @author Matheus Januario
#'
#' @export
#'

extract.translate.block.nexus <- function(text = NULL, file = NULL) {
  if (is.null(text)) {
    # read file
    lines <- readLines(file)
  } else {
    lines <- text
  }
  
  # find start and end of Translate block
  # using Position to avoid scanning entire file
  start <- Position(function(ln) grepl("^\\s*Translate\\b", ln,
                                       ignore.case = TRUE, perl = TRUE),
                    lines)
  if (is.na(start)) return(NULL)
  rel <- Position(function(ln) grepl("^\\s*;", ln, perl = TRUE),
                  lines[(start + 1):length(lines)])
  end <- if (!is.null(rel)) start + rel else integer(0)

  # extract block (remove "Translate" line)
  block <- lines[(start + 1):(end - 1)]
  
  # remove commas and trim whitespace
  block <- gsub(",", "", block)
  block <- trimws(block)
  
  # remove empty lines
  block <- block[nzchar(block)]
  
  # extract number and name
  num <- as.integer(sub("^([0-9]+)\\s+.*", "\\1", block))
  name <- sub("^[0-9]+\\s+", "", block)
  
  # return data frame
  res <- data.frame(id = num, taxon = name, stringsAsFactors = FALSE)
  
  return(res)
}

#' Parse and organize annotations from a newick tree
#'
#' This function finds and sorts any node-specific annotation in a nexus file.
#'
#' @param newick A character containing an annotated tree in newick
#' (parenthetic) format (e.g., as output from \code{parse.nexus}).
#' 
#' @param taxonDiction A data frame similar to the output of 
#' \code{extract.translate.block.nexus}, containing the node \code{id} and 
#' the \code{taxon} associated with it.
#'
#' @return A data frame containing node information on lineage (taxon +
#' first/last occurrence), taxon (taxonomic name associated with the branch),
#' orientation, branch length, node name, and any other annotations
#' associated with every node in a tree. Note the function does not return a
#' \code{buddPhylo} object.
#'
#' @author Matheus Januario
#'
#' @export
#'

parse.newick.annotations <- function(newick, taxonDiction) {
  # extract tokens
  pattern <- "([A-Za-z0-9_]+)(\\[&[^]]*\\])?(?::([0-9\\.eE+-]+))?"
  matches <- gregexpr(pattern, newick, perl = TRUE)
  tokens <- regmatches(newick, matches)[[1]]
  
  # get labels
  labels <- sub("\\[&.*|:.*", "", tokens)
  
  # get annotations, when present
  has_annot <- grepl("[&", tokens, fixed = TRUE)
  annots <- rep(NA_character_, length(tokens))
  annots[has_annot] <- sub(".*\\[&([^]]*)\\].*", "\\1", tokens[has_annot])
  
  # get lengths, when present
  has_length <- grepl(":", tokens, fixed = TRUE)
  lengths <- rep(NA_real_, length(tokens))
  lengths[has_length] <- as.numeric(
    sub(".*:([0-9\\.eE+-]+).*", "\\1", tokens[has_length])
  )
  
  # remove labels and lengths for NA annots
  labels <- labels[!is.na(annots)]
  lengths <- lengths[!is.na(annots)]
  annots <- annots[!is.na(annots)]
  
  # split all annotations into each annotation
  split_list <- split_list <- strsplit(
    annots,
    "\\s*,\\s*(?=(?:[^{}]*\\{[^{}]*\\})*[^{}]*$)",
    perl = TRUE
  )
  
  
  # number of total nodes
  N <- length(annots)
  
  # number of constituent annotations per node
  parts_per_annot <- lengths(split_list)
  
  # the index of the rows for each annotation
  row_idx <- rep(seq_len(N), parts_per_annot)
  
  # flatten split_list to get all individual annotations
  flat <- unlist(split_list, use.names = FALSE)
  
  # keys and values
  key_vec <- sub("=.*",  "", flat)
  val_vec <- sub(".*=",  "", flat)
  
  # which annotations, if any, are ranges
  is_interval <- grepl("{", val_vec, fixed = TRUE)
  
  # get values per interval annotations
  int_raw <- gsub("[{}]", "", val_vec[is_interval])
  int_pair <- strsplit(int_raw, ",", fixed = TRUE)
  int_min <- as.numeric(vapply(int_pair, `[`, character(1), 1))
  int_max <- as.numeric(vapply(int_pair, `[`, character(1), 2))
  
  # start list of columns
  cols <- list()
  
  # scalar rows
  s_rows <- row_idx[!is_interval]
  s_keys <- key_vec[!is_interval]
  s_vals <- val_vec[!is_interval]
  
  # iterate through scalar keys
  for (k in unique(s_keys)) {
    # make the entire column NAs
    col <- rep(NA_character_, N)
    
    # choose which keys are this iter
    sel <- s_keys == k
    
    # add values to this column
    col[s_rows[sel]] <- s_vals[sel]
    
    # add col to cols
    cols[[k]] <- col
  }
  
  # interval rows
  i_rows <- row_idx[is_interval]
  i_keys <- key_vec[is_interval]
  
  
  # iterate through interval keys
  for (k in unique(i_keys)) {
    # choose which keys are this iter
    sel <- i_keys == k
    
    # start columns as NA
    col_min <- rep(NA_real_, N)
    col_max <- rep(NA_real_, N)
    
    # get min and max columns
    col_min[i_rows[sel]] <- int_min[sel]
    col_max[i_rows[sel]] <- int_max[sel]
    
    # put them into cols
    cols[[paste0(k, "_min")]] <- col_min
    cols[[paste0(k, "_max")]] <- col_max
  }
  
  # assemble data frame
  df <- data.frame(
    name   = labels,
    length = lengths,
    cols,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # extract and store lineage info:
  df$lineage <- 
    if (!is.null(taxonDiction)) taxonDiction$taxon[match(df$name, taxonDiction$id)]else df$name
  
  # If there is no "taxon" column, make it
  if (!"taxon" %in% colnames(df)) {
    df$taxon <- sub("_first|_last", "", df$lineage)
  }
  
  # last tidying:
  rownames(df) <- NULL
  
  # re ordering columns:
  neworder <- match(c("lineage", "taxon", "orientation", "length", "name"), colnames(df))
  df <- df[, c(neworder, (1:ncol(df))[-neworder])]
  
  return(df)
}

#' Build a buddPhylo object
#'
#' This function builds a \code{buddPhylo} object from the data contained in
#' a newick tree.
#'
#' @param newick A character containing an annotated tree in newick
#' (parenthetic) format (e.g., as output from \code{parse.nexus}).
#'
#' @param parsedTxt A data frame similar to the one generated by
#' \code{parse.newick.annotations}. It must contain node information on lineage
#' (taxon + first/last occurrence), taxon (taxonomic name associated with the
#' branch), orientation, branch length, node name, and any other annotations
#' associated with every node in a tree.
#'
#' @param tolExtant Rounding error tolerance for assigning a tip as extant,
#' given the distance in time from the tip that is further from the root.
#' Default is \code{10E-10}.
#'
#' @return A \code{buddPhylo} object.
#'
#' @author Matheus Januario
#'
#' @export
#'

build.buddPhylo <- function(newick, parsedTxt, tolExtant = 10E-10) {
  
  branches <- parsedTxt
  
  # 1. Read input as a phylo object
  phyl <- ape::read.tree(text = newick)
  # plot(phyl, cex=.5, show.node.label = T)
  tip_names <- phyl$tip.label
  if ("node.label" %in% names(phyl)) {
    node_names <- phyl$node.label
    assignNodeNamesManually <- FALSE
  } else {
    assignNodeNamesManually <- TRUE
  }
  
  ################################
  # Add parent branch information
  
  # Build full vector of node names (tips + internal nodes)
  Ntip <- length(phyl$tip.label)
  Nnode <- phyl$Nnode
  
  # If node names were missing, you already created them
  node_names <- phyl$node.label
  tip_names <- phyl$tip.label
  
  all_names <- c(phyl$tip.label, node_names)
  
  # edge matrix: parent -> child (in numeric indices)
  edges <- phyl$edge
  
  # Match each branch name to its parent
  parent_of <- setNames(all_names[edges[, 1]], all_names[edges[, 2]])
  branches$parent <- parent_of[branches$name]
  
  # The root (unique ancestor) will naturally get NA
  branches$orientation[is.na(branches$parent)] <- "uca"
  branches$parent[branches$parent == ""] <- NA # branches$name[branches$parent==""]
  
  # now mark which are sampled ancestors:
  sampAncs_names <- tip_names
  sampAncs_names <- sampAncs_names[sampAncs_names %in% branches$name[branches$length == 0]]
  
  branches$type <- "node"
  branches$type[branches$name %in% tip_names] <- "tip"
  branches$type[branches$name %in% sampAncs_names] <- "sampAnc"
  
  phylo <- ape::keep.tip(phyl, branches$name[branches$type == "tip"])
  
  # if phylo only has one tip, we need to use the whole tree to find the coords
  if (length(phylo$tip.label) == 1) phylo <- phyl
  
  branches <- fix.coords(branches, phylo)
  
  # Flagging extant species:
  branches$extant <- FALSE
  branches$extant[max(branches$x_coord, na.rm = T) - branches$x_coord < tolExtant] <- TRUE
  
  # correct lineage and taxon to NA for nodes
  branches$lineage[branches$type == "node"] <- 
    branches$taxon[branches$type == "node"] <- NA
  
  # set some helper variables for reordering
  N_row <- nrow(branches)
  par_idx <- match(branches$parent, branches$name)
  ch_by_par <- split(seq_len(N_row), par_idx)
  
  # reorder with a BFS depth from the root
  depth <- integer(N_row)
  frontier <- which(is.na(par_idx))
  while (length(frontier)) {
    ch <- unlist(ch_by_par[as.character(frontier)], use.names = FALSE)
    if (!length(ch)) break
    depth[ch] <- depth[par_idx[ch]] + 1
    frontier <- ch
  }
  
  # stable sort by ascending depth — root first, deepest last
  branches <- branches[order(depth), , drop = FALSE]
  rownames(branches) <- NULL
  
  class(branches) <- c("buddPhylo", "data.frame")
  
  return(branches)
}

#' Read a buddPhylo tree from a file in nexus format
#'
#' This function reads one or several trees in a NEXUS file.
#'
#' @param file 	A character containing a nexus file name.
#'
#' @param treeNumbers A numeric vector containing which trees should be read.
#' If \code{NULL} (default), the function reads all trees within \code{file}.
#'
#' @param timeFromRoot A logical indicating whether time should be expressed as
#' time since the root (\code{timeFromRoot=TRUE}). If set to \code{FALSE}
#' (default), time is expressed as millions of years ago (Mya), taking the tip extant tip that is further from the root as a reference for "present".
#'
#' @return A \code{buddPhylo} with the time adjusted accordingly.
#'
#' @author Matheus Januario
#' 
#' @export
#'

read.nexus.buddPhylo <- function(file, treeNumbers = NULL, timeFromRoot = FALSE) {
  full_text <- readLines(file)
  
  # parse file info & build buddPhylo:
  txt <- parse.nexus(text = full_text, treeNumbers = treeNumbers)
  
  translate_block <- extract.translate.block.nexus(text = full_text)
  branches <- lapply(txt, parse.newick.annotations, 
                     taxonDiction = translate_block)
  
  res <- list()
  for (i in 1:length(txt)) {
    phy <- build.buddPhylo(newick = txt[[i]], parsedTxt = branches[[i]])
    res[[i]] <- adjust.timescale.buddPhylo(phy)
  }
  
  if (length(res) == 1) {
    res <- res[[1]]
  } else {
    class(res) <- "multi.buddPhylo"
  }
  
  return(res)
}

#' Read a buddPhylo tree from a file in Newick format
#'
#' This function reads one or several trees in Newick file.
#'
#' @param file 	A character containing a nexus file name.
#'
#' @param treeNumbers A numeric vector containing which trees should be read.
#' If \code{NULL} (default), the function reads all trees within \code{file}.
#'
#' @param timeFromRoot A logical indicating whether time should be expressed as
#' time since the root (\code{timeFromRoot=TRUE}). If set to \code{FALSE}
#' (default), time is expressed as millions of years ago (Mya), taking the tip extant tip that is further from the root as a reference for "present".
#'
#' @return A \code{buddPhylo} with the time adjusted accordingly.
#'
#' @author Bruno do Rosario Petrucci
#' 
#' @export
#'

read.trees.buddPhylo <- function(file, treeNumbers = NULL, timeFromRoot = FALSE) {
  txt <- readLines(file)
  
  # parse file info & build buddPhylo:
  branches <- lapply(txt, parse.newick.annotations, 
                     taxonDiction = NULL)
  
  res <- list()
  for (i in 1:length(txt)) {
    phy <- build.buddPhylo(newick = txt[[i]], parsedTxt = branches[[i]])
    res[[i]] <- adjust.timescale.buddPhylo(phy)
  }
  
  if (length(res) == 1) {
    res <- res[[1]]
  } else {
    class(res) <- "multi.buddPhylo"
  }
  
  return(res)
}