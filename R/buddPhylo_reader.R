#' Subset and parse a nexus file
#'
#' This function parses one or several budding phylogenetic trees in a Nexus 
#' file.
#'
#' @param text A string with the tree in a Newick format. Optional, but at least
#' one of \code{text} or \code{file} must be not \code{NULL}.
#' 
#' @param file A character containing a Nexus filename. Optional, but at least
#' one of \code{text} or \code{file} must be not \code{NULL}.
#'
#' @param treeNumbers A numeric vector containing which trees should be read.
#' If \code{NULL} (default), the function reads all trees within \code{file}.
#'
#' @return A character with the parenthetic format within the nexus file.
#' If more than one tree is read, a list of characters is returned in the
#' order given by \code{treeNumbers}.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'
#' @details The does not take note on  whichever name was assigned for any tree
#' in the nexus file.
#'

parse.nexus <- function(text = NULL, file = NULL, treeNumbers = NULL) {
  # check if text is null
  if (is.null(text)) {
    # check if file is also null
    if (is.null(file)) {
      stop("Either file or text must not be null.")
    }
    
    # otherwise, read file
    lines <- readLines(file)
  } else {
    # otherwise, read text
    lines <- text
  }
  
  # find the lines with the tree-indicating pattern
  target_line <- grep("tree ", lines, value = TRUE, ignore.case = TRUE)
  
  # restrict to the treeNumbers lines if it is not NULL
  if (!is.null(treeNumbers)) target_line <- target_line[treeNumbers]
  trees <- sub('^\\s*tree\\s+\\S+\\s*=\\s*"?', "", target_line, 
               perl = TRUE, ignore.case = TRUE)
  trees <- sub('"\\s*$', "", trees, perl = TRUE)
  
  return(trees)
}

#' Extract the "translation" block in a nexus file
#'
#' This function extracts the tip translation block (the block relating numbers
#' to taxon names) in a Nexus file.
#'
#' @inheritParams parse.nexus
#'
#' @return A data frame containing the number \code{id} and the \code{taxon}
#' associated with it.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'

extract.translate.block.nexus <- function(text = NULL, file = NULL) {
  # check if text is null
  if (is.null(text)) {
    # check if file is also null
    if (is.null(file)) {
      stop("Either file or text must not be null.")
    }
    
    # otherwise, read file
    lines <- readLines(file)
  } else {
    # otherwise, read text
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
#' @return A data frame containing node information on orientation, branch 
#' length, node name, and any other annotations associated with every node in 
#' the tree. Note the function does not return a \code{buddPhylo} object.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci
#'

parse.newick.annotations <- function(newick, taxonDiction = NULL) {
  # ensure we are only using this if it's a buddPhylo
  if (nrow(df) == 0) {
    stop("No annotated nodes found. buddPhylo trees require [&...] ",
         "annotations (e.g. orientation) on each node; a plain Newick ",
         "tree can be read with ape::read.tree instead.")
  }
  
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
    name = labels,
    length = lengths,
    cols,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # translate name via taxonDiction if provided 
  if (!is.null(taxonDiction)) {
    m <- match(df$name, taxonDiction$id)
    df$name[!is.na(m)] <- taxonDiction$taxon[m[!is.na(m)]]
  }
  
  # tidy
  rownames(df) <- NULL
  
  # put name, length, orientation up front; leave annotation cols after
  front <- intersect(c("orientation", "length", "name"), colnames(df))
  df <- df[, c(front, setdiff(colnames(df), front))]
  
  return(df)
}

#' Build a buddPhylo object
#'
#' This function builds a \code{buddPhylo} object from the data contained in
#' a newick tree.
#'
#' @inheritParams parse.newick.annotations
#'
#' @param parsedTxt A data frame similar to the one generated by
#' \code{parse.newick.annotations}. It must contain node information on 
#' orientation, branch length, node name, and any other annotations associated 
#' with every node in a tree.
#'
#' @param tolExtant Rounding error tolerance for assigning a tip as extant,
#' given the distance in time from the tip that is further from the root.
#' Default is \code{10E-10}.
#'
#' @return A \code{buddPhylo} object.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci.
#'

build.buddPhylo <- function(newick, parsedTxt, 
                            taxonDiction = NULL, tolExtant = 10E-10) {
  # start the branches data frame with the parsed annotations
  branches <- parsedTxt
  
  # get the tree as a phylo object
  phyl <- ape::read.tree(text = newick)
  
  # get tip names
  tip_names <- phyl$tip.label
  
  # match tip names with taxon dictionary, if it exists
  if (!is.null(taxonDiction)) {
    m <- match(tip_names, taxonDiction$id)
    tip_names[!is.na(m)] <- taxonDiction$taxon[m[!is.na(m)]]
    phyl$tip.label <- tip_names
  }
  
  # build full vector of node names (tips + internal nodes)
  Ntip <- length(phyl$tip.label)
  Nnode <- phyl$Nnode
  
  # collect them together
  all_names <- c(tip_names, phyl$node.label)
  
  # edge matrix: parent -> child (in numeric indices)
  edges <- phyl$edge
  
  # match each branch name to its parent
  parent_of <- setNames(all_names[edges[, 1]], all_names[edges[, 2]])
  branches$parent <- parent_of[branches$name]
  
  # the root (unique ancestor) will get its own orientation
  branches$orientation[is.na(branches$parent)] <- "uca"

  # now mark which are sampled ancestors
  sampAncs_names <- tip_names
  sampAncs_names <- sampAncs_names[sampAncs_names %in% 
                                     branches$name[branches$length == 0]]
  
  # and assign the type of node for each branch
  branches$type <- "node"
  branches$type[branches$name %in% tip_names] <- "tip"
  branches$type[branches$name %in% sampAncs_names] <- "sampAnc"
  
  # get a phylo object with no sampled ancestors to find coords
  phylo <- ape::keep.tip(phyl, branches$name[branches$type == "tip"])
  
  # if phylo only has one tip, we need to use the whole tree to find the coords
  if (length(phylo$tip.label) == 1) phylo <- phyl
  
  # find the x and y coords of branches
  branches <- fix.coords(branches, phylo)
  
  # flag extant species
  branches$extant <- FALSE
  branches$extant[max(branches$x_coord, na.rm = TRUE) - 
                    branches$x_coord < tolExtant] <- TRUE
  
  # derive lineage and taxon from name + type
  tip_like <- branches$type %in% c("tip", "sampAnc")
  branches$lineage <- NA_character_
  branches$lineage[tip_like] <- branches$name[tip_like]
  branches$taxon <- NA_character_
  branches$taxon[tip_like] <- 
    sub("_(first|last)$", "", branches$name[tip_like])
  
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
  
  # final reorder of columns
  front <- c("lineage", "taxon", "orientation", "length", "name")
  branches <- branches[, c(front, setdiff(colnames(branches), front))]
  
  # set classes
  class(branches) <- c("buddPhylo", "data.frame")
  
  return(branches)
}

#' Read a buddPhylo tree from a file in Nexus format
#'
#' This function reads one or several trees in a Nexus file.
#'
#' @inheritParams parse.nexus
#'
#' @param timeFromRoot A logical indicating whether time should be expressed as
#' time since the root (\code{TRUE}). If set to \code{FALSE} (default), time is 
#' expressed as millions of years ago (Mya), taking the tip extant tip that is 
#' further from the root as a reference for "present".
#'
#' @return A \code{buddPhylo} with the time adjusted accordingly.
#'
#' @author Matheus Januario and Bruno do Rosario Petrucci
#' 
#' @export
#'

read.nexus.buddPhylo <- function(file, treeNumbers = NULL, 
                                 timeFromRoot = FALSE) {
  # read the entire file
  full_text <- readLines(file)
  
  # parse file
  txt <- parse.nexus(text = full_text, treeNumbers = treeNumbers)
  
  # get translate block
  translate_block <- extract.translate.block.nexus(text = full_text)
  
  # parse newick annotations
  branches <- lapply(txt, parse.newick.annotations, 
                     taxonDiction = translate_block)
  
  # start list of trees
  res <- list()
  
  # iterate through trees in the file
  for (i in 1:length(txt)) {
    # build this tree
    phy <- build.buddPhylo(newick = txt[[i]], parsedTxt = branches[[i]],
                           taxonDiction = translate_block)

    # adjust the timescale
    res[[i]] <- adjust.timescale.buddPhylo(phy, timeFromRoot = timeFromRoot)
  }
  
  # if we only have one tree, just take the first element
  if (length(res) == 1) {
    res <- res[[1]]
  } else {
    # otherwise, set class
    class(res) <- "multi.buddPhylo"
  }
  
  return(res)
}

#' Read a buddPhylo tree from a file in Newick format
#'
#' This function reads one or several trees in Newick file.
#'
#' @param file A character containing a Newick file name.
#' 
#' @param text A string with one or several trees in Newick format, separated
#' by newlines (\code{\\n}). \code{NULL} by default. If not \code{NULL}, the 
#' argument \code{file} is ignored.
#'
#' @inheritParams read.nexus.buddPhylo
#'
#' @return A \code{buddPhylo} with the time adjusted accordingly.
#'
#' @author Bruno do Rosario Petrucci.
#' 
#' @export
#'

read.trees.buddPhylo <- function(file, text = NULL, treeNumbers = NULL, 
                                 timeFromRoot = FALSE) {
  # check if text is null
  if (is.null(text)) {
    # if so, read from file
    txt <- readLines(file)
  } else {
    # otherwise, just use text
    txt <- text
  }
  
  # take just the treeNumbers supplied
  if (!is.null(treeNumbers)) txt <- txt[treeNumbers]
  
  # parse newick annotations
  branches <- lapply(txt, parse.newick.annotations, 
                     taxonDiction = NULL)
  
  # start list of trees
  res <- list()
  
  # iterate through trees in the file
  for (i in 1:length(txt)) {
    # build this tree
    phy <- build.buddPhylo(newick = txt[[i]], parsedTxt = branches[[i]])
    
    # adjust the timescale
    res[[i]] <- adjust.timescale.buddPhylo(phy, timeFromRoot = timeFromRoot)
  }
  
  # if we only have one tree, just take the first element
  if (length(res) == 1) {
    res <- res[[1]]
  } else {
    # otherwise, set class
    class(res) <- "multi.buddPhylo"
  }
  
  return(res)
}