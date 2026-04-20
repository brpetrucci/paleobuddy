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
  
  # Find the line with the tree-indicating pattern
  target_line <- grep("tree ", lines, value = TRUE)
  
  if (length(target_line) > 1) {
    lines <- as.list(target_line)
    foo <- function(x) {
      sub('.*tree MCC_tree =\\s*"([^"]*)".*', "\\1", x)
    }
    
    if (!is.null(treeNumbers)) {
      trees <- lapply(lines[treeNumbers], foo)
    } else {
      trees <- lapply(lines, foo)
    }
    return(trees)
  } else {
    # Extract text inside quotes after the pattern
    mcc_tree <- sub('.*tree MCC_tree =\\s*"([^"]*)".*', "\\1", target_line)
    return(mcc_tree)
  }
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
  start <- grep("^\\s*Translate\\b", lines, ignore.case = TRUE)
  end <- grep("^\\s*;", lines)
  
  # keep the first ";" after Translate
  end <- end[end > start][1]
  
  if (length(start) == 0 || length(end) == 0) {
    return(NULL)
    #stop("Translate block not found.")
  }
  
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
  # Step 1: Extract tokens (label + annotation + length)
  pattern <- "([A-Za-z0-9_]+)(\\[&[^]]*\\])?(?::([0-9\\.eE+-]+))?"
  matches <- gregexpr(pattern, newick, perl = TRUE)
  tokens <- regmatches(newick, matches)[[1]]
  
  # Step 2: Extract components
  get_label <- function(x) sub("\\[&.*|:.*", "", x)
  
  get_annot <- function(x) {
    if (grepl("\\[&", x)) {
      sub(".*\\[&([^]]*)\\].*", "\\1", x)
    } else {
      NA
    }
  }
  
  get_length <- function(x) {
    if (grepl(":", x)) {
      sub(".*:([0-9\\.eE+-]+).*", "\\1", x)
    } else {
      NA
    }
  }
  
  labels <- sapply(tokens, get_label)
  annots <- sapply(tokens, get_annot)
  lengths <- as.numeric(sapply(tokens, get_length))
  
  #  Step 3: split annotation safely (ignores commas inside {})
  split_annot <- function(a) {
    parts <- strsplit(a, ",(?=(?:[^{}]*\\{[^{}]*\\})*[^{}]*$)", perl = TRUE)[[1]]
    trimws(parts)
  }
  
  #  Step 4: Collect all keys (including interval-related keys)
  all_keys <- c()
  
  for (a in annots[!is.na(annots)]) {
    parts <- split_annot(a)
    for (p in parts) {
      key <- sub("=.*", "", p)
      
      if (grepl("\\{", p)) {
        all_keys <- c(all_keys, paste0(key, "_min"), paste0(key, "_max"))
      } else {
        all_keys <- c(all_keys, key)
      }
    }
  }
  
  all_keys <- unique(all_keys)
  
  #  Step 5:  Initialize dataframe
  df <- data.frame(
    name = labels,
    length = lengths,
    stringsAsFactors = FALSE
  )
  
  for (k in all_keys) df[[k]] <- NA
  
  #  Step 6:  Fill data.frame
  for (i in seq_along(annots)) {
    if (!is.na(annots[i])) {
      parts <- split_annot(annots[i])
      
      for (p in parts) {
        key <- sub("=.*", "", p)
        val <- sub(".*=", "", p)
        
        # interval case {a,b}
        if (grepl("^\\{.*\\}$", val)) {
          nums <- gsub("[\\{\\}]", "", val)
          nums <- strsplit(nums, ",")[[1]]
          
          df[i, paste0(key, "_min")] <- as.numeric(nums[1])
          df[i, paste0(key, "_max")] <- as.numeric(nums[2])
        } else {
          # normal scalar value
          df[i, key] <- val
        }
      }
    }
  }
  
  # extract and sotre lineage info:
  df$lineage <- taxonDiction$taxon[match(df$name, taxonDiction$id)]
  
  # If there is no "taxon" column, make it
  if (!"taxon" %in% colnames(df)) {
    df$taxon <- sub("_first", "", df$lineage)
    df$taxon <- sub("_last", "", df$taxon)
  }
  
  # last tidying:
  rownames(df) <- NULL
  df=df[-grep("tree|STATE_", df$name),]
  
  
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

build.buddPhylo <- function(newick, parsedTxt, tolExtant = 10E-10, old_way = FALSE) {
  
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
  
  # Convert to data.frame with names instead of indices
  edge_df <- data.frame(
    parent = all_names[edges[, 1]],
    child = all_names[edges[, 2]],
    stringsAsFactors = FALSE
  )
  
  # Match each branch name to its parent
  branches$parent <- edge_df$parent[match(branches$name, edge_df$child)]
  
  # The root (unique ancestor) will naturally get NA
  branches$orientation[is.na(branches$parent)] <- "uca"
  branches$parent[branches$parent == ""] <- NA # branches$name[branches$parent==""]
  
  # now mark which are sampled ancestors:
  sampAncs_names <- tip_names
  sampAncs_names <- sampAncs_names[sampAncs_names %in% branches$name[branches$length == 0]]
  
  branches$type <- "node"
  branches$type[branches$name %in% tip_names] <- "tip"
  branches$type[branches$name %in% sampAncs_names] <- "sampAnc"
  
  dropSampAncs <- TRUE
  if (dropSampAncs) {
    treeWithSampAncs <- phyl
    phylo <- ape::drop.tip(phyl, sampAncs_names)
  }
  
  phylo <- ape::keep.tip(phylo, branches$name[branches$type %in% c("tip")])
  phylo <- ape::ladderize(phylo)
  
  if (old_way) {
    obj <- ape::plot.phylo(phylo, plot = FALSE)
    
    if (is.null(obj$yy)) {
      grDevices::pdf(NULL)  # open null device (no file created)
      on.exit(grDevices::dev.off(), add = TRUE)
      
      obj <- ape::plot.phylo(phylo, plot = FALSE)
      obj <- get("last_plot.phylo", envir = utils::getFromNamespace(".PlotPhyloEnv", "ape"))
    }
    
    yy <- obj$yy
    xx <- obj$xx
    
    if ("node.label" %in% names(phylo)) {
      names(yy) <- c(phylo$tip.label, phylo$node.label)
      names(xx) <- names(yy)
    } else {
      names(yy)[1:Ntip(phylo)] <- phylo$tip.label
      names(xx) <- names(yy)
    }
    
    # drop yy thata re irrelevant to budding:
    valid_refs <- which(branches$type == "tip" | branches$orientation == "none")
    yy <- yy[names(yy) %in% branches$name[valid_refs]]
    xx <- xx[names(xx) %in% branches$name[valid_refs]]
    
    branches$y_coord <- NA
    branches$x_coord <- NA
    
    matching_ids <- match(names(yy), branches$name)
    
    branches$y_coord[matching_ids] <- yy
    branches$x_coord[matching_ids] <- xx
    
    # now re-assign y values recursively from tip to root:
    tipNames <- unique(branches$name[branches$type == "tip"])
    branches$x_par <- NA
    branches$y_par <- NA
    tt <- tip_names[1]
    for (tt in tipNames) {
      todo <- TRUE
      while (todo) {
        ttID <- which(branches$name == tt)
        Parental <- getParents.buddPhylo(branches, tt)[1]
        ParID <- which(branches$name == Parental)
        
        # update the refences coordinates:
        if (!is.na(branches$x_coord[ParID])) {
          branches$x_par[ttID] <- branches$x_coord[ttID] - branches$length[ttID]
          todo <- FALSE
        } else {
          branches$y_coord[ParID] <- branches$y_coord[ttID]
          branches$x_coord[ParID] <- branches$x_coord[ttID] - branches$length[ttID]
          branches$x_par[ttID] <- branches$x_coord[ttID] - branches$length[ttID]
        }
        
        if (branches$orientation[ParID] == "uca") {
          branches$x_par[ParID] <- branches$x_coord[ParID] - branches$length[ParID]
          branches$y_par[ParID] <- branches$y_coord[ParID]
          todo <- FALSE
        } else {
          tt <- Parental
        }
      }
    }
    
    # Now do the sampling ancestors:
    SAncNames <- unique(branches$name[branches$type == "sampAnc"])
    ss <- SAncNames[1]
    for (ss in SAncNames) {
      SAID <- which(branches$name == ss)
      SAParID <- which(branches$name == branches$parent[SAID])
      
      branches$y_coord[SAID] <- branches$y_coord[SAParID]
      branches$x_coord[SAID] <- branches$x_coord[SAParID]
      branches$x_par[SAID] <- branches$x_coord[SAParID]
    }
    
    min_x <- min(c(branches$x_coord, branches$x_par), na.rm = T)
    
    branches$x_coord <- branches$x_coord - min_x
    branches$x_par <- branches$x_par - min_x
    
    # Now for descendants
    branches$y_par <- branches$y_coord[match(branches$parent, branches$name)]
    branches$y_par[branches$orientation == "uca"] <- branches$y_coord[branches$orientation == "uca"]
    # branches$x_par = branches$x_coord[match(branches$parent, branches$name)]
  } else {
    branches <- fix.coords(branches, phylo)
  }
  
  # Flagging extant species:
  branches$extant <- FALSE
  branches$extant[max(branches$x_coord, na.rm = T) - branches$x_coord < tolExtant] <- TRUE
  
  class(branches) <- c("buddPhylo", "data.frame")
  
  return(branches)
}

#' Parse and organize annotations from a newick tree
#'
#' This function finds and sorts any node-specific annotation in a nexus file.
#'
#' @param newick A character containing an annotated tree in newick 
#' (parenthetic) format (e.g., as output from \code{parse.nexus}).
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

parse.newick.annotations.notrans <- function(newick) {
  
  # Step 1: Extract tokens (label + annotation + length)
  pattern <- "([A-Za-z0-9_]+)(\\[&[^]]*\\])?(?::([0-9\\.eE+-]+))?"
  matches <- gregexpr(pattern, newick, perl = TRUE)
  tokens <- regmatches(newick, matches)[[1]]
  
  # Step 2: Extract components
  get_label <- function(x) sub("\\[&.*|:.*", "", x)
  
  get_annot <- function(x) {
    if (grepl("\\[&", x)) {
      sub(".*\\[&([^]]*)\\].*", "\\1", x)
    } else {
      NA
    }
  }
  
  get_length <- function(x) {
    if (grepl(":", x)) {
      sub(".*:([0-9\\.eE+-]+).*", "\\1", x)
    } else {
      NA
    }
  }
  
  labels  <- sapply(tokens, get_label)
  annots  <- sapply(tokens, get_annot)
  lengths <- as.numeric(sapply(tokens, get_length))
  
  #  Step 3: split annotation safely (ignores commas inside {})
  split_annot <- function(a) {
    parts <- strsplit(a, ",(?=(?:[^{}]*\\{[^{}]*\\})*[^{}]*$)", perl = TRUE)[[1]]
    trimws(parts)
  }
  
  #  Step 4: Collect all keys (including interval-related keys)
  all_keys <- c()
  
  for (a in annots[!is.na(annots)]) {
    parts <- split_annot(a)
    for (p in parts) {
      key <- sub("=.*", "", p)
      
      if (grepl("\\{", p)) {
        all_keys <- c(all_keys, paste0(key, "_min"), paste0(key, "_max"))
      } else {
        all_keys <- c(all_keys, key)
      }
    }
  }
  
  all_keys <- unique(all_keys)
  
  #  Step 5:  Initialize dataframe
  df <- data.frame(
    name   = labels,
    length = lengths,
    stringsAsFactors = FALSE
  )
  
  for (k in all_keys) df[[k]] <- NA
  
  #  Step 6:  Fill data.frame
  for (i in seq_along(annots)) {
    if (!is.na(annots[i])) {
      parts <- split_annot(annots[i])
      
      for (p in parts) {
        key <- sub("=.*", "", p)
        val <- sub(".*=", "", p)
        
        # interval case {a,b}
        if (grepl("^\\{.*\\}$", val)) {
          nums <- gsub("[\\{\\}]", "", val)
          nums <- strsplit(nums, ",")[[1]]
          
          df[i, paste0(key, "_min")] <- as.numeric(nums[1])
          df[i, paste0(key, "_max")] <- as.numeric(nums[2])
          
        } else {
          # normal scalar value
          df[i, key] <- val
        }
      }
    }
  }
  
  # If there is no "taxon" column, make it
  if(!"taxon" %in% colnames(df)){
    #taxonDiction = extract.translate.block.nexus(file)
    df$lineage= df$name#taxonDiction$taxon[match(df$name, taxonDiction$id)]
    df$taxon = sub("_first", "", df$lineage)
    df$taxon = sub("_last", "", df$taxon)
  }
  
  # last tidying:
  rownames(df)=NULL
  df=df[-grep("tree|STATE_", df$name),]

  # re ordering columns:
  neworder = match(c("lineage", "taxon", "orientation", "length", "name"), colnames(df))
  df = df[, c(neworder, (1:ncol(df))[-neworder])]
  
  return(df)
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

read.nexus <- function(file, treeNumbers = NULL, timeFromRoot = FALSE, old_way = FALSE) {
  full_text <- readLines(file)
  
  # parse file info & build buddPhylo:
  txt <- parse.nexus(text = full_text, treeNumbers = treeNumbers)
  if (length(txt) > 1) {
    translate_block <- extract.translate.block.nexus(text = full_text)
    if (is.null(translate_block)) {
      branches <- lapply(txt, function(x) parse.newick.annotations.notrans(x))
    } else {
      branches <- lapply(txt, function(x) 
        parse.newick.annotations(x, taxonDiction = translate_block))
    }
    res <- list()
    for (i in 1:length(txt)) {
      phy <- build.buddPhylo(newick = txt[[i]], parsedTxt = as.data.frame(branches[i]), old_way = old_way)
      res[[i]] <- adjust.timescale.buddPhylo(phy)
    }
    class(res) <- "multi.buddPhylo"
  } else {
    txt <- txt[[1]]
    
    translate_block <- extract.translate.block.nexus(text = full_text)
    if (is.null(translate_block)) {
      branches <- parse.newick.annotations.notrans(txt)
    } else {
      branches <- parse.newick.annotations(txt, taxonDiction = translate_block)
    }
    
    phy <- build.buddPhylo(newick = txt, parsedTxt = branches, old_way = old_way)
    res <- adjust.timescale.buddPhylo(phy)
  }
  
  return(res)
}
