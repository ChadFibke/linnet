#' Find the most recent common ancestor (MRCA)
#'
#' Internal function used by `detect_groups()`.
#'
#' @keywords internal
.mrca <- function(graph_obj, lineages) {
  
  # Find location of lineages
  lineage_indices <- match(lineages, igraph::V(graph_obj)$name)
  
  # Collect the pathway between each lineage to root (inclusive)
  anscestor_list <- lapply( lineage_indices,
                            function(lineage)  c(names( unlist(
                              igraph::all_simple_paths( graph_obj, from = lineage,
                                                to = 1, mode = "in")
                            )), igraph::V(graph_obj)$name[lineage] ) )
  
  # Remove null values
  anscestor_list <- Filter(Negate(is.null), anscestor_list)
  
  # Find common lineages found across all pathways
  common_anscestors <- Reduce(intersect, anscestor_list)
  
  # Return the most recent common lineage
  mrca <- igraph::V(graph_obj)$name[ max(match(common_anscestors, igraph::V(graph_obj)$name) ) ]
  
  return( mrca )
}


#' Subsets network to observed lineages
#'
#' Internal function used by `detect_groups()`.
#'
#' @keywords internal
.subset_graph <- function(graph_obj, root, lineages){

  # Find location of lineages in graph
  lineage_indices <- match(lineages, igraph::V(graph_obj)$name )
  
  # Find new root location
  root_index <- match(root, igraph::V(graph_obj)$name )
  
  # Return singleton
  if ( all(lineage_indices == root_index) ) {
    
    return( igraph::induced_subgraph(graph_obj, root_index )  )
    
  }
  
  # Collect all pathways between root and lineages of interest
  sub_list <- anscestor_list <- lapply( lineage_indices,
                                        function(lineage)  names( unlist(
                                          igraph::all_simple_paths( graph_obj,
                                                                    from = lineage,
                                                                    to = root_index,
                                                                    mode = "in")
                                          )
                                          )
                                        )
  
  # Find all unique lineages in all pathways
  subset_members <- Reduce(union, sub_list)
  
  # Extract tree of interest
  subtree <- igraph::induced_subgraph(graph_obj, subset_members)
  
  return(subtree)
}


#' Assign a numeric attribute to vertices in network
#'
#' Internal function used by `detect_groups()`.
#'
#' @keywords internal
.annotate <- function(graph_obj, lineages, values, default = 0 ){
  
  # Collect list of lineages in a tree
  all_lineages <- igraph::V(graph_obj)$name
  
  # Create a vector of attributes for each lineage
  base_attributes <- rep( default, (length(all_lineages)) )
  
  # Identify location of lineages within provided graph
  matching_indices <- match(lineages, all_lineages)
  
  # Assign attributes to matching lineages
  base_attributes[matching_indices] <- values
  
  annotated_graph <- igraph::set_vertex_attr(graph = graph_obj,
                                             name = "value", 
                                             value = base_attributes)
  
  return(annotated_graph)
}


# Clip and append recombinant clades to root
#'
#' Internal function used by `detect_groups()`.
#'
#' @keywords internal
.split_recombinants <- function(graph_obj){
  
  # Collect all recombinants
  graph_degrees <- igraph::degree(graph_obj, mode = "in")
  
  rcbnt <- names(graph_degrees)[ graph_degrees > 1 ]
  
  rcbnt_indices <- match(rcbnt, igraph::V(graph_obj)$name)
  
  # Find root name
  root <- igraph::V(graph_obj)$name[1]
  
  # For each recombinant clade, attach clade to root of tree
  for (rcbnt_start in rcbnt_indices) {
    
    # Cut of the recombinant clade
    headss <- igraph::E(graph_obj)[ unlist(
      igraph::incident_edges( graph_obj, rcbnt_start, mode = "in")
      ) ]
    
    graph_obj <- igraph::delete_edges(graph_obj, headss)
    
    
    # Stitch recombinant head to the root node using names
    rcbnt_init <- igraph::V(graph_obj)$name[rcbnt_start]
    
    graph_obj <- igraph::add_edges(graph_obj, c(root, rcbnt_init))
    
  }
  
  return(graph_obj)
}


# Cumulatively sums with reset above threshold from leaves to root using a breadth first search
#'
#' Internal function used by `detect_groups()`.
#'
#' @keywords internal
.threshold_cumsum <- function(graph_obj, threshold){
  
  # Set up
  size <- igraph::vcount(graph_obj)
  
  cumsum_values <- numeric(size)
  capped  <- logical(size)
  
  segment_id <- integer(size)
  current_segment <- 1
  
  # Collect bfs orders to properly sum while traversing upward
  topo_order <- as.integer(igraph::topo_sort(graph_obj, mode = "in") )
  
  # Collect direct children per lineage
  children_list <- lapply( V(graph_obj),
                           function(LINEAGE) as.integer(igraph::neighbors(graph_obj, LINEAGE, mode="out")))
  
  # Propagate values as cumsum with reset
  for (INDEX in topo_order) {
    
    children <- children_list[[INDEX]]
    
    # Sum uncapped children
    if (length(children) > 0) {
      
      contribution <- sum(cumsum_values[children] * !capped[children] )
      
    } else {
      
      contribution <- 0
      
    }
    
    cumulative <- contribution + igraph::V(graph_obj)$value[INDEX]
    
    cumsum_values[INDEX] <- cumulative
    
    
    # Cap if threshold reached
    if ( cumulative >= threshold ) {
    
      capped[INDEX] <- TRUE
      
      segment_id[INDEX] <- current_segment
      
      current_segment <- current_segment + 1
       
      
    } else {
      
      segment_id[INDEX] <- current_segment
      
    }
  }
  
  # Assign cumulative value
  igraph::V(graph_obj)$value <- cumsum_values
  igraph::V(graph_obj)$segment <- segment_id
  
  return(graph_obj)
}


# Iteratively prune tree from bottom up until N non-zero lineages remain
#'
#' Internal function used by `detect_groups()`.
#'
#' @keywords internal
.topn <- function(graph_obj, n){
  
  lineage_tracker <- igraph::V(graph_obj)$value
  names(lineage_tracker) <- igraph::V(graph_obj)$name
  
  
  # Identify distance between all lineages ordered by placement in tree
  distances_matrix <- igraph::distances(graph_obj,
                                v = names(lineage_tracker),
                                to = names(lineage_tracker),
                                mode = "all")
  
  distances_matrix[is.infinite(distances_matrix) | distances_matrix == 0 ] <- NA
  
  # Child tracker
  children_list <- lapply( V(graph_obj),
                           function(LINEAGE) igraph::as_ids(igraph::neighbors(graph_obj, LINEAGE, mode="out")))
  
  
  while (length( lineage_tracker[ lineage_tracker > 0]) > n) {
    
    # Collect leaf lineage
    leaves <- lengths( children_list[ names(lineage_tracker[ lineage_tracker > 0])] ) 
    
    # order by value, while prioritizing leaves
    mate1 <- names( which.min(lineage_tracker[lineage_tracker > 0][ order( lineage_tracker[lineage_tracker > 0], leaves ) ] ))
    
    
    # Collect the closest lineage(s)
    mate2 <- names(which( distances_matrix[mate1, names(lineage_tracker) ] == min(distances_matrix[mate1, names(lineage_tracker)], na.rm = TRUE)) )
    
    # If there are multiple lineage collect the lineage with most similar values
    # or smaller
    if ( length(mate2) > 1) {
      
      mate2 <-  names( which.min( lineage_tracker[mate2] - lineage_tracker[mate1] ) )
    }
    
    collapsed_group = .mrca(graph_obj, lineages = c(mate1, mate2) )
    
    # Sum and assign values 
    lineage_tracker[ collapsed_group ] <- sum( lineage_tracker[ c(mate1, mate2)], na.rm = TRUE )
    
    to_remove <- setdiff( c(mate1, mate2), collapsed_group )
    
    # Update tracker
    lineage_tracker <- lineage_tracker[ !( names(lineage_tracker) %chin% to_remove ) ]
    
    # Update children list
    children_list[[to_remove]] <- NULL
    children_list <- lapply(children_list, function(x) { x[!x %in% to_remove] })
    
  }
  
  
  return(names( lineage_tracker[lineage_tracker > 0 ]) )
}


#' Empirically Determine Predominant Lineage Groups
#' 
#' `detect_groups()` returns SARS-CoV-2 lineages which can be used to as groups
#'  for `categorize()`.
#'  
#' @importFrom igraph V
#' 
#' @param graph_obj an igraph object constructed from `read_lineages()`
#'     depicting all SARS-CoV-2 lineage parent/child relationships.
#' 
#' @param lineages  A vector of observed SARS-CoV-2 lineages.
#' 
#' @param values  A vector of observed SARS-CoV-2 lineage quantities (numeric).
#'  
#' @param threshold Minimum numeric quantity required to identify lineages (numeric).
#' Must be greater than 0. Predominant lineage groups are identified via cumulative
#' sum of observed lineages with reset after a threshold is achieved. 
#' 
#' @param n Number of predominant lineages which capture all observed sub-lineage frequencies.
#' Top n lineage groups are identified via bottom up pruning.
#' 
#' @returns A character vector of SARS-CoV-2 lineages 
#' 
#' @examples
#' # Identify predominant lineages which, when pooled with their sub-lineages, meet a count of 10
#' lineage_graph <- read_lineages()
#' lineages <- c("KP.3", "KP.3.1", "KP.3.1.1", "KP.3.1.4", "XEF" )
#' values = c( 1, 1, 5, 5, 10) 
#' 
#' detect_groups(lineage_graph, lineages, values, threshold = 10)
#' 
#' #Identify the top 2 lineages which capture all observed sub-lineages
#' detect_groups(lineage_graph, lineages, values, n = 2)
#' 
#' 
#' @export
detect_groups <- function(graph_obj, lineages, values, threshold = NULL, n = NULL ){
  
  if (!igraph::is.igraph(graph_obj)) {
    stop("Must supply an igraph object for 'graph_obj'")
  }
  
  if (!is.character(lineages)) {
    stop("Must supply a character vector for 'lineages'")
  }
  
  if ( any(duplicated(lineages)) ) {
    stop("Lineage entries must be unique")
  }
  
  if (!is.numeric(values)) {
    stop("Must supply a numeric vector for 'values'")
  }
  
  if (length(lineages) != length(values)) {
    stop("Must supply equal length vectors for 'lineages' and 'values'")
  }
  
  # Check if both arguments are provided
  if ( (is.null(threshold) && is.null(n)) || (!is.null(threshold) && !is.null(n)) ) {
    stop("Must provide 'threshold' or 'n' arguments.")
  }

  if (!is.null(threshold)) {
    if ( threshold <= 0) {
      stop("Must supply a threshold > 0")
    }
  }
  
  
  if (!is.null(n)) {
    if ( n <= 0) {
      stop("Must supply a n > 0")
    }
  }
  
  # Find a simplified root/MRCA 
  new_root <- .mrca(graph_obj, lineages)
  
  # Re-root the tree 
  graph_subset <- .subset_graph(graph_obj, new_root, lineages)
  
  # Assign counts
  graph_annotated <- .annotate(graph_subset, lineages, values, default = 0)
  
  # make separate graphs for proper count
  graph_split <- .split_recombinants(graph_annotated)
  
  
  if ( !is.null(threshold) ) {
    
    # Propagate values up graphs with cumsum, reset count as threshold is met
    graph_cummulative <- .threshold_cumsum(graph_split, threshold)
    
    # Identify groups which meet threshold 
    identified <- V(graph_cummulative)$name[ V(graph_cummulative)$value >= threshold ]
    
  }
  
  
  if ( !is.null(n) ) {
    
    # Iterativly prune tree until n amount of lineages remain  
    identified <- .topn(graph_split, n)
    
  }

  return(identified)
}
