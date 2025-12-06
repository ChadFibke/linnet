#' Identify closest parental group for a given lineage of interest
#'
#' Internal function used by `categorize()`.
#'
#' @keywords internal
.closest_parent <- function(group_indices, lineage_indices, graph_obj) {
  
  # Calculate distances from each lineage node to each parental node 
  distance_matrix <- igraph::distances(graph_obj,
                                       v = lineage_indices,
                                       to = group_indices,
                                       mode = "in")
  
  # Standardize inf values
  distance_matrix[ is.infinite(distance_matrix) ] <- Inf
  
  # Collect the closest parental group
  ancestor_indices <- max.col(-distance_matrix, ties.method = "first")
  
  # Account for those who fall in no group
  distance_matrix[ is.infinite(distance_matrix) ] <- NA
  
  ancestor_indices[ is.na(rowMeans(distance_matrix, na.rm = TRUE)) ] <- NA
  
  # Assign name of closest ancestor
  closest_parent <- colnames(distance_matrix)[ancestor_indices]
  
  
  return(closest_parent)
}


#' Assign recombination status
#'
#' Internal function used by `categorize()`.
#'
#' @keywords internal
.assign_rcbnt <- function(graph_obj, parental_group, closest_parents, lineage_of_interest) {
  
  rcbnts <- names( which( igraph::degree( graph_obj, mode = "in" ) > 1 ) )
  
  # Collect all children of recombinants
  rcbnt_direct_children <- unique(
    unlist( 
      lapply( 
        rcbnts, function(name) igraph::as_ids( igraph::subcomponent(graph_obj,
                                                                    name,
                                                                    mode = "out")
        )
      )
    )
  )
  
  # Collect parent sub clades nested within recombinants
  adjust_parents <- parental_group[ parental_group %chin% rcbnt_direct_children ]
  
  spawned <- unique(
    unlist(
      lapply(
        adjust_parents,
        function(name) igraph::as_ids( igraph::subcomponent(graph_obj,
                                                            name,
                                                            mode = "out")
        )
      )
    )
  )
  
  # Identify recombinants downstream of parent sub clades
  nested_recombinants <- names(
    which( igraph::degree(graph_obj,
                          v = setdiff(spawned, adjust_parents),
                          mode = "in") == 2) )
  
  # Expand all nested recombinants
  nested_recombinants <- unique(
    unlist(
      lapply(
        nested_recombinants,
        function(name) igraph::as_ids( igraph::subcomponent(graph_obj,
                                                            name,
                                                            mode = "out")
        )
      )
    )
  )
  
  to_remove <- setdiff( spawned, nested_recombinants )
  
  
  # Remove nested lineages of interest from recombinant
  rcbnt_direct_children <- setdiff( rcbnt_direct_children, to_remove )
  
  closest_parents[ lineage_of_interest %chin% rcbnt_direct_children ] <- "Recombinants"
  
  
  return(closest_parents)
}


#' Categorize a Lineage into Nearest Supplied Parental Lineages
#' 
#' `categorize()` returns the nearest SARS-CoV-2 lineage relative to lineage of
#'  interest. 
#'  
#' @importFrom igraph V
#' 
#' @importFrom data.table `%chin%`
#' 
#' @param graph_obj an igraph object constructed from `read_lineages()`
#'     depicting all SARS-CoV-2 lineage parent/child relationships.
#' 
#' @param lineage_of_interest  A character string depicting a SARS-CoV-2 lineage
#'  to be collapsed.
#' 
#' @param parental_group A character vector of a lineage(s) which we would like
#'  to categorize into.
#'  
#' @param default_lineage If the lineage is not grouped, what value would you 
#'  like to assign to the lineage (default: NA)?
#' 
#' @returns A character vector of SARS-CoV-2 lineages 
#' 
#' @examples
#' #An example of identifying the nearest parental lineage to "KP.3.1.1"
#' lineages <- read_lineages()
#' 
#' categorize(graph_obj = lineages, lineage_of_interest = "KP.3.1.1",
#'  parental_group =  c("KP.1", "KP.2", "KP.3") )
#'
#' @export
categorize <- function( graph_obj = NULL,
                        lineage_of_interest = NULL,
                        parental_group = NULL,
                        default_lineage = NA_character_){
  
  if ( !igraph::is.igraph(graph_obj) ) {
    stop("Must supply an igraph object for 'graph_obj'")
  }
  
  if ( !is.character(lineage_of_interest) ) {
    stop("Must supply a character vector for 'lineage_of_interest'")
  }
  
  if ( !all( lineage_of_interest %chin% V(graph_obj)$name) ) {
    stop( 
      sprintf("All lineages of interest must be in tree, the following are not: %s",
              lineage_of_interest[ !(lineage_of_interest %chin% V(graph_obj)$name) ] )
    )
  }
  
  if ( !is.character(parental_group) ) {
    stop("Must supply a character vector for 'parental_group'")
  }
  
  if ( !all( parental_group %chin% V(graph_obj)$name) ) {
    stop( 
      sprintf("All parents must be in tree, the following are not: %s",
              parental_group[ !(parental_group %chin% V(graph_obj)$name) ] )
    )
  }
  
  # Identify location of each lineage and parental group in network
  
  lineage_indices <- match( lineage_of_interest, V(graph_obj)$name )
  group_indices <- match( parental_group, V(graph_obj)$name )
  
  # Identify crude closest parents
  category <- .closest_parent(group_indices, lineage_indices, graph_obj)
  
  # Assign recombinant status to lineages which have 2 direct parents,
  # and descend this recombinant lineage. However, exclude lineages which
  # are not a parental group/or direct descendant of a recombinant parental group.
  
  category <- .assign_rcbnt(graph_obj, parental_group, category, lineage_of_interest)
  
  # Fill remaining values with default
  category[ is.na( category ) ] <- default_lineage
  
  
  return( category )
}
