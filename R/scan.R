#' Lists All Sub-lineages of a Given Lineage
#' 
#' `decedents()` returns all SARS-CoV-2 lineages descending the lineage of
#'  interest. 
#' 
#' @param graph_obj an igraph object constructed from `read_lineages()`
#'     depicting all SARS-CoV-2 lineage parent/child relationships.
#' 
#' @param lineage_of_interest  A character string depicting a SARS-CoV-2 lineage. 
#' 
#' @param remove_recombinants  Logical value to specify removal of recombination 
#' 
#' @returns A character vector of SARS-CoV-2 lineages 
#' 
#' @examples
#' #An example of pulling all SARS-CoV-2 sub-lineages of KP.3.1.1, while removing
#' #all recombinants and their decedents.
#' lineages <- read_lineages()
#' decedents(graph_obj = lineages, lineage_of_interest = "KP.3.1.1", remove_recombinants = TRUE)
#'
#' #An example of pulling all sub-lineages of XFG.
#' decedents(graph_obj = lineages, lineage_of_interest = "XFG")
#' @export
decedents <- function(graph_obj = NULL,
                 lineage_of_interest = NULL,
                 remove_recombinants = FALSE) {
  
  if (!igraph::is.igraph(graph_obj)) {
    stop("Must supply an igraph object for 'graph_obj'")
  }
  
  if (!is.character(lineage_of_interest)) {
    stop("Must supply a character vector for 'lineage_of_interest'")
  }
  
  
  # Provide NA if lineage is not found in tree 
  if ( is.na(match(lineage_of_interest,  igraph::V(graph_obj)$name )) ) {
    return( NA_character_ )
  }
  

  # Find all descendants
  children <- igraph::as_ids( igraph::subcomponent(graph_obj,
                                                   lineage_of_interest,
                                                   mode = "out") ) 
  
  
  
  if ( remove_recombinants ) {
    
    # Identify and remove nested recombinants under lineages of interest
    nested_recombinants <- names(which(igraph::degree(graph_obj,
                                                       v = setdiff(children,
                                                                   lineage_of_interest),
                                                       mode = "in") == 2) )
    
    # Expand all nested recombinants
    nested_recombinants <- unique(
      unlist(lapply(nested_recombinants,
                    function(name) igraph::as_ids(igraph::subcomponent(graph_obj,
                                                                       name,
                                                                       mode = "out"))
                    )
             )
      )
    
    # Remove nested recombinants
    return( setdiff( children, c(nested_recombinants, lineage_of_interest) ) )
    
    
  } else {
    
    return( setdiff( children, lineage_of_interest) )
    
  }

}


#' Lists All Direct Parents to a Given Lineage
#' 
#' `parents()` returns all direct SARS-CoV-2 parental lineages to the lineage of
#'  interest. 
#' 
#' @param graph_obj an igraph object constructed from `read_lineages()`
#'     depicting all SARS-CoV-2 lineage parent/child relationships.
#' 
#' @param lineage_of_interest  A character string depicting a SARS-CoV-2 lineage. 
#' 
#' @returns A character vector of SARS-CoV-2 lineages 
#' 
#' @examples
#' #An example of pulling all direct parent(s) of KP.3.1.1.
#' lineages <- read_lineages()
#' parents(graph_obj = lineages, lineage_of_interest = "KP.3.1.1")
#'
#' #An example of pulling all direct parents of XFG.
#' parents(graph_obj = lineages, lineage_of_interest = "XFG")
#' 
#' @export
parents <- function(graph_obj = NULL, lineage_of_interest = NULL) {
  
  if (!igraph::is.igraph(graph_obj)) {
    stop("Must supply an igraph object for 'graph_obj'")
  }
  
  if (!is.character(lineage_of_interest)) {
    stop("Must supply a character vector for 'lineage_of_interest'")
  }
  
  # Provide NA if lineage is not found in tree 
  if ( is.na(match(lineage_of_interest, V(graph_obj)$name)) ) {
    return( NA_character_ )
  }
  
  
  parents <- unlist(igraph::as_ids( igraph::neighbors(graph_obj,
                                                     lineage_of_interest,
                                                     mode = "in")))
  
  return( parents )
}
