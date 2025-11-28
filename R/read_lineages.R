#' Fetch current outbreak.info yaml
#'
#' Internal function used by `read_lineages()`.
#'
#' @keywords internal
.fetch_lineages <- function(url) {
  lineage_yaml <- yaml::yaml.load( httr::content( httr::GET(url), as = "text" ) )
  
  return(lineage_yaml)
}


#' Convert YAML to a tabular lineage summary
#'
#' Internal function used by `read_lineages()`.
#'
#' @keywords internal
#' @importFrom data.table `:=`
.yaml_to_dt <- function(yaml) {
  
  # Scrape corresponding information from yaml.
  # Collect unique parents where recombination exists
  dt_lineage <- data.table::rbindlist(lapply(yaml, function(entry) {
    list( 
      name = entry$name,
      alias = entry$alias,
      parent = if( length(entry$parent) == 0 ) NA_character_ else entry$parent,
      recombinant_parent = ifelse(
        !is.null(entry$recombinant_parent),
        paste( unique(unlist(strsplit( entry$recombinant_parent, split = ","))),collapse = ","),
        NA
        )
      ) 
    }),
    fill = TRUE )

  # Collect the maximum number of parents needed to split up the recombinant parents
  n_splits <- max(lengths(strsplit( dt_lineage$recombinant_parent, split = ",")))
  
  # Clean parent-child data and reformat to all parent child relations
  dt_tree <- dt_lineage[ is.na(parent) & is.na(recombinant_parent), parent := "root" ][
    is.na(parent) & !is.na(recombinant_parent), parent := recombinant_parent ][
      , paste0("parent_", 1:n_splits) := data.table::tstrsplit( parent, ",", fixed = TRUE )][
        , c("parent", "recombinant_parent") := NULL ][
          , data.table::melt( .SD, id.vars = c("name"),
                              measure.vars = patterns("^parent"),
                              value.name = "parent") ][
                                !is.na(parent), .( parent, child = name) ][
                                  , parent := gsub(pattern = "\\*$", replacement = "",x = parent) ]
  
  return(dt_tree)
}


#' Reads and Converts SARS-CoV-2 Lineages into an igraph object 
#' 
#' `read_lineages()` downloads the current SARS-CoV-2 PANGO lineage data from
#'  outbreak.info, converts it to a tabular structure, and returns a directed
#'  network
#' 
#' @returns An igraph object.
#' 
#' @examples
#' lineage_network <- read_lineages()
#' 
#' @export
read_lineages <- function() {
  
  # Download current data
  lineage_url <- "https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/curated_reports_prep/lineages.yml"
  
  # Collect lineage YAML 
  yaml_lineage <- .fetch_lineages(lineage_url)
  
  # Reformat to dt and clean parent-child relations
  dt_tree <- .yaml_to_dt(yaml_lineage)
  
  # Convert to directed network
  network_lineage <- igraph::graph_from_data_frame(d = dt_tree, directed = TRUE)

  return( network_lineage )
  
}