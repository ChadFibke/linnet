#' Reads in SARS-CoV-2 Lineage Data 
#' 
#' `read_lineages()` reads in the current SARS-CoV-2 PANGO lineage data from
#'  outbreak.info and represents this information in a directed network
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

  lineage_yaml <- yaml::yaml.load( httr::content( httr::GET(lineage_url), as = "text" ) )

  # Convert YAML to a tabular lineage summary
  suppressWarnings(
  lineage_dt <- data.table::rbindlist(lapply(lineage_yaml, function(entry) {
    list( name = entry$name,
          alias = entry$alias,
          parent = entry$parent,
          recombinant_parent = ifelse( !is.null(entry$recombinant_parent),
                                       paste( unique(unlist(strsplit( entry$recombinant_parent, split = ","))),
                                              collapse = ","),
                                       NA) )
  }), fill = TRUE)
  )
  
  
  # Collect the maximum number of parents needed to split up the recombinant parents
  n_splits <- max(lengths(strsplit( lineage_dt$recombinant_parent, split = ",")))
  
  # Clean parent-child data and reformat to all parent child relations
  tree_dt <- lineage_dt[  is.na(parent) & is.na(recombinant_parent), parent := "root" ][
    is.na(parent) & !is.na(recombinant_parent), parent := recombinant_parent ][
      , paste0("parent_", 1:n_splits) := data.table::tstrsplit( parent, ",", fixed = TRUE )][
        , c("parent", "recombinant_parent") := NULL ][
          , data.table::melt( .SD, id.vars = c("name"),
                  measure.vars = patterns("^parent"),
                  value.name = "parent") ][
            !is.na(parent), .( child = name, parent) ][
              , parent := gsub(pattern = "\\*$", replacement = "",x = parent)]
  
  # Return network from parent-child relations
  return( igraph::graph_from_data_frame(d = tree_data[, .(parent, child)],
                                                 directed = TRUE)
  )
  
}