
test_that(".mrca identifies appropriate lineages", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "D", "C", "C"), 
               to  = c("B", "C", "D", "E", "F", "G") ),
    directed = TRUE)
  
  # Provides character lineage in tree (singleton)
  expect_type( .mrca( test_graph, c("F")), "character" )
  expect_equal( .mrca( test_graph, c("F")), "F" )
  
  # Test self and decedent
  expect_equal( .mrca( test_graph, c("C", "F")), "C" )
  
  # All leaves
  expect_equal( .mrca( test_graph, c("E", "F", "G")), "A" )
  
})


test_that(".subset_graph reduces graph to expected lineages", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "D", "C", "C"), 
               to  = c("B", "C", "D", "E", "F", "G") ),
    directed = TRUE)
  
  # Structure
  expect_s3_class( .subset_graph( test_graph, "A", "A"), "igraph" )
  expect_true( igraph::is.directed( .subset_graph( test_graph, "A", "F") ) )
  
  # Only one MRCA exists
  expect_equal( length( which( igraph::degree(.subset_graph(test_graph, "A", c( "D", "F")),
                                      mode = "in") == 0 ) ),
                           1 )
            
  # Singleton case
  expect_equal( igraph::V(.subset_graph( test_graph, "A", "A"))$name, "A" )
  
  # Are intermediates present
  expect_true( any( igraph::V(.subset_graph( test_graph, "A", "F"))$name == "C" ) )
  
  # Proper subset
  expect_true( igraph::subgraph_isomorphic( .subset_graph( test_graph, "A", c( "D", "F")),
                                            test_graph) )
  
  
})



test_that(".annotate assigns correct value to vertices", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "D", "C", "C"), 
               to  = c("B", "C", "D", "E", "F", "G") ),
    directed = TRUE)
  
  test_lineages <- c("D", "E", "C", "F", "G")
  
  test_counts <- c( 1, 4, 1, 2, 3)
  
  graph_anno <- .annotate(test_graph, test_lineages, test_counts)
  
  # Structure
  expect_s3_class( graph_anno, "igraph" )
  expect_true( igraph::is.directed(graph_anno) )
  expect_equal( igraph::vertex_attr_names( graph_anno ), c("name", "value") )
  expect_type( igraph::vertex_attr(graph_anno, "value") , "double" )
  expect_true( all(!is.na(igraph::vertex_attr(graph_anno, "value"))) )
  
  
  # Default used
  expect_equal( igraph::V(graph_anno)$name[which(igraph::V(graph_anno)$value == 0)],
                c("A", "B") )
  
  # Expected values received
  expect_equal( igraph::V(graph_anno)$name[which(igraph::V(graph_anno)$value == 1)],
                c("D", "C") )
  
  
  
})


test_that(".split_recombinants attaches recombinant clades to root", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "C", "X", "X.1", "G", "B" ), 
               to  = c("B", "C", "X", "X", "X.1", "G", "XX", "XX") ),
    directed = TRUE)
  
  
  graph_split <- .split_recombinants( test_graph )
  
  # Structure
  expect_s3_class( graph_split, "igraph" )
  expect_true( igraph::is.directed(graph_split) )
  
  # All recombinants have root as direct parent
  expect_equal( igraph::neighbors(graph_split, c("X", "XX"), mode = "in")$name,
                "A" )
  
  # Only one MRCA exists
  expect_equal( length( which( igraph::degree(graph_split, mode = "in") == 0 ) ),
                1 )
  
  # Nested recombinant has expected children
  expect_equal( igraph::subcomponent(graph_split, "X", mode = "out")$name,
                c("X", "X.1", "G"))

})



test_that(".threshold_cumsum returns expeted lineages", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "A", "A", "X", "X.1" ), 
               to  = c("B", "C", "X", "XX", "X.1", "G")),
    vertices =  data.frame(name = c("A", "B", "C", "X", "X.1", "G", "XX"), 
                           value = c(0, 1, 1, 0, 1, 5, 2) ),
    directed = TRUE)
  
  
  graph_sum <- .threshold_cumsum( test_graph, 2 )
  
  # Structure
  expect_s3_class( graph_sum, "igraph" )
  expect_true( igraph::is.directed(graph_sum) )
  expect_equal( igraph::vertex_attr_names( graph_sum ), c("name", "value", "segment") )
  expect_type( igraph::vertex_attr(graph_sum, "value") , "double" )
  expect_true( all(!is.na(igraph::vertex_attr(graph_sum, "value"))) )
  
  # Threshold is met before restarting counts
  vertex_order <- as.integer(igraph::topo_sort(graph_sum, mode = "in"))
  values <- igraph::vertex_attr(graph_sum, "value")[vertex_order] 
  segment_id <- igraph::vertex_attr(graph_sum, "segment")[vertex_order] 
  
  expect_true( all( values[ !duplicated( segment_id, fromLast = TRUE) ] >= 2 ) )
  expect_true( all( values[ duplicated( segment_id, fromLast = TRUE) ] < 2 ) )
  
  # Total unique counts remain in graph
  expect_equal( sum( values[ !duplicated( segment_id, fromLast = TRUE) ] ) ,
                sum( igraph::vertex_attr(test_graph, "value") ) )
  
  # Expected values for chain from X, X.1, G
  expect_equal( igraph::get.vertex.attribute(graph_sum, name = "value", index = c(4,5,6)),
                c(1, 1, 5) )
})


test_that("cumsum_binning returns expeted errors", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "A", "A", "X", "X.1" ), 
               to  = c("B", "C", "X", "XX", "X.1", "G")),
    directed = TRUE)
  
    expect_error( cumsum_binning( data.frame(),"X", 1, 1),
                  regexp = "Must supply an igraph object")
    
    expect_error( cumsum_binning( test_graph, 1, 1, 1),
                  regexp = "Must supply a character vector")
    
    expect_error( cumsum_binning( test_graph, c("X", "X"), c(1,1), 1),
                  regexp = "Lineage entries must be unique")
    
    expect_error( cumsum_binning( test_graph, c("X", "Y"), c(TRUE,FALSE), 1),
                  regexp = "Must supply a numeric vector")
    
    expect_error( cumsum_binning( test_graph, c("X", "Y"), c(1), 1),
                  regexp = "Must supply equal length vectors")
    
    expect_error( cumsum_binning( test_graph, c("X", "Y"), c(1, 1), -1),
                  regexp = "Must supply a threshold > 0")
})


test_that("cumsum_binning returns expeted values", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "A", "A", "X", "X.1" ), 
               to  = c("B", "C", "X", "XX", "X.1", "G")),
    directed = TRUE)
  
  lineages <- c("A", "B", "C", "X", "X.1", "G", "XX")
  values <- c(0, 1, 1, 0, 1, 5, 2)
  
  expect_equal( cumsum_binning( test_graph, lineages, values, 2), c("A", "XX", "G") )
  expect_equal( cumsum_binning( test_graph, lineages, values, 3), c("A", "G") )
  expect_equal( cumsum_binning( test_graph, lineages, values, 6), "X.1" )
  expect_equal( cumsum_binning( test_graph, lineages, values, 10), "A" )

})
