# Internal tests

# Unit test distance based grouping
test_that(".closest_parent identifies appropriate group", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "D", "C", "C"), 
               to  = c("B", "C", "D", "E", "F", "G") ),
    directed = TRUE)
  
  # No present
  expect_true( is.na( .closest_parent("B", "G", test_graph)) )
  
  # Nested group
  expect_equal( .closest_parent( c("A", "D"), "E", test_graph), "D" )
  
  # Lineage equal to group
  expect_equal( .closest_parent( c("A", "D"), "D", test_graph), "D" )
  
  # Vectorized
  expect_equal( .closest_parent( c("A", "D"), c("D", "F"), test_graph), c("D", "A") )
})



# Unit test recombinant assignment
test_that(".assign_rcbnt identifies recombinants outside of parental groups", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "C", "X", "X.1", "G", "B" ), 
               to  = c("B", "C", "X", "X", "X.1", "G", "XX", "XX") ),
    directed = TRUE)
  
  # No recombinants
  expect_equal( .assign_rcbnt(test_graph, "A", "A", "B"), "A" )
  
  # Lineage is recombinant
  expect_equal( .assign_rcbnt(test_graph, "A", "A", "X.1"), "Recombinants" )
  
  # Parent is recombinant
  expect_equal( .assign_rcbnt(test_graph, "X", "X", "X.1"), "X" )
  
  # Assign recombinant to nested recombinant
  expect_equal( .assign_rcbnt(test_graph, "X", "X", "XX"), "Recombinants" )

})


test_that("categorize reacts to parents/children not in graph", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "D", "C", "C"), 
               to  = c("B", "C", "D", "E", "F", "G") ),
    directed = TRUE)
  
  expect_error( categorize(test_graph,"X","A"), regexp = "lineages of interest must be in tree" )
  expect_error( categorize(test_graph,"C","X"), regexp = "parents must be in tree" )
  
})


test_that("categorize defualt catches all lineages in the graph", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "B", "D", "C", "C"), 
               to  = c("B", "C", "D", "E", "F", "G") ),
    directed = TRUE)
  
  expect_true( all(categorize(test_graph, c("C", "A"), "B", "Other") == "Other") )
  
})


test_that("categorize groups appropriatly", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("A", "A", "A", "B", "C", "X", "X.1", "G", "J" ), 
               to  = c("B", "C", "J",  "X", "X", "X.1", "G", "XX", "XX") ),
    directed = TRUE)
  
  # No recombinant involvement
  expect_true( all(categorize(test_graph, c("C", "B", "J"), "A") == "A") )
  
  # Recombinant nesting
  expect_equal( categorize(test_graph, c("X", "G", "XX"), "X"), c("X", "X", "Recombinants") )
  
  # All categories observed
  expect_equal( categorize(test_graph, c("C", "G", "J", "XX"), c("C", "X"), "Other") ,
                c("C", "X", "Other", "Recombinants") )
  
  
})
