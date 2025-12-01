# Decedents tests
test_that("decedents returns all non recombinant from graph", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("root","A","A", "A", "B", "X", "X"), 
               to  = c("A","B", "C", "X", "X", "X.1", "X.2") ),
                                directed = TRUE)
  
  test_children <- decedents(graph_obj = test_graph, lineage_of_interest = "A",
            remove_recombinants = TRUE )
  
  
  expect_equal( class(test_children), "character")
  expect_false( any(is.na(match(test_children, igraph::V(test_graph)$name ))) )
  expect_false( any(grepl("X", test_children)) )
  expect_equal( sort(test_children), c("B", "C") )
  

})


test_that("decedents returns all descending lineages from graph", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("root","A","A", "A", "B", "X", "X"), 
               to  = c("A","B", "C", "X", "X", "X.1", "X.2") ),
    directed = TRUE)
  
  test_children <- decedents(graph_obj = test_graph, lineage_of_interest = "A",
                             remove_recombinants = FALSE )
  
  
  expect_equal( class(test_children), "character")
  expect_false( any(is.na(match(test_children, igraph::V(test_graph)$name ))) )
  expect_true( any(grepl("X", test_children)) )
  expect_equal( sort(test_children), c("B", "C", "X", "X.1", "X.2") )
  
})

test_that("decedents returns NA for lineages not in tree", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("root","A","A", "A", "B", "X", "X"), 
               to  = c("A","B", "C", "X", "X", "X.1", "X.2") ),
    directed = TRUE)
  
  test_children <- decedents(graph_obj = test_graph, lineage_of_interest = "J",
                             remove_recombinants = FALSE )
  
  
  expect_true( any(is.na(test_children)) )
  expect_equal( class(test_children), "character")
  
})

# Parent tests

test_that("parents returns both parents for recombinants from graph", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("root","A","A", "A", "B", "X", "X"), 
               to  = c("A","B", "C", "X", "X", "X.1", "X.2") ),
    directed = TRUE)
  
  test_parents <- parents(graph_obj = test_graph, lineage_of_interest = "X" )
  
  
  expect_equal( class(test_parents), "character")
  expect_false( any(is.na(match(test_parents, igraph::V(test_graph)$name ))) )
  expect_equal( sort(test_parents), c("A", "B") )
  
})

test_that("parents returns single parent from evolving lineage from graph", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("root","A","A", "B", "B", "X", "X"), 
               to  = c("A","B", "C", "X", "X", "X.1", "X.2") ),
    directed = TRUE)
  
  test_parents <- parents(graph_obj = test_graph, lineage_of_interest = "C" )
  
  
  expect_equal( class(test_parents), "character")
  expect_false( any(is.na(match(test_parents, igraph::V(test_graph)$name ))) )
  expect_equal( test_parents, "A" )
  
})


test_that("parents returns NA for lineages not in tree", {
  
  test_graph <- igraph::graph_from_data_frame( 
    data.frame(from = c("root","A","A", "B", "B", "X", "X"), 
               to  = c("A","B", "C", "X", "X", "X.1", "X.2") ),
    directed = TRUE)
  
  test_parents <- parents(graph_obj = test_graph, lineage_of_interest = "J" )
  
  
  expect_true( any(is.na(test_parents)) )
  expect_equal( class(test_parents), "character")
  
})

