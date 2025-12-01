# Testing lineage collection offline

# Unit test fetching
test_that(".fetch_lineages loads and parses YAML correctly", {
  
  # Offline lineages
  test_yaml <- list(
    list(name = "A", alias = "A", parent = NULL),
    list(name = "B", alias = "B", parent = "A"),
    list(name = "C", alias = "C", parent = "A")
  )
  
  # Create mocks for each step in collection
  mock_get <- mockery::mock("test_http")
  
  mock_content  <- mockery::mock("
  - name: A
  alias: A
  parent:
- name: B
  alias: B
  parent: A
- name: C
  alias: C
  parent: A
  ")
  
  mock_yaml <- mockery::mock(test_yaml)

  
  # Insert mocks
  mockery::stub(.fetch_lineages, "httr::GET", mock_get)
  mockery::stub(.fetch_lineages, "httr::content", mock_content)
  mockery::stub(.fetch_lineages, "yaml::yaml.load", mock_yaml)
  
  # Chaining get, content, and load work
  expect_equal(.fetch_lineages("test_url"), test_yaml)
  
})


# Unit test conversion
test_that(".yaml_to_dt converts YAML to parent–child table", {
  
  test_yaml <- list(
    list(name = "A", alias = "A", parent = NULL),
    list(name = "B", alias = "B", parent = "A"),
    list(name = "X", alias = "X", recombinant_parents = "A,A*,B,B*"),
    list(name = "X", alias = "X", parent = "X") )
  
  dt <- .yaml_to_dt(test_yaml)
  
  test_dt <- data.table::data.table(
    parent = c("root","A","A", "B"),
    child  = c("A","B","X", "X") )
  
  expect_equal( dt[order(parent, child)], test_dt[order(parent, child)] )
  expect_false( any(grepl("\\*$", dt$parent)) )
  expect_false( any( dt$parent == dt$child ) )
  expect_true( all( complete.cases(dt) ) )
})


# Integration test on full read
test_that("read_lineages returns correct graph for recombinants", {
  
  
  # Mock yaml  
  mock_fetch <- mockery::mock(
    list(
      list(name = "A", alias = "A", parent = NULL),
      list(name = "B", alias = "B", parent = "A"),
      list(name = "X", alias = "X", recombinant_parents = "A,A*,B,B*"),
      list(name = "X", alias = "X", parent = "X") )
  )
  
  mockery::stub(read_lineages, ".fetch_lineages", mock_fetch)
  
  test_dt <- data.frame(
    from = c("root","A","A", "B"),
    to  = c("A","B","X", "X") )
  
  
  test_g <- read_lineages()
  
  expect_s3_class(test_g, "igraph")
  expect_equal( igraph::as_data_frame(test_g), test_dt )
  expect_equal( igraph::vcount(test_g), 4 )
  expect_equal( igraph::ecount(test_g), 4 )
  expect_equal( min(igraph::degree(test_g, mode = "out")), 0 )
  expect_equal( max(igraph::degree(test_g, mode = "out")), 2 )
  expect_true( igraph::is.named(test_g) )

})


