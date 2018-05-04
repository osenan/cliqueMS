context("Network creation")

netlist = createNetwork(exmsSet, exmsSet@peaks, filter = TRUE)

test_that("Network is of class igraph", {
    expect_identical(class(netlist$network), "igraph")
})

test_that("Network has not any edge with weigth zero", {
    expect_equal(sum(igraph::E(netlist$network)$weight == 0), 0)
})
