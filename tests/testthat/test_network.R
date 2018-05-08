context("Network creation")

mzfile <- system.file("standards.mzXML", package = "cliqueMS")
msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
                        ppm = 15, peakwidth = c(5,20), snthresh = 10)
netlist = createNetwork(msSet, msSet@peaks, filter = TRUE)

test_that("Network is of class igraph", {
    expect_identical(class(netlist$network), "igraph")
})

test_that("Network has not any edge with weigth zero", {
    expect_equal(sum(igraph::E(netlist$network)$weight == 0), 0)
})
