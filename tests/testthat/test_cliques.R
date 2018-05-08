context("Clique creation")

mzfile <- system.file("standards.mzXML", package = "cliqueMS")
msSet <- xcms::xcmsSet(files = mzfile, method = "centWave",
                        ppm = 15, peakwidth = c(5,20), snthresh = 10)

test_that("model likelihood grows after obtaining cliques", {
    out <- capture.output(getCliques(msSet))
    likelihood.start <- out[grep("Beggining", out, fixed = T)]
    likelihood.start <- as.numeric(
        gsub("[[:blank:]]","",
             gsub("[[:alpha:]]", "", likelihood.start)))
    likelihood.end <- out[grep("Finishing", out, fixed = T)]
    likelihood.end <- as.numeric(
        gsub("[[:blank:]]","",
             gsub("[[:alpha:]]", "", likelihood.end)))
    expect_true(likelihood.end > likelihood.start)
})

cliques <- getCliques(msSet)

test_that("two peaks are in separate cliques ", {
    expect_false(cliques$peaklist$cliqueGroup[1] ==
                 cliques$peaklist$cliqueGroup[50])
})

test_that("two peaks are in the same clique ", {
    expect_true(cliques$peaklist$cliqueGroup[205] ==
                 cliques$peaklist$cliqueGroup[197])
})
