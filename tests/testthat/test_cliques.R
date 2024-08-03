context("Clique creation")

mzfile <- system.file("standards.mzXML", package = "cliqueMS")
library(xcms)
mzraw <- MSnbase::readMSData(files = mzfile, mode = "onDisk")
cpw <- CentWaveParam(ppm = 15, peakwidth = c(5,20), snthresh = 10)
mzData <- findChromPeaks(object = mzraw, param = cpw)


test_that("model likelihood grows after obtaining cliques", {
    out <- capture.output(getCliques(mzData))
    likelihood.start <- out[grep("Beggining", out, fixed = TRUE)]
    likelihood.start <- as.numeric(
        gsub("[[:blank:]]","",
             gsub("[[:alpha:]]", "", likelihood.start)))
    likelihood.end <- out[grep("Finishing", out, fixed = TRUE)]
    likelihood.end <- as.numeric(
        gsub("[[:blank:]]","",
             gsub("[[:alpha:]]", "", likelihood.end)))
    expect_true(likelihood.end > likelihood.start)
})

cliques <- getCliques(mzData)

test_that("two peaks are in separate cliques ", {
    expect_false(getPeaklistanClique(cliques)$cliqueGroup[1] ==
                 getPeaklistanClique(cliques)$cliqueGroup[50])
})

test_that("two peaks are in the same clique ", {
    expect_true(getPeaklistanClique(cliques)$cliqueGroup[38] ==
                 getPeaklistanClique(cliques)$cliqueGroup[65])
})
