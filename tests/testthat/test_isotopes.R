context("Find isotopes")

data(ex.cliqueGroups)
isotopes <- getIsotopes(ex.cliqueGroups)

test_that("Filter isotopes is working ", {
    expect_equal(length(unique(getIsolistanClique(isotopes)[,1])),
                 length(getIsolistanClique(isotopes)[,1]))
})

test_that("This features are correct isotopes", {
    iso0 = getIsolistanClique(isotopes)[
        getIsolistanClique(isotopes)$feature == 122,]
    iso1 = getIsolistanClique(isotopes)[
        getIsolistanClique(isotopes)$feature == 115,]
    expect_identical(iso0$cluster, iso1$cluster)
    expect_equal(iso0$grade, 0)
    expect_equal(iso1$grade, 1)
})
