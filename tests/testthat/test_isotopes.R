context("Find isotopes")

cliques <- getCliques(exmsSet)
isotopes <- getIsotopes(cliques)

test_that("Filter isotopes is working ", {
    expect_equal(length(unique(isotopes$isotopes[,1])),
                 length(isotopes$isotopes[,1]))
})

test_that("This features are correct isotopes", {
    iso0 = isotopes$isotopes[isotopes$isotopes$feature == 205,]
    iso1 = isotopes$isotopes[isotopes$isotopes$feature == 197,]
    expect_identical(iso0$cluster, iso1$cluster)
    expect_equal(iso0$grade, 0)
    expect_equal(iso1$grade, 1)
})
