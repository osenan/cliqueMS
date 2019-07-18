context("Find annotation")

isotopes <- getIsotopes(ex.cliqueGroups)
adducts <- getAnnotation(isotopes, positive.adinfo, "positive")

test_that("This features are correct adducts", {
    expect_identical(adducts$peaklist$an1[69], "[M+H-NH3]+")
    expect_identical(adducts$peaklist$an1[73], "[M+NH4]+")
    expect_identical(adducts$peaklist$an2[93], "[M+H]+")
})
