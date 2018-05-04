context("Find annotation")

cliques <- getCliques(exmsSet, filter = TRUE)
isotopes <- getIsotopes(cliques)
adducts <- getAnnotation(isotopes, positive.adinfo, "positive")

test_that("This features are correct adducts", {
    expect_identical(adducts$peaklist$an1[132], "[M+H]+")
    expect_identical(adducts$peaklist$an1[143], "[M+Na]+")
    expect_identical(adducts$peaklist$an1[241], "[M+K]+")
})
