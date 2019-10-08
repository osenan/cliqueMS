context("Find annotation")
set.seed(2)

data(ex.cliqueGroups)
data(positive.adinfo)
print(positive.adinfo)
print(positive.adinfo$log10freq+2.2)
isotopes <- getIsotopes(ex.cliqueGroups)
adducts <- getAnnotation(isotopes, positive.adinfo, "positive")
CP069 <- getPeaklistanClique(adducts)["75.1",]
anCP069 <- which(c(CP069$an1,CP069$an2) != "")
CP073 <- getPeaklistanClique(adducts)["75.3",]
anCP073 <- which(c(CP073$an1,CP073$an2) == "[M+NH4]+")
       
test_that("This features are correct adducts", {
    expect_identical(getPeaklistanClique(adducts)["75.1","score1"],
        getPeaklistanClique(adducts)["75.1","score2"])
    expect_identical(getPeaklistanClique(adducts)["75.1",
        paste("an",anCP069, sep = "")], "[M+H-NH3]+")
    expect_identical(getPeaklistanClique(adducts)["75.3","score1"],
        getPeaklistanClique(adducts)["75.3","score2"])
    expect_identical(getPeaklistanClique(adducts)["75.3",
       paste("an",anCP073, sep = "")], "[M+NH4]+")
    expect_identical(getPeaklistanClique(adducts)["91.1","an2"], "[M+H]+")
})
