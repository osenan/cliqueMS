context("Find annotation")
set.seed(2)

data(ex.cliqueGroups)
data(positive.adinfo)
print(positive.adinfo)
print(positive.adinfo$log10freq+2.2)
isotopes <- getIsotopes(ex.cliqueGroups)
adducts <- getAnnotation(isotopes, positive.adinfo, "positive")
CP069 <- getPeaklistanClique(adducts)["CP069",]
anCP069 <- which(c(CP069$an1,CP069$an2) != "")
CP073 <- getPeaklistanClique(adducts)["CP073",]
anCP073 <- which(c(CP073$an1,CP073$an2) == "[M+NH4]+")
       
test_that("This features are correct adducts", {
    expect_identical(getPeaklistanClique(adducts)["CP069","score1"],
        getPeaklistanClique(adducts)["CP069","score2"])
    expect_identical(getPeaklistanClique(adducts)["CP069",
        paste("an",anCP069, sep = "")], "[M+H-NH3]+")
    expect_identical(getPeaklistanClique(adducts)["CP073","score1"],
        getPeaklistanClique(adducts)["CP073","score2"])
    expect_identical(getPeaklistanClique(adducts)["CP073",
       paste("an",anCP073, sep = "")], "[M+NH4]+")
    expect_identical(getPeaklistanClique(adducts)["CP093","an2"], "[M+H]+")
})
