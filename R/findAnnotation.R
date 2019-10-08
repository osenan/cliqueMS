isotopeAnnotation <- function(df.annotation, anclique) {
    ## Isotopes of grade > 0 are excluded from annotation
    ## This function adds the annotation, in case it exists,
    ## to all isotopes of grade > 0 as they should have the same
    ## annotation than its isotope of grade 0
    isofeatures <- anclique@isotopes$feature[
        which(anclique@isotopes$grade == 0)]
    extraAnList <- lapply(isofeatures, function(x) {
        cluster <- anclique@isotopes$cluster[
            which(anclique@isotopes$feature == x)]
        extraf <- anclique@isotopes$feature[
            which(anclique@isotopes$cluster == cluster)]
        ## now drop the isotope of grade0 because it already has annotation
        extraf <- extraf[-1*which(extraf == x)]
        an <- df.annotation[which(df.annotation$feature == x),]
        ## annotation for all isotopes should be the same
        extraan <- lapply(extraf, function(x) {
            an$feature = x
            an
        })
        extraan <- do.call(rbind, extraan)
    })
    extraAn <- do.call(rbind, extraAnList)
    return (extraAn)
}


checkadinfo <- function(adinfo, polarity) {
    sameCols = match(colnames(adinfo),
        c("adduct","log10freq","massdiff","nummol","charge"))
    if(sum(!is.na(sameCols)) < 5) {
        stop(paste("some missing colnames, check that column",
            "'adduct','log10freq','massdiff','nummol' and 'charge'",
            "are included"))
    }
    if(polarity == "positive") {
        ## Separate them in adducts with charge >1, 
        ## adducts with nmol >1 && charge ==1, 
        ## and adducts with charge and nmol = 1
        chargead = adinfo[which(adinfo$charge > 1),]
        chargead = chargead[order(chargead$massdiff),]

        nummolad = adinfo[which(adinfo$nummol > 1),]
        nummolad = nummolad[order(nummolad$massdiff),]

        normalad = adinfo[which(adinfo$charge == 1),]
        normalad = normalad[which(normalad$nummol == 1),]

        # Sort them from adducts with charge >1, 
        # normal and adducts with nummol > 1, 
        # trying to virtually sort all adducts from smaller massDiff to bigger
        returnAdinfo = rbind(chargead, normalad, nummolad)
    } else {
        ## Separate them in adducts with charge < -1, 
        ## adducts with nmol >1 && charge == -1
        ## and adducts with charge and nmol = 1
        chargead = adinfo[which(adinfo$charge < -1),]
        chargead = chargead[order(chargead$massdiff),]

        nummolad = adinfo[which(adinfo$nummol > 1),]
        nummolad = nummolad[order(nummolad$massdiff),]

        normalad = adinfo[which(adinfo$charge == -1),]
        normalad = normalad[which(normalad$nummol == 1),]

        ## Sort them from adducts with charge >1
        ## normal and adducts with nummol > 1
        ## trying to virtually sort all adducts from smaller massDiff to bigger
        returnAdinfo = rbind(chargead, normalad, nummolad)
    }
    return(returnAdinfo)
}


getIsoCharge <- function(df.clique, iso) {
    ## this function sets charge to features that are isotopes 
    ## because with isotope annotation they already have charge
    ## and therefore, adduct annotation is limited to that charge
    chargeV <- vapply(seq_len(nrow(df.clique)), function(x) {
        charge <- 0
        if( length(grep("[",df.clique[x,"isotope"],
            fixed = TRUE)) != 0 ) {
            feat <- df.clique[x,"feature"]
            charge <- iso[which(iso$feature == feat),"charge"]
        }
        charge
    }, numeric(1))
    df.ret <- cbind(df.clique, chargeV)
    colnames(df.ret) <- c("mz","isotope","feature","charge")
    return(df.ret)
}

addIsoAnnotation <- function(isotopes, annotation, anclique) {
    newannotation <- annotation
    if( sum(is.na(unlist(isotopes))) !=
        length(unlist(isotopes)) ) {
        isoAn <- isotopeAnnotation(newannotation, anclique)
        newannotation <- rbind(newannotation, isoAn)
    }
    return(newannotation)
}


transformAnnotation <- function(peaklist, annotation) {
    ## function to write the annotation into the peaklist
    ## of the anclique object, and transforming into a
    ## character
    newpeaklist <- peaklist
    newannotation <- annotation
    newannotation <- newannotation[order(newannotation$feature),]
    newannotation <- newannotation[,-1]
    newpeaklist <- cbind(newpeaklist, newannotation)
    ## transform annotation in character

    newpeaklist$an1 = as.character(newpeaklist$an1)
    newpeaklist$an1[newpeaklist$an1 == "NA"] <- NA
    newpeaklist$an2 = as.character(newpeaklist$an2)
    newpeaklist$an2[newpeaklist$an2 == "NA"] <- NA
    newpeaklist$an3 = as.character(newpeaklist$an3)
    newpeaklist$an3[newpeaklist$an3 == "NA"] <- NA
    newpeaklist$an4 = as.character(newpeaklist$an4)
    newpeaklist$an4[newpeaklist$an4 == "NA"] <- NA
    newpeaklist$an5 = as.character(newpeaklist$an5)
    newpeaklist$an5[newpeaklist$an5 == "NA"] <- NA
    ## transform masses that are 0 in na
    newpeaklist$mass1[newpeaklist$mass1 == 0] <- NA
    newpeaklist$mass2[newpeaklist$mass2 == 0] <- NA
    newpeaklist$mass3[newpeaklist$mass3 == 0] <- NA
    newpeaklist$mass4[newpeaklist$mass4 == 0] <- NA
    newpeaklist$mass5[newpeaklist$mass5 == 0] <- NA
    return(newpeaklist)
}


#' @export
#' @title Annotate adducts and fragments
#'
#' @description This function annotates adducts after
#' isotope annotation. For each clique group,
#' it searches for combinations of two or more features 
#' compatible with the same neutral mass and two or more
#' adducts in 'adinfo' list. For clique groups than have more
#' than one annotation solution, it scores all
#' possibilities and returns the top five solutions.
#' @param anclique Object of class 'anClique' with isotope
#' annotation
#' @param adinfo data.frame with columns 'adduct' with
#' adduct name, column 'log10freq' with the log10 frequency of each
#' adduct in the list, column 'massdiff' with the adduct mass
#' diff, column 'nummol' with the number of metabolite's
#' molecule necessary for that adduct and column 'charge' with the
#' charge of that adduct.
#' @param polarity Polarity of the adducts, choose between
#' 'positive' or 'negative'
#' @param topmasstotal All neutral masses in the group are 
#' ordered based on their adduct log-frequencies and their 
#' number of adducts. From that list, a number of "topmasstotal" 
#' masses are selected for the final annotation.
#' @param topmassf In addition to 'topmasstotal', for each feature
#' the list of ordered neutral masses is subsetted to the masses with
#' an adduct in that particular feature. For each sublist, a number 
#' 'topmassf' neutral masses are also selected for the final annotation.
#' @param sizeanG After neutral mass selection, if a clique group 
#' has a number of monoisotopic features bigger than 'sizeanG', 
#' the annotation group is divided into non-overlapping annotation groups.
#' Each subdivision is annotated independently.
#' @param ppm Relative error in ppm in which we consider two or 
#' more features compatible with a neutral mass and two or more
#' adducts in 'adinfo'.
#' @param filter This parameter removes redundant annotations.
#' If two neutral masses in the same annotation group 
#' have a relative mass difference smaller than 'filter' and the same
#' features and adducts, drop the neutral mass with less adducts
#' @param emptyS Score given to non annotated features. If you use your own
#' 'adinfo', do not set 'emptyS' bigger than any adduct log frequency in
#' your list.
#' @param normalizeScore If 'TRUE', the reported score is normalized and scaled.
#' Normalized score goes from 0, when it means that the raw score
#' is close to the minimum score (all features with
#' empty annotations), up to 100, which is the score value of the theoretical
#' maximum annotation (all the adducts of the list with the
#' minimum number of neutral masses).
#' @details The default 'adinfo' lists are 'positive.adinfo' and
#' 'negative.adinfo'. For use load them with 'data(positive.adinfo)' or
#' data(negative.adinfo) commands.
#' Reported scores do not always refer to the entire clique group.
#' There might be features whose annotation is independent
#' from other features of the clique group. This occurs when there are 
#' no neutral masses with adducts in both groups of features.
#' Therefore, the clique
#' group is divided in non overlapping regions, called annotation groups.
#' Scores report for these annotation groups. To compare scores between
#' different groups use 'normalizeScore' = TRUE.
#'
#' If clique groups have a lot of features, there are
#' many combinations of neutral masses and adducts. This could
#' lead to long running times to score the top annotations.
#' Parameters 'topmassf' and 'topmasstotal' are relevant in
#' those cases to drop the less likely neutral masses to
#' speed up the time of computation and still obtain
#' the most plausible annotation. If the clique group is small
#' usually no neutral masses are discarded for the scoring.
#' @return An 'anClique' object with annotation columns added
#' to the peaklist
#' @examples
#' data(ex.cliqueGroups)
#' show(ex.cliqueGroups)
#' ex.isoAn <- getIsotopes(ex.cliqueGroups)
#' show(ex.isoAn)
#' data(positive.adinfo)
#' ex.adductAn <- getAnnotation(ex.isoAn, positive.adinfo, 'positive')
getAnnotation <- function(anclique, adinfo, polarity, topmasstotal = 10, 
    topmassf = 1, sizeanG = 20, ppm = 10, filter = 1e-4, emptyS = -6,
    normalizeScore = TRUE) {
    if( (polarity != "positive")&&(polarity != "negative") ) {
        stop("Polarity has to be 'positive' or 'negative'")
    }
    orderadinfo = checkadinfo(adinfo, polarity)
    if(anclique@anFound == TRUE) {
        warning("Annotation has already been computed for this object")
    }
    if(anclique@isoFound == FALSE) {
        warning("Isotopes have not been annotated.
            This could lead to some errors in adduct annotation\n")
    }
    if(anclique@cliquesFound == FALSE) {
        warning("Cliques have not been computed.
            This could lead to long computing times
            for adduct annotation")
    }
    message("Computing annotation")
    ## Compute annotation for each clique
    ppm = 1e-6*ppm
    anList <- lapply(anclique@cliques, function(x) {
        df.clique <- as.data.frame(
            cbind(anclique@peaklist[x, c("mz","isotope")],x)
        )
        colnames(df.clique) <- c("mz","isotope","feature")
        df.clique <- df.clique[grep("^M0", df.clique$isotope),]
        ## set charge to isotopic features
        df.clique <- getIsoCharge(df.clique, anclique@isotopes)
        df.clique <- df.clique[order(df.clique$mz, decreasing = FALSE),]
        annotation <- returnAnnotation(df.clique, orderadinfo,
            topmassf, topmasstotal,
            sizeanG, ppm,
            filter, emptyS, normalizeScore)
    })
    df.annotation <- do.call(rbind, anList)
    message("Annotation computed, updating peaklist")
    df.annotation <- addIsoAnnotation(anclique@isotopes, df.annotation,
        anclique)
    anclique@peaklist <- transformAnnotation(anclique@peaklist, df.annotation)
    ## update object information
    anclique@anFound <- TRUE
    return(anclique)
}
