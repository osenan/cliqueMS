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
