load_echinoDB <- function(file) {
    echinoDB <- read.csv(file,
                         stringsAsFactors=FALSE)
    invisible(echinoDB)
}

load_cukeDB_noLabels <- function(echinoDB=load_echinoDB()) {

    cukeDB <- subset(echinoDB, class_ == "Holothuroidea")
    cukeDB <- subset(cukeDB, pass.seq != "GenBank")
    cukeDB <- subset(cukeDB, pass.seq != "fix")
    cukeDB <- subset(cukeDB, pass.seq != "no_seq_yet")
    cukeDB <- subset(cukeDB, pass.seq != "no")
    cukeDB <- subset(cukeDB, Notes != "MH sequence")
    cukeDB <- subset(cukeDB, pass.seq != "duplicate")

    ## Taxonomic check
    testGenera <- as.matrix(xtabs(~ genusorhigher + family, data=cukeDB,
                                  subset=family != "Uncertain"))
    resGenera <- apply(testGenera, 1, function(x) sum(x != 0))
    stopifnot(all(resGenera == 1))

    testFamily <- as.matrix(xtabs(~ family + order, data=cukeDB,
                                  subset=family != "Uncertain"))
    resFamily <- apply(testFamily, 1, function(x) sum(x != 0))
    stopifnot(all(resFamily == 1))

    ## only non-ambiguous bp
    lSeq <- sapply(cukeDB$Sequence, function(x)
        length(gregexpr("[actgACTG]", x)[[1]]))

    ## all bp
    lAmb <- sapply(cukeDB$Sequence, function(x)
        length(gregexpr("[^-]", x)[[1]]))

    ## this also takes care of empty sequences (only -)
    cukeDB <- cukeDB[lAmb >= 500, ]

    ## check for duplicated samples
    dup <- cukeDB[duplicated(cukeDB$Sample), "Sample"]
    stopifnot(length(dup) == 0)

    ## These 3 sequences are not represented by other representative

    ##  it might be worth trying to figure out if we can clean up the
    ##  sequences to deal with the issues
    ##  - FRM-194
    ##  - NMV F112128
    ##  - NIWA 38032

    ## fix GPS coordinates
    cukeDB[cukeDB$decimalLatitude==19.95 & cukeDB$Loc == "MexicoPac",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(20.3, -105.5)
    cukeDB[cukeDB$Loc == "Tanzania", "decimalLatitude"] <- -cukeDB[cukeDB$Loc == "Tanzania", "decimalLatitude"]
    cukeDB[cukeDB$Loc == "Eparses", "decimalLatitude"] <- -cukeDB[cukeDB$Loc == "Eparses", "decimalLatitude"]
    cukeDB[cukeDB$Sample == "MOLAF_0139",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-9.536, 147.289)
    cukeDB[cukeDB$Sample == "8928",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
    cukeDB[cukeDB$Sample == "FRM069",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(NA, NA) ## prob cont.
    cukeDB[cukeDB$Sample == "8919",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.125, 39.191)
    cukeDB[cukeDB$Sample == "RUMF-ZE-00072",  ## not from okinawa but Xmas Island
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-10.501823, 105.685488)
    cukeDB[cukeDB$Sample == "8858F",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(NA, NA) ## prob cont.
    cukeDB[cukeDB$Sample == "9166",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-22.33917, 40.3388) ## prob cont.
    cukeDB[cukeDB$Sample == "MOLAF_0108",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-9.536, 147.289)
    cukeDB[cukeDB$Sample == "8931",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
    cukeDB[cukeDB$Sample == "8932",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.14623, 39.13786)
    cukeDB[cukeDB$Sample == "9190",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-22.34657, 40.33203)
    cukeDB[cukeDB$Sample == "8937",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
    cukeDB[cukeDB$Sample == "Hickman_needed3",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-0.41, -91.48)
    cukeDB[cukeDB$Sample == "8923",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.146230, 39.13786)
    cukeDB[cukeDB$Sample == "8933",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
    cukeDB[cukeDB$Sample == "8930",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.146230, 39.13786)
    cukeDB[cukeDB$Sample == "8938",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-6.38733, 39.28727)
    cukeDB[cukeDB$Sample == "6355",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-21.1008, 55.2437)
    cukeDB[cukeDB$Sample == "6322",
           c("decimalLatitude", "decimalLongitude")] <- data.frame(-21.0, 55.2437)

    ## check taxonomy
    testGenera <-  as.matrix(xtabs(~ genusorhigher + family, data=cukeDB, subset=family != "Uncertain"))
    resGenera <- apply(testGenera, 1, function(x) sum(x != 0))
    stopifnot(all(resGenera == 1))
    testFamily <- as.matrix(xtabs(~ family + order, data=cukeDB, subset=family != "Uncertain"))
    resFamily <- apply(testFamily, 1, function(x) sum(x != 0))
    stopifnot(all(resFamily == 1))

    invisible(cukeDB)
}

load_cukeDB <- function(cukeDB=load_cukeDB_noLabels()) {

    ## add labels
    dataLbls <- character(nrow(cukeDB))
    for (i in 1:nrow(cukeDB)) {
        dataLbls[i] <- genLabel(cukeDB[i, ])
    }

    ## using the same pattern as in seqManagement::cleanSeqLabels
    cukeDB$Labels <- gsub(":|,|\\(|\\)|;|\\[|\\]|\\'|\\s|\t", "", dataLbls)

    ## This should be something separate and called "check_labels"
    ## or something like that

    ## treeH <- load_cukeTree_k2p_phylo4(overwrite)
    ## treeTips <- data.frame(Labels_withAmb = tipLabels(treeH), Labels =
    ## gsub("_\\d+amb$", "", tipLabels(treeH)),
    ## stringsAsFactors=FALSE)

    ## stopifnot(all(treeTips$treeLabels %in% cukeDB$Labels))
    ## cukeDB_lbls <- merge(cukeDB, treeTips, by="Labels")

    invisible(cukeDB)
}
