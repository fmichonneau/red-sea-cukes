source("R/packages.R")
source("R/genFasta.R")
source("R/findGroups.R")

source("R/load_cukeDB.R")
source("R/load_cukeAlg.R")
source("R/load_cukeDist.R")
source("R/load_cukeTree.R")
#source("R/load_taxonomyDf.R")
source("R/load_clusterGrps.R")

source("R/pairwise-groups-functions.R")
source("R/test-allopatry-functions.R")

phylobase.options(allow.duplicated.labels="ok")

load_labelsFromTaxa <- function(taxa="all") {
    taxonomyDf <- load_taxonomyDf()
    taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
    cukeDB <- load_cukeDB()

    if (taxa == "all") {
        lbls <- cukeDB[, "Labels_withAmb"]
    } else if (taxonomyDf[taxonomyDf$taxa == taxa, "rank"] == "Order") {
        lbls <- cukeDB[cukeDB$order == taxa, "Labels_withAmb"]
    } else if (taxonomyDf[taxonomyDf$taxa == taxa, "rank"] == "Family") {
        lbls <- cukeDB[cukeDB$family == taxa, "Labels_withAmb"]
    } else {
        stop("something is wrong with ", taxa)
    }
    lbls
}

load_thresholdPairwise <- function() {
    c(seq(1, 5, by=.5), 6:8)/100
}

load_thresholdClusters <- function() {
    load_thresholdPairwise()/2
}

load_tree_raxml <- function(overwrite=FALSE) {
    fnm <- "data/cukeTree-raxml.rds"
    origTreeNm <- "data/raxml/RAxML_bipartitions.cukeBarcodes"
    if (file.exists(fnm) && !overwrite) {
        tree <- readRDS(file=fnm)
    }
    else {
        if (file.exists(origTreeNm) && !overwrite) {
            tree <- ape::read.tree(file=origTreeNm)
            saveRDS(tree, fnm)
        }
        else {
            build_raxml_tree()
            tree <- load_tree_raxml(FALSE)
        }
    }
    tree
}

load_tree_phylo4 <- function(distance="raw", taxa="all") {
    distance <- match.arg(distance, c("raw", "K80"))
    if(identical(distance, "raw")) {
        tree <- load_cukeTree_raw_phylo4()
    }
    else {
        tree <- load_cukeTree_k2p_phylo4()
    }

    if (taxa == "all") {
        tree
    }
    else {
        toKeep <- load_labelsFromTaxa(taxa)
        tree <- subset(tree, tips.include=toKeep)
        stopifnot(all(toKeep %in% tipLabels(tree)) ||
                  all(!is.na(toKeep)))
        tree
    }
}

load_tree_raxml_phylo4 <- function(taxa="all", overwrite=FALSE) {
    fnm <- "data/cukeTree-raxml-phylo4.rds"
    if (file.exists(fnm) && !overwrite) {
        tr <- readRDS(fnm)
    }
    else {
        tr <- load_tree_raxml()
        ## TODO -- double check that using 1 blindly doesn't cause
        ## any issues
        ## TODO -- use mid-point rooting, or rooting on longest branch
        tr <- ape::root(tr, 1, resolve.root=TRUE)
        tr <- as(tr, "phylo4")
        bs <- nodeLabels(tr)
        bs <- data.frame(bs, stringsAsFactors=FALSE)
        bs$bs <- as.numeric(bs$bs)
        tr <- phylo4d(tr, node.data=bs)
        tr <- removeNodeLabels(tr)
        saveRDS(tr, file=fnm)
    }

    if (taxa == "all") {
        invisible(tr)
    }
    else {
        toKeep <- load_labelsFromTaxa(taxa)
        tr <- subset(tr, tips.include=toKeep)
        stopifnot(all(toKeep %in% tipLabels(tr)) ||
                  all(!is.na(toKeep)))
        invisible(tr)
    }
}


load_species_pairwiseGrps <- function(distance="raw", taxa="all",
                                      threshold=0.03) {

    taxonomyDf <- load_taxonomyDf()
    thres <- load_thresholdPairwise()

    taxa <- match.arg(as.character(taxa), taxonomyDf$taxa)
    distance <- match.arg(distance, c("raw", "K80"))
    stopifnot(length(threshold) == 1 && threshold %in% thres)

    pairwiseGrpRes <- load_pairwiseGrpRes()
    nmRes <- paste(distance, taxa, sep="-")
    res <- pairwiseGrpRes[[match(nmRes, names(pairwiseGrpRes))]][[which(thres == threshold)]]
    stopifnot(! is.null(res))
    res
}

load_tree_pairwiseGrps <- function(distance="raw", taxa="all",
                                   threshold=0.03) {
    spp <- load_species_pairwiseGrps(distance=distance, taxa=taxa,
                                     threshold=threshold)
    lSpp <- sapply(spp, length)
    tmpGrps <- mapply(rep, 1:length(lSpp), lSpp)
    tmpGrps <- data.frame(Groups=unlist(tmpGrps))
    rownames(tmpGrps) <- unlist(spp)

    tree <- load_tree_phylo4(distance=distance, taxa=taxa)
    addData(tree, tip.data=tmpGrps)
}



load_manESU <- function(taxa="Holothuriidae") {
    manESU <- read.csv(file="data/raw/manualESUs.csv", stringsAsFactors=FALSE)
    manESU$ESU_noGeo <- gsub("_[A-Z]{2}$", "", manESU$ESU_genetic)
    manESU <- manESU[-grep("\\d+amb$", manESU$Labels), ]
    manESU
}

load_species_manualGrps <- function(taxa="Holothuriidae") {
    manESU <- load_manESU()
    split(manESU$Labels, manESU$ESU_noGeo)
}

load_tree_manualGrps <- function(taxa="Holothuriidae", overwrite=FALSE) {
    fnm <- "data/cukeTree-manualESUs.rds"
    ##distance <- match.arg(distance, c("raw", "K80"))
    taxa <- match.arg(taxa) ## only Holothuriidae for now
    if (file.exists(fnm) && !overwrite) {
        manESU <- readRDS(file=fnm)
    }
    else {
        ## holTree <- load_tree_clusterGrps("raw", "Holothuriidae", threshold=0.02)
        ## write.csv(tdata(holTree, "tip"), file="data/raw_manualESUs.csv")
        ## ESU coding
        ## - for ESU_genetic: name of the species, (genus_species), followed
        ##   by alphanum to distinguish ESUs, (_1, _1a, _2, _ESU1), followed
        ##   by geography as needed (_IO, _PO):
        ##   - If same complex but different ESUs (e.g., the ind from a region
        ##     form a rec. monophyletic clade): hol_ver_1a_IO and hol_ver_1b_PO
        ##   - If same complex and same ESUs: hol_ver_1_IO and hol_ver_1_PO
        ##     (IO or PO not rec. monophyletic, or sing. ind. w/ low div.)
        ##   - If different ESUs but not sisters: hol_ver_1 and hol_ver_2
        manESU <- load_manESU()
        tree <- load_tree_raxml_phylo4(taxa)
        tDat <- data.frame(as.numeric(factor(manESU$ESU_noGeo)))
        names(tDat) <- "Groups"
        rownames(tDat) <- manESU$Labels
        manESU <- addData(tree, tip.data=tDat)
        saveRDS(manESU, file=fnm)
    }
    invisible(manESU)
}

load_localGap <- function(taxa="Holothuriidae", overwrite=FALSE) {

    fnm <- "data/localGap-manualESUs.rds"
    taxa <- match.arg(taxa) # only Holothuriidae for now
    if (file.exists(fnm) && !overwrite) {
        localGap <- readRDS(file=fnm)
    } else {
        summaryInterDist <- function(listSpecies, cukeAlg) {
            lapply(listSpecies, function(x) {
                lapply(listSpecies, function(y) {
                    interESUDist(x, y, cukeDistRaw)
                })
            })
        }
        manESU <- load_manESU()
        cukeDistRaw <- load_cukeDist_raw()
        cukeDB <- load_cukeDB()

        noGeoGrps <- load_species_manualGrps()

        summInterDist <- summaryInterDist(noGeoGrps, load_cukeAlg())

        minInter <- mapply(function(summ, nm) {
            minDist <- sapply(summ, function(x) x$min)
            minDist <- minDist[-match(nm, names(minDist))]
            minDist[which.min(minDist)]
        }, summInterDist, names(summInterDist))

        intraDist <- lapply(noGeoGrps, intraESUDist, cukeDistRaw)
        maxIntra <- sapply(intraDist, function(x) x$max)

        esuSpatial <- spatialFromSpecies(noGeoGrps, cukeDB)
        nmEsuSpatial <- names(esuSpatial[[1]])

        rgType <- lapply(names(minInter), function(x) {
            spp <- unlist(strsplit(x, "\\."))
            i <- grep(paste0("^", spp[1], "-"), nmEsuSpatial)
            j <- grep(paste0("^", spp[2], "-"), nmEsuSpatial)
            rangeType(i, j, esuSpatial[[2]])
        })

        localGap <- data.frame(minInter, maxIntra, row.names=names(minInter))

        localGap$rangeType <- sapply(rgType, function(x) x$rangeType)
        localGap$rangeType[is.na(localGap$rangeType)] <- "unknown"
        localGap <- localGap[complete.cases(localGap), ]

        localGap$species <- NA
        localGap$species[localGap$maxIntra > localGap$minInter] <-
            rownames(localGap)[localGap$maxIntra > localGap$minInter]

        saveRDS(localGap, file=fnm)
    }
    localGap
}
