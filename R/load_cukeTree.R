removeNodeLabels <- function(phy) {
    intNd <- nodeId(phy, "internal")
    is.na(phy@label[intNd]) <- TRUE
    phy
}

build_cukeTree <- function(alg, model, Nrep) {
    model <- match.arg(model, c("raw", "K80"))
    if (identical(model, "raw")) {
        dMat <- load_cukeDist_raw()
    }
    else if (identical(model, "K80")) {
        dMat <- load_cukeDist_k2p()
    } else stop("Houston, we have a problem")

    treH <- ape::nj(dMat)
    bootH <- ape::boot.phylo(treH, alg, function(xx) {
        ape::nj(ape::dist.dna(xx, model=model))
    }, B=Nrep)
    treH$node.label <- bootH
    treH
}

build_cukeTree_phylo4 <- function(tree) {
    tree$edge.length[tree$edge.length < 0] <- 1e-6
    ## TODO -- double check that using 1 blindly doesn't cause
    ## any issues
    ## TODO -- use mid-point rooting, or rooting on longest branch
    treeRooted <- ape::root(tree, 1, resolve.root=TRUE)
    treeP4 <- as(treeRooted, "phylo4")
    bs <- nodeLabels(treeP4)
    is.na(bs[as.character(rootNode(treeP4))]) <- TRUE
    bs <- data.frame(bs, stringsAsFactors=FALSE)
    bs$bs <- as.numeric(bs$bs)
    treeP4 <- phylo4d(treeP4, node.data=bs)
    treeP4 <- removeNodeLabels(treeP4)
    treeP4
}


load_cukeTree_k2p <- function(overwrite=FALSE, Nrep, ...) {
    fnm <- "data/cukeTree-k2p.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree <- readRDS(file=fnm)
        cukeTree$tip.label <- gsub("\\\"", "", cukeTree$tip.label)
    }
    else {
        if(missing(Nrep)) Nrep <- 200
        cukeAlg <- load_cukeAlg()
        cukeTree <- build_cukeTree(cukeAlg, model="K80", Nrep=Nrep, ...)
        saveRDS(cukeTree, file=fnm)
    }
    invisible(cukeTree)
}

load_cukeTree_k2p_phylo4 <- function(overwrite=FALSE) {
    fnm <- "data/cukeTree-k2p-phylo4.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeTree4 <- readRDS(fnm)
    }
    else {
        tmpTr <- load_cukeTree_k2p()
        cukeTree4 <- build_cukeTree_phylo4(tmpTr)
        saveRDS(cukeTree4, file=fnm)
    }
    invisible(cukeTree4)
}
