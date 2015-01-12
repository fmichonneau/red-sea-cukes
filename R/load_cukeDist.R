load_cukeDist_raw <- function(overwrite=FALSE, ...) {
    fnm <- "data/cukeDist-raw.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeDist <- readRDS(file=fnm)
    } else {
        cukeDist <- ape::dist.dna(load_cukeAlg(), as.matrix=TRUE, model="raw",
                                  pairwise.deletion=FALSE)
        saveRDS(cukeDist, file=fnm)
    }
    invisible(cukeDist)
}

load_cukeDist_k2p <- function(overwrite=FALSE, ...) {
    fnm <- "data/cukeDist-k2p.rds"
    if (file.exists(fnm) && !overwrite) {
        cukeDist <- readRDS(file=fnm)
    } else {
        cukeDist <- ape::dist.dna(load_cukeAlg(), as.matrix=TRUE, model="K80",
                                  pairwise.deletion=FALSE)
        saveRDS(cukeDist, file=fnm)
    }
    invisible(cukeDist)
}
