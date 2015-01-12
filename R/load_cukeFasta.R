

load_cukeFasta <- function(overwrite=FALSE) {
    fnm <- "data/cukeAlg-cleaned.fas"
    if(file.exists(fnm) && !overwrite) {
        cukeFasta <- readRDS("data/cukeAlg-cleaned.fas.rds")
    } else {

        tmpDir <- tempdir()
        unaligned <- "cukeAlg-unaligned.fas"
        aligned <- "cukeAlg-aligned.fas"
        cleaned <- "cukeAlg-cleaned.fas"

        algFile <- file.path("data", aligned)
        cleanFile <- file.path("data", cleaned)

        ## Generate FASTA file
        if (genFasta(load_cukeDB_noLabels(overwrite), out=file.path(tmpDir, unaligned))) {
            message("Unaligned fasta file generated from CSV")
        }

        mafftOut <- file.path("tmp", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),
                                                  "-mafft.out"))

        mafftCmd <- paste("mafft --auto --op 10 --thread -1",
                          file.path(tmpDir, unaligned), ">", file.path(tmpDir, aligned),
                          "2>", mafftOut)
        message(mafftCmd)
        message("mafft output is written to ", maffOut)
        system(mafftCmd)

        if (file.copy(file.path(tmpDir, aligned),
                      algFile, overwrite=TRUE)) {
            message("Alignment generated")
        }

        ## Identify sequences with internal gaps
        seqHolC <- ape::read.dna(file=algFile, format="fasta", as.character=TRUE)
        seqHolC <- apply(seqHolC, 1, function(x) paste(x, sep="", collapse=""))
        intGap <- sapply(seqHolC, function(x)
            gregexpr("[actgnrmsykw]-+[actgnrmsykw]", x)[[1]][1] != -1)
        seqWithIntGap <- names(intGap[intGap])

        ## Identify sequences with stop codons
        seqHol <- read.dna(file=algFile, format="fasta")
        tranE <- foreach (i = 1:nrow(seqHol)) %dopar% {
            seqinr::translate(as.character(seqHol[i, ]), frame=1, numcode=9)
        }
        seqWithStop <- dimnames(seqHol)[[1]][grep("\\*", tranE)]

        ## Remove sequences with internal gaps and stop codons
        toRm <- union(seqWithStop, seqWithIntGap)
        ## dimnames(seqHol)[[1]][match(toRm, dimnames(seqHol)[[1]])] <- paste("stop-intgap", toRm, sep="_")
        toRmInd <- match(toRm, dimnames(seqHol)[[1]])
        seqHol <- seqHol[-toRmInd, ]

        ## Write working copy of fasta file
        ape::write.dna(seqHol, file=cleanFile, format="fasta", colsep="")
        cukeFasta <- seqHol
    }
    cukeFasta
}
