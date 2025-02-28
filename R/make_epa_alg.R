##' Function to create the alignment to be used with the Evolutionary
##' Placement Algorithm in RAxML (places sequence in Tree). It needs a
##' tree and an alignment that contains the sequences used for the
##' tree + the sequences to be assigned.  The original alignment and
##' tree come from the barcoding chapter in dissertation. The file
##' with the sequences to be assigned were sent to me by Gustav in
##' fasta format (and I added on 2016-02-10 the latest sequences from
##' the Red Sea/Oman).
##'
##'
##' This function (1) takes both files, (2) adds prefix to sequences
##' to assign to find them more easily a posteriori, (3) align the
##' files
##' @param orig_alg the original alignment used to build the ML tree
##'     used
##' @param assign_alg the sequences to be assigned in the tree
make_epa_alg <- function(orig_alg = "data/raxml-epa-input/cukeBarcodes-flagAmb.phy",
                         assign_alg = "data/raxml-epa-input/Cukes_Arabia.fas",
                         output = "data/raxml-epa-input/epa-input.phy") {

    orig_ape <- ape::read.dna(file=orig_alg, format="sequential")
    assign_ape <- ape::read.dna(file=assign_alg, format="fasta", as.matrix = TRUE)

    dimnames(assign_ape)[[1]] <- paste0("RS_", dimnames(assign_ape)[[1]])

    tmp_unalg <- tempfile()
    tmp_alg <- tempfile()

    message("unaligned file: ", tmp_unalg)
    ape::write.dna(orig_ape, file=tmp_unalg, format="fasta", colw=10000, colsep="")
    ape::write.dna(assign_ape, file=tmp_unalg, append=TRUE, format="fasta", colw=10000, colsep="")

    message("aligned file: ", tmp_alg)
    system(paste("mafft --auto --op 10 --thread -1", tmp_unalg, ">", tmp_alg))

    tmp_ape <- ape::read.dna(tmp_alg, format="fasta")
    ape::write.dna(tmp_ape, file=output, colw=10000, colsep="")

    invisible(file.exists(output))
}

## This function is needed to remove one of the taxa (Stichopus aff
## horrens from Galapagos that is not spelled correctly).
normalize_tree <- function(alg = "data/raxml-epa-input/epa-input.phy",
                           tree = "data/raxml-epa-input/RAxML_bestTree.cukeBarcodes",
                           outtree = "data/raxml-epa-input/RAxML_bestTree.cukeBarcodes.normalized") {
    tr <- ape::read.tree(file=tree)
    alg <- ape::read.dna(file=alg, format="sequential")
    toDrop <- tr$tip.label[!tr$tip.label %in% dimnames(alg)[[1]]]
    message("Dropping: ", paste(toDrop, collapse = ", "))
    tr <- ape::drop.tip(tr, toDrop)
    ape::write.tree(tr, file = outtree)
    invisible(file.exists(outtree))
}


run_epa <- function(alg = "epa-input.phy",
                    tree = "RAxML_bestTree.cukeBarcodes.normalized",
                    model="GTRCAT", prefix="epaTest") {

    system(paste("cd data/raxml-epa-input/; /home/francois/Software/RAxML-8.1.15/./raxmlHPC-PTHREADS-SSE3 -f v -s", alg, "-t", tree,
                 "-m", model, "-n", prefix, "-T7"))

}
