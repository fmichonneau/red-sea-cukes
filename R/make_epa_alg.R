##' function to create the alignemnt to be used with the Evolutionary
##' Placement Algorithm in RAxML (places sequence in Tree). It needs a
##' tree and an alignment that contains the sequences used for the
##' tree + the sequences to be assigned.  The original alignment and
##' tree come from the barcoding chapter in dissertation. The file
##' with the sequences to be assigned were sent to me by Gustav in
##' fasta format.  This function (1) takes both files, (2) adds prefix
##' to sequences to assign to find them more easily a posteriori, (3)
##' align the files
##' @param orig_alg the original alignment used to build the ML tree
##' used
##' @param assign_alg the sequences to be assigned in the tree
make_epa_alg <- function(orig_alg = "data/raxml-epa-input/cukeBarcodes-flagAmb.phy",
                         assign_alg = "data/raxml-epa-input/Cukes_Arabia.fas",
                         output = "data/raxml-epa-input/epa-input.phy") {

    orig_ape <- ape::read.dna(file=orig_alg, format="sequential")
    assign_ape <- ape::read.dna(file=assign_alg, format="fasta")

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

## changed manually "immature" in tree and "immature in alg to _immature

normalize_tree <- function(alg = "data/raxml-epa-input/epa-input.phy.reduced",
                           tree = "data/raxml-epa-input/RAxML_bestTree.cukeBarcodes") {
    tr <- ape::read.tree(file=tree)
    alg <- ape::read.dna(file=alg, format="sequential")
    toDrop <- tr$tip.label[!tr$tip.label %in% dimnames(alg)[[1]]]
    tr <- ape::drop.tip(tr, toDrop)
    tr
}

reduced_tree <- normalize_tree()
ape::write.tree(reduced_tree, file="data/raxml-epa-input/reduced_tree.phy")

run_epa <- function(alg = "epa-input.phy",
                    tree = "RAxML_bestTree.cukeBarcodes",
                    model="GTRCAT", prefix="epaTest") {

    system(paste("cd data/raxml-epa-input/; /home/francois/Software/RAxML-8.1.15/./raxmlHPC-PTHREADS-SSE3 -f v -s", alg, "-t", tree,
                 "-m", model, "-n", prefix, "-T7"))

}


find_sister <- function(tree = "data/raxml-epa-input/RAxML_labelledTree.epaTest") {

    tr <- ape::read.tree(file = tree)
    tr$edge.length <- NULL
    tr <- as(tr, "phylo4")

    qry <- grep("^QUERY", tipLabels(tr))

    res <- vector("list", length(qry))
    names(res) <- tipLabels(tr)[qry]

    for (i in seq_along(qry)) {
        cur_node <- qry[i]
        while(TRUE) {
            anc_node <- ancestor(tr, cur_node)
            chi_node <- children(tr, anc_node)
            check_nds <- descendants(tr, chi_node[chi_node != cur_node])
            if (length(check_nds) > 1 &&
                length(grep("^QUERY", names(check_nds))) == 0) {
                break
            } else {
                cur_node <- anc_node
            }
        }
        res[[i]] <- descendants(tr, anc_node)
    }
    res[!duplicated(res)]
}

make_find_sister <- function(overwrite=FALSE) {
    fnm <- "data/find_sister.rds"
    if (file.exists(fnm) && !overwrite) {
        res <- readRDS(fnm)
    } else {
        res <- find_sister()
        saveRDS(res, file=fnm)
    }
    res
}

rm_duplicate_seq <- function(sis_seq) {

    ## only the QUERY sequences should be duplicated, the other ones
    ## were checked during the creation of the alignment
    res <- lapply(sis_seq, function(x) {
        ufnb <- regmatches(names(x), gregexpr("UF_?[0-9]+", names(x)))
        ufnb <- gsub("_", "", ufnb)
        x[!duplicated(ufnb)]
    })
    res
}

is_arabia <- function(sis_seq=rm_duplicate_seq(make_find_sister())) {

    cukeDB <- load_cukeDB()
    res <- lapply(sis_seq, function(x) {
        latlong <- cukeDB[match(gsub("_[0-9]+amb$", "", names(x)),
                                cukeDB$Labels), c("decimalLatitude", "decimalLongitude")]
        ## Red Sea + Gulf of Aden polygon:
        ## lat: 30 -- long: 30
        ## lat: 3  -- long: 62
        isara <- ifelse(latlong$decimalLatitude > 3 & latlong$decimalLatitude < 30 &
                            latlong$decimalLongitude > 30 & latlong$decimalLongitude < 100,
                        "arabia", "outside")
        isara[grep("^QUERY", names(x))] <- "arabia"
        isara
    })

}
