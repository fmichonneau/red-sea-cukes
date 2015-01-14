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

is_arabia <- function(sis_seq, cukeDB) {

    res <- lapply(sis_seq, function(x) {
        uniqStr <- gsub("_t_", "_", names(x))
        mtch <- match(gsub("_[0-9]+amb$", "", names(x)), cukeDB$Labels)
        latlong <- cukeDB[mtch, c("decimalLatitude", "decimalLongitude")]
        ## Red Sea + Gulf of Aden polygon corners:
        ## lat: 30 -- long: 30
        ## lat: 3  -- long: 62
        isara <- ifelse(latlong$decimalLatitude > 3 & latlong$decimalLatitude < 30 &
                            latlong$decimalLongitude > 30 & latlong$decimalLongitude < 100,
                        "arabia", "outside")
        isara[grep("^QUERY", names(x))] <- "arabia"
        isara
    })
    res
}
