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

make_im_file <- function(seqs, pop, output_dir="data/im_files",
                         alg_file, cukeDB) {

    alg <- ape::read.dna(file=alg_file, format="sequential", as.character=TRUE,
                         as.matrix=TRUE)

    if (! file.exists(output_dir)) {
        message("creating output_dir", output_dir)
        dir.create(output_dir)
    }

    stopifnot(length(seqs) == length(pop))

    ## filter out species for which only 1 population is represented and species without enough individuals
    ## role of gap only sites?

    for (i in seq_along(seqs)) {

        if (length(seqs[[i]]) < 2) {
            message("not enough sequences for ", i)
            next
        }

        pop_names <- na.omit(unique(pop[[i]]))
        if (length(pop_names) > 2) {
            stop("there should only be 2 populations")
        } else if (length(pop_names) == 1) {
            message("not enough populations for ", i)
            next
        }

        pop_sizes <- as.numeric(table(pop[[i]])[pop_names])
        loc_length <- dim(alg)[2]

        ## order sequences according to their respective populations
        alg_names <- gsub("^(RS_)", "QUERY___\\1", dimnames(alg)[[1]])
        sub_seq <- alg[match(names(seqs[[i]]), alg_names),, drop=FALSE ]
        sub_seq <- sub_seq[c(which(pop[[i]] == pop_names[1]),   # to reorder the sequences
                             which(pop[[i]] == pop_names[2])),, drop=FALSE ]

        ## extract species name for 1 line of file and file name
        spp_name <- strsplit(gsub("^QUERY___", "", dimnames(sub_seq)[[1]]), "_")
        spp_name <- sapply(spp_name, function(x) {
            paste(x[2:3], collapse="_")
        })
        spp_name <- names(table(spp_name)[which.max(table(spp_name))])
        spp_name <- gsub("~", "", spp_name)
        spp_name <- gsub("_+$", "", spp_name)

        ## converts sequence names to something short
        rs_names <- grep("^RS_", dimnames(sub_seq)[[1]])
        match_names <- match(dimnames(sub_seq)[[1]], cukeDB$Labels)
        new_names <- cukeDB[match_names, "Sample"]
        new_names[rs_names] <- gsub("(.+)_+(UF_?[0-9]+)(_.+)?$", "\\2",
                                    dimnames(sub_seq)[[1]][rs_names])
        new_names[is.na(new_names)] <- gsub("(.+)_+(UF_?[^_]+)(_.+)?$", "\\2",
                                            dimnames(sub_seq)[[1]][is.na(new_names)])

        if (any(is.na(new_names))) {
            stop("don't know what to do about:",
                 paste(dimnames(sub_seq)[[1]][is.na(new_names)], collapse=", "))
        }

        if (length(new_names) == length(pop[[i]][!is.na(pop[[i]])])) {
            new_names <- paste0(new_names, abbreviate(pop[[i]][!is.na(pop[[i]])],
                                                      minlength=1))
        } else {
            stop("i = ", i, ". problem with NAs")
        }

        ## clean up the sequence names
        new_names <- gsub("UF_?", "", new_names)
        new_names <- gsub("^(.+)-", "", new_names)
        new_names <- gsub("_", "", new_names)

        add_spc <- sapply(new_names, function(x) {
            nbSpc <- 10 - nchar(x)
            if (nbSpc < 0) {
                message("sequence name too long for ", x, ". Using first 10 char.")
                x <- substr(x, 1, 10)
                nbSpc <- 0
            }
            paste(rep.int(" ", nbSpc), collapse="")
        })

        sub_seq <- apply(sub_seq, 1, function(x) paste0(x, collapse=""))
        sub_seq <- paste0(new_names, add_spc, sub_seq)

        output <- file.path(output_dir, paste0(spp_name, ".im"))

        ## line 1: arbitrary text -- getting species names
        cat(spp_name, "\n", file=output)

        ## line 2: population names
        cat(paste(pop_names, collapse=" "), "\n", file=output, append=TRUE)

        ## line 3: an integer representing the number of loci
        ## here only mtDNA, so 1
        cat(1, "\n", file=output, append=TRUE)

        ## line 4: information about the locus
        ## name - number of ind for pop1 - number of ind for pop2
        ## length of sequence - letter for model of molecular evolution (H for HKY)
        ## inheritence scalar (1 for nuc, .75 for X linked, 0.25 for Y or mtDNA)
        ## (mutation rate, apparently not needed for msbayes)
        loc_name <- "COI"
        cat(loc_name, pop_sizes, loc_length, "H", 0.25, "\n",
            file=output, append=TRUE)

        ## line 5: sequence data
             cat(sub_seq, sep="\n", file=output, append=TRUE)
    }
    TRUE
}
