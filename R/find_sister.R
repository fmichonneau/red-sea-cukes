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
            check_nds <- phylobase::descendants(tr, chi_node[chi_node != cur_node])
            if (length(check_nds) > 1 &&
                length(grep("^QUERY", names(check_nds))) == 0) {
                break
            } else {
                cur_node <- anc_node
            }
        }
        res[[i]] <- phylobase::descendants(tr, anc_node)
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

    uniqStr <- gsub("_t_", "_", names(sis_seq))
    mtch <- match(gsub("_[0-9]+amb$", "", names(sis_seq)), cukeDB$Labels)
    latlong <- cukeDB[mtch, c("decimalLatitude", "decimalLongitude")]
    ## Red Sea + Gulf of Aden polygon corners:
    ## lat: 30 -- long: 30
    ## lat: 3  -- long: 62
    isara <- ifelse(latlong$decimalLatitude > 3 & latlong$decimalLatitude < 30 &
                        latlong$decimalLongitude > 30 & latlong$decimalLongitude < 100,
                    "arabia", "outside")
    isara[grep("^QUERY", names(sis_seq))] <- "arabia"
    isara
}

read_alg_file <- function(alg_file) {
    ## Get complete alignment used to build tree
    alg <- ape::read.dna(file=alg_file, format="sequential", as.character=TRUE,
                         as.matrix=TRUE)
    alg
}

get_pop_names <- function(pop) {
    sort(unique(pop))
}

check_sequences <- function(seqs) {

    to_keep <- sapply(seqs, function(seq) {

        tmp_seqs <- seq[["sequences"]]
        tmp_pop <- seq[["populations"]]
        pop_names <- get_pop_names(tmp_pop)

        ## At least 2 sequences for each species
        if (length(tmp_seqs) < 2) {
            message("not enough sequences for ", dimnames(tmp_seqs)[[1]][1])
            return(FALSE)
        }

        ## At least 2 populations represented
        if (length(pop_names) > 2) {
            stop("there should only be 2 populations")
        } else if (length(pop_names) == 1) {
            message("not enough populations for ", dimnames(tmp_seqs)[[1]][1])
            return(FALSE)
        }

        if (any(is.na(tmp_pop)))
            stop("There are NAs in populations")

        if (nrow(tmp_seqs) != length(tmp_pop))
            stop("not the same number of sequences and populations")

        return(TRUE)
    })
    to_keep
}

extract_sequences <- function(seqs, alg, cukeDB, to_ignore_file, to_add_file) {

    to_ignore <- scan(file=to_ignore_file, what="character")
    to_ignore <- paste0(to_ignore, "_?")
    to_ignore <- paste0(to_ignore, collapse="|")
    to_ignore <- sapply(seqs, function(x) length(grep(to_ignore, names(x))) > 0)
    seqs <- seqs[! to_ignore]

    to_add <- scan(file=to_add_file, what="character", sep="\n")
    to_add <- strsplit(to_add, " ")
    to_add <- lapply(to_add, function(x) paste0(x, "_?"))
    pos_add <- lapply(to_add, function(x) sapply(x, function(y) {
        grep(y, dimnames(alg)[[1]])
    }))

    alg_names <- gsub("^(RS_)", "QUERY___\\1", dimnames(alg)[[1]])
    pos_seq <- lapply(seqs, function(x) match(names(x), alg_names))
    pos_seq <- c(pos_seq, pos_add)
    pos_seq <- lapply(pos_seq, function(x) {
        names(x) <- alg_names[x]
        x
    })

    res_seq <- lapply(pos_seq, function(x) {
        tmp_seq <- alg[x, ]
        tmp_pop <- is_arabia(x, cukeDB=cukeDB)

        ## Remove sequences for which the population can't be identified
        if (any(is.na(tmp_pop))) {
            tmp_seq <- tmp_seq[which(!is.na(tmp_pop)), ]
            tmp_pop <- tmp_pop[which(!is.na(tmp_pop))]
        }

        list(sequences=tmp_seq,
             populations=tmp_pop)
    })

    to_keep <- check_sequences(res_seq)
    res_seq <- res_seq[to_keep]
    res_seq
}

get_family <- function(seqs) {
    lbls <- dimnames(seqs[["sequences"]])[[1]]
    to_ignore <- grep("^RS_", lbls)
    lbls <- lbls[-to_ignore]

    fams <- strsplit(lbls, "_")
    fams <- unique(sapply(fams, function(x) x[1]))
    if (length(fams) != 1)
        stop("More than 1 family for", lbls[1])
    else
        fams
}

make_nexus_files <- function(seqs, cukeDB, output_dir="data/nexus_files") {

    ## Create directory if non-existing
    if (! file.exists(output_dir)) {
        message("creating output_dir", output_dir)
        dir.create(output_dir)
    }

    fams <- sapply(seqs, get_family)

    uniq_fam <- unique(fams)

    out_files <- character(length(uniq_fam))

    for (i in seq_along(uniq_fam)) {
        idx_fam <- which(fams == uniq_fam[i])
        tmp_seq <- lapply(seqs[idx_fam], function(x) x[["sequences"]])
        tmp_pop <- lapply(seqs[idx_fam], function(x) x[["populations"]])
        tmp_seq <- mapply(function(x, y, z) {
            dimnames(x)[[1]] <- paste(dimnames(x)[[1]], y,
                                      paste0("spp", z), sep="_")
            x
        }, tmp_seq, tmp_pop, seq_along(tmp_seq))
        tmp_alg <- do.call("rbind", tmp_seq)
        out_files[i] <- file.path(output_dir, paste0(uniq_fam[i], ".nex"))
        ape::write.nexus.data(tmp_alg, file=out_files[i], interleaved=FALSE)
    }
    out_files
}

make_im_file <- function(seqs, output_dir="data/im_files", cukeDB) {

    ## Create directory if non-existing
    if (! file.exists(output_dir)) {
        message("creating output_dir", output_dir)
        dir.create(output_dir)
    }

    list_files <- character(length(seqs))
    tree_labels <- vector("list", length(seqs))

    for (i in seq_along(seqs)) {

        sub_seq <- seqs[[i]][["sequences"]]
        tree_labels[[i]] <- dimnames(sub_seq)[[1]]

        sub_pop <- seqs[[i]][["populations"]]
        pop_names <- get_pop_names(sub_pop)
        pop_sizes <- as.numeric(table(sub_pop)[pop_names])

        ## reorder the sequences
        sub_seq <- sub_seq[c(which(sub_pop == pop_names[1]), # to reorder the sequences
                             which(sub_pop == pop_names[2])), ]
        sub_pop <- sub_pop[c(which(sub_pop == pop_names[1]), # to reorder the sequences
                             which(sub_pop == pop_names[2]))]
        loc_length <- dim(sub_seq)[2]

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
        match_names <- match(dimnames(sub_seq)[[1]], cukeDB[["Labels"]])
        new_names <- cukeDB[match_names, "Sample"]
        new_names[rs_names] <- gsub("(.+)_+(UF_?[0-9]+)(_.+)?$", "\\2",
                                    dimnames(sub_seq)[[1]][rs_names])
        new_names[is.na(new_names)] <- gsub("(.+)_+(UF_?[^_]+)(_.+)?$", "\\2",
                                            dimnames(sub_seq)[[1]][is.na(new_names)])

        if (any(is.na(new_names))) {
            stop("don't know what to do about:",
                 paste(dimnames(sub_seq)[[1]][is.na(new_names)], collapse=", "))
        }

        new_names <- paste0(new_names, abbreviate(sub_pop, minlength=1))

        ## clean up the sequence names
        new_names <- gsub("UF_?", "", new_names)
        new_names <- gsub("^(.+)-", "", new_names)
        new_names <- gsub("_", "", new_names)

        add_spc <- sapply(new_names, function(x) {
            nbSpc <- 10 - nchar(x)
            if (nbSpc < 0) {
                message("sequence name too long for ", x, ". Using first 10 char.")
                nm <- substr(x, 1, 9)
                x <- paste0(nm, substr(x, nchar(x), nchar(x)))
                nbSpc <- 0
            }
            paste0(x, paste(rep.int(" ", nbSpc), collapse=""))
        })

        sub_seq <- apply(sub_seq, 1, function(x) paste0(x, collapse=""))
        sub_seq <- paste0(add_spc, sub_seq)

        ## adds the i^th in file name to make unique
        output <- file.path(output_dir, paste0(spp_name, "_", i, ".im"))
        list_files[i] <- output

        ## line 1: arbitrary text -- getting species names
        cat(spp_name, "\n", file=output)

        ## line 2: population names
        cat(paste(get_pop_names(sub_pop), collapse=" "), "\n", file=output, append=TRUE)

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
    res <- list_files[nzchar(list_files)]
    attr(res, "tree_labels") <- tree_labels
    res
}

get_npop <- function(file) {
    fcon <- readLines(file)
    fcon <- fcon[-(1:4)]
    fcon <- sapply(fcon,  function(x) substr(x, 1, 10))
    fcon <- unname(gsub("\\s+$", "",  fcon))
    table(sapply(fcon,  function(x) substr(x, nchar(x), nchar(x))))
}

make_infile <- function(im_files, dest, min_n_ind) {
    if (! missing(min_n_ind)) {
        too_few <- lapply(im_files, function(x) {
                              np <- get_npop(x)
                              all(np >= min_n_ind)
                          })
        too_few <- unlist(too_few)
        im_files <- im_files[too_few]
    }
    cat(basename(im_files), sep="\n", file=dest)
}

check_tree_labels <- function(im_files, tree_file, dest) {
    tree <- rncl::read_newick_phylo(file=tree_file)
    tree_labels <- attr(im_files, "tree_labels")
    tree_labels <- lapply(tree_labels, function(x) {
        gsub("^(RS_)", "QUERY___\\1", x)
    })


    pdf(file=dest, width=15)

    for (i in seq_along(tree_labels)) {
        if (length(tree_labels[[i]]) > 1) {
            to_rm <- tree$tip.label[! tree$tip.label %in% tree_labels[[i]]]
            if (length(to_rm) == length(tree$tip.label)) browser()
            sub_tree <- ape::drop.tip(tree, to_rm)
            plot(sub_tree)
            add.scale.bar(x=0, y=length(sub_tree$tip.label))
        } else next
    }
    dev.off()
    TRUE
}

run_convertIM <- function(infile, outfile, wd="data/im_files") {
    outfile <- basename(outfile)
    infile <- basename(infile)
    cmd <- paste("cd", wd, ";", "convertIM.pl", "-o", outfile, infile)
    message(cmd)
    system(cmd)
}

make_obsSS <- function(batch_file,  output="obsSS.txt", wd="data/im_files") {
    system(paste("cd", wd, ";",
                 "obsSumStats.pl -T obsSS.table", batch_file,
                 "> ", output))
}
