combine_trees <- function(list_files,
                          path = "../beast-analyses",  pattern = "_run[0-9]",
                          tree_ext = "\\.trees", log_combiner = "-b 40",
                          logcombiner_path = "~/Software/BEASTv2.1.3/bin/./logcombiner",
                          output_ext = "_combined.trees",  quiet = FALSE) {

    if (missing(list_files)) {
        lst_tr <- list.files(path = path,  pattern = paste0(pattern, tree_ext),
                             recursive = TRUE, full.names = TRUE)
    } else {
        stopifnot(all(file.exists(list_files)))
        lst_tr <- list_files
    }

    runs <- basename(lst_tr) %>% gsub(paste0(pattern, tree_ext), "", .)
    uniq_runs <- unique(runs)

    res <- character(length(uniq_runs))
    for (i in seq_along(uniq_runs)) {
        ind_runs <- lst_tr[uniq_runs[i] == runs]
        output <- gsub(paste0(pattern, tree_ext), output_ext, ind_runs[1])
        cmd <- paste(logcombiner_path, paste("-log", ind_runs, collapse = " "),
                     log_combiner, "-o", output)
        if (!quiet) message(cmd)
        res[i] <- output
        system(cmd)
    }
    if (!all(file.exists(res))) {
        warning(res[!file.exists(res)], "not created.")
    }
    res[file.exists(res)]
}

summarize_trees <- function(list_files,
                            path = "../beast-analyses", tree_ext = "_combined.trees",
                            treeannotator_path = "~/Software/BEASTv2.1.3/bin/./treeannotator",
                            treeannotator =  "-burnin 0 -heights ca",
                            output_ext = ".tree.nex", quiet = FALSE) {
    if (missing(list_files)) {
        lst_tr <- list.files(path = path, pattern = tree_ext, recursive = TRUE, full.names = TRUE)
    } else {
        stopifnot(all(file.exists(list_files)))
        lst_tr <- list_files
    }

    res <- character(length(lst_tr))
    for (i in seq_along(lst_tr)) {
        output <- gsub(tree_ext, output_ext, lst_tr[i])
        cmd <- paste(treeannotator_path, treeannotator, lst_tr[i], output)
        res[i] <- output
        if (!quiet) message(cmd)
        system(cmd)
    }
    if (!all(file.exists(res))) {
        warning(res[!file.exists(res)], " not created.")
    }
    res[file.exists(res)]
}

move_trees <- function(list_files, dest = "data/beast_trees") {
    stopifnot(all(file.exists(list_files)))
    stopifnot(file.exists(dest))
    new_files <- file.path(dest, basename(list_files))
    fnm <- cbind(list_files, new_files)
    apply(fnm, 1, function(x) file.copy(x[1], x[2], overwrite = TRUE))
}

combined_trees <- combine_trees()
summarized_trees <- summarize_trees(list_files = combined_trees)
move_trees(summarized_trees)
