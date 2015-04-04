read_beast_trees <- function() {
    taxa <- c("holothuriidae/Holothuriidae.tree.nex",
              "chiridotidae/Chiridotidae.tree.nex",
              "sclerodactylidae/Sclerodactylidae.tree.nex",
              "stichopodidae/Stichopodidae.tree.nex",
              "synaptidae/Synaptidae.tree.nex")
    tr_files <- file.path("..", "beast-analyses", taxa)

    trees <- lapply(tr_files, read.beast)
    names(trees) <- gsub("\\..+$", "", basename(taxa))
    trees
}

id_nodes <- function(tree) {
    lbl <- tree$tip.label
    spp <- strsplit(lbl, "_")
    spp <- sapply(spp,  function(x) x[length(x)])
    spt <- split(lbl,  spp)
    tr4 <- as(tree, "phylo4")
    nds <- sapply(spt,  function(x) {
                      MRCA(tr4, x)
                  })
    rec_mono <- sapply(nds, function(x) {
                           ch <- phylobase::children(tr4, x)
                           dd_ch <- lapply(ch, function(x) { names(phylobase::descendants(tr4, x)) })
                           le_ar <- sapply(dd_ch, function(x) grep("_arabia_", x))
                           le_ot <- sapply(dd_ch, function(x) grep("_outside_", x))
                           wh_ar <- which(sapply(le_ar, length) > 0)
                           wh_ot <- which(sapply(le_ot, length) > 0)
                           if (length(wh_ar) != 1 || length(wh_ot) != 1) {
                               res <- FALSE
                           } else {
                               if (length(dd_ch[[wh_ot]]) == length(le_ot[[wh_ot]]) &&
                                   length(dd_ch[[wh_ar]]) == length(le_ar[[wh_ar]])) {
                                   res <- TRUE
                               } else {
                                   res <- FALSE
                               }
                           }
                           res
                       })
    list(nodes = nds, rec_mono = rec_mono)
}

summary_nodes <- function(trees) {

    all_nds <- lapply(trees, id_nodes)
    tt <- mapply(function(tr, tr_nm, full_nds) {
                     nds_orig <- full_nds$nodes
                     recmono <- full_nds$rec_mono
                     n_desc <- lapply(nds_orig,  function(nd) {
                                          tr4 <- as(tr,  "phylo4")
                                          dd <- phylobase::descendants(tr4, nd, "tip")
                                          ot <- length(grep("_outside_", names(dd)))
                                          ar <- length(grep("_arabia_",  names(dd)))
                                          c("outside" = ot, "arabia" = ar)
                                      })
                     n_ot <- sapply(n_desc,  function(x) x[["outside"]])
                     n_ar <- sapply(n_desc,  function(x) x[["arabia"]])
                     nds <- nds_orig - length(tr$tip.label)
                     h_hpd_min <- tr$'height_95%_HPD_MIN'[nds]
                     h_hpd_max <- tr$'height_95%_HPD_MAX'[nds]
                     h_hpd_med <- tr$'height_median'[nds]
                     cbind(family = rep(tr_nm, length(nds)),
                           rec_mono = recmono,
                           node = nds_orig,
                           h_hpd_min, h_hpd_max, h_hpd_med, n_ot, n_ar)
                 }, trees, names(trees), all_nds,  SIMPLIFY = FALSE)
    tt <- do.call("rbind",  tt)
    tt <- data.frame(tt, row.names = NULL, stringsAsFactors = FALSE) %>%
        arrange(desc(h_hpd_med))
    tt
}

if (FALSE) {

tt <- fetch("summary_beast_trees")

tt %<>% mutate(n_ot = as.numeric(n_ot),
               n_ar = as.numeric(n_ar)) %>%
    mutate(node = as.numeric(node)) %>%
    mutate(h_hpd_min = as.numeric(h_hpd_min)) %>%
    mutate(h_hpd_max = as.numeric(h_hpd_max)) %>%
    mutate(h_hpd_med = as.numeric(h_hpd_med)) %>%
    mutate(order_nds = order(h_hpd_med, decreasing = TRUE)) %>%
    filter(!is.na(h_hpd_med) &  n_ot >= 3 & n_ar >= 3)

pdf(file = "node_ages.pdf")

ggplot(tt, aes(x = factor(order_nds), y = h_hpd_med, colour = family)) +
  geom_pointrange(aes(ymin = h_hpd_min, ymax = h_hpd_max)) +
  xlab(element_blank()) + ylab("Divergence time (My)")

dev.off()

}

## pdf(file = "test_plot.pdf")#,  height = 70, width = 12)
## for (i in seq_along(taxa)) {
##     plot.phylo(trees[[i]], show.tip.label = FALSE)
##     nodelabels(text = rep(NA, length(all_nds[[i]])), node = all_nds[[i]],
##                frame = "circ", col = "red", bg = "red", cex = 10)
##     HPDbars(trees[[i]], nodes = all_nds[[i]], lwd = 2, col = "red")
##     axisPhylo()
## }
## dev.off()
