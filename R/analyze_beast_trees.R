read_beast_trees <- function() {
    tr_files <- list.files(path = "data/beast_trees", pattern = "\\.tree\\.nex$",
                           full.names = TRUE)

    trees <- lapply(tr_files, phyloch::read.beast)
    names(trees) <- gsub("\\..+$", "", basename(tr_files))
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
    nd_data <- mapply(function(tr, tr_nm, full_nds) {
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
    nd_data <- do.call("rbind",  nd_data)
    nd_data <- data.frame(nd_data, row.names = NULL, stringsAsFactors = FALSE) %>%
      arrange(desc(h_hpd_med))

    sub_trees <- trees[names(trees) %in% nd_data$family]
    sub_nodes <- lapply(split(nd_data, nd_data$family),
                        function(x) x$node)
    sub_tips <- mapply(function(tr, nodes) {
                           tr <- as(tr, "phylo4")
                           res <- lapply(nodes, function(x)
                               phylobase::descendants(tr, as.numeric(x), "tips"))
                           res
                       }, sub_trees, sub_nodes)

    spp <- lapply(sub_tips,  function(y) {
                      lapply(y, function(x) {
                                 x <- gsub("_~|_sp|_aff|_cf", "", names(x))
                                 x <- gsub("truncata", "gracilis", x)
                                 x <- gsub("impatiens", "aff. impatiens", x)
                                 tt <- sapply(strsplit(x, "_"), function(x)
                                     paste0(x[2:3], collapse = " "))
                                 res <- table(tt)
                                 paste(names(res[res > 2]), collapse = "/")
                             })
                  })
    nd_data$spp <- unlist(spp)
    nd_data
}

if (FALSE) {

node_data <- fetch("summary_beast_trees") %>%
    mutate(n_ot = as.numeric(n_ot),
           n_ar = as.numeric(n_ar),
           rec_mono = as.logical(rec_mono)) %>%
    mutate(node = as.numeric(node)) %>%
    mutate(h_hpd_min = as.numeric(h_hpd_min)) %>%
    mutate(h_hpd_max = as.numeric(h_hpd_max)) %>%
    mutate(h_hpd_med = as.numeric(h_hpd_med)) %>%
    mutate(order_nds = order(h_hpd_med, decreasing = TRUE)) %>%
    filter(!is.na(h_hpd_med) & rec_mono & n_ot > 1 & n_ar > 1) %>%
    mutate(spp = paste0("(", order_nds, ") ", spp)) %>%
    mutate(Species = factor(spp, levels = as.character(spp)))

pdf(file = "tmp/node_ages.pdf", width = 10, height = 4)
ggplot(node_data, aes(x = factor(order_nds), y = h_hpd_med, colour = Species)) +
  geom_pointrange(aes(ymin = h_hpd_min, ymax = h_hpd_max)) +
  xlab("Node ID") + ylab("Divergence time (My)") +
  theme_bw()
dev.off()

rees <- fetch("get_beast_trees")


pdf(file = "tmp/family_trees.pdf", paper = "US")#,  height = 70, width = 12)
par(mar = c(2, 0, 1, 0))
for (i in seq_len(length(trees) - 1)) {
    cur_node_data <- node_data[node_data$family == names(trees)[i], ]
    nds <-  cur_node_data$"node"
    nd_order <-  cur_node_data$order_nds
    if (length(nds) < 1) {
        message("no node for: ", names(trees)[i])
        next
    }
    tip_lbl <- trees[[i]]$tip.label
    ed <- trees[[i]]$edge
    col_tips <- rep("black", nrow(ed))
    ed_ar <- match(grep("arabia",  tip_lbl),  ed[, 2])
    ed_ot <- match(grep("outside", tip_lbl),  ed[, 2])
    col_tips[ed_ot] <- "steelblue"
    col_tips[ed_ar] <- "darkgreen"
    plot.phylo(trees[[i]], show.tip.label = FALSE, main = names(trees)[i],
               edge.color = col_tips, no.margin = FALSE)
    ## nodelabels(text = rep(NA, length(nds)), node = nds,
    ##           frame = "circ", col = "red", bg = "red", width = 3)
    HPDbars_fixed(trees[[i]], nodes = nds, lwd = 2, col = "red")
    nodelabels(text = nd_order, node = nds, frame = "circ", bg = "white")
    axisPhylo()
}
dev.off()


}
