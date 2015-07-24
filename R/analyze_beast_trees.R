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

    nd_data <- data.frame(nd_data, row.names = NULL, stringsAsFactors = FALSE)

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
                      each_nm <- lapply(y, function(x) {
                                            x <- gsub("_~|_aff|_cf", "", names(x))
                                            x <- gsub("sp_([0-9]+)", "sp.\\1", x)
                                            x <- gsub("truncata", "gracilis", x)
                                            x <- gsub("pentard", "flavomaculata", x)
                                            x <- gsub("impatiens", "aff. impatiens", x)
                                            x <- gsub("Labidodemas_NEEDS_CHECKING",
                                                      "Holothuria_hartmeyeri", x)
                                            x <- gsub("Labidodemas__Saudi_Arabia",
                                                      "Labidodemas_sp.1_Saudi_Arabia", x)
                                            x <- gsub("_Labidodemas_Moorea_UF5127",
                                                      "_Labidodemas_semperianum_Moorea_UF5127", x)
                                            tt <- sapply(strsplit(x, "_"), function(x)
                                                paste0(x[2:3], collapse = " "))
                                            res <- table(tt)
                                            paste(names(res[res >=  2]), collapse = "/")
                                        })
                      unlist(each_nm)
                  })
    nd_data$spp <- unlist(spp)

    nd_data <- arrange(nd_data, desc(h_hpd_med))
    attr(nd_data, "spp") <- sub_tips
    nd_data
}

if (FALSE) {

node_data_common <- fetch("summary_beast_trees") %>%
    mutate(n_ot = as.numeric(n_ot),
           n_ar = as.numeric(n_ar),
           rec_mono = as.logical(rec_mono)) %>%
    mutate(node = as.numeric(node)) %>%
    mutate(h_hpd_min = as.numeric(h_hpd_min)) %>%
    mutate(h_hpd_max = as.numeric(h_hpd_max)) %>%
    mutate(h_hpd_med = as.numeric(h_hpd_med)) %>%
    filter(!is.na(h_hpd_med)) %>%
    arrange(desc(h_hpd_med)) %>% ## not sure why needed... bug in dplyr?
    mutate(order_nds = order(h_hpd_med, decreasing = TRUE))

node_data <- node_data_common %>%
    filter(rec_mono & n_ot > 1 & n_ar > 1) %>%
    filter(grepl("cousteaui|aphanes|cinerascens|flavomaculata|hartmeyeri|impatiens|olivacea|immobilis|Ohshimella",
            .$spp)) %>%
    mutate(spp = paste0("(", order_nds, ") ", spp)) %>%
    mutate(Species = factor(spp, levels = as.character(spp)))

levels(node_data$Species) <- gsub("(aff\\.\\simpatiens)/", "\\1 (ESU Red Sea)/", levels(node_data$Species))
levels(node_data$Species) <- gsub("olivacea", "aff. olivacea", levels(node_data$Species))
levels(node_data$Species) <- gsub("(aff\\.\\simpatiens$)", "\\1 (ESU tiger)", levels(node_data$Species))

pdf(file = "tmp/node_ages.pdf", width = 10, height = 4)
ggplot(node_data, aes(x = factor(order_nds), y = h_hpd_med, colour = Species)) +
  geom_pointrange(aes(ymin = h_hpd_min, ymax = h_hpd_max)) +
  xlab("Node ID") + ylab("Divergence time (My)") +
  theme_bw()
dev.off()

node_data_supp <- node_data_common %>%
  ##filter(! node %in% node_data$node) %>%
  ##filter(!is.na(h_hpd_med) & rec_mono & n_ot <=  2 |  n_ar <= 2) %>%
  filter(grepl("stuhl|parva|polyplectana|hawaii?ensis", .$spp, ignore.case = T)) %>%
  mutate(spp = paste0("(", order_nds, ") ", spp)) %>%
  mutate(Species = factor(spp, levels = as.character(spp)))

levels(node_data_supp$Species) <- gsub("(stuhlmanni)", "\\1/Chiridota sp.5", levels(node_data_supp$Species))
levels(node_data_supp$Species) <- gsub("(kefersteinii)", "\\1/Polyplectana sp.4", levels(node_data_supp$Species))
levels(node_data_supp$Species) <- gsub("hawaiiensis", "aff. hawaiiensis", levels(node_data_supp$Species))



pdf(file = "tmp/node_ages_supp.pdf", width = 10, height = 4)
ggplot(node_data_supp, aes(x = factor(order_nds), y = h_hpd_med, colour = Species)) +
  geom_pointrange(aes(ymin = h_hpd_min, ymax = h_hpd_max)) +
  xlab("Node ID") + ylab("Divergence time (My)") +
  theme_bw()
dev.off()

####

trees <- fetch("get_beast_trees")

pdf(file = "tmp/family_trees.pdf", paper = "US")#,  height = 70, width = 12)
par(mar = c(2, 0, 1, 0))
for (i in seq_len(length(trees))) {
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
    plot.phylo(trees[[i]], show.tip.label = T, main = names(trees)[i],
               edge.color = col_tips, no.margin = FALSE, cex = .3)
    ## nodelabels(text = rep(NA, length(nds)), node = nds,
    ##           frame = "circ", col = "red", bg = "red", width = 3)
    HPDbars_fixed(trees[[i]], nodes = nds, lwd = 2, col = "red")
    ape::nodelabels(text = nd_order, node = nds, frame = "circ", bg = "white")
    axisPhylo()
}
dev.off()


library(tidyr)

mat <- read.csv(file = "data/matrix.csv")
mat[["RedSea-Aden"]] <- ifelse(mat$Red.Sea == 1 | mat$Aden == 1, 1, 0)

dt_mat <- mat %>% select(order,  family, `RedSea-Aden`, Oman) %>%
  group_by(order) %>%
  summarize("RedSea+Aden" = sum(`RedSea-Aden`),
            "Oman" = sum(Oman)) %>%
  mutate("Red Sea - Aden" = `RedSea+Aden`/sum(`RedSea+Aden`),
         "Oman" = Oman/sum(Oman)) %>%
  select(order, `Red Sea - Aden`, Oman) %>%
  gather(region, proportion, -order) %>%
  filter(order != "Molpadida")

dt_mat$order <- factor(dt_mat$order, levels = c("Apodida", "Dendrochirotida", "Aspidochirotida"))

pdf(file = "tmp/compare_fauna_1.pdf", height = 4)
ggplot(dt_mat) + geom_bar(aes(x = order, y = proportion, fill = region),
                          stat = "identity", position = "dodge") +
  coord_flip() +
  theme(legend.position = "top")
dev.off()

pdf(file = "tmp/compare_fauna_2.pdf", height = 4)
ggplot(dt_mat) + geom_bar(aes(x = region, y = proportion, fill = order),
                          stat = "identity") +
  coord_flip() +
  theme(legend.position = "top")
dev.off()



}
