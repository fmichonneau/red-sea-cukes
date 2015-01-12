
distFromTip <- function(tr, node, trPost, parallel=TRUE, cores=8) {
### trying to use Rcpp to do the intersect doesn't make it faster
    if (missing(trPost)) {
        trPost <- reorder(tr, order="postorder")
    }
    allDesc <- descendants(tr, node, type="all")
    lDesc <- allDesc[allDesc < nTips(tr)]
    if (parallel) {
        tmpD <- mclapply(lDesc, function(x) {
            pth <- intersect(allDesc, ancestors(trPost, x, type="ALL"))
            sumEdgeLength(tr, pth)
        }, mc.cores=cores)
    }
    else {
        tmpD <- lapply(lDesc, function(x) {
            pth <- intersect(allDesc, ancestors(trPost, x, type="ALL"))
            sumEdgeLength(tr, pth)
        })
    }
    unlist(tmpD[which.max(tmpD)])
}

findGroups <- function(tr, threshold=.015, experimental=FALSE, parallel=TRUE) {
### idea for optimization see if it makes a difference to supply the tree
### both in pre- and post-order when computing the distance to tip
  stopifnot(inherits(tr, "phylo4"))
  stopifnot(!hasDuplicatedLabels(tr))
  if (! require(igraph)) {
      stop("Install the igraph package")
  }

  grp <- vector("list", nTips(tr))

  ## find the distance to the tips for each internal node and select the nodes
  ## below the threshold

  if (FALSE) {
      getAncDist <- function(nd, tr, trPost) {
          tmpD <- distFromTip(tr, nd, trPost)
          if (tmpD < threshold) {
              distTip[as.character(nd)] <<- tmpD
              ancNd <- ancestors(trPost, nd, "parent")
              getAncDist(ancNd, tr, trPost)
          }
          else {
              allAnc <- ancestors(trPost, nd, "ALL")
              distTip[as.character(allAnc)] <<- threshold + 1
          }
      }

      trP <- reorder(tr, order="postorder")
      intNodes <- unique(trP@edge[, 1])
      intNodes <- intNodes[intNodes != 0]
      distTip <- numeric(length(intNodes))
      is.na(distTip) <- TRUE
      names(distTip) <- intNodes

      if (parallel) {
          xx <- mclapply(intNodes, function(nd) {
              if (is.na(distTip[as.character(nd)])) {
                  getAncDist(nd, tr, trP)
              }
          }, mc.cores=7)
          browser()
      }
      else {
          ## need to parallelize this!
          for (nd in intNodes) {
              if (is.na(distTip[as.character(nd)])) {
                  getAncDist(nd, tr, trP)
              }
          }
          lGrp <- distTip
      }
  }
  if (experimental) {
      trP <- reorder(tr, "postorder")
      intNodes <- (nTips(tr)+1):(nEdges(tr))
      lGrp <- sapply(intNodes, function(x) distFromTip(tr, x, trPost=trP))
  }
  else {
      intNodes <- (nTips(tr)+1):(nEdges(tr))
      if (parallel) {
          lGrp <- foreach(i=intNodes, .final=unlist) %dopar% {
              distFromTip(tr, i)
          }
          names(lGrp) <- intNodes
      }
      else {
          lGrp <- sapply(intNodes, function(x) distFromTip(tr, x))
      }
  }
  lGrp <- lGrp[lGrp <= threshold]

  ## find all the descendants for the nodes below the threshold
  descGrp <- sapply(as.numeric(names(lGrp)), function(x) descendants(tr, x))

  ## remove overlapping sets
  snglGrp <- sapply(descGrp, function(x) length(x) == 1)
  edgeGrp <- do.call("rbind", lapply(descGrp, function(x) {
      if(length(x) > 1) cbind(head(x, -1), tail(x, -1)) else NULL
  }))
  if (!is.null(edgeGrp)) {
      graphGrp <- graph.data.frame(edgeGrp)
      descGrp <- c(split(V(graphGrp)$name, clusters(graphGrp)$membership),
                   descGrp[snglGrp])

      ## add singletons
      missingGrps <- setdiff(nodeId(tr, "tip"), as.numeric(unlist(descGrp)))
      descGrp <- c(descGrp, as.list(missingGrps))
  } else {
      descGrp <- nodeId(tr, "tip")
  }

  ## Return tip labels
  grp <- sapply(descGrp, function(x) tipLabels(tr)[as.numeric(x)])

  ## build a phylo4d object for the results
  dTip <- data.frame(Groups=rep(1:length(grp), sapply(grp, length)))
  rownames(dTip) <- unlist(grp)
  if (inherits(tr, "phylo4d"))
    return(addData(tr, tip.data=dTip))
  else if (inherits(tr, "phylo4"))
    return(phylo4d(tr, tip.data=dTip))
}


### Example code
## spGrp <- findGroups(trP4, threshold=.025)
## spGrpCopy <- spGrp
## tipLabels(spGrp) <- paste(tipData(spGrp)$Group, tipLabels(spGrp), sep="_")

## grpLbl <- paste("^", 1:max(tipData(spGrp)), "_", sep="")

## pdf(file="treeWithBars-025.pdf", height=100, width=10)
## par(mai=c(0.5,0,2,0), xpd=T)
## plot(as(spGrp, "phylo"), cex=.5, show.tip.label=T, no.margin=F, label.offset=0)
## barMonophyletic(grpLbl, as(spGrp, "phylo"), extra.space=0.01, cex.plot=.5, cex.text=.5,
##                 bar.at.tips=TRUE, include.tip.label=TRUE)
## add.scale.bar()
## dev.off()
