genLabel <- function(dbTmp) {
  seqNm <- paste(dbTmp$family, dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, dbTmp$Loc,
                 paste(dbTmp$"Collection.Code", dbTmp$Catalog_number, sep=""),
                 dbTmp$Sample, dbTmp$Type, sep="_")
  if ( dbTmp$Type != "")
      seqNm <- paste(seqNm, dbTmp$Sample, sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
  seqNm <- gsub("_$", "", seqNm)
  seqNm <- gsub("\\s+$", "", seqNm)
  seqNm
}

genSp <- function(dbTmp) {
  seqNm <- paste(dbTmp$genusorhigher, dbTmp$modifier, dbTmp$species, sep="_")
  seqNm <- gsub("_{2,}", "_", seqNm)
  seqNm
}


genFasta <- function(db, out=file.path(tempdir(), paste(format(Sys.time(), "%Y%m%d-%H%M%S"), "seq.fas", sep="_"))) {
### db -- database in which the data is stored
### out -- file name of the fasta that will be generated
  for (i in 1:nrow(db)) {
    dbTmp <- db[i, ]
    seqNm <- genLabel(dbTmp)
    seqNm <- paste(">", seqNm, sep="")
    seqNm <- gsub("\\s+", "", seqNm)
    seqTmp <- dbTmp$Sequence
    cat(seqNm, "\n", dbTmp$Sequence, "\n", file=out,
        append=ifelse(i == 1, FALSE, TRUE), sep="")
  }
  TRUE
}
