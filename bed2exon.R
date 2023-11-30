bed2exon <- function(bed) {
    nexon <- as.integer(bed[10])
    exon <- data.frame(chrom = rep(as.character(bed[1]), nexon), start = rep(0, nexon), end = rep(0, nexon),
                    name = rep("E", nexon), score = rep(0, nexon), strand = rep(as.character(bed[6]), nexon),
                    cds_start = rep(0, nexon), cds_end = rep(0, nexon), rgb = rep("0,0,0", nexon),
                    nexon = rep(1, nexon), length_exon = rep(0, nexon), pos_exon = rep(0, nexon))
    lexon <- strsplit(as.character(bed[11]), ",")[[1]]
    pexon <- strsplit(as.character(bed[12]), ",")[[1]]
    for (i in 1:nexon) {
        exon[i, "start"] <- as.integer(pexon[i]) + as.integer(bed[2])
        exon[i, "end"] <- as.integer(pexon[i]) + as.integer(bed[2]) + as.integer(lexon[i])
        exon[i, "cds_start"] <- as.integer(pexon[i]) + as.integer(bed[2])
        exon[i, "cds_end"] <- as.integer(pexon[i]) + as.integer(bed[2]) + as.integer(lexon[i])
        exon[i, "name"] <- paste("E", i, sep = "")
        exon[i, "length_exon"] <- as.integer(lexon[i])
        if (i == 1) {
            exon[i, "pos_exon"] <- 0
        } else {
            exon[i, "pos_exon"] <- as.integer(lexon[i - 1]) + exon[i - 1, "pos_exon"]
        }
    }
    return(exon)
}
