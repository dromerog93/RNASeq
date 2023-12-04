snp2gene <- function(exon_bed, snp,
                    gene_name = "UNK", utr = 0, filter = FALSE) {
    nexon <- nrow(exon_bed)
    snp <- snp[snp$CHROM == exon_bed$chrom[1] &
                snp$POS >= min(exon_bed$start) - as.integer(utr) &
                snp$POS <= max(exon_bed$end) + as.integer(utr), ]
    nsnp <- nrow(snp)
    if (nsnp == 0) {
        message(sprintf("There are no mutations for %s", gene_name))
        return()
    }
    for (j in 1:nsnp) {
        check_snp <- 0
        for (i in 1:nexon) {
            if (snp$POS[j] > exon_bed[i, "start"] &&
            snp$POS[j] <= exon_bed[i, "end"]) {
                snp$genepos[j] <- snp$POS[j] -
                                    exon_bed[i, "start"] +
                                    exon_bed[i, "pos_exon"]
                snp$exon[j] <- exon_bed[i, "name"]
                snp$gene[j] <- gene_name
                check_snp <- 1
                break
            }
        }
        if (check_snp == 0) {
            snp$genepos[j] <- -1
            snp$exon[j] <- "out"
            snp$gene[j] <- gene_name
        }
    }
    if (filter) {
        return(snp[snp$exon != "out", ])
    } else {
        return(snp)
    }
}
