seqmut <- function(seq, snp) {
    nsnp <- nrow(snp)
    seq <- tolower(strsplit(seq, "")[[1]])
    for (i in 1:nsnp) {
        if (seq[snp$genepos[i]] == tolower(snp$REF[i])) {
            seq[snp$genepos[i]] <- snp$ALT[i]
        } else {
            message(sprintf("SNP number %s (POS: %s) 
                    does not match the reference.\n", i, snp$POS[i]))
        }
    }
    return(paste(seq, collapse = ""))
}
