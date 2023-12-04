SNPanalysis <- function(bed, snp, seq, type = "seq") {
    source("bed2exon.R")
    source("snp2gene.R")
    source("seqmut.R")

    exon <- bed2exon(bed)
    snp_gene <- snp2gene(exon, snp,
                        gene_name = as.character(bed[4]),
                        utr = 0, filter = TRUE)

    if (type == "indel") {
        if (class(snp_gene) == "NULL") {
            return()
        }
        indel_gene <- snp_gene[snp_gene$VARIANT == "INDEL", ]
        n_indel <- nrow(indel_gene)
        if (n_indel != 0) {
            j <- 0
            for (i in 1:n_indel) {
                if (grepl(",", indel_gene$ALT[i])) {
                    indel_gene[2 * j + 1, ] <- indel_gene[i, ]
                    indel_gene[2 * j + 2, ] <- indel_gene[i, ]
                    indel_gene[2 * j + 1, "ALT"] <- sub(",.*", "", indel_gene[i, "ALT"])
                    indel_gene[2 * j + 2, "ALT"] <- sub(".*,", "", indel_gene[i, "ALT"])
                    j <- j + 1
                }
            }
            indel_gene <- indel_gene[!grepl(",", indel_gene$ALT), ]
            
            n_indel <- nrow(indel_gene)
            for (i in 1:n_indel) {
                indel_gene$diff[i] <- (nchar(indel_gene$REF[i]) - nchar(indel_gene$ALT[i])) %% 3
            }
            return(indel_gene[indel_gene$diff != 0, ])
        } else {
            message(sprintf("There are no INDELs for %s", as.character(bed[4])))
            return()
        }
    } else if (type == "seq") {
        snp_gene <- snp_gene[snp_gene$VARIANT == "SNP", ]
        if (class(snp_gene) == "NULL") {
            return(seq)
        } else {
            snp_gene[, "ALT"] <- sub(",.*", "", snp_gene$ALT)
            return(seqmut(seq, snp_gene))
        }
    }
}
