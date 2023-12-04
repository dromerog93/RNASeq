GeneSNP <- function(bed, vcf, fasta, filename, type = "orf", aa_matrix = "aa.dat", n_threads = 1, only_file = TRUE) {
    require(seqinr)
    require(foreach)
    require(parallel)
    source("snp2gene.R")
    source("bed2exon.R")
    source("SNPanalysis.R")
    source("translate_CDS.R")
    source("reverse_DNA.R")
    aa_df <- read.delim(aa_matrix, sep = "\t")

    cat("Creating cluster...\n")
    ngenes <- nrow(bed)
    cluster <- parallel::makeCluster(as.integer(n_threads), type = "FORK")
    doParallel::registerDoParallel(cl = cluster)

    if (type == "snp" || type == "all") {
        cat("Getting SNP in genes...\n")
        vcf_snp <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
            snp2gene(bed2exon(bed[i, ]), vcf,
                    gene_name = bed[i, 4], filter = TRUE)
        }
    }

    if (type == "orf" || type == "cds" || type == "all") {
        cat("Getting the mutated sequences...\n")
        if (type == "all") {
            vcf_genes <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
                c(bed[i, 4],
                paste(seqinr::getSequence(fasta[bed[i, 4]])[[1]], collapse = ""),
                SNPanalysis(bed[i, ], vcf_snp[vcf_snp$gene == bed[i, 4], ],
                        paste(seqinr::getSequence(fasta[bed[i, 4]])[[1]], collapse = ""),
                        filte = FALSE))
            }
        } else {
            vcf_genes <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
                c(bed[i, 4],
                paste(seqinr::getSequence(fasta[bed[i, 4]])[[1]], collapse = ""),
                SNPanalysis(bed[i, ], vcf, paste(seqinr::getSequence(fasta[bed[i, 4]])[[1]], collapse = "")))
            }
        }
        vcf_genes <- as.data.frame(vcf_genes)
        colnames(vcf_genes) <- c("gene", "seq", "seqmut")
    }

    if (type == "orf" || type == "all") {
        cat("Translating sequences to amino acids...\n")
        orfs <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
            if (bed[bed$gene == vcf_genes$gene[i], 6] == "-") {
                translate_CDS(reverse_dna(vcf_genes[i, "seq"]))
            } else {
                translate_CDS(vcf_genes[i, "seq"])
            }
        }
        orfsmut <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
            if (bed[bed$gene == vcf_genes$gene[i], 6] == "-") {
                translate_CDS(reverse_dna(vcf_genes[i, "seqmut"]))
            } else {
                translate_CDS(vcf_genes[i, "seqmut"])
            }
        }
        vcf_genes$orf <- orfs
        vcf_genes$orfmut <- orfsmut
    }

    if (type == "indel" || type == "all") {
        cat("Getting significant INDELs in genes...\n")
        if (type == "all"){
            vcf_indel <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
                SNPanalysis(bed[i, ], vcf_snp[vcf_snp$gene == bed[i, 4], ],
                        type = "indel", filter = FALSE)
            }
        } else {
            vcf_indel <- foreach(i = 1:ngenes, .combine = rbind) %dopar% {
                SNPanalysis(bed[i, ], vcf, type = "indel")
            }
        }
    }

    cat("Creating file(s)...\n")
    if (type == "cds" || type == "all") {
        for (i in 1:nrow(vcf_genes)) {
            if (!(identical(vcf_genes$seq[i], vcf_genes$seqmut[i]))) {
                if (file.exists(paste(filename, "fasta", sep = "."))) {
                    write(c(paste(">", vcf_genes$gene[i], sep = ""), vcf_genes$seq[i]), file = paste(filename, "fasta", sep = "."), append = TRUE)
                    write(c(paste(paste(">", vcf_genes$gene[i], sep = ""), "_mut", sep = ""), vcf_genes$seqmut[i]), file = paste(filename, "fasta", sep = "."), append = TRUE)
                } else {
                    write(c(paste(">",vcf_genes$gene[i],sep=""), vcf_genes$seq[i]), paste(filename, "fasta", sep = "."))
                    write(c(paste(paste(">", vcf_genes$gene[i], sep = ""), "_mut", sep = ""), vcf_genes$seqmut[i]), file = paste(filename, "fasta", sep = "."), append = TRUE)
                }
            }
        }
    }

    if (type == "orf" || type == "all") {
        for (i in 1:nrow(vcf_genes)) {
            if (!(identical(vcf_genes$orf[i],vcf_genes$orfmut[i]))) {
                if (file.exists(paste(filename, "pep", sep = "."))) {
                    write(c(paste(">", vcf_genes$gene[i], sep = ""), vcf_genes$orf[i]), file = paste(filename, "pep", sep = "."), append = TRUE)
                    write(c(paste(paste(">", vcf_genes$gene[i], sep = ""), "_mut", sep = ""), vcf_genes$orfmut[i]), file = paste(filename, "pep", sep = "."), append = TRUE)
                } else {
                    write(c(paste(">", vcf_genes$gene[i], sep = ""), vcf_genes$orf[i]), file = paste(filename, "pep", sep = "."))
                    write(c(paste(paste(">",vcf_genes$gene[i],sep=""),"_mut",sep=""), vcf_genes$orfmut[i]), file = paste(filename, "pep", sep = "."), append = TRUE)
                }
            }
        }

        for (i in 1:nrow(vcf_genes)) {
            vcf_genes$length[i] <- nchar(vcf_genes$orf[i])
            orf_s <- str_split(vcf_genes$orf[i],"")[[1]]
            orfmut_s <- str_split(vcf_genes$orfmut[i],"")[[1]]
            vcf_genes$diff[i] <- sum(!(orf_s == orfmut_s))
            score_aa <- 0
            for (j in 1:length(orf_s)) {
                score_aa <- score_aa + as.integer(aa_df[(aa_df$Var1 == orf_s[j] & aa_df$Var2 == orfmut_s[j]), "Value"])
            }
            vcf_genes$score[i] <- score_aa
        }
        write.table(vcf_genes[, c(1, 6:8)], file = paste(filename, "dat", sep = "."), sep = "\t", row.names = FALSE, quote = FALSE)
    }

    if (type == "snp" || type == "all") {
        write.table(vcf_snp, file = paste(filename, "_SNP.tsv", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    }

    if (type == "indel" || type == "all") {
        write.table(vcf_indel, file = paste(filename, "_Indel.tsv", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
    }

    stopCluster(cluster)
    if (!only_file) {
        if (type == "orf" || type == "cds" || type == "all") {
            return(vcf_genes)
        } else if (type == "snp") {
            return(vcf_snp)
        } else if (type == "indel") {
            return(vcf_indel)
        }
    }
}
