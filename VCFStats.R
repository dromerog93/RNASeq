#-----------------------------------------------------------------------------
#VCFStats obtains basic information (Number of variants, snp and indel) from
#a VCF data.frame obtained with read_VCF.
#The function can apply different filters according to the quality value (QUAL),
# depth (AD) and SNP index for each position.
#-----------------------------------------------------------------------------
VCFStats <- function(vcf, standard = TRUE, chrom = FALSE, index = 0.70,
                    QUAL = "", AD = "", sample = "") {
    #vcf: VCF data.frame.
    #standard: Default threshold values for QUAL and AD: Q0-AD0, Q20-AD0, 
    #Q30-AD0, Q0-AD10 and Q0-AD20. By default: TRUE.
    #chrom: Analyze data for each chromosome. By default: FALSE.
    #QUAL: Threshold for quality value (>=). By default: empty.
    #AD: Threshold for total read depth [AD_REF + AD_ALT] (>=). By default: empty.
    #index: Threshold for SNPIndex [AD_ALT / (AD_REF + AD_ALT)] (>=). By default: 0.70
    #sample: Name of the sample(s). By default: Sample_#.
    #-------------------------------------------------------------------------
    
    cat("Initializing data frame and variables...\n")
    #Output data.frame
    if (chrom) {
        region <- unique(vcf$CHROM)
        region_chr <- region[grep("chrom|CHR|chr|CHROM",region)]
        n_chr <- length(region_chr)
        vcf_stats <- as.data.frame(matrix(rep(0, 9 + n_chr + 1), nrow = 1))
        names(vcf_stats) <- c("Sample", "QUAL", "AD",
                            "Variant", "SNP", "INDEL",
                            paste("Variant_I", index * 100, sep = ""),
                            paste("SNP_I", index * 100, sep = ""),
                            paste("INDEL_I", index * 100, sep = ""),
                            region_chr, "not_chr")
    } else {
        vcf_stats <- as.data.frame(matrix(rep(0, 9), nrow = 1))
        names(vcf_stats) <- c("Sample", "QUAL", "AD",
                            "Variant", "SNP", "INDEL",
                            paste("Variant_I", index * 100, sep = ""),
                            paste("SNP_I", index * 100, sep = ""),
                            paste("INDEL_I", index * 100, sep = ""))
    }
    
    #Threshold values
    if (standard){
        fQUAL <- c(0, 20, 30, 0, 0)
        fAD <- c(0, 0, 0, 10, 20)
    }
    if (QUAL != "" & AD != "") {
        if (length(QUAL) != length(AD)) {
            stop("QUAL and AD must have the same length.")
        } else {
            if (standard) {
                fQUAL <- c(c(0, 20, 30, 0, 0), QUAL)
                fAD <- c(c(0, 0, 0, 10, 20), fAD)
            } else {
                fQUAL <- QUAL
                fAD <- fAD
            }
        }
    } else {
        if (standard) {
            fQUAL <- c(0, 20, 30, 0, 0)
            fAD <- c(0, 0, 0, 10, 20)
        } else {
            message("No filter has been chosen.")
            fQUAL <- 0
            fAD <- 0
        }
    }
    n_filter <- length(fQUAL)

    #Sample names
    if (sample == ""){
        n_sample <- as.integer(length(grep("AD",colnames(vcf))) / 2)
        sample <- rep("Sample", n_sample)
        for (i in 1:n_sample) {
            sample[i] <- paste(sample[i], i, sep = "_")
        }
    } else {
        n_sample <- length(sample)
    }

    #Calculate statistics
    cat("Calculating statistics...\n")
    for (i in 1:n_sample) {
        adref <- grep("AD_REF", colnames(vcf))
        adalt <- grep("AD_ALT", colnames(vcf))
        for (j in 1:n_filter) {
            vcf_stats[(i - 1) * n_filter + j, 1] <- sample[i]
            vcf_stats[(i - 1) * n_filter + j, 2] <- fQUAL[j]
            vcf_stats[(i - 1) * n_filter + j, 3] <- fAD[j]
            vcf_stats[(i - 1) * n_filter + j, 4] <- nrow(vcf[vcf$QUAL >= fQUAL[j] &
                                                        vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j], ])
            vcf_stats[(i - 1) * n_filter + j, 5] <- nrow(vcf[vcf$VARIANT == "SNP" &
                                                        vcf$QUAL >= fQUAL[j] &
                                                        vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j], ])
            vcf_stats[(i - 1) * n_filter + j, 6] <- nrow(vcf[vcf$VARIANT == "INDEL" &
                                                        vcf$QUAL >= fQUAL[j] &
                                                        vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j], ])
            vcf_stats[(i - 1) * n_filter + j, 7] <- nrow(vcf[vcf$QUAL >= fQUAL[j] &
                                                        vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j] &
                                                        vcf[, colnames(vcf)[adalt[i]]] / (vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]]) >= Index, ])
            vcf_stats[(i - 1) * n_filter + j, 8] <- nrow(vcf[vcf$VARIANT == "SNP" &
                                                        vcf$QUAL >= fQUAL[j] &
                                                        vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j] &
                                                        vcf[, colnames(vcf)[adalt[i]]] / (vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]]) >= Index, ])
            vcf_stats[(i - 1) * n_filter + j, 9] <- nrow(vcf[vcf$VARIANT == "INDEL" &
                                                        vcf$QUAL >= fQUAL[j] &
                                                        vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j] &
                                                        vcf[, colnames(vcf)[adalt[i]]] / (vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]]) >= Index, ])
            if (chrom) {
                for (k in 1:n_chr) {
                    vcf_stats[(i - 1) * n_filter + j, 9 + k] <- nrow(vcf[vcf$QUAL >= fQUAL[j] &
                                                                    vcf[, colnames(vcf)[adref[i]]] + vcf[, colnames(vcf)[adalt[i]]] >= fAD[j] &
                                                                    vcf$CHROM == region_chr[k], ])
                }
            }
            vcf_stats[(i - 1) * n_filter + j, 9 + k + 1] <- vcf_stats[(i - 1) * n_filter + j, "Variant"] - sum(vcf_stats[(i - 1) * n_filter + j, c(10:(9 + k))])
        }
    }

    return(vcf_stats)
}
