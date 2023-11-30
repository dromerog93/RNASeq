#-----------------------------------------------------------------------------
#read_VCF read standard VCF file with "--annotate FORMAT/AD" option and 
#converts them into data frame with:
# - CHROM: The name of the sequence or genome region.
# - POS: The 1-based position of the variation on CHROM.
# - REF: The reference base(s) in the reference sequence.
# - ALT: The list of alternative variations/alleles in CHROM:POS.
# - QUAL: The quality score of the variation.
# - VARIANT: Variation type: SNP or INDEL.
# - AD_REF_XX: Read depth for the reference base.
# - AD_ALT_XX: Read depth for the alternative variation(s).
#-----------------------------------------------------------------------------
read_VCF <- function(file, f_QUAL = 0, f_AD = 0, output_file = "", 
                     sample_name = "", stats = TRUE){
  #file: Pathname of the VCF file.
  #f_QUAL: Threshold for quality value (>=). By default: 0.
  #f_AD: Threshold for total read depth [AD_REF + AD_ALT] (>=). By default: 0.
  #output_file: Pathname of the VCF-table file. By default: Output file is not 
  #generated.
  #sample_name: Name of the sample(s). By default: Named according to the file.
  #stats: Basic statistics screen printout. By default: TRUE.
  #---------------------------------------------------------------------------
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  if (!require("readr", quietly = TRUE)){install.packages("readr")}
  if (!require("tidyr", quietly = TRUE)){install.packages("tidyr")}
  if (!require("stringr", quietly = TRUE)){install.packages("stringr")}
  library(dplyr)
  library(readr, include.only = 'read_lines')
  library(tidyr, include.only = 'separate')
  library(stringr, include.only = c('str_split', 'str_detect'))
  defaultW <- getOption("warn")
  
  cat("Reading VCF file\n")
  
  #Skip VCF file header
  for (i in 0:10) {
    vcf <- readr::read_lines(file, n_max = 5000, skip = i * 5000)
    for (n in 1:5000) {
      if(substr(vcf[n], 1, 2) != "##"){
        break
      }
    }
    if (n != 5000){
      n <- n + (i * 5000)
      break
    }
  }
  rm(i, vcf)
  
  cat(sprintf("Header lines: %i\n", n - 1))
  #Obtain number and names of samples
  namescol_df <- read.delim(file = file,
                            header = FALSE, sep = "\t",
                            skip = n - 1,
                            nrows = 1)[1, ] %>% unlist
  if(length(namescol_df) < 10){
    stop("There are not enough columns in the VCF file.")
  }
  
  namescol_df[1] <- "CHROM"
  n_col_read <- length(namescol_df)
  n_sample <- n_col_read - 9
  cat(sprintf("Samples: %i\n", n_sample))
  if (sample_name[1] != "" & length(sample_name) != n_sample){
    stop("Incorrect sample names. Please indicate the name of the samples 
         with "sample_name" option.")
  }

  for(i in 10:n_col_read){
    if (sample_name == ""){
        if(substr(namescol_df[i], nchar(namescol_df[i]) - 3, nchar(namescol_df[i])) == ".bam"){
            namescol_df[i] <- substr(namescol_df[i], 1, nchar(namescol_df[i]) - 4)
        }
    } else {
        namescol_df[i] <- sample_name[i - 9]
    }
    cat(sprintf(" - %s\n", namescol_df[i]))
  }


  cat("VCF file check completed. Generating VCF table...\n")
  steps <- n_sample + 3
  cat(sprintf("[Step 1/%i] Data reading...\n",steps))

  #Reading data from the VCF file
  vcf_df <- read.delim(file = file,
                       header = FALSE, sep = "\t",
                       col.names = namescol_df,
                       skip = n)
  rm(i, n, n_col_read)
  cat(sprintf("Number of variations: %i\n", nrow(vcf_df)))
  
  cat(sprintf("[Step 2/%i] Data filtering by QUAL\n", steps))
  #Create the data frame with the filtered data
  vcf_df <- vcf_df %>% filter(QUAL >= f_QUAL)
  cat(sprintf("Number of variations after quality filtering: %i\n", nrow(vcf_df)))
  
  #Reading the allelic depth
  if(sum(stringr::str_detect(vcf_df$FORMAT, "AD")) == 0){
    stop("The vcf file does not contain allelic depth data (annotate FORMAT/AD).")
  }
  
  formato <- vcf_df$FORMAT[1] %>% stringr::str_split(pattern = ":") %>% unlist
  
  vcf_df_output <- vcf_df %>% select(CHROM, POS, REF, ALT, QUAL) %>%
    mutate(VARIANT = ifelse(stringr::str_detect(vcf_df$INFO, "INDEL") == TRUE, "INDEL", "SNP"))

  options(warn = -1)  
  for(i in 1:n_sample){
    cat(sprintf("[Step %i/%i] Allelic depth reading of %s...\n", i + 2, steps, sample_name[i]))
    #Obtaining AD data (REF and ALT).
    #NOTE: The depths of different variants for the same position will be added together.

    df_test <- data.frame(INFO = vcf_df[, i + 9]) %>%
      tidyr::separate(col = INFO, into = formato, sep = ":") %>%
      select(AD) %>%
      tidyr::separate(col = AD, into = c("REF", "ALT1", "ALT2"), sep = ",")
    df_test[] <- lapply(df_test, function(x) as.numeric(as.character(x)))
    df_test$ALT <- rowSums(df_test[, c(2:3)], na.rm = TRUE)
    
    vcf_df_output[, (7 + (i - 1) * 2)] <- df_test$REF
    vcf_df_output[, (8 + (i - 1) * 2)] <- df_test$ALT
    names(vcf_df_output)[(7 + (i - 1) * 2)] <- sprintf("AD_REF_%s", namescol_df[9 + i])
    names(vcf_df_output)[(8 + (i - 1) * 2)] <- sprintf("AD_ALT_%s", namescol_df[9 + i])
    if(sum(is.na(vcf_df_output[, (7 + (i - 1) * 2)])) != 0 | sum(is.na(vcf_df_output[, (8 + (i - 1) * 2)])) != 0){
      stop(sprintf("Some variant of %s does not have a valid depth value. Please, check VCF file.", sample_name[i]))
    }
  }
  options(warn = defaultW)
  
  cat(sprintf("[Step %i/%i] Data filtering by number of depth per position.\n",steps,steps))
  if(f_AD != 0){
    vcf_df_output <- subset(vcf_df_output, rowSums(vcf_df_output[, c(7:ncol(vcf_df_output))]) >= f_AD)
  }
  cat(sprintf("Number of variations after depth filtering: %i\n", nrow(vcf_df_output)))

  if(stats){
    cat(sprintf("SNP: %i\n",sum(grepl("SNP", vcf_df_output$VARIANT))))
    cat(sprintf("InDel: %i\n",sum(grepl("INDEL", vcf_df_output$VARIANT))))
    for(i in 1:n_muestras){
        cat(sprintf("Sample %s:\n", namescol_df[9 + i]))
        vcf_df_output[, paste("Index", sample_name[i], sep = "_")] <- vcf_df_output[, 6 + (2 * i)]/(vcf_df_output[, 5 + (2 * i)] + vcf_df_output[, 6 + (2 * i)])
        cat(sprintf("Variations (Index > 0.6): %i\n", sum(vcf_df_output[, paste("Index", sample_name[i], sep = "_")] > 0.6)))
        cat(sprintf("Heterozygosity (Index between 0.3 and 0.6): %i\n", sum(vcf_df_output[, paste("Index", sample_name[i], sep = "_")] <= 0.6 &
                                                                            vcf_df_output[, paste("Index", sample_name[i], sep = "_")] >= 0.3)))
        cat(sprintf("Reference (Index < 0.3): %i\n", sum(vcf_df_output[, paste("Index", sample_name[i], sep = "_")] < 0.3)))
    }
  }
  

  if(output_file != ""){
    cat(sprintf("Saving in file %s.\n",output_file))
    write.table(vcf_df_output, file = output_file, sep="\t", row.names = FALSE, quote = FALSE)
  }
  cat("Done!\n")
  
  return(vcf_df_output)
}
