translate_CDS <- function(seq){
  #Librerias
  if (!require("dplyr", quietly = TRUE)){install.packages("dplyr")}
  if (!require("seqinr", quietly = TRUE)){install.packages("seqinr")}
  library(dplyr)
  library(seqinr, include.only = "s2c")
  
  #Data frame de codon-aminoacido
  seq <- tolower(seq)
  df_codons <- data.frame(codon = c("ttt","ttc","tta","ttg",
                                    "tct","tcc","tca","tcg",
                                    "tat","tac","taa","tag",
                                    "tgt","tgc","tga","tgg",
                                    "ctt","ctc","cta","ctg",
                                    "cct","ccc","cca","ccg",
                                    "cat","cac","caa","cag",
                                    "cgt","cgc","cga","cgg",
                                    "att","atc","ata","atg",
                                    "act","acc","aca","acg",
                                    "aat","aac","aaa","aag",
                                    "agt","agc","aga","agg",
                                    "gtt","gtc","gta","gtg",
                                    "gct","gcc","gca","gcg",
                                    "gat","gac","gaa","gag",
                                    "ggt","ggc","gga","ggg"),
                          AA = c("F","F","L","L",
                                 "S","S","S","S",
                                 "Y","Y","*","*",
                                 "C","C","*","W",
                                 "L","L","L","L",
                                 "P","P","P","P",
                                 "H","H","Q","Q",
                                 "R","R","R","R",
                                 "I","I","I","M",
                                 "T","T","T","T",
                                 "N","N","K","K",
                                 "S","S","R","R",
                                 "V","V","V","V",
                                 "A","A","A","A",
                                 "D","D","E","E",
                                 "G","G","G","G"))
  
  #Division de la secuencia en codones
  seq_cds <- strsplit(seq,"(?<=.{3})",perl=T)[[1]]
  
  #Traduccion
  seq_aa <- ""
  for (i in 1:length(seq_cds)){
    if(nchar(seq_cds[i]) == 3){
      seq_aa[i] <- df_codons[df_codons$codon == seq_cds[i],2]
    } else {
      seq_aa[i] <- "-"
    }
  }
  
  return(paste(seq_aa,collapse=""))
}
