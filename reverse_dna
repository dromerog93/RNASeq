reverse_dna <- function(seq){
  seq_tmp <- strsplit(tolower(seq), "")[[1]]
  seq_reverse <- rev(seq_tmp)
  seq_reverse_tmp <- gsub("a","1",
                          gsub("t","2",
                               gsub("c","3",
                                    gsub("g","4",seq_reverse))))
  seq_reverse  <- gsub("1","t",
                       gsub("2","a",
                            gsub("3","g",
                                 gsub("4","c",seq_reverse_tmp))))
  
  return(paste(seq_reverse,collapse=""))
}
