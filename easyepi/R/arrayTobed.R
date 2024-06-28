library(data.table)
library(dplyr)
# array to bed-----------------------------------------------------------------
#' Convert a data frame to a bed format
#'
#' This function takes a data frame and the name of a methylation chip, and returns a new data frame
#' with the data formatted for the specified chip. The function reads in chip-specific probe information
#' and maps the data in the input data frame to the chip format.
#'
#' @param df A data frame containing the data to be formatted.
#' @param chip A character string indicating the chip to be used. Supported chips are '40K', '450K', and '850K'.
#' @return A data frame with the data formatted for the specified chip.
arrayTobed <- function(df, chip){
  if (chip =='40K'){
    chip_path <- '../Illumina HorvathMammalianMethylChip40 BeadChip/GPL28271-57075.txt'
    chip <- read.csv(chip_path, comment.char = '#', fill = T, header = T, sep = '\t')
    map <- data.frame(ID = chip$ID, chr = chip$Human.CHR, start = chip$Human.MAPINFO-1, end = chip$Human.MAPINFO+1,
                      site = paste(chip$Human.CHR, chip$Human.MAPINFO, sep = '.'))
    map <- map[!grepl("rs|ch", map$ID), ]
    map <- map[!duplicated(map$site), ]
  } else if (chip =='450K'){
    chip_path <- '../Illumina HumanMethylation450 BeadChip/GPL13534-11288.txt'
    chip <- read.csv(chip_path, comment.char = '#', fill = T, header = T, sep = '\t')
    map <- data.frame(ID = chip$ID, chrom = paste0('chr', chip$CHR), start = chip$MAPINFO-1, end = chip$MAPINFO+1,
                      site = paste(paste0('chr', chip$CHR), chip$MAPINFO, sep = '.'))
    map <- map[!grepl("rs|ch", map$ID), ]
    map <- map[!duplicated(map$site), ]
  } else if (chip =='850K'){
    chip_path <- '../Infinium MethylationEPIC/GPL21145-48548.txt'
    chip <- read.csv(chip_path, comment.char = '#', fill = T, header = T, sep = '\t')
    map <- data.frame(ID = chip$ID, chrom = paste0('chr', chip$CHR), start = chip$MAPINFO-1, end = chip$MAPINFO+1,
                      site = paste(paste0('chr', chip$CHR), chip$MAPINFO, sep = '.'))
    map <- map[!grepl("rs|ch", map$ID), ]
    map <- map[!duplicated(map$site), ]
  }
  df$ID <- rownames(df)
  result <- inner_join(map, df, by = "ID")
  rownames(result) <- result$ID
  result$ID <- NULL
  result$site <- NULL
  return(result)
}
