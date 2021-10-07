#' Function to create a list of exons from bed file of capture regions from exon custom kit
#'
#' @param pathToBed Path to bed file (tsv format assumed)
#' @param hg19_or_38 hg19 or hg38 version of genome
#' @return data frame of TCRA exons as provided by capture kit
#' @name createExonDFBed
#' @export
createExonDFBed <- function(pathToBed, hg19or38 = 'hg19'){
  warning('Note: exome capture kits not officially supported by TcellExTRECT may have unknown biases or lack necessary coverage within TCRA locus')
  input_bed <- readr::read_tsv(pathToBed, col_names = FALSE)

  tcra_chr <- c('chr14','14')

  if(hg19or38 == 'hg19'){
    tcra_start <- 22090057
    tcra_end <- 23221076
  }
  if(hg19or38 == 'hg38'){
    tcra_start <- 21621904
    tcra_end <- 22752132
  }

  output.df <- input_bed %>%
    as.data.frame() %>%
    dplyr::filter(X1 %in% tcra_chr) %>%
    dplyr::mutate(X2 = as.numeric(X2)) %>%
    dplyr::filter(X2 > tcra_start) %>%
    dplyr::mutate(X3 = as.numeric(X3)) %>%
    dplyr::filter(X3 < tcra_end) %>%
    dplyr::select(X1, X2, X3) %>%
    dplyr::mutate(X1 = 'chr14')

  return(output.df)
}
