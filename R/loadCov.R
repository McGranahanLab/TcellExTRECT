#' Function to load coverage from getCovFromBam
#'
#' @param covFile Path to coverage file from getCovFromBam output
#' @name loadCov

loadCov <- function(covFile){
  cov_df <- readr::read_tsv(covFile, col_names = FALSE)
  cov_df <- cov_df[,c(2,3)]
  colnames(cov_df) <- c('pos','reads')
  return(cov_df)
}
