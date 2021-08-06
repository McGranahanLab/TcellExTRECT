#' Method for calculating log ratio from single region
#'
#' @param region.df coverage values with a specific region
#' @param segs Location of segments used for normalisation and focal region
#' @param minCov Minimum GC coverage required
#' @return Calculated single log ratio based on segments
#' @importFrom stats median
#' @name getLogRdf

getLogRdf <-  function(region.df, segs, minCov = 0){

  pos <- reads <- NULL

  col_input <- 'reads'
  col.sym <- rlang::sym(col_input)
  # For random locations with no VDJ effect use beginning and end of TCRA
  tumour.random.covs1 <- region.df %>%
    dplyr::filter(pos <= segs[3,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs2 <- region.df %>%
    dplyr::filter(pos >= segs[4,]$start) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs <- c(tumour.random.covs1, tumour.random.covs2)
  # Use median of these values for normalisation to get "logR"
  n1 <- median(tumour.random.covs)

  # Select coverage values across TCRA
  tumour.test.covs1 <- region.df %>%
    dplyr::filter(pos >= segs[1,]$start & pos <= segs[1,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)

  # Select positions
  tumour.test.covs.pos <- region.df %>%
    dplyr::filter(pos >= segs[1,]$start & pos <= segs[1,]$end) %>%
    dplyr::filter(reads >=  minCov) %>% dplyr::select(pos) %>% `[[`(1)

  # Adjust for bin width
  # tumour.test.covs.pos <- tumour.test.covs.pos + bin.width/2

  # Get LogR df
  tumour.logR <- data.frame(pos = tumour.test.covs.pos,
                            Ratio = log2(tumour.test.covs1/n1))

  return(tumour.logR)

}
