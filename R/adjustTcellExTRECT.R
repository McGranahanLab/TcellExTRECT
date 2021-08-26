#' Adjust TCRA T cell fraction based on tumour purity and copy number
#'
#' @param TCRA.out Output of runTcellExTRECT
#' @param purity Tumour purity, between 0 and 1
#' @param TCRA.cn Tumour copynumber at TCRA locus
#' @param trustPurity If lower CI of TCRA score > 1 - purity assume tumour copy number is wrong and instead assume 2
#'
#' @return data frame of TCRA scores with purity/copy number adjusted outputs
#' @name adjustTcellExTRECT
#' @export

adjustTcellExTRECT <- function(TCRA.out, purity, TCRA.cn, trustPurity = TRUE){
  # solve visible binding issue
  TCRA.tcell.fraction <- TCRA.tcell.fraction.lwr <- TCRA.tcell.fraction.upr <- NULL
  rawRatio <- rawRatio.lwr <- rawRatio.upr <- maxPossible <- highTcellFlag <- NULL
  TCRA.tcell.fraction.adj <- TCRA.tcell.fraction.adj.lwr <-TCRA.tcell.fraction.adj.upr <- NULL

  # check purity is correct
  if(purity > 1 | purity < 0) stop('purity needs to be between 0 and 1')

  TCRA.out$purity <- purity
  TCRA.out$TCRA.cn <- TCRA.cn

  TCRA.out <- TCRA.out %>%
    dplyr::mutate(rawRatio = 1-TCRA.tcell.fraction) %>%
    dplyr::mutate(rawRatio.lwr = 1-TCRA.tcell.fraction.lwr) %>%
    dplyr::mutate(rawRatio.upr = 1-TCRA.tcell.fraction.upr) %>%
    dplyr::mutate(TCRA.tcell.fraction.adj= 1 - ((1-purity+(purity*TCRA.cn)/2)*rawRatio) - purity + ((purity*TCRA.cn)/2)) %>%
    dplyr::mutate(TCRA.tcell.fraction.adj.lwr= 1 - ((1-purity+(purity*TCRA.cn)/2)*rawRatio.lwr) - purity + ((purity*TCRA.cn)/2)) %>%
    dplyr::mutate(TCRA.tcell.fraction.adj.upr= 1 - ((1-purity+(purity*TCRA.cn)/2)*rawRatio.upr) - purity + ((purity*TCRA.cn)/2)) %>%
    dplyr::mutate(maxPossible = 1- purity) %>%
    dplyr::mutate(highTcellFlag = TCRA.tcell.fraction.lwr > maxPossible)

  if(trustPurity){
    TCRA.out <- TCRA.out %>%
      dplyr::mutate(TCRA.tcell.fraction.adj = ifelse(highTcellFlag, TCRA.tcell.fraction,TCRA.tcell.fraction.adj)) %>%
      dplyr::mutate(TCRA.tcell.fraction.adj.lwr = ifelse(highTcellFlag, TCRA.tcell.fraction.lwr,TCRA.tcell.fraction.lwr.adj)) %>%
      dplyr::mutate(TCRA.tcell.fraction.adj.upr = ifelse(highTcellFlag, TCRA.tcell.fraction.upr, TCRA.tcell.fraction.lwr.adj))
  }
  TCRA.out <- TCRA.out %>%
    dplyr::select(sample, purity, TCRA.cn,
           TCRA.tcell.fraction, TCRA.tcell.fraction.lwr, TCRA.tcell.fraction.upr,
           TCRA.tcell.fraction.adj, TCRA.tcell.fraction.adj.lwr, TCRA.tcell.fraction.adj.upr,
           highTcellFlag)
  return(TCRA.out)

}
