#' Adjust TCRA T cell fraction based on tumour purity and copy number
#'
#' @param TCRA.out Output of runTcellExTRECT
#' @param purity Tumour purity
#' @param TCRA.cn Tumour copynumber at TCRA locus
#' @param trustPurity If lower CI of TCRA score > 1 - purity assume tumour copy number is wrong and instead assume 2
#'
#' @return data frame of TCRA scores with purity/copy number adjusted outputs
#' @name adjustTcellExTRECT
#' @export

adjustTcellExTRECT <- function(TCRA.out, purity, TCRA.cn, trustPurity = TRUE){
  TCRA.out$purity <- purity
  TCRA.out$TCRA.cn <- TCRA.cn

  TCRA.out <- TCRA.out %>%
    mutate(rawRatio = 1-TCRA.tcell.fraction) %>%
    mutate(rawRatio.lwr = 1-TCRA.tcell.fraction.lwr) %>%
    mutate(rawRatio.upr = 1-TCRA.tcell.fraction.upr) %>%
    mutate(TCRA.tcell.fraction.adj= 1 - ((1-purity+(purity*TCRA.cn)/2)*rawRatio) - purity + ((purity*TCRA.cn)/2)) %>%
    mutate(TCRA.tcell.fraction.adj.lwr= 1 - ((1-purity+(purity*TCRA.cn)/2)*rawRatio.lwr) - purity + ((purity*TCRA.cn)/2)) %>%
    mutate(TCRA.tcell.fraction.adj.upr= 1 - ((1-purity+(purity*TCRA.cn)/2)*rawRatio.upr) - purity + ((purity*TCRA.cn)/2)) %>%
    mutate(maxPossible = 1- purity) %>%
    mutate(highTcellFlag = TCRA.tcell.fraction.lwr > maxPossible)

  if(trustPurity){
    TCRA.out <- TCRA.out %>%
      mutate(TCRA.tcell.fraction.adj = ifelse(highTcellFlag, maxPossible, TCRA.tcell.fraction)) %>%
      mutate(TCRA.tcell.fraction.adj.lwr = ifelse(highTcellFlag, maxPossible, TCRA.tcell.fraction.lwr)) %>%
      mutate(TCRA.tcell.fraction.adj.upr = ifelse(highTcellFlag, maxPossible, TCRA.tcell.fraction.upr))
  }
  TCRA.out <- TCRA.out %>%
    select(sample, purity, TCRA.cn,
           TCRA.tcell.fraction, TCRA.tcell.fraction.lwr, TCRA.tcell.fraction.upr,
           TCRA.tcell.fraction.adj, TCRA.tcell.fraction.adj.lwr, TCRA.tcell.fraction.adj.upr,
           highTcellFlag)
  return(TCRA.out)

}
