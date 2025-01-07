#' Method for calculating final VDJ fraction
#'
#' @param tumour.logR log ratio in region
#' @param segs Location of segments used for normalisation and focal region
#' @param norm.ci.95 Adjusted value for 95\% CI
#' @param GC.correct Whether to use GC corrected value or not
#' @param ci_type Method to use for calculation of the confidence interval
#' @return VDJ fraction and upper and lower CI values
#' @importFrom stats confint
#' @import gratia
#' @name getVDJFraction
#' @export getVDJFraction

getVDJFraction <- function(tumour.logR, segs, norm.ci.95,
                                GC.correct = FALSE, ci_type = 'simultaneous'){
  if(GC.correct){
    tumour.genomic.region.model <- mgcv::gam(Ratio.gc.correct~s(pos, bs = 'cs'), data = tumour.logR)
  }else{
    tumour.genomic.region.model <- mgcv::gam(Ratio~s(pos, bs = 'cs'), data = tumour.logR)
  }

  focal.start <- segs[2,2]
  focal.end <- segs[2,3]

  fit.model <- confint(tumour.genomic.region.model, parm = "s(pos)",
                                   partial_match = TRUE, type = ci_type,
                       data = data.frame(pos = seq(segs[2,2], segs[2,3],by=100)),
                       shift = TRUE)
  
  
  if ("est" %in% colnames(fit.model)) {
    estimate_col <- "est"
    upper_col <- "upper"
    lower_col <- "lower"
  } else {
    estimate_col <- ".estimate"
    upper_col <- ".upper_ci"
    lower_col <- ".lower_ci"
  }
  
  
  fit.loc <- which(fit.model$pos > focal.start & fit.model$pos < focal.end)
  if(length(fit.loc) == 0){
    fit.loc <- c(which(fit.model$pos > focal.start)[1]-1,
                 which(fit.model$pos > focal.start)[1])}

  tumour.log2.score <- mean(fit.model[[estimate_col]][fit.loc])
  tumour.log2.score.min <- mean(fit.model[[lower_col]][fit.loc]) - norm.ci.95
  tumour.log2.score.max <- mean(fit.model[[upper_col]][fit.loc]) + norm.ci.95

  tumour.tcell.purity <- 1 - (2^tumour.log2.score)
  tumour.tcell.purity.max <- 1 - (2^tumour.log2.score.min)
  tumour.tcell.purity.min <- 1 - (2^tumour.log2.score.max)
  return(c(tumour.tcell.purity, tumour.tcell.purity.min, tumour.tcell.purity.max))
}
