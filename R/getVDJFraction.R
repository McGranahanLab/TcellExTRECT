#' Original method for calculating final VDJ fraction
#'
#' @param tumour.logR log ratio in region
#' @param segs Location of segments used for normalisation and focal region
#' @param GC.correct Whether to use GC corrected value or not
#' @return VDJ fraction and upper and lower CI values
#' @name getVDJFraction
#'
getVDJFraction <- function(tumour.logR, segs, GC.correct = FALSE){
  # Create model in ggplot - REPLACE with ASPCF
  if(GC.correct){
    tumour.genomic.region.p1 <- tumour.logR %>%
      ggplot2::ggplot(aes(pos, Ratio.gc.correct)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth()
  }else{
    tumour.genomic.region.p1 <- tumour.logR %>%
      ggplot2::ggplot(aes(pos, Ratio)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth()
  }


  # Extract mean logR at location of predicted focal VDJ recombination from LOESS model
  tumour.log2.score <- coverageScoreFun_log2_vdj(tumour.genomic.region.p1, segs)

  # Calculate
  tumour.tcell.purity <- 1 - (2^tumour.log2.score)
  return(tumour.tcell.purity)
}


# Function to calculate logR from GAM (generalised additive model) if observations > 1000, loess if not
# The fitted model is used in the calculation of the raw VDJ fraction
coverageScoreFun_log2_vdj <- function(plot.object, seg.df){
  focal.start <- seg.df[2,2]
  focal.end <- seg.df[2,3]
  fit.model <- ggplot2::ggplot_build(plot.object)$data[[2]]

  fit.loc <- which(fit.model$x > focal.start & fit.model$x < focal.end)
  if(length(fit.loc) == 0){
    fit.loc <- c(which(fit.model$x > focal.start)[1]-1,
                 which(fit.model$x > focal.start)[1])}

  vdj.y <- mean(fit.model$y[fit.loc])
  return(vdj.y)
}
