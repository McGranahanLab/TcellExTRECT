#' Calculate median values of coverage within exons in rolling windows
#'
#' @param vdj.region.df data frame of positions and coverage values
#' @param exons.selected Locations of exons
#' @param median.k rolling median window
#' @param median.thresh threshold to remove exons with low coverage
#' @param exons.to.use option to manually select exons
#' @return Data frame of positions and rolling median of coverage values within chosen exons
#' @importFrom stats median
#' @name medianExonCoverage

medianExonCoverage <- function(vdj.region.df, exons.selected, median.k = 50, median.thresh = 15, exons.to.use){

  pos <- NULL

  # Filter for positions within exons as expected and apply median filter and remove any exons required
  vdj.region.df.filt.exons.median <- lapply(seq_len(dim(exons.selected)[1]), function(x){
    tmp <- vdj.region.df %>%
      dplyr::filter(pos >= exons.selected$X2[x] & pos <= exons.selected$X3[x])
    tmp$reads <- medianFilter(tmp$reads, median.k)
    return(tmp)})[exons.to.use]

  # Remove any exons with very low values (suspected exon failure - possible 100% T cell)
  # Threshold of < 15 may not be appropriate on low coverage data sets/genomic regions
  median.values.exons <- sapply(vdj.region.df.filt.exons.median,
                                function(x) median(x$reads, na.rm = TRUE))
  exon.remove <- which(median.values.exons < median.thresh | is.na(median.values.exons))
  if(length(exon.remove) > 0){
    vdj.region.df.filt.exons.median <- vdj.region.df.filt.exons.median[-exon.remove]
  }
  vdj.region.df.filt.exons.median <- Reduce(rbind, vdj.region.df.filt.exons.median)

  # Get rid of any repeated rows that might exist - in theory I don't think they will
  vdj.region.df.filt.exons.median <- dplyr::distinct(vdj.region.df.filt.exons.median)
  return(list(vdj.region.df.filt.exons.median, exon.remove))
}
