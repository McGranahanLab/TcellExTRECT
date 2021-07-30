#' Method for calculating final VDJ fraction
#'
#' @param tumour.logR log ratio in region
#' @param segs Location of segments used for normalisation and focal region
#' @param norm.ci.95 Adjusted value for 95% CI
#' @param GC.correct Whether to use GC corrected value or not
#' @param ci_type Method to use for calculation of the confidence interval
#' @return VDJ fraction and upper and lower CI values
#' @name getVDJFraction_upd2
#' @export getVDJFraction_upd2
medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  # Maybe better than median?
  runMedian <- stats::runmed(x,k=filtWidth,endrule="constant")
  return(runMedian)
}
