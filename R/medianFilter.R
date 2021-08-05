#' Method for creating running medians
#'
#' @param x Vector
#' @param k Location of segments used for normalisation and focal region
#' @return vector with running median values
#' @name medianFilter

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
