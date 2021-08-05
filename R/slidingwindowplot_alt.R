#' Method for getting GC content in different window sizes
#'
#' @param windowsize Length of window to calculate GC content in
#' @param inputseq Sequence to calculate GC content for
#' @return data frame of GC content
#' @name slidingwindowplot_alt

slidingwindowplot_alt <- function(windowsize, inputseq){
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkGCs <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- seqinr::GC(chunk)
    chunkGCs[i] <- chunkGC
  }
  return(data.frame(loc = starts,GC = chunkGCs))
}
