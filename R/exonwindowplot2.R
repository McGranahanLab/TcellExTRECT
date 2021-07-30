#' Method for getting GC content in selected exons
#'
#' @param exon.loc Location of exons
#' @param inputseq Sequence to calculate GC content for
#' @param extend1 extend exon sequence by specific number of bases
#' @return data frame of GC content for exons
#' @name exonwindowplot2
#' @export

exonwindowplot2 <- function(exon.loc, inputseq, extend1){
  n <- length(exon.loc)
  chunkGCs <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[(exon.loc[[i]][1]-extend1):c(exon.loc[[i]][2] + extend1)]
    chunkGC <- GC(chunk)
    chunkGCs[i] <- chunkGC
  }
  plot(seq(n),chunkGCs,type="b",xlab="Exon",ylab="GC content")
  return(data.frame(exon = seq(n), GC = chunkGCs))
}
