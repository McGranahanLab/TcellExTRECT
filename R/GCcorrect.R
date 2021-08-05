#' Method for correcting for GC biases
#'
#' @param tumour.logR LogR values
#' @param exons Location of exons
#' @param exonList GC content of exons
#' @param gene.fasta FASTA file for VDJ gene e.g TCRA
#' @param gene.fasta offset for FASTA file, e.g 21999999 for TCRA in hg19
#' @param hg19_or38 What genome version to use, must be 'hg19' or 'hg38'
#' @param sliding number of bp for gc windows
#' @return data frame of GC correct logR values
#' @name GCcorrect

GCcorrect <- function(tumour.logR, exons , exonList, gene.fasta, hg19_or38 = 'hg19',sliding = 1000){

  gene.fasta.start <- ifelse(hg19_or38 == 'hg19', 21999999, 21531846) # Get number for hg38)

  TCRA.gc.df <- slidingwindowplot_alt(sliding, gene.fasta[[1]])
  TCRA.gc.df$pos <- TCRA.gc.df$loc + gene.fasta.start
  gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = TCRA.gc.df)

  # This is the function to get the smoothed GC content at any position
  get_gc_prediction <- function(x){mgcv::predict.gam(gam.model, newdata = data.frame(pos = x))}

  tumour.logR <- tumour.logR %>%
    dplyr::mutate(exon = exonPosFun_v(pos, exons)) %>%
    dplyr::left_join(exonList, 'exon') %>%
    dplyr::rename(exon.gc = GC) %>%
    dplyr::mutate(exon.gc2 = exon.gc^2) %>%
    dplyr::mutate(smooth.gc = get_gc_prediction(pos))  %>%
    dplyr::mutate(smooth.gc2 = smooth.gc^2)


  gc.lm = lm(Ratio ~ exon.gc + exon.gc2 + smooth.gc + smooth.gc2,
             y = TRUE, data = tumour.logR)

  # Look at if interested
  # summary(gc.lm)

  tumour.logR$Ratio.gc.correct <- gc.lm$residuals
  return(tumour.logR)

}

exonPosFun <- function(x, exons = TCRA.exons){
  rev(which(exons$X2 <= x))[1]
}
exonPosFun_v <- Vectorize(exonPosFun, vectorize.args = 'x')

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




