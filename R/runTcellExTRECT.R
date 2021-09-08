#' Function to run T cell ExTRECT
#'
#' @param vdj.region.df data frame containing coverage values by position
#' @param exons.selected list of exon positions based on exome capture kit used
#' @param vdj.seg locations of segments used for calculation of TCRA score
#' @param hg19_or_38 hg19 or hg38 version of genome
#' @param exons.to.use option to manually select which exons to use (defaults to all)
#' @param GC_correct whether to use GC correction or not for output
#' @param median.k rolling median window
#' @param median.thresh threshold to remove exons with low coverage
#' @param sample_name name of sample run
#' @return data frame of TCRA T cell fractions with 95\% CI and QC fit value for quality of the solution
#' @name runTcellExTRECT
#' @export

runTcellExTRECT <- function(vdj.region.df, exons.selected,
                            vdj.seg, hg19_or_38 = 'hg19',
                            exons.to.use = NULL, GC_correct = TRUE,
                            median.k = 50, median.thresh = 15,
                            sample_name = 'test'){

  # Make sure colnames are correct
  colnames(exons.selected) <- c('X1','X2','X3')

  # By default use all exons
  if(is.null(exons.to.use)){
    exons.to.use <-  seq_len(dim(exons.selected)[1])
  }

  # Get GC content for exons
  exon.adjust.loc <- ifelse(hg19_or_38 == 'hg19',21999999,21531846)

  TCRA.exons.loc <- list()
  for(i in seq_len(dim(exons.selected)[1])){
    TCRA.exons.loc[[i]] <- c(exons.selected$X2[i] - exon.adjust.loc,
                             exons.selected$X3[i] - exon.adjust.loc)

  }
  # These are the GC content within the exons
  exons.gc.content <- exonwindowplot2(TCRA.exons.loc, TCRA_fasta[[1]],0)


  # 1. Check VDJ file:
  if(dim(vdj.region.df)[1] == 0){
    return(NULL)
  }

  # median filter coverage for normalisation
  median.exon.output <- medianExonCoverage(vdj.region.df, exons.selected,
                                           median.k, median.thresh,
                                           exons.to.use)
  vdj.region.df.filt.exons.median <- median.exon.output[[1]]
  exon.remove <- median.exon.output[[2]]

  # Calculate log ratio
  vdj.logR.df <- getLogRdf(vdj.region.df.filt.exons.median, vdj.seg, minCov = 0)

  # VDJ.QC check
  vdj.logR.df <- vdj.logR.df[!is.infinite(vdj.logR.df$Ratio), ]
  vdj.logR.df <- vdj.logR.df[!is.na(vdj.logR.df$Ratio), ]

  if(dim(vdj.logR.df)[1] == 0 | length(exon.remove) > 30){
    return(NA)}

  if(GC_correct){
    vdj.logR.df <- GCcorrect(vdj.logR.df,
                             exons = exons.selected,
                             exonList = exons.gc.content,
                             gene.fasta = TCRA_fasta)
    baselineAdj.out <- baselineAdj(vdj.logR.df,
                                   vdj.seg, GCcorrect = TRUE)
    vdj.logR.df <-baselineAdj.out[[1]]
    ci.95.value <- baselineAdj.out[[2]]

    vdj.fraction.output  <- getVDJFraction(vdj.logR.df, vdj.seg,ci.95.value, TRUE)

    qc.value <- calcQCvalue(vdj.logR.df, vdj.seg, GC_correct = TRUE)
    if(qc.value[2] > 4){
      warning('Fitted GAM model is very noisy and results may be less accurate. Run plotTcellExTRECT() to visually see fit.')
    }
  }else{
    baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = FALSE)
    vdj.logR.df <-baselineAdj.out[[1]]
    ci.95.value <- baselineAdj.out[[2]]

    vdj.fraction.output  <- getVDJFraction(vdj.logR.df, vdj.seg,ci.95.value, FALSE)
    qc.value <- calcQCvalue(vdj.logR.df, vdj.seg, GC_correct = FALSE)
    if(qc.value[2] > 4){
      warning('Fitted GAM model is very noisy and results may be less accurate. Run plotTcellExTRECT() to see fit.')
    }
  }

  out.df <- data.frame(sample = sample_name,
                       TCRA.tcell.fraction = vdj.fraction.output[1],
                       TCRA.tcell.fraction.lwr = vdj.fraction.output[2],
                       TCRA.tcell.fraction.upr = vdj.fraction.output[3],
                       qcFit = qc.value[2])
  return(out.df)
}
