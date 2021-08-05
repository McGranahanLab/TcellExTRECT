#' Function to plot output from T cell ExTRECT
#'
#' @param vdj.region.df data frame containing coverage values by position
#' @param exons.selected list of exon positions based on exome capture kit used
#' @param exons.gc.content GC content of exons
#' @param vdj.seg locations of segments used for calculation of TCRA score
#' @param hg19_or_38 hg19 or hg38 version of genome
#' @param exons.to.use option to manually select which exons to use (defaults to all)
#' @param sample_name name of sample run
#' @param median.k rolling median window
#' @param median.thresh threshold to remove exons with low coverage
#' @param txt.size size of annotation text
#' @param txt.height location on y axis of annotation txt
#' @return data frame of TCRA T cell fractions with 95% CI
#' @name plotTcellExTRECT
#' @export

plotTcellExTRECT <- function(vdj.region.df, exons.selected,
                             vdj.seg, hg19_or_38, exons.to.use = NULL,
                             median.k = 50, median.thresh = 15,
                             sample_name = '',txt.size = 4, txt.height = 1.5){

  # Make sure colnames are correct
  colnames(exons.selected) <- c('X1','X2','X3')

  data("TCRA_fasta")
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
    return(NULL)}

  data("TCRA_fasta")
  vdj.logR.df <- GCcorrect(vdj.logR.df,
                           exons = exons.selected,
                           exonList = exons.gc.content,
                           gene.fasta = TCRA_fasta)
  baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = TRUE)
  vdj.logR.df <-baselineAdj.out[[1]]
  ci.95.value <- baselineAdj.out[[2]]


  baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = FALSE)
  vdj.logR.df <-baselineAdj.out[[1]]
  ci.95.value <- baselineAdj.out[[2]]

  # Also add in segment info
  data(TCRA_segments)
  vdj.logR.df  <-vdj.logR.df  %>%
    dplyr::mutate(TRA_segment = TRA_segment_find_v(pos, hg19_or_38))

  vdj.logR.df$TRA_segment_short <- gsub('-','',gsub('[0-9]*','',
                                                    gsub('TR','',
                                                         vdj.logR.df$TRA_segment)))


  vdj.name.convert <- data.frame(TRA_segment_short = c('AC','AV','AJ', 'AVDV','DC','DD','DJ','DV','None'),
                                 TRA_segment_goodname = c('TCRA-C','TCRA-V','TCRA-J','TCRA-V/TCRD-V','TCRD-C','TCRD-D',
                                                          'TCRD-J','TCRD-V','None'))

  custom.cols <- c("TCRA-C" = "#E41A1C", "TCRA-V" = "#377EB8","TCRA-J" = "#4DAF4A",
                   "TCRA-V/TCRD-V"= "#984EA3",
                   "TCRD-C"=  "#FF7F00","TCRD-D"= "#FFFF33","TCRD-J"= "#A65628",
                   "TCRD-V"= "#F781BF",
                   "None" = 'lightgrey')

  p1 <- vdj.logR.df %>%
    dplyr::left_join(vdj.name.convert, 'TRA_segment_short') %>%
    ggplot2::ggplot(ggplot2::aes(pos, Ratio)) +
    ggplot2::geom_point(ggplot2::aes(col = TRA_segment_goodname)) + ggplot2::geom_smooth() +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_manual(name = 'VDJ segment', values =custom.cols  ) +
    ggplot2::ggtitle(paste0(sample_name,': TCRA loci (pre GC correction)'))+
    ggplot2::ylab('Log2 Read Depth Ratio') + ggplot2::xlab('Chr14 position') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'focal'), 'start']),
             col = 'red', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'focal'), 'end']),
               col = 'red', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local1'), 'start']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local1'), 'end']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local2'), 'start']),
                      col = 'black', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local2'), 'end']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::annotate("text", x = mean(unlist(vdj.seg[which(vdj.seg$segName == 'local1'), c('start','end')])),
                      y =  txt.height, label = "Norm region \n start", size = txt.size) +
    ggplot2::annotate("text", x = mean(unlist(vdj.seg[which(vdj.seg$segName == 'focal'), c('start','end')])),
                      y =  txt.height, label = "Focal region", size = txt.size) +
    ggplot2::annotate("text", x = mean(unlist(vdj.seg[which(vdj.seg$segName == 'local2'), c('start','end')])),
                      y =  txt.height, label = "Norm region \n end", size = txt.size)


  p2 <- vdj.logR.df %>%
    dplyr::left_join(vdj.name.convert, 'TRA_segment_short') %>%
    ggplot2::ggplot(ggplot2::aes(pos, Ratio.gc.correct)) +
    ggplot2::geom_point(ggplot2::aes(col = TRA_segment_goodname)) + ggplot2::geom_smooth() +
    ggplot2::scale_colour_manual(name = 'VDJ segment', values = custom.cols  ) +
    ggplot2::ggtitle(paste0(sample_name,': TCRA loci (post GC correction)'))+
    ggplot2::ylab('Log2 Read Depth Ratio (GC corrected)') + ggplot2::xlab('Chr14 position') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'focal'), 'start']),
               col = 'red', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'focal'), 'end']),
               col = 'red', linetype = 'dashed') +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local1'), 'start']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local1'), 'end']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local2'), 'start']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::geom_vline(ggplot2::aes(xintercept =vdj.seg[which(vdj.seg$segName == 'local2'), 'end']),
                        col = 'black', linetype = 'dashed') +
    ggplot2::annotate("text", x = mean(unlist(vdj.seg[which(vdj.seg$segName == 'local1'), c('start','end')])),
             y =  txt.height, label = "Norm region \n start", size = txt.size) +
    ggplot2::annotate("text", x = mean(unlist(vdj.seg[which(vdj.seg$segName == 'focal'), c('start','end')])),
             y =  txt.height, label = "Focal region", size = txt.size) +
    ggplot2::annotate("text", x = mean(unlist(vdj.seg[which(vdj.seg$segName == 'local2'), c('start','end')])),
             y =  txt.height, label = "Norm region \n end", size = txt.size)


  print(p1)
  print(p2)

}

TRA_segment_find <- function(x,hg19_or_38 = 'hg19'){
  if(hg19_or_38 == 'hg19'){
    gene.loc <- setdiff(intersect(which(TCRA_segments$start_hg19 <= x),
                                  which(x <= TCRA_segments$end_hg19)), c(1,57))
  }else{
    gene.loc <- setdiff(intersect(which(TCRA_segments$start_hg38 <= x),
                                  which(x <= TCRA_segments$end_hg38)), c(1,57))
  }

  return(ifelse(length(gene.loc) > 0, TCRA_segments$Gene[gene.loc], 'None'))
}

TRA_segment_find_v <- Vectorize(TRA_segment_find)
