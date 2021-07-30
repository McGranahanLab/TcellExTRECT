# Tcell ExTRECT new template

# To add to data:
# TCRA exons - agilent/nimblegen
# TCRA fasta
# Seg files defining TCRA segments - hg19 and hg38
# Exons not to use in agilent v2.

library(TcellExTRECT)


# Things to put into a data file ----
seg.list <- list()
seg.list[['TCRA']]  <- data.frame(segName = c('all','focal', 'local1','local2'),
                                  start = c(22090057, 22800000, 22090057, 23016447),
                                  end = c(23221076, 22880000, 22298223, 23221076))

TCRA.exons <- readRDS('data/sureselect_TCRA_exons_hg19.RDS')
TCRA.exons <- TCRA.exons %>%
  dplyr::select(X1=chr,X2=start,X3=stop) %>% distinct()

vdj.chr.df <- data.frame(gene = c('TCRA','TCRB','TCRG','IGH','IGL','IGK'),
                         chr = c('chr14','chr7','chr7','chr14','chr22','chr2'))

# FASTA file of TCRA gene used for gc correction
TCRA.fasta <- read.fasta(file = "data/TCRA.fasta")

# Main script ----

# Need input (skip generation of cov.output.file for now)
cov.output.file  <- 'data/example_TCRA.txt'

# This should be region of TCRA gene
vdj.region.df <- read_tsv(cov.output.file, col_names = FALSE)

# This should be just what the input looks like!
vdj.region.df <- vdj.region.df[,c(2,3)]
colnames(vdj.region.df) <- c('pos','reads')

# input:
# vdj.regions
# TCRA.exons
# TCRA.exon.gc2 (can I generate this within function?)
# vdj.seg - have preset for TCRA etc
# Need to make hg19 and hg38 version
runTcellExTRECT <- function(vdj.region.df, exons.selected, exons.gc.content, vdj.seg, hg19_or_38){
  # 1. Check VDJ file:
  if(dim(vdj.region.df)[1] == 0){
    return(NULL)
  }

  # median filter coverage for normalisation
  median.exon.output <- medianExonCoverage(vdj.region.df, exons.selected)
  vdj.region.df.filt.exons.median <- median.exon.output[[1]]
  exon.remove <- median.exon.output[[2]]

  # Calculate log ratio
  vdj.logR.df <- getLogRdf(vdj.region.df.filt.exons.median, vdj.seg, minCov = 0)

  # VDJ.QC check
  vdj.logR.df <- vdj.logR.df[!is.infinite(vdj.logR.df$Ratio), ]
  vdj.logR.df <- vdj.logR.df[!is.na(vdj.logR.df$Ratio), ]

  if(dim(vdj.logR.df)[1] == 0 | length(exon.remove) > 30){
    return(NA)}

  vdj.logR.df <- GCcorrect(vdj.logR.df, exons = exons.selected, exonList = exons.gc.content)
  baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = TRUE)
  vdj.logR.df <-baselineAdj.out[[1]]
  ci.95.value <- baselineAdj.out[[2]]

  vdj.fraction.gc  <- getVDJFraction_upd2(vdj.logR.df, vdj.seg,ci.95.value, TRUE)

  baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = FALSE)
  vdj.logR.df <-baselineAdj.out[[1]]
  ci.95.value <- baselineAdj.out[[2]]

  vdj.fraction  <- getVDJFraction_upd2(vdj.logR.df, vdj.seg,ci.95.value, FALSE)

  return(c(vdj.fraction, vdj.fraction.gc))

}


plotTcellExTRECT <- function(vdj.region.df, exons.selected, exons.gc.content, vdj.seg, hg19_or_38){
  # 1. Check VDJ file:
  if(dim(vdj.region.df)[1] == 0){
    return(NULL)
  }

  # median filter coverage for normalisation
  median.exon.output <- medianExonCoverage(vdj.region.df, exons.selected)
  vdj.region.df.filt.exons.median <- median.exon.output[[1]]
  exon.remove <- median.exon.output[[2]]

  # Calculate log ratio
  vdj.logR.df <- getLogRdf(vdj.region.df.filt.exons.median, vdj.seg, minCov = 0)

  # VDJ.QC check
  vdj.logR.df <- vdj.logR.df[!is.infinite(vdj.logR.df$Ratio), ]
  vdj.logR.df <- vdj.logR.df[!is.na(vdj.logR.df$Ratio), ]

  if(dim(vdj.logR.df)[1] == 0 | length(exon.remove) > 30){
    return(NULL)}

  vdj.logR.df <- GCcorrect(vdj.logR.df, exons = exons.selected, exonList = exons.gc.content)
  baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = TRUE)
  vdj.logR.df <-baselineAdj.out[[1]]
  ci.95.value <- baselineAdj.out[[2]]

  #vdj.fraction.gc  <- getVDJFraction_upd2(vdj.logR.df, vdj.seg,ci.95.value, TRUE)

  baselineAdj.out <- baselineAdj(vdj.logR.df, vdj.seg, GCcorrect = FALSE)
  vdj.logR.df <-baselineAdj.out[[1]]
  ci.95.value <- baselineAdj.out[[2]]

  p1 <- vdj.logR.df %>%
    ggplot(aes(pos, Ratio)) +
      geom_point() + geom_smooth() +
      ggtitle('Pre GC correction')

  p2 <- vdj.logR.df %>%
    ggplot(aes(pos, Ratio.gc.correct)) +
    geom_point() + geom_smooth() +
    ggtitle('Post GC correction')

  #vdj.fraction  <- getVDJFraction_upd2(vdj.logR.df, vdj.seg,ci.95.value, FALSE)

  # return(c(vdj.fraction, vdj.fraction.gc))

}



medianExonCoverage <- function(vdj.region.df, exons.selected, median.k = 50, median.thresh = 15){
  # Filter for positions within exons as expected and apply median filter and remove any exons required
  vdj.region.df.filt.exons.median <- lapply(seq_len(dim(exons.selected)[1]), function(x){
    tmp <- vdj.region.df %>%
      filter(pos >= exons.selected$X2[x] & pos <= exons.selected$X3[x])
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
  vdj.region.df.filt.exons.median <- vdj.region.df.filt.exons.median %>%
    Reduce(rbind, .)

  # Get rid of any repeated rows that might exist - in theory I don't think they will
  vdj.region.df.filt.exons.median <- distinct(vdj.region.df.filt.exons.median)
  return(list(vdj.region.df.filt.exons.median, exon.remove))
}

baselineAdj <- function(vdj_logR_input, vdj_seg, GCcorrect = TRUE){

  ratio.col <- ifelse(GCcorrect, 'Ratio.gc.correct','Ratio')
  ratio.col <- rlang::sym(ratio.col)

  adjust.baseline.value <- vdj_logR_input %>%
  filter((pos >= vdj_seg[3,'start'] & pos <= vdj_seg[3,'end']) |
           (pos >= vdj_seg[4,'start'] & pos <= vdj_seg[4,'end'])) %>%
  summarise(gc.adjust = mean(!!ratio.col),
            CI.95.range = 1.96*sd(!!ratio.col)/sqrt(length(!!ratio.col)))

  # Look at focal too in GAM model
  adjust.model <- mgcv::gam(!!ratio.col~s(pos, bs = 'cs'), data = vdj_logR_input)
  adjust.fit.model <- confint(adjust.model, parm = "s(pos)",
                            partial_match = TRUE, type = 'simultaneous',
                            newdata = seq(vdj_seg[2,2], vdj_seg[2,3],by=100),
                            shift = TRUE)
  fit.loc <- which(adjust.fit.model$pos > vdj_seg[2,2] & adjust.fit.model$pos < vdj_seg[2,3])
  adjust.baseline.value2 <- list(mean(adjust.fit.model$est[fit.loc]),
                               mean(adjust.fit.model$upper[fit.loc]) - mean(adjust.fit.model$est[fit.loc]))

  adjust.value <- max(adjust.baseline.value[[1]], adjust.baseline.value2[[1]])
  ci.95.value <- ifelse(adjust.baseline.value2[[1]] > adjust.baseline.value[[1]],
                      adjust.baseline.value2[[2]], adjust.baseline.value[[2]])


  vdj_logR_input[!!ratio.col] <- vdj_logR_input[!!ratio.col] - adjust.value
return(list(vdj_logR_input, ci.95.value))
}



