###############################################################################
### TIL-ExTRECT TCRA score generation - GC correction                       ###
### Template with annotation                                                ###
###                                                                         ###
### February 18th 2021 - Bobby Bentham                                      ###
###############################################################################

# Libraries needed
library(tidyverse)
library(rlang)
library(readxl)
library(seqinr)

# TIL ExTRECT setup ----

# seg.list is a list of genomic segments
# all defines the gene of interest (TCRA), focal is the region of VDJ deletion, local1/2
# are the regions at the beginning/end of the gene used for normalisation
# set up as a list for future support of other VDJ recombination genes.
# Code is set up so that the order all, focal, local1, local2 must be maintained
seg.list <- list()
seg.list[['TCRA']]  <- data.frame(segName = c('all','focal', 'local1','local2'),
                                  start = c(22090057, 22800000, 22090057, 23016447),
                                  end = c(23221076, 22880000, 22298223, 23221076))

# You only want to look at regions within exons captured by the kit used.
# Here has support for sureselect and nimblegen kir.
# For other regions generate from bams as appropriate.
TCRA.exons.nimblegen <- readRDS('data/nimblegen_TCRA_exons_hg19.RDS')
TCRA.exons <- readRDS('data/sureselect_TCRA_exons_hg19.RDS')
TCRA.exons <- TCRA.exons %>%
  dplyr::select(X1=chr,X2=start,X3=stop) %>% distinct()

# Data frame with chromosome of genes, always using TCRA now so not really used.
vdj.chr.df <- data.frame(gene = c('TCRA','TCRB','TCRG','IGH','IGL','IGK'),
                         chr = c('chr14','chr7','chr7','chr14','chr22','chr2'))

# FASTA file of TCRA gene used for gc correction
TCRA.fasta <- read.fasta(file = "data/TCRA.fasta")

# Functions  ----

# From a data frame and defined segs (focal etc) calculate "log ratio"
# Calculates the overage in the two local regionsm and calculates the median value
# which is then used for normalisation
getLogRdf <-  function(region.df, segs, minCov = 0){
  col_input <- 'reads'
  col.sym <- rlang::sym(col_input)
  # For random locations with no VDJ effect use beginning and end of TCRA
  tumour.random.covs1 <- region.df %>%
    filter(pos <= segs[3,]$end) %>%
    filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs2 <- region.df %>%
    filter(pos >= segs[4,]$start) %>%
    filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)
  tumour.random.covs <- c(tumour.random.covs1, tumour.random.covs2)
  # Use median of these values for normalisation to get "logR"
  n1 <- median(tumour.random.covs)

  # Select coverage values across TCRA
  tumour.test.covs1 <- region.df %>%
    filter(pos >= segs[1,]$start & pos <= segs[1,]$end) %>%
    filter(reads >=  minCov) %>% dplyr::select(!!col.sym) %>% `[[`(1)

  # Select positions
  tumour.test.covs.pos <- region.df %>%
    filter(pos >= segs[1,]$start & pos <= segs[1,]$end) %>%
    filter(reads >=  minCov) %>% dplyr::select(pos) %>% `[[`(1)

  # Adjust for bin width
  # tumour.test.covs.pos <- tumour.test.covs.pos + bin.width/2

  # Get LogR df
  tumour.logR <- data.frame(pos = tumour.test.covs.pos,
                            Ratio = log2(tumour.test.covs1/n1))

  return(tumour.logR)

}

# From a 'logR' coverage data frame generated with getLogRdf calculate
# the raw VDJ fraction pre adjustment of purity/tumour ploidy.
getVDJFraction <- function(tumour.logR, segs, GC.correct = FALSE){
  # Create model in ggplot - REPLACE with ASPCF
  if(GC.correct){
    tumour.genomic.region.p1 <- tumour.logR %>%
      ggplot(aes(pos, Ratio.gc.correct)) +
      geom_point() +
      geom_smooth()
  }else{
    tumour.genomic.region.p1 <- tumour.logR %>%
      ggplot(aes(pos, Ratio)) +
      geom_point() +
      geom_smooth()
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
  fit.model <- ggplot_build(plot.object)$data[[2]]

  fit.loc <- which(fit.model$x > focal.start & fit.model$x < focal.end)
  if(length(fit.loc) == 0){
    fit.loc <- c(which(fit.model$x > focal.start)[1]-1,
                 which(fit.model$x > focal.start)[1])}

  vdj.y <- mean(fit.model$y[fit.loc])
  return(vdj.y)
}

# Function used in normalisation and smoothing of exon coverage values using a rolling median
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
  runMedian <- runmed(x,k=filtWidth,endrule="constant")
  return(runMedian)
}

# GC correction functions:
# Functions to calculate both the GC content within each exon and smoothed across region.
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

slidingwindowplot_alt <- function(windowsize, inputseq){
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkGCs <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    chunkGCs[i] <- chunkGC
  }
  return(data.frame(loc = starts,GC = chunkGCs))
}

exonPosFun <- function(x, exons = TCRA.exons){
  rev(which(exons$X2 <= x))[1]
}
exonPosFun_v <- Vectorize(exonPosFun, vectorize.args = 'x')

# Main GC correct method, input is logR data frame and uses linear model
# to get the residules post GC correction (similar to ASCAT method)
GCcorrect <- function(tumour.logR, exons =  TCRA.exons, exonList = TCRA.exon.gc2){
  tumour.logR <- tumour.logR %>%
    mutate(exon = exonPosFun_v(pos, exons)) %>%
    left_join(exonList, 'exon') %>%
    dplyr::rename(exon.gc = GC) %>%
    mutate(exon.gc2 = exon.gc^2) %>%
    mutate(smooth.gc = get_gc_prediction(pos))  %>%
    mutate(smooth.gc2 = smooth.gc^2)


  gc.lm = lm(Ratio ~ exon.gc + exon.gc2 + smooth.gc + smooth.gc2,
             y = TRUE, data = tumour.logR)

  # Look at if interested
  # summary(gc.lm)

  tumour.logR$Ratio.gc.correct <- gc.lm$residuals
  return(tumour.logR)

}

# Calculating GC values ----
# For TCRA the GC values are calculated as follows
TCRA.gc.df <- slidingwindowplot_alt(1000,TCRA.fasta[[1]])
# plus 21999999 as TCRA fasta begins at position 21999999 and needs to be corrected
TCRA.gc.df$pos <- TCRA.gc.df$loc + 21999999
gam.model <- mgcv::gam(GC~s(pos, bs = 'cs'), data = TCRA.gc.df)

# This is the function to get the smoothed GC content at any position
get_gc_prediction <- function(x){mgcv::predict.gam(gam.model, newdata = data.frame(pos = x))}

# Again subtract to make compatible with FASTA file
TCRA.exons.loc <- list()
for(i in 1:192){
  TCRA.exons.loc[[i]] <- c(TCRA.exons$X2[i] - 21999999,
                           TCRA.exons$X3[i] - 21999999)

}
# These are the GC content within the exons
TCRA.exon.gc2 <- exonwindowplot2(TCRA.exons.loc, TCRA.fasta[[1]],0)

# Repeat as necessary for nimblegen exons
TCRA.exons.ng.loc <- list()
for(i in 1:length(TCRA.exons.nimblegen$X2)){
  TCRA.exons.ng.loc[[i]] <- c(TCRA.exons.nimblegen$X2[i] - 21999999,
                              TCRA.exons.nimblegen$X3[i] - 21999999)

}
TCRA.exon.ng.gc2 <- exonwindowplot2(TCRA.exons.ng.loc, TCRA.fasta[[1]],0)


# Updated functions:
# Simultaneous - 1000 points altogehter
getVDJFraction_upd2 <- function(tumour.logR, segs,norm.ci.95, GC.correct = FALSE, ci_type = 'simultaneous'){
  # uses con
  if(GC.correct){
    tumour.genomic.region.model <- mgcv::gam(Ratio.gc.correct~s(pos, bs = 'cs'), data = tumour.logR)
  }else{
    tumour.genomic.region.model <- mgcv::gam(Ratio~s(pos, bs = 'cs'), data = tumour.logR)
  }

  focal.start <- segs[2,2]
  focal.end <- segs[2,3]

  fit.model <- confint(tumour.genomic.region.model, parm = "s(pos)", partial_match = TRUE, type = ci_type,
                       newdata = seq(segs[2,2], segs[2,3],by=100), shift = TRUE)
  fit.loc <- which(fit.model$pos > focal.start & fit.model$pos < focal.end)
  if(length(fit.loc) == 0){
    fit.loc <- c(which(fit.model$pos > focal.start)[1]-1,
                 which(fit.model$pos > focal.start)[1])}

  tumour.log2.score <- mean(fit.model$est[fit.loc])
  tumour.log2.score.min <- mean(fit.model$lower[fit.loc]) - norm.ci.95
  tumour.log2.score.max <- mean(fit.model$upper[fit.loc]) + norm.ci.95

  tumour.tcell.purity <- 1 - (2^tumour.log2.score)
  tumour.tcell.purity.max <- 1 - (2^tumour.log2.score.min)
  tumour.tcell.purity.min <- 1 - (2^tumour.log2.score.max)
  return(c(tumour.tcell.purity, tumour.tcell.purity.min, tumour.tcell.purity.max))
}




# Running TIL ExTRECT -----
# To extend for non TCRA genes/locations, the segments need to be in the seg list
# and the exons as well as their GC values used and calculated.

# Needed on CAMP:
samtools.load <- 'module load SAMtools/1.3.1-foss-2016b'

# Data set specific input
all.bam.files <- c('test1.bam','test2.bam') # a vector of your bam files
regions <- c('test1','test2') # names of each case
outDir <- '/DIR/TO/OUT/' # an output directory to save coverage files

# For loop for different VDJ genes but currently only TCRA
# The next section at some point will be simplified with more functions
for(i in 1:length(seg.list)){
  vdj.gene <- names(seg.list)[i]
  vdj.seg <- seg.list[[i]]

  print(paste0('Running for VDJ gene ', vdj.gene))
  # The cov.output file is where a temp file of the coverage within TCRA is saved
  cov.output.files <- paste0(outDir, regions, '_', vdj.gene, '.txt')

  # Extract chromosome + info for samtools
  vdj.chr <- vdj.chr.df$chr[which(vdj.chr.df$gene == vdj.gene)]
  vdj.start <- vdj.seg[1,2]
  vdj.end <- vdj.seg[1,3]

  # Extract coverage from bams and write to file
  samtools.cmds <- paste0(samtools.load,'; samtools depth ',all.bam.files,
                          ' -q 20 -Q 20 -r ',vdj.chr,':', vdj.start,'-',vdj.end,' > ',
                          cov.output.files)

  sapply(samtools.cmds, system)

  # Create two lists for output and loop through all coverage files created
  vdj.fraction <- list()
  vdj.fraction.gc <- list()
  for(j in c(1:length(cov.output.files))){
    print(j)
    print(cov.output.files[j])

    # Read in the coverage of TCRA from one bam file
    vdj.region.df <- read_tsv(cov.output.files[j], col_names = FALSE)

    # Manually you can discard exons e.g. low coverage, process is now automated so should not be needed
    # Often used for TCGA runs where I have a list of exons which have a low coverage bias
    # set currently for use on the 192 TCRA sureselect exons
    exons.to.use <- seq(192)

    if(dim(vdj.region.df)[1] == 0){
      # If there is no coverage for whatever reason
      vdj.fraction[[j]] <- NA
      vdj.fraction.gc[[j]]  < NA
    }else{

      vdj.region.df <- vdj.region.df[,c(2,3)]
      colnames(vdj.region.df) <- c('pos','reads')

      # Use standard agilent tracerx exons - remember to change if needed
      exons.selected <- TCRA.exons
      # Load exon content
      exons.gc.content <- TCRA.exon.gc2

      # Filter for positions within exons as expected and apply median filter and remove any exons required
      vdj.region.df.filt.exons.median <- lapply(seq_len(dim(exons.selected)[1]), function(x){
        tmp <- vdj.region.df %>%
          filter(pos >= exons.selected$X2[x] & pos <= exons.selected$X3[x])
        tmp$reads <- medianFilter(tmp$reads, 50)
        return(tmp)})[exons.to.use]

      # Remove any exons with very low values (suspected exon failure - possible 100% T cell)
      # Threshold of < 15 may not be appropriate on low coverage data sets/genomic regions
      median.values.exons <- sapply(vdj.region.df.filt.exons.median,
                                    function(x) median(x$reads, na.rm = TRUE))
      exon.remove <- which(median.values.exons < 15 | is.na(median.values.exons))
      if(length(exon.remove) > 0){
        vdj.region.df.filt.exons.median <- vdj.region.df.filt.exons.median[-exon.remove]
      }
      vdj.region.df.filt.exons.median <- vdj.region.df.filt.exons.median %>%
        Reduce(rbind, .)

      # Get rid of any repeated rows that might exist - in theory I don't think they will
      vdj.region.df.filt.exons.median <- distinct(vdj.region.df.filt.exons.median)

      # Calculate logR based on the normalisation regions
      vdj.logR.df <- getLogRdf(vdj.region.df.filt.exons.median, vdj.seg, minCov = 0)

      #Remove any problematic values (NAs or infinite) - again I think these shouldn't exist
      vdj.logR.df <- vdj.logR.df[!is.infinite(vdj.logR.df$Ratio), ]
      vdj.logR.df <- vdj.logR.df[!is.na(vdj.logR.df$Ratio), ]

      # If lots of exons are being removed something is probably wrong and return NA
      if(dim(vdj.logR.df)[1] == 0 | length(exon.remove) > 30){
        # QC return NA if many exons fail or now there is absolutely no coverage
        vdj.fraction[[j]] <- NA
        vdj.fraction.gc[[j]] <- NA
      }else{
        # Run GC correction on logR
        vdj.logR.df <- GCcorrect(vdj.logR.df, exons = exons.selected, exonList = exons.gc.content)

        # GC correction sometimes moves the center away from 0 (needed for VDJ calculation)
        # adjust baseline by subtracting mean at start/end (works better than median)
        # Adjust baseline
        adjust.baseline.value <- vdj.logR.df %>%
          filter((pos >= vdj.seg[3,'start'] & pos <= vdj.seg[3,'end']) |
                   (pos >= vdj.seg[4,'start'] & pos <= vdj.seg[4,'end'])) %>%
          summarise(gc.adjust = mean(Ratio.gc.correct),
                    CI.95.range = 1.96*sd(Ratio.gc.correct)/sqrt(length(Ratio.gc.correct)))

        # Look at focal too in GAM model
        adjust.model <- mgcv::gam(Ratio.gc.correct~s(pos, bs = 'cs'), data = vdj.logR.df)
        adjust.fit.model <- confint(adjust.model, parm = "s(pos)",
                                    partial_match = TRUE, type = 'simultaneous',
                                    newdata = seq(vdj.seg[2,2], vdj.seg[2,3],by=100),
                                    shift = TRUE)
        fit.loc <- which(adjust.fit.model$pos > vdj.seg[2,2] & adjust.fit.model$pos < vdj.seg[2,3])
        adjust.baseline.value2 <- list(mean(adjust.fit.model$est[fit.loc]),
                                       mean(adjust.fit.model$upper[fit.loc]) - mean(adjust.fit.model$est[fit.loc]))

        adjust.value <- max(adjust.baseline.value[[1]], adjust.baseline.value2[[1]])
        ci.95.value <- ifelse(adjust.baseline.value2[[1]] > adjust.baseline.value[[1]],
                              adjust.baseline.value2[[2]], adjust.baseline.value[[2]])


        vdj.logR.df$Ratio.gc.correct <- vdj.logR.df$Ratio.gc.correct - adjust.value
        vdj.fraction.gc[[j]]  <- getVDJFraction_upd2(vdj.logR.df, vdj.seg,ci.95.value, TRUE)


        # Adjust baseline
        adjust.baseline.value <- vdj.logR.df %>%
          filter((pos >= vdj.seg[3,'start'] & pos <= vdj.seg[3,'end']) |
                   (pos >= vdj.seg[4,'start'] & pos <= vdj.seg[4,'end'])) %>%
          summarise(gc.adjust = mean(Ratio),
                    CI.95.range = 1.96*sd(Ratio)/sqrt(length(Ratio)))

        adjust.model <- mgcv::gam(Ratio~s(pos, bs = 'cs'), data = vdj.logR.df)
        adjust.fit.model <- confint(adjust.model, parm = "s(pos)",
                                    partial_match = TRUE, type = 'simultaneous',
                                    newdata = seq(vdj.seg[2,2], vdj.seg[2,3],by=100),
                                    shift = TRUE)
        fit.loc <- which(adjust.fit.model$pos > vdj.seg[2,2] & adjust.fit.model$pos < vdj.seg[2,3])
        adjust.baseline.value2 <- list(mean(adjust.fit.model$est[fit.loc]),
                                       mean(adjust.fit.model$upper[fit.loc]) - mean(adjust.fit.model$est[fit.loc]))

        adjust.value <- max(adjust.baseline.value[[1]], adjust.baseline.value2[[1]])
        ci.95.value <- ifelse(adjust.baseline.value2[[1]] > adjust.baseline.value[[1]],
                              adjust.baseline.value2[[2]], adjust.baseline.value[[2]])


        vdj.logR.df$Ratio<- vdj.logR.df$Ratio- adjust.value

        vdj.fraction[[j]]  <- getVDJFraction_upd2(vdj.logR.df, vdj.seg,ci.95.value, FALSE)

      }
    }
  }
}



TCRA.out.df <- data.frame(sample = regions,
                          T.cell.fraction.raw = sapply(vdj.fraction,`[`,1),
                          T.cell.fraction.raw.lwr = sapply(vdj.fraction,`[`,2),
                          T.cell.fraction.raw.upr = sapply(vdj.fraction,`[`,3),
                          T.cell.fraction.raw.gc = sapply(vdj.fraction.gc,`[`,1),
                          T.cell.fraction.raw.gc.lwr = sapply(vdj.fraction.gc,`[`,2),
                          T.cell.fraction.raw.gc.upr = sapply(vdj.fraction.gc,`[`,3))

# Adjusting for tumour purity/local TCRA copy number ----



# The following code:
# 1. joins purity and copy number data
# 2. Marks if it is a germline sample
# 3. If it is a germline sample, set purity = 0 and cnTotal = 2
# 4. If raw value is < 0 set to 0
# 5. Calculate the raw log ratio (used in calculation later)
# 6. Calculate the adjusted T cell fraction (see paper)
# 7. If > 1 - purity set it to this max value
# 8. Repeat 4-7 for GC corrected version.



TCRA.out.df   <- TCRA.out.df   %>%
  left_join(purity.df, 'sample') %>%
  left_join(seg.df, 'sample') %>%
  mutate(germline = grepl('GL',sample)) %>%
  mutate(purity = ifelse(germline, 0, purity)) %>%
  mutate(cnTotal = ifelse(germline, 2, cnTotal)) %>%
  mutate(rawRatio = 1-T.cell.fraction.raw) %>%
  mutate(rawRatio.lwr = 1-T.cell.fraction.raw.lwr) %>%
  mutate(rawRatio.upr = 1-T.cell.fraction.raw.upr) %>%
  mutate(TCRA.tcell.fraction= 1 - ((1-purity+(purity*cnTotal)/2)*rawRatio) - purity + ((purity*cnTotal)/2)) %>%
  mutate(TCRA.tcell.fraction,lwr= 1 - ((1-purity+(purity*cnTotal)/2)*rawRatio.lwr) - purity + ((purity*cnTotal)/2)) %>%
  mutate(TCRA.tcell.fraction.upr= 1 - ((1-purity+(purity*cnTotal)/2)*rawRatio.upr) - purity + ((purity*cnTotal)/2)) %>%
  mutate(maxPossible = 1- purity) %>%
  mutate(highTcellFlag = TCRA.tcell.fraction.lwr > maxPossible) %>%
  mutate(TCRA.tcell.fraction = ifelse(highTcellFlag.gc, maxPossible, T.cell.fraction.raw)) %>%
  mutate(TCRA.tcell.fraction.lwr = ifelse(highTcellFlag.gc, maxPossible, T.cell.fraction.raw.lwr)) %>%
  mutate(TCRA.tcell.fraction.upr = ifelse(highTcellFlag.gc, maxPossible, T.cell.fraction.raw.upr)) %>%
  mutate(rawRatio.gc = 1-T.cell.fraction.raw) %>%
  mutate(rawRatio.gc.lwr = 1-T.cell.fraction.raw.gc.lwr) %>%
  mutate(rawRatio.gc.upr = 1-T.cell.fraction.raw.gc.upr) %>%
  mutate(TCRA.tcell.fraction.gc= 1 - ((1-purity+(purity*cnTotal)/2)*rawRatio.gc) - purity + ((purity*cnTotal)/2)) %>%
  mutate(TCRA.tcell.fraction.gc.lwr= 1 - ((1-purity+(purity*cnTotal)/2)*rawRatio.gc.lwr) - purity + ((purity*cnTotal)/2)) %>%
  mutate(TCRA.tcell.fraction.gc.upr= 1 - ((1-purity+(purity*cnTotal)/2)*rawRatio.gc.upr) - purity + ((purity*cnTotal)/2)) %>%
  mutate(highTcellFlag.gc = TCRA.tcell.fraction.gc.lwr > maxPossible) %>%
  mutate(TCRA.tcell.fraction.gc = ifelse(highTcellFlag.gc, maxPossible, T.cell.fraction.raw.gc)) %>%
  mutate(TCRA.tcell.fraction.gc.lwr = ifelse(highTcellFlag.gc, maxPossible, T.cell.fraction.raw.gc.lwr)) %>%
  mutate(TCRA.tcell.fraction.gc.upr = ifelse(highTcellFlag.gc, maxPossible, T.cell.fraction.raw.gc.upr))


