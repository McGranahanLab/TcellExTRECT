#' TCRA fasta file
#'
#' Sequence around the TCRA locus from fasta file
#' chr14:22000000-23110000 (hg19)
#' chr14:???(hg38)
#'
#' @format A list containing a SeqFastadna object
#' @return NA
"TCRA_fasta"


#' TCRA segments file (hg19)
#'
#' Segments used for running T cell ExTRECT, defining TCRA loci,
#' local regions used for normalisation and region of maximum VDJ
#' recombination (hg19 version)
#'
#' @format A list containing regions used to run T cell ExTRECT
#' @return NA
"tcra_seg_hg19"

#' TCRA segments file (hg38)
#'
#' Segments used for running T cell ExTRECT, defining TCRA loci,
#' local regions used for normalisation and region of maximum VDJ
#' recombination (hg38 version)
#'
#' @format A list containing regions used to run T cell ExTRECT
#' @return NA
"tcra_seg_hg38"


#' TCRA exons (hg19)
#'
#' Exons used in Agilent v4/5 capture kit within TCRA loci, extracted
#' from bed file (hg19 version)
#'
#' @format A data frame of exon locations
#' @return NA
"TCRA_exons_hg19"


#' TCRA exons (hg38)
#'
#' Exons used in Agilent v4/5 capture kit within TCRA loci, extracted
#' from bed file (hg38 version)
#'
#' @format A data frame of exon locations
#' @return NA
"TCRA_exons_hg38"

#' TCRA segments
#'
#' Data frame of the VDJ gene segments within the TCRA locu
#'
#' @format A data frame of TCRA segments
#' @return NA
"TCRA_segments"

#' TCRA exons nimblegen (hg19)
#'
#' Exons used in Nimblegen capture kit within TCRA loci, extracted
#' from bed file (hg19 version)
#'
#' @format A data frame of exon locations
#' @return NA
"TCRA_exons_nimblegen_hg19"

#' Example coverage file
#'
#' Coverage values of the TCRA gene from the TRACERx100 CRUK0023 R1 sample (hg19 format)
#'
#' @format A data frame of coverage values within TCRA gene
#' @return NA
"cov_example"
