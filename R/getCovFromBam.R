#' Function to create TCRA coverage file from bams
#'
#' @param bamPath Path to bam file
#' @param outPath Path to output directory
#' @param vdj.seg Location of TCRA file
#' @name getCovFromBam
#' @export


getCovFromBam <- function(bamPath, outPath, vdj.seg){
  # Requires samtools to be installed and working!

  vdj.start <- vdj.seg[vdj.seg$segName == 'all','start']
  vdj.end <- vdj.seg[vdj.seg$segName == 'all','end']

  cov.name <- gsub('.bam','',basename(bamPath))
  cov.output.files <- paste0(outPath,cov.name, '_TCRA.txt')


  samtools.cmds <- paste0('samtools depth ',bamPath,
                          ' -q 20 -Q 20 -r ','chr14:', vdj.start,'-',vdj.end,' > ',
                          cov.output.files)

  # Replace with processx!
  sapply(samtools.cmds, system)

  return(cov.output.files)
}
