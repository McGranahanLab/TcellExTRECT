#' Function to create TCRA coverage file from bams
#'
#' @param bamPath Path to bam file
#' @param outPath Path to output directory
#' @param vdj.seg Location of TCRA file
#' @name getCovFromBam
#' @export


getCovFromBam <- function(bamPath, outPath, vdj.seg){
  # Requires samtools to be installed and working!
  if(!(file.exists(bamPath))){
    stop('Can not find bam file')
  }

  vdj.start <- vdj.seg[vdj.seg$segName == 'all','start']
  vdj.end <- vdj.seg[vdj.seg$segName == 'all','end']

  cov.name <- gsub('.bam','',basename(bamPath))
  cov.output.files <- paste0(outPath,cov.name, '_TCRA.txt')

  # Check does bai file exist
  bai.path <- paste0(bamPath,'.bai')
  bai.path2 <- paste0(gsub('.bam','',bamPath),'.bai')
  if(!(file.exists(bai.path) | file.exists(bai.path2))){
    stop('No index bai file found for bam, please index first before proceeding')
  }

  # Check if bam has chr or not before
  samtools.chr.check <- paste0('samtools idxstats ', bamPath,' | head -n 2')
  chr.check.output <- system(samtools.chr.check, intern = TRUE, ignore.stderr = TRUE)[2]
  chr.check.output2 <- strsplit(chr.check.output,'\t')[[1]][1]
  chr.present <- grepl('chr',chr.check.output2)
  chr14 <- ifelse(chr.present, 'chr14:','14:')


  samtools.cmds <- paste0('samtools depth ',bamPath,
                          ' -q 20 -Q 20 -r ',chr14, vdj.start,'-',vdj.end,' > ',
                          cov.output.files)

  # Replace with processx!
  sapply(samtools.cmds, system)

  return(cov.output.files)
}
