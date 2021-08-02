# Tcell ExTRECT new template

# To add to data:
# Exons not to use in agilent v2.

library(TcellExTRECT)

data(tcra_seg_hg19)
data(TCRA_exons_hg19)
data(cov_example)

example_bam <- '/path/to/file.bam'

cov.file <- getCovFromBam(bamPath = example_bam,
                          outPath = '',
                          vdj.seg = tcra_seg_hg19)

cov_df <- loadCov(cov.file)

TCRA.out <- runTcellExTRECT(cov_df, TCRA_exons_hg19, tcra_seg_hg19, 'hg19')


plotTcellExTRECT(cov_example, TCRA_exons_hg19,
                tcra_seg_hg19,'hg19', sample_name = 'TEST')


adjustTcellExTRECT(TCRA.out, purity = 0.5, TCRA.cn = 3)



