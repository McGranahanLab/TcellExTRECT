# T Cell ExTRECT version 1.0.0

T Cell Exome TREC Tool (T Cell ExTRECT) is an R package to calculate T cell fractions from WES data from hg19 or hg38 aligned genomes.

## License

T Cell ExTRECT is free for academic use.

INSERT LICENSE

## Instalation guide

T cell ExTRECT can either be installed from github using the `install_github()` function in the package devtools or directly from source.

### Method 1 - Using devtools

```r
# install.packages('devtools')
library(devtools)
install_github("McGranahanLab/T-Cell-ExTRECT/")

```

### Method 2 - Downloaded from source


```r
install.packages('PATH/To/T-Cell=ExTRECT/', repos=NULL, type ='source')
```

## Requirements

Samtools (>v0)

## Example use
Running T cell ExTRECT on your data is both fast and easy!

```r
library(TcellExTRECT)
```

First take an aligned bam file (hg38 or hg19) of your choice

```r
example_bam <- '/path/to/file.bam'
```

Then use the pre-defined *TCRA* gene segments in the `tcra_seg_hg19` or `tcra_seg_hg38` data file to extract the coverage values. 

```r
data("tcra_seg_hg19")
```

```r
cov.file <- getCovFromBam(bamPath = example_bam,
                          outPath = '',
                          vdj.seg = tcra_seg_hg19)
cov_df <- loadCov(cov.file)
```

The `cov.df` object should be a data frame with two columns names `pos` (position on chr14) and `reads` (number of coverage reads). The function getCovFromBam calls samtools so that will need to be installed! Alternatively extract the coverage in the TCRA region using your own method and load before preceding.


Once the coverage values have been loaded you can run TCellExTRECT with the following function:

```r
data(TCRA_exons_hg19)
data(cov_example)

TCRA.out <- runTcellExTRECT(cov_example, TCRA_exons_hg19, tcra_seg_hg19, 'hg19')
TCRA.out
```

The `cov_example` data frame is an inbuilt example of coverage reads from the *TCRA* locus. Different capture kits have different exon positions and not all capture kits cover the *TCRA* locus. The `TCRA_exons_hg19` data contains exon locations from the Agilent v4/5 exome capture kits for genomes aligned to hg38. The TcellExTRECT package also comes with data for hg38 Agilent v4/5 genomes (`TCRA_exons_hg38`) and an exon set for Nimblegen kits `TCRA_exons_nimblegen_hg19`, more will be added soon!

We can also visualise the calculate log ratio from T Cell ExTRECT, this can be very useful to check that everything is working.

```r
plotTcellExTRECT(cov_example, TCRA_exons_hg19,
                tcra_seg_hg19,'hg19', sample_name = 'TEST')
```

This will plot both the pre and post GC corrected versions. 

Finally the TCRA T cell fractions can be adjusted for tumour copy number and purity with the following function:

```r
TCRA.out <- adjustTcellExTRECT(TCRA.out, purity = 0.5, TCRA.cn = 3)
TCRA.out
```


