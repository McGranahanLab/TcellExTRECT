
This software and associated documentation files (the “Software”) are protected by copyright. This Software is provided “as is” (at your own risk) for internal non-commercial academic research purposes only. Please read the Non-Commercial Academic License in detail before downloading a copy. By installing or using this Software, you agree to be bound by the terms and conditions of the Non-Commercial Academic License (included in the “non-commercial-academic-license.exe” file, available in the main directory of this software repository). 

All commercial use of the Software or any modification, manipulation or derivative of the Software, including but not limited to transfer, sale or licence to a commercial third party or use on behalf of a commercial third party (including but not limited to use as part of a service supplied to any third party for financial reward) is strictly prohibited and requires a commercial use licence. For further information please email commercial@cancer.org.uk

# T Cell ExTRECT version 1.0.0

T Cell Exome TREC Tool (T Cell ExTRECT) is an R package to calculate T cell fractions from WES data from hg19 or hg38 aligned genomes.


## Instalation guide

T cell ExTRECT can either be installed from github using the `install_github()` function in the package devtools or directly from source.

### Method 1 - Using devtools

```r
# install.packages('devtools')
library(devtools)
install_github("McGranahanLab/TcellExTRECT/")

```

### Method 2 - Downloaded from source


```r
install.packages('PATH/To/TcellExTRECT/', repos=NULL, type ='source')
```

## Requirements

Samtools (>v1.3.1)

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

The `cov_example` data frame is an inbuilt example of coverage reads from the *TCRA* locus. Different capture kits have different exon positions and not all capture kits cover the *TCRA* locus. The `TCRA_exons_hg19` data contains exon locations from the Agilent v4/5 exome capture kits for genomes aligned to hg10. The TcellExTRECT package also comes with data for hg38 Agilent v4/5 genomes (`TCRA_exons_hg38`) and an exon set for Nimblegen kits `TCRA_exons_nimblegen_hg19`, more will be added soon or can be requested!

The output of `runTcellExTRECT` is a data frame containing 5 columns:

* `sample` - the name of the sample run
* `TCRA.tcell.fraction` - The estimated fraction of the sample that are T cells
* `TCRA.tcell.fraction.lwr` and `TCRA.tcell.fraction.upr` - The 95\% confindence interval values of the estimated T cell fraction
* `qcFit` - A QC value describing the noise in the fit of the GAM model used to estimate the T cell fraction. Minimum value is `1` which represents a sample with a clear signal and little noise, if this value is `>4` a warning message will be output and viewing of the GAM model using `plotTcellExTRECT` (see below) is recommended.

We can also visualise the calculated log ratio from T Cell ExTRECT, this can be very useful to check that everything is working. The following produces the pre and post GC corrected versions of the log ratio within the *TCRA* loci, reads are coloured by the class of VDJ segment they are, e.g. *TCRA-V* segments are blue. Note that the TCRA loci also includes segments related to TCR delta (*TCRD* or *TRD*), e.g. *TCRD-V*. Vertical dotted lines represent the regions of the genome used for the normalised baseline (Norm region start and Norm region end) as well as the 'Focal region' that is used in the calculation of the *TCRA* fraction as the location we expect to see maximum signal.

If there is a signal coming from T cells you expect to see a dip in the Log2 Read Depth Ratio around the focal region. A quick visual QC check is to see if there are any exons or reads that are outliers and potentially interfering with the calculation of the T cell fraction.

```
plotTcellExTRECT(cov_example, TCRA_exons_hg19,
                tcra_seg_hg19,'hg19', sample_name = 'TEST')
```

This will plot both the pre and post GC corrected versions. 

As a final step for cancer samples with known purities and copy number states around the TCRA loci, the TCRA T cell fractions can be adjusted with the following function:

```
TCRA.out <- adjustTcellExTRECT(TCRA.out, purity = 0.5, TCRA.cn = 3)
TCRA.out
```

Finally, if the raw data used to produce the scores or the plots is required it can be extracted with the following function:

```
TCRA.df <- exonsTcellExTRECT(cov_example, TCRA_exons_hg19, tcra_seg_hg19, 'hg19',GC_correct = TRUE,
                             summariseExons = FALSE)

```

The option `summariseExons = TRUE` will return a smaller data frame summarised using the median value over the input exons. This function is intended to be used for additional QC checks on the exons used but could also be used for manual plotting.


