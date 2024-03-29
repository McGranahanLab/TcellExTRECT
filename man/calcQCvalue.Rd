% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcQCvalue.R
\name{calcQCvalue}
\alias{calcQCvalue}
\title{Method for calculating QC values to judge GAM fit}
\usage{
calcQCvalue(tumour.logR, segs, GC_correct = TRUE, ci_type = "simultaneous")
}
\arguments{
\item{tumour.logR}{LogR values}

\item{segs}{GC content of exons}

\item{GC_correct}{whether to use GC correction or not for output}

\item{ci_type}{FASTA file for VDJ gene e.g TCRA}
}
\value{
QC values for judging fit of GAM model
}
\description{
Method for calculating QC values to judge GAM fit
}
