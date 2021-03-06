% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/RGTvcftools.R
\name{readVariantsFromFile}
\alias{readVariantsFromFile}
\title{Read variants from disk}
\usage{
readVariantsFromFile(ff, n = 4096, ignorelines = c("thesaurushard",
  "thesaurusmany"), getcolumns = 10, withindels = FALSE)
}
\arguments{
\item{ff}{filename}

\item{n}{integer, determines how many lines are read from the input to scan for the
vcf header. Set this to a moderatly large number, i.e. larger than the expected number
of header lines in the vcf file.}

\item{ignorelines}{vector of filter codes. Lines in the vcf that contain one
of these filter codes will be eliminated, i.e. not reported in the output}

\item{getcolumns}{number of columns to read/return. Use this to select only a few of the columns of
the vcf table (ignoring the info or format fields can save memory)}

\item{withindels}{logical. Set to TRUE to return all variants in input vcf, or
set to FALSE to obtain only single-nucleotide substitutions,}
}
\description{
Read variants from disk and perform some pre-processing. In contrast to
obtaining a data.frame using read.table(), this function does not crash
when a vcf file does not have any valid rows (i.e. no variants detected).
This function also automatically allows the user to discard some columns
of the input table.
}

