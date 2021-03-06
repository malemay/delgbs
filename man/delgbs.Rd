% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/delgbs.R
\docType{package}
\name{delgbs}
\alias{delgbs}
\alias{delgbs-package}
\title{delgbs: Detection of deletions and duplications from GBS data}
\description{
This package was designed to allow the detection of copy number
variation (CNV) from standard genotyping-by-sequencing (GBS) data. Developed
with populations of fast neutron mutants in mind, it compares the depth of
sequencing in discrete genomic windows between individuals, segments the
sequencing depth profiles to identify homozygous deletions, hemizygous
deletions, and duplications. The functions included in this package provide
the utilities needed for going from a set of BAM alignment files to a set
of CNVs as well as high-level plotting functions to allow visualizing the
output.
}
