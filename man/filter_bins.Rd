% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_bins.R
\name{filter_bins}
\alias{filter_bins}
\title{Filtering bins according to user-defined thresholds}
\usage{
filter_bins(read_counts, max_reads, min_reads, odisp_filter)
}
\arguments{
\item{read_counts}{a \code{data.frame} of read counts as generated from
function \code{link{read_counts}}.}

\item{max_reads}{a single numeric value. The maximum mean number of reads per
individual that a bin can have to be kept for further processing.}

\item{min_reads}{a single numeric value. The minimum mean number of reads per
individual that a bin must have to be kept for further processing.}

\item{odisp_filter}{a single numeric value. The maximum ratio of the variance
of the number of reads per individual in a bin to the mean number of reads
per individual in this bin. The purpose of this filter is to remove bins
that display overdispersion. A value of 3 has commonly been used for this
filter.}
}
\value{
a \code{data.frame} of read counts similar to that given as input,
  but with bins removed according to the specified filters
}
\description{
This function takes a \code{data.frame} of normalized read counts per sample
and per bin and removes the bins considered uninformative for the purposes
of CNV calling according to three different criteria: the minimum mean number
of reads, the maximum mean number of reads, and the ratio of the variance
in the number of reads per sample to the mean number of reads per sample (
overdispersion filter)
}
\examples{
NULL
}
