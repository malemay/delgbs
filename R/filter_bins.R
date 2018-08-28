#' Filtering bins according to user-defined thresholds
#'
#' This function takes a \code{data.frame} of normalized read counts per sample
#' and per bin and removes the bins considered uninformative for the purposes
#' of CNV calling according to three different criteria: the minimum mean number
#' of reads, the maximum mean number of reads, and the ratio of the variance
#' in the number of reads per sample to the mean number of reads per sample (
#' overdispersion filter)
#'
#' @param read_counts a \code{data.frame} of read counts as generated from
#'   function \code{link{read_counts}}.
#' @param max_reads a single numeric value. The maximum mean number of reads per
#'   individual that a bin can have to be kept for further processing.
#' @param min_reads a single numeric value. The minimum mean number of reads per
#'   individual that a bin must have to be kept for further processing.
#' @param odisp_filter a single numeric value. The maximum ratio of the variance
#'   of the number of reads per individual in a bin to the mean number of reads
#'   per individual in this bin. The purpose of this filter is to remove bins
#'   that display overdispersion. A value of 3 has commonly been used for this
#'   filter.
#'
#' @return a \code{data.frame} of read counts similar to that given as input,
#'   but with bins removed according to the specified filters
#' @export
#'
#' @examples
#' NULL
filter_bins <- function(read_counts, max_reads, min_reads, odisp_filter) {
  means <- apply(read_counts[, -c(1:3)], 1, mean)
  vars <-  apply(read_counts[, -c(1:3)], 1, stats::var)
  read_counts[(means <= max_reads) & (means >= min_reads) & ((vars / means) <= odisp_filter), ]
}
