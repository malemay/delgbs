#' Count the number of reads falling into genomic bins
#'
#' This function finds the number of reads falling into each bin from a
#' specified set. The aligned reads are read directly from .bam files matching
#' a pattern.
#'
#' @param patterns a character vector indicating patterns which will be used
#'   to look ofr .bam files.
#' @param binlist a named list of numeric vectors of breakpoints for the bins
#'   that will be used to tally the reads. The names of the list elements should
#'   be chromosome names matching the names used in the .bam file, whereas the
#'   numeric vectors should start at position 0 and identify the limits of the
#'   bins.
#' @param minq an integer or numeric of length 1. The minimum mapping quality
#'   for a read to be used for binning.
#' @param binwidth a integer or numeric of length 1. The width of the bins used.
#'   This value is not dynamically obtained from the
#' @param strand_diff a logical of length 1 (\code{TRUE} or \code{FALSE}).
#'   If \code{TRUE}, reads mapping to the "+" strand are tallied according to
#'   their 3'-most position, whereas reads mapping to the "-" strand are tallied
#'   according to their 5'-most position. If \code{FALSE}, all reads are tallied
#'   according to their 5'-most position. This option should be chosen according
#'   to the expected mean size of the restriction fragments, so as to maximize
#'   the probability that two reads originating from the same restriction
#'   fragment but with differing orientation are counted in the same bin.
#'
#' @return a \code{data.frame} with the first three columns giving the position
#'   of the bin and the remaining columns corresponding to the
#' @export
#'
#' @examples
#' NULL
count_reads <- function(patterns, binlist, minq, binwidth, strand_diff = FALSE) {

  ### +++ Initializing the bin data.frame
  # Initializing a data frame which will hold one kbp bin per row
  # The first element of each breakpoints vector is removed ; this
  #  represents the upper bound of each bin
  bin_df <- data.frame(upper = unlist(lapply(binlist, `[`, -1)))
  # The lower bound is computed by subtracting the binwidth from the upper bound
  bin_df[["lower"]] <- bin_df[["upper"]] - binwidth
  # Add chromosome names by repeating chr names as many times as there
  #  are bins in every chromosome (number of breakpoints minus 1)
  bin_df[["chr"]] <- rep(names(binlist), lengths(binlist) - 1)
  # Reordering the column names of the bin data frame
  bin_df <- bin_df[ , c("chr", "lower", "upper")]
  # Chromosome name is to be considered as a factor
  bin_df[["chr"]] <- as.factor(bin_df[["chr"]])
  # Row names are removed as they are irrelevant
  rownames(bin_df) <- NULL
  ### ---

  # Reading the data from all the individuals one after the other
  for (i in patterns) {
    # Keeping only the .bam file (not .bai) that matches the individual i
    filename <- dir(pattern = paste0(i, ".*bam$"))

    # Functionality should be added where several files belong to the same
    #  individual (see call_deletions.R for a solution to this)
    # For the moment let's just throw an error
    stopifnot(length(filename) == 1)

    # Reading the bam file with the specified mapping quality filter
    alignments <- get_pos(filename, minq = minq, strand_diff = strand_diff)

    # Binning the reads
    bin_df[[i]] <- bin_reads(alignments, binlist)
    cat(i, "processed.\n")
  }

  # Removing bins without any reads
  empty_bins <- apply(bin_df[ , -c(1:3)], 1, sum) == 0
  bin_df <- bin_df[!empty_bins, ]

  return(bin_df)
}
