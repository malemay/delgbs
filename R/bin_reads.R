#' Count reads in genomic bins
#'
#' A function to put every the position of every read listed in a
#' data.frame into its corresponding bin.
#'
#' @param reads a data.frame of read positions for a given individual
#' @param bins a list (length of which = number of chr) of breakpoints
#'
#' @return To be completed
#' @export
#'
#' @examples
#' NULL
#'
bin_reads <- function(reads, bins) {
  # Initializing a list of counts per bin per chromosome
  counts <- list()

  # Bin the reads chromosome-wise
  for(i in 1:length(bins)) {
    # Extract the read positions corresponding to chromosome i
    chr_reads <- reads[as.character(reads$chr) == names(bins)[i], ]
    # Counting the number of reads falling into each bin
    counts[[i]] <- table(cut(chr_reads[["pos"]], bins[[i]],
                             labels = as.character(bins[[i]][-1])))
  }

  #unlisting the results allows putting them as a column into a DF
  unlist(counts)
}
