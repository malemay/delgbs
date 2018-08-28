#' Merge two read count datasets
#'
#' This function takes two read count datasets as output by function
#' \code{link{count_reads}} and merges them, adding new rows and computing
#' the sum of read counts in shared bins as needed. The output is sorted
#' according to the genomic position of the bins, as the inputs were.
#'
#' @param x,y the two read count \code{data.frame}'s to merge.
#'
#' @return a \code{data.frame} of read counts similar to those given as input,
#'   but that combines the data of the inputs.
#' @export
#'
#' @examples
#' NULL
merge_counts <- function(x, y) {
  # Checking that the column names are identical
  stopifnot(identical(names(x), names(y)))

  # A first raw combination of both and reordering thereof
  combined <- rbind(x, y)
  combined <- combined[order(combined$chr, combined$lower), ]

  # Identifying which are duplicated (need columns chr [[1]] and pos [[2]])
  which_duplicates <- which(duplicated(combined[, 1:2]))
  # We want a vector to extract ALL rows involved in duplicates
  which_duplicates <- c(which_duplicates, which_duplicates - 1)
  which_duplicates <- sort(which_duplicates)

  # Extracting the data.frame corresponding to these data
  #  and separating the (unique) metadata columns from the counts
  duplicate_counts <- combined[which_duplicates, ]
  duplicate_df <- duplicate_counts[seq(1, nrow(duplicate_counts), by = 2), 1:3]

  # Safety checks
  length(rle(as.character(duplicate_counts[["chr"]]))$lengths) == 20
  all(rle(duplicate_counts[["lower"]])$lengths == 2)
  all(rle(duplicate_counts[["upper"]])$lengths == 2)

  # We can use two matrices of summed counts (odd and even indices) and add them
  duplicate_counts <-
    as.matrix(duplicate_counts[seq(1, nrow(duplicate_counts), by = 2), -c(1:3)]) +
    as.matrix(duplicate_counts[seq(2, nrow(duplicate_counts), by = 2), -c(1:3)])

  # Merging the metadata and the counts together
  duplicate_df <- cbind(duplicate_df, duplicate_counts)

  # Combining the counts from duplicate and non-duplicated regions
  combined <- rbind(combined[-which_duplicates, ], duplicate_df)
  # Reordering the counts
  combined <- combined[order(combined$chr, combined$lower), ]
  # Output
  combined
}
