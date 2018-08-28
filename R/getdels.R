#' Extract putative CNVs from a set of segments
#'
#' This function identifies a set of segments representing putative CNVs (
#' homozygous deletions, hemizygous deletions or duplications) from a set of
#' segments and associated mean log2 ratio values, and outputs them as a
#' \code{data.frame} along with some metadata describing them.
#'
#' @param data a data.frame representing a set of segments.
#' @param ind a character vector of length one. The sample (individual) for
#'   which to extract the putative CNVs.
#' @param all_segments logical of length one (\code{TRUE} or \code{FALSE}).
#'   Should all segments be extracted and output? Defaults to \code{FALSE},
#'   meaning that only segments representing putative CNVs are output.
#' @param homdel a single numeric value. The threshold to be used for calling
#'   homozygous deletions; all segments with a mean log2 ratio less than this
#'   values is considered called as a homozygous deletion.
#' @param hetdel a single numeric value. The threshold to be used for calling
#'   hemizygous deletions; all segments with a mean log2 ratio less than this
#'   value and great than or equal to \code{homdel} is called as a hemizygous
#'   deletion.
#' @param dup a single numeric value. The threshold to be used for calling
#'   duplications; all segments with a mean log2 ratio larger than this value
#'   is called as a duplication.
#'
#' @return a data.frame of the segments and associated metadata.
#' @export
#'
#' @examples
#' NULL
getdels <- function(data, ind, all_segments = FALSE,
                    homdel = -2, hetdel = -0.2, dup = 0.2) {

  # Extracting the different segment values from the data
  segments <- getseg(data, ind)

  # A different procedure is used whether all segments are seeked or only those
  #  with values warranting the call of a CNV
  if(!all_segments) {
    # Keeping only segments that are at least hemizygous deletions or duplications
    which_signif <- which(segments$mean_log2 < hetdel | segments$mean_log2 > dup)
  } else {
    # Otherwise all segments are kept (from the first to the last)
    which_signif <- 1:nrow(segments)
  }

  # An empty data.frame is returned if not segment is significant
  if(!length(which_signif)) {
    return(data.frame(ind = character(),
                      chr = character(),
                      start = numeric(),
                      end = numeric(),
                      kbp = numeric(),
                      length = numeric(),
                      mean_log2 = numeric(),
                      type = character(),
                      stringsAsFactors = FALSE))
  }

  # The output is initialized from the significant segments
  output <- segments[which_signif, ]

  # The individual will be the same for all rows
  output[["ind"]]   <- ind
  # Initializing columns with coherent classes
  output[["chr"]]   <- character(nrow(output))
  output[["start"]] <- numeric(nrow(output))
  output[["end"]]   <- numeric(nrow(output))

  # Iterating over all the significant segments
  for(i in 1:length(which_signif)) {
    # Extracting the data for this particular segment
    i_segment <- which_signif[i]
    # Retrieving the indices of the event boundaries in the log2 ratio dataset
    index_1 <- sum(segments[1:i_segment, "length"]) - segments[i_segment, "length"] + 1
    index_2 <- index_1 + segments[i_segment, "length"] - 1
    # A subset of the log2 ratio data containing the extent of the event
    subframe <- data[index_1:index_2, ]

    # A safety check to ensure that all datapoints are on the same chromosome
    stopifnot(length(unique(subframe$chr)) == 1)


    output[i, "chr"] <- subframe[["chr"]][1]

    # This defines a lower boundary on the size of the deletion
    # This boundary is sufficient for verifying the overlap between CGH and GBS
    # Upper boundaries will also be computed for the publication of the results
    output[i, "start"] <- subframe[["upper"]][1]
    output[i, "end"] <- subframe[["lower"]][nrow(subframe)]
  }

  output[["kbp"]]  <- (output[["end"]] - output[["start"]]) / 1000

  # Initializing a "type" column indicating the event type
  output[["type"]] <- "none"
  # Finding those that are duplications
  output[["type"]][output[["mean_log2"]] > dup] <- "dup"
  # Hemizygous deletions are first called, then homozygous deletions
  output[["type"]][output[["mean_log2"]] < hetdel] <- "hetdel"
  output[["type"]][output[["mean_log2"]] < homdel] <- "homdel"
  # Reordering and returning the output
  output[ , c("ind", "chr", "start", "end", "kbp", "length", "mean_log2", "type")]
}
