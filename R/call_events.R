#' Call CNV events from a set of read counts
#'
#' This function performs all the steps needed to go from a \code{data.frame}
#' read counts per bin and per individual to a set of CNVs. The steps involved
#' are the normalization of the read counts per sample, the removal of bins
#' deemed uninformative for the purposes of CNV calling, the computation of
#' log2 ratios and segmentation based on these ratios, and finally CNV calling
#' based on user-defined log2 ratio thresholds.
#'
#' @param read_counts a \code{data.frame} of read counts as generated from
#'   function \code{link{read_counts}}.
#' @param min_reads a single numeric value. The minimum mean number of reads per
#'   individual that a bin must have to be kept for further processing.
#' @param max_reads a single numeric value. The maximum mean number of reads per
#'   individual that a bin can have to be kept for further processing.
#' @param seg_threshold a single numeric value. The threshold parameter used
#'   for segmentation. See function \code{\link{segment}} for more details.
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
#' @return a list of two elements :
#' - \code{segdels}: a \code{data.frame} containing the CNV events called and
#'   some metadata about them.
#' - \code{log2_ratios}: a \code{data.frame} of log2 ratios per bin and per
#'   sample. The first three columns indicate the genomic position of the
#'   bin, whereas each remaining column represents the data for a single
#'   sample (column names correspond to sample names).
#' @export
#'
#' @examples
#' NULL
call_events <- function(read_counts, min_reads, max_reads, seg_threshold,
                        homdel = -2.5, hetdel = -0.5, dup = 0.2) {

  # STEP 1: NORMALIZING
  # Normalizing the number of reads per individual
  ## Computing the overall mean first
  message("Normalizing the read counts...")
  mean_count <- mean(colSums(read_counts[, -c(1:3)]))
  ## Then normalizing per individual
  read_counts <- cbind(read_counts[ , 1:3],
                       apply(read_counts[ , -c(1:3)], 2,
                             function(x) x / sum(x) * mean_count))
  cat("Mean number of reads per individual:", mean_count, "\n")

  # STEP 2: FILTERING AND COMPUTING LOG2 RATIO
  # Filtering steps
  ## Filtering to remove bins with a very high mean number of reads
  message("Filtering the bins...")
  read_counts <- read_counts[apply(read_counts[ , -c(1:3)], 1, mean) <= max_reads, ]
  # Now applying the minimum read count filter
  read_counts <- read_counts[apply(read_counts[ , -c(1:3)], 1, mean) >= min_reads, ]
  ## Filtering to remove overdispersion
  means <- apply(read_counts[ , -c(1:3)], 1, mean)
  vars  <- apply(read_counts[ , -c(1:3)], 1, stats::var)
  read_counts <- read_counts[vars / means <= 3, ]
  # Computing the log2 ratio of these data
  log2_ratios <- cbind(read_counts[ , 1:3],
                       t(apply(read_counts[ , -c(1:3)], 1,
                               function(x) log2((x + 1) / mean(x + 1)))))
  log2_ratios$chr <- as.character(log2_ratios$chr)
  cat("Bins kept:", nrow(log2_ratios), "\n")

  # STEP 3: PERFORMING SEGMENTATION ON THE LOG2 RATIOS
  message("Segmenting the log2 ratios...")
  segments <- segment(log2_ratios, info_indices = 1:3,
                      chrom_col = "chr", pos_col = "lower",
                      threshold_param = seg_threshold, lambda_param = 2,
                      verbose = FALSE)
  segments$chr <- as.character(segments$chr)

  # STEP 4: CALLING CNVS FROM THE SEGMENTS
  message("Calling CNV events...")
  segdels <- list()

  # CNV calling is done separately for each individual
  for(i in names(segments)[-c(1:3)]) {
    segdels[[i]] <- delgbs::getdels(segments, ind = i, homdel = homdel,
                                    hetdel = hetdel, dup = dup)
  }

  segdels <- do.call("rbind", segdels)

  # Returning the raw deletions (still have to be filtered)
  cat("Number of homozygous deletions:", sum(segdels$type == "homdel"),"\n")
  cat("Number of hemizygous deletions:", sum(segdels$type == "hetdel"),"\n")
  cat("Number of duplications:", sum(segdels$type == "dup"), "\n")

  list(segdels, log2_ratios)
}
