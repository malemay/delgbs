#' Segment the log2 ratios for a set of individuals
#'
#' @param data_matrix
#' @param info_indices
#' @param chrom_col
#' @param pos_col
#' @param threshold_param
#' @param lambda_param
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
segment <- function(data_matrix, info_indices, chrom_col, pos_col,
                    threshold_param, lambda_param = 2, verbose = TRUE) {

  # Creating a data.frame that will hold the segment values at the end
  segment_frame <- data_matrix
  segment_frame[ , -info_indices] <- NA

  # Checking that all chromosomes are grouped
  stopifnot(length(rle(segment_frame[[chrom_col]])$lengths) ==
              length(unique(segment_frame[[chrom_col]])))
  # Checking that all positions are in increasing order
  stopifnot(all(tapply(data_matrix[[pos_col]],
                       data_matrix[[chrom_col]],
                       function(x) all(diff(x) >= 0))))

  # We are going to perform the segmentation for each individual
  for(ind in names(data_matrix)[-info_indices]) {

    # The standard deviation estimation is done for all chromosomes at once
    ## This differs from what was originally done by Guillem
    ## He would compute the standard deviation for each chromosome
    ## Doing it individual-wise makes more sense to me
    sd_estimate <- var_diff(data_matrix[[ind]])
    # The threshold for segmentation depends on this estimate
    threshold <- threshold_param * sd_estimate

    # Creating a list holding the values for every chromosome
    segment_list <- list()

    # Looping over all the chromosomes
    for(chr in unique(data_matrix[[chrom_col]])) {

      # Extracting the log2 ratio values for this chromosome
      log2_values <- data_matrix[[ind]][data_matrix[[chrom_col]] == chr]
      # Setting the segmentation lamba parameter (depends on chr size)
      lambda <- lambda_param * log(length(log2_values)) * sd_estimate
      # Actually performing the segmentation and storing the values in the list
      segment_list[[chr]] <- robseg::Rob_seg.std(x = log2_values,
                                                 loss = "Outlier",
                                                 lambda = lambda,
                                                 lthreshold = threshold)$smt
    }

    # Values are inserted in the data.frame once all chromosomes are processed
    segment_list <- unlist(segment_list)
    stopifnot(length(segment_frame[[ind]]) == length(segment_list))
    segment_frame[[ind]] <- segment_list

    # Informing the user that this individual has been processed
    if(verbose) message("Individual ", ind, " segmented.")
  }


  # Return a data.frame with the segment value associated with each data point
  #  for each individual
  stopifnot(all(stats::complete.cases(segment_frame))) # no NA values supported
  return(segment_frame)
}

#' Difference based estimation of the variance
#'
#' The value obtained from this function is used to set the parameter of the
#' segmentation (penalty and threshold).
#'
#' @param x a numeric vector of log2 ratio values.
#'
#' @return a numeric vector of length one. The estimated variance of the log2
#'   ratios.
#'
#' @examples
#' NULL
var_diff <- function(x) {
  n <- length(x)

  wei <- c(0.1942, 0.2809, 0.3832, -0.8582)

  mat <- wei %*% t(x)
  mat[2, -n] <- mat[2, -1]
  mat[3, -c(n-1, n)] <- mat[3, -c(1, 2)]
  mat[4, -c(n-2, n-1, n)] <- mat[4, -c(1, 2, 3)]

  sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3))
}
