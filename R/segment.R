#' Segment the log2 ratios for a set of individuals
#'
#' This function is a wrapper around the function \code{\link[robseg]{Rob_seg.std}}
#' of package \code{robseg} which segments all the log2 ratio profiles in the
#' \code{data.frame} provided as input. This wrapper implements a set of
#' parameters is is very narrow compared to the complete functionality of
#' function  \code{\link[robseg]{Rob_seg.std}}. Notably, the loss function used
#' here is the one that is robust to outliers as the cost function is bounded
#' at a certain threshold. More information can be found on the \code{robseg}
#' package's Github page and in the research paper describing the method (see
#' the "see also")
#'
#' @param data_matrix a \code{data.frame} containing optional metadata columns
#'   and one or several columns of log2 ratios corresponding to individuals
#'   after which the columns are named.
#' @param info_indices an integer vector. The indices of the metadata columns,
#'   such that these columns are not segmented. Typically, this will be
#'   \code{1:3} as the first three columns will be the chromosome and the lower
#'   and upper bounds of the genomic bins.
#' @param chrom_col a character vector of length one. The name of the column
#'   containing the names of the chromosomes.
#' @param pos_col a character vector of length one. The name of a column giving
#'   the genomic position of the bin (either lower or upper). This argument is
#'   only used to check that the bins are in increasing order in the
#'   \code{data.frame}, as this is needed for segmentation. It therefore does
#'   not matter whether this position is the lower or upper breakpoint of the
#'   bin.
#' @param threshold_param a single numeric value. The multiplier that is applied
#'   to the standard deviation estimate to get the threshold that bounds the
#'   increase of the cost function. Higher values make the segmentation more
#'   sensitive to outliers, whereas lower values make it more robust to outliers.
#' @param lambda_param a single numeric value. The penalty value used by the
#'   segmentation. Defaults to two.
#' @param verbose a single logical value (\code{TRUE} or \code{FALSE}). Whether
#'   the progress of the segmentation should be printed. Defaults to \code{TRUE}.
#'
#' @return To be completed.
#' @export
#'
#' @section Source:
#' Link to \code{robseg}'s Github page: \url{https://github.com/guillemr/robust-fpop}
#'
#' Original description of the segmentation approach:
#' Fearnhead, P., & G. Rigaill (2018) \emph{Changepoint Detection in the
#' Presence of Outliers}. Journal of the American Statistical Association.
#' DOI: 10.1080/01621459.2017.1385466
#'
#' @examples
#' NULL
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
