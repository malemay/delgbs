#' Extract the segments from a given individual
#'
#' @param data a \code{data.frame} of mean log2 ratio values per bin, following
#'   segmentation with function \code{\link{segment}}.
#' @param ind a character vector of length one. The individual for which the
#'   data should be extracted (the name must correspond to one of the columns
#'   of \code{data}).
#'
#' @return a \code{data.frame} with the segments of a given individual.
#'
#' @examples
#' NULL
getseg <- function(data, ind) {
  rle_obj <- rle(data[[ind]])
  data.frame(length = rle_obj$lengths, mean_log2 = rle_obj$values)
}
