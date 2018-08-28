#' Length of the chromosomes of genome assembly 2 of soybean
#'
#' A dataset listing the lengths (in base pairs) of the 20 chromosomes of the
#' soybean (\emph{Glycine max}) genome assembly 2 (Glycine_max_v2.0)
#'
#' @format A \code{data.frame} with 20 rows and 2 columns:
#' \describe{
#'   \item{chr}{character; the name of the chromosome in the assembly}
#'   \item{length}{integer; the number of base pairs of the chromosome}
#' }
#'
#' @source \url{http://plants.ensembl.org/Glycine_max/Info/Annotation/#assembly}
"gmax_chr_sizes"

#' List of 1-kb bins for tallying reads in the soybean genome
#'
#' This list can be used as input to the function \code{\link{count_reads}} to
#' tally reads from a .bam file into 1-kb bins covering the whole assembly 2 of
#' the soybean (\emph{Glycine max}) genome.
#'
#' @format A list of 20 numeric vectors. Elements of the list are named after
#'   the 20 chromosomes of the genome assembly, while the numeric vectors give
#'   the breakpoints of the 1-kb bins of each chromosome.
#'
"gmax_binlist"
