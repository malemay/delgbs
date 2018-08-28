#' A function to get the positions of the reads given a BAM file
#'
#' Outputs a data.frame containing the chromosome, initial base
#' position, and mapping quality for every read in the data set. Filters
#' based on mapping quality (default = 30, but different thresholds
#' should be tested).
#'
#' @param filename character, the name of the bam file
#' @param minq integer, the minimum mapping quality for a read to be kept
#' @param strand_diff a logical of length 1 (\code{TRUE} or \code{FALSE}).
#'   If \code{TRUE}, reads mapping to the "+" strand are tallied according to
#'   their 3'-most position, whereas reads mapping to the "-" strand are tallied
#'   according to their 5'-most position. If \code{FALSE}, all reads are tallied
#'   according to their 5'-most position. This option should be chosen according
#'   to the expected mean size of the restriction fragments, so as to maximize
#'   the probability that two reads originating from the same restriction
#'   fragment but with differing orientation are counted in the same bin.
#'
#' @return To be completed
#' @export
#'
#' @examples
#' NULL
get_pos <- function(filename, minq = 30, strand_diff = FALSE) {
  # Initializing the scanning parameters
  if(strand_diff) {
    scanpar <- Rsamtools::ScanBamParam(what = c("rname", "pos",
                                                "qwidth", "strand"),
                                       mapqFilter = minq)
  } else {
    scanpar <- Rsamtools::ScanBamParam(what = c("rname", "pos"),
                                       mapqFilter = minq)
  }

  # Reading the .bam file and converting the object to a DataFrame
  bam <- S4Vectors::DataFrame(Rsamtools::scanBam(filename, param = scanpar)[[1]])
  # Removing reads that have not been assigned a position
  bam <- bam[!is.na(bam$pos), ]
  # Keeping only reads that mapped to a chromosome (not scaffolds)
  bam <- bam[grep("Chr", bam$rname),]

  # Reads on the "+" strand are tallied according to their end position if
  #  requested (this is to allow a more consistent talling when a larger size
  #  selection window is applied
  if(strand_diff) {
    bam[["pos"]][bam$strand == "+"] <-
      bam[["pos"]][bam$strand == "+"] + bam[["qwidth"]][bam$strand == "+"]
    bam <- bam[ , c("rname", "pos")]
  }

  # Converting to a data.frame and changing chromosome name and levels
  read_pos <- as.data.frame(bam)
  names(read_pos) <- c("chr", "pos")
  read_pos$chr <- droplevels(read_pos$chr)
  #returning the data frame with read positions
  read_pos
}
