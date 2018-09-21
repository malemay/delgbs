#' A function to visualize log2 ratios and events for all chromosomes of a given
#' individual
#'
#' @param log2_data a data.frame, the log2 ratio data to be plotted. The first
#'   column must be named "chr" and hold chromosome names. The second column
#'   must be named "lower" and hold the starting position of the probe. The
#'   second column must be named "upper" and hold the end position of the probe.
#'   The fourth column must be named "ind" and hold the name of the individual
#'   about which a particular data record is. The alst column must be named
#'   "log2_ratio" and gives the log2 ratio observed at a particular location
#'   for a given individual.
#' @param ind a character vector of length one. The individual for which the
#'   data is the be plotted.
#' @param event_data an optional data.frame holding the positions of CNV events
#'   (homozygous deletions, heterozygous deletions or duplications) to be
#'   plotted alongside the log2 ratio data. The data.frame can contain an
#'   arbitrary number of columns, but at least the 6 following columns must be
#'   present: ind, chr, start, end, mean_log2, type. "ind" corresponds to the
#'   individual in which the event is found. "chr" corresponds to the
#'   chromosome on which the event is located. "start" gives the starting
#'   position of the vent. "end" gives the end position of the event.
#'   "mean_log2" corresponds to the mean log2 ratio over the extent of the
#'   event. "type" can take one of three mutually exclusive values: "homdel"
#'   stands for homozygous deletions, "hetdel" stands for heterozygous deletions,
#'   and "dup" stands for duplication.
#' @param het_sites a data.frame giving the positions of heterogeneity regions
#'   among the chromosomes to be plotted. These regions will be indicated by
#'   a red-colored semi-transparent rectangle for easier interpretation of the
#'   results. The data.frame must contain three columns, "chr", "start", "stop",
#'   collectively indicating the positions of the heterogeneity regions.
#' @param event_colors a character vector of length three indicating the colors
#'   to use for plotting, respectively, homozygous deletions, heterozygous
#'   deletions, and duplications.
#'
#' @return a ggplot graph
#' @export
#'
#' @examples
#' NULL
segplot <- function(log2_data, ind, event_data = NULL, het_sites = NULL,
                    event_colors = c("red", "orange", "forestgreen")) {

  # Testing the inputs
  stopifnot(all(c("chr", "lower", "upper", "ind", "log2_ratio") %in% names(log2_data)))
  stopifnot(length(ind) == 1 && ind %in% log2_data$ind)

  # Subsetting the input data
  log2_data <- log2_data[log2_data$ind == ind, ]

  # Preparing data columns needed for plotting
  log2_data$POSITION <- (log2_data$lower + log2_data$upper) / (2*10^6)

  # Creating the main plot, first without CNV and heterogeneity information
  splot <-
    ggplot2::ggplot(log2_data) +
    ggplot2::geom_point(ggplot2::aes_string(x = "POSITION", y = "log2_ratio"),
               size = 0.08, col = "blue") +
    ggplot2::facet_wrap(~ chr) +
    ggplot2::xlab("Position (Mb)") + ggplot2::ylab("log2 ratio") +
    ggplot2::theme_bw()

  # Event (CNV) data are plotted if provided
  if(!is.null(event_data) && any(event_data$ind == ind)) {
    # Input testing
    stopifnot(all(c("ind", "chr", "start", "end", "mean_log2", "type") %in% names(event_data)))
    stopifnot(is.character(event_colors) && length(event_colors) == 3)

    # Extracting the data corresponding to the individual to plot
    events <- event_data[event_data$ind == ind, ]

    # These are the y-positions of the arrows that will highlight the events
    #  the arrows should occupy one fifth of the plotting height
    arrow_y1 <- max(log2_data$log2_ratio, na.rm = TRUE)
    arrow_y2 <- (arrow_y1 - min(log2_data$log2_ratio, na.rm = TRUE)) / 5

    # Preparing data columns needed for plotting
    events$START    <- events$start / 10^6
    events$END      <- events$end / 10^6

    # Adding the geom data to the plot
    splot <- splot +
      # Adding thick horizontal lines to indicate the extent and location
      ggplot2::geom_segment(data = events,
                   ggplot2::aes_string(x = "START", xend = "END",
                       y = "mean_log2", yend = "mean_log2", col = "type"),
                   size = 0.5) +
      # Adding arrows pointing down to highlight the location of the events
      ggplot2::geom_segment(data = events,
                   ggplot2::aes_string(x = "START", xend = "START", col = "type"),
                   y = arrow_y1, yend = arrow_y2,
                   arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "npc")),
                   size = 0.8) +
      # Plotting different event types with different colors
      #  Colors can be edited through the argument event_colors
      ggplot2::scale_color_manual(name = "Type",
                         values = c("homdel" = event_colors[1],
                                    "hetdel" = event_colors[2],
                                    "dup" = event_colors[3]),
                         labels = c("homdel" = "homdel",
                                    "hetdel" = "hetdel",
                                    "dup" = "dup")) +
      # The legend is not shown in order to spare plotting space
      ggplot2::guides(color = "none")
  }

  # Heterogeneity regions are plotted if provided
  if(!is.null(het_sites)) {
    # Input testing
    stopifnot(all(c("chr", "start", "stop") %in% names(het_sites)))
    # Finding the minimum and maximum y positions of the rectangles
    ymin <- min(log2_data$log2_ratio, na.rm = TRUE)
    ymax <- max(log2_data$log2_ratio, na.rm = TRUE)

    # Preparing data columns needed for plotting
    het_sites$START    <- het_sites$start / 10^6
    het_sites$STOP      <- het_sites$stop / 10^6

    # Adding the data to the plot
    splot <- splot +
      ggplot2::geom_rect(data = het_sites, ggplot2::aes_string(xmin = "START", xmax = "STOP"),
                ymin = ymin, ymax = ymax, col = "transparent", fill = "red", alpha = 0.1)
  }

  # Outputting the final graph
  splot
}
