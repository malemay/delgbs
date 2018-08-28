#' A function to visualize log2 ratios and events for all individuals
#'  on a given chromosome
#'
#' @param log2_data a data.frame, the log2 ratio data to be plotted. The first
#'   column must be named "chr" and hold chromosome names. The second column
#'   must be named "lower" and hold the starting position of the probe. The
#'   second column must be named "upper" and hold the end position of the probe.
#'   The fourth column must be named "ind" and hold the name of the individual
#'   about which a particular data record is. The alst column must be named
#'   "log2_ratio" and gives the log2 ratio observed at a particular location
#'   for a given individual.
#' @param chrom a character vector of length one. The chromosome for which the
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
#'   among the individuals to be plotted. These regions will be indicated by
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
chromplot <- function(log2_data, chrom, event_data = NULL, het_sites = NULL,
                      event_colors = c("red", "orange", "forestgreen")) {

  # Testing the inputs
  stopifnot(all(c("chr", "lower", "upper", "ind", "log2_ratio") %in% names(log2_data)))
  stopifnot(length(chrom) == 1 && chrom %in% log2_data$chr)

  # Generating subsets of the data according the the chromosome to be plotted
  log2_data <- log2_data[log2_data$chr == chrom, ]

  # This will ensure that the labels on the plot match their individual
  log2_data$ind <- as.factor(log2_data$ind)
  ind_labels <- levels(log2_data$ind)

  # Creating the main plot, first without CNV and heterogeneity information
  cplot <-
    # Initializing the plot
    ggplot2::ggplot(data = log2_data) +
    # Plotting the log2 ratio data
    ggplot2::geom_point(ggplot2::aes_string(x = ("lower" + "upper") / (2*10^6),
                                            y = "log2_ratio"),
                        col = "blue", size = 0.1) +
    # Plotting each individual on its own row
    ggplot2::facet_wrap(~ ind, ncol = 1) +
    # Setting the y-axis so as to maximize the plotting area
    ggplot2::scale_y_continuous(limits = c(min(log2_data$log2_ratio),
                                  max(log2_data$log2_ratio)),
                                expand = c(0, 0)) +
    # Same here with the x-axis, leaving space for the labels
    ggplot2::scale_x_continuous(limits = c(0, max(log2_data$upper) / 10^6 + 3),
                                expand = c(0,0),
                                breaks = seq(0, max(log2_data$upper) / 10^6, 5)) +
    # Adding the ID of the individual to the right of each log2 profile
    ggplot2::annotate("text", label = ind_labels, hjust = 0,
                      x = max(log2_data$upper) / 10^6 + 0.1, y = 0) +
    # Theme elements mostly aim at maximizing the area for plotting
    ggplot2::theme_bw() +
    ggplot2::theme(panel.spacing = grid::unit(0, "lines"),
                   strip.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())

  # Event (CNV) data are plotted if provided
  if(!is.null(event_data) && any(event_data$chr == chrom)) {
    # Input testing
    stopifnot(all(c("ind", "chr", "start", "end", "mean_log2", "type") %in% names(event_data)))
    stopifnot(is.character(event_colors) && length(event_colors) == 3)

    # Extracting the data corresponding to the chromosome to plot
    events <- event_data[event_data$chr == chrom, ]

    # These are the y-positions of the arrows that will highlight the events
    #  the arrows shoudl occupy one fifth of the plotting height
    arrow_y1 <- max(log2_data$log2_ratio, na.rm = TRUE)
    arrow_y2 <- (arrow_y1 - min(log2_data$log2_ratio, na.rm = TRUE)) / 5

    # Adding the CNV information to the plot
    cplot <- cplot +
      # Adding thick horizontal lines to indicate the extent and location
      ggplot2::geom_segment(data = events,
                            ggplot2::aes_string(x = "start" / 10^6,
                                                xend = "end" / 10^6,
                                                y = "mean_log2",
                                                yend = "mean_log2",
                                                col = "type"),
                            size = 0.5) +
      # Adding arrows pointing down to highlight the location of the events
      ggplot2::geom_segment(data = events,
                            ggplot2::aes_string(x = "start" / 10^6,
                                                xend = "start" / 10^6,
                                                col = "type"),
                            y = arrow_y1, yend = arrow_y2,
                            arrow = grid::arrow(length = grid::unit(0.05, "npc")),
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
  if(!is.null(het_sites) && any(het_sites$chr == chrom)) {
    # Input testing
    stopifnot(all(c("chr", "start", "stop") %in% names(het_sites)))
    # Subsetting the data
    het_sites <- het_sites[het_sites$chr == chrom, ]
    # Finding the minimum and maximum y positions of the rectangles
    ymin <- min(log2_data$log2_ratio, na.rm = TRUE)
    ymax <- max(log2_data$log2_ratio, na.rm = TRUE)
    # Adding the data to the plot
    cplot <- cplot +
      ggplot2::geom_rect(data = het_sites,
                         ggplot2::aes_string(xmin = "start" / 10^6, xmax = "stop" / 10^6),
                         ymin = ymin, ymax = ymax, col = "transparent",
                         fill = "red", alpha = 0.1)

  }

  # Outputting the final plot
  cplot
}
