#' Copy number plot for the whole genome
#'
#' Whole genome copy number plot with partitioning lines between chromosomes.
#'
#' \code{copy_number} values can be on the relative or absolute scales or can
#' be log2 ratios but the same scale should be used consistently in each of the
#' \code{copy_number}, \code{segments} and \code{copy_number_steps} data frames.
#'
#' @param copy_number a data frame containing \code{chromosome}, \code{start},
#' \code{end} and \code{copy_number} columns where there is a row for each copy
#' number bin; may optionally contain a \code{sample} column if the data frame
#' contains data for multiple samples.
#' @param segments a data frame containing \code{chromosome}, \code{start},
#' \code{end} and \code{copy_number} columns; may optionally contain a
#' \code{sample} column if the data frame contains data for multiple samples.
#' @param sample the sample (required if the \code{copy_number} and
#' \code{segments} contain data for multiple samples)
#' @param chromosome_lengths a data frame containing \code{chromosome} and
#' \code{length} columns (optional).
#' @param copy_number_steps a data frame containing \code{absolute_copy_number}
#' and \code{copy_number} columns.
#' @param max_points_to_display maximum number of copy number points to display
#' (downsampling may be carried out if there are more copy number values than
#' this number).
#' @param min_copy_number the minimum \code{copy_number} to display.
#' @param max_copy_number the maximum \code{copy_number} to display.
#' @param copy_number_breaks breaks at which grid lines will be displayed.
#' @param point_colour the colour of the copy number points.
#' @param point_alpha the transparency of the copy number points.
#' @param point_size = the size of the copy number points.
#' @param segment_colour the colour of the copy number segments.
#' @param segment_alpha the transparency of the copy number segments.
#' @param segment_line_size the size of the lines for copy number segments.
#' @param copy_number_step_colour the colour of the copy number step lines.
#' @param copy_number_step_alpha the transparency of the copy number step lines.
#' @param copy_number_step_line_size the size of the lines for the copy number
#' steps.
#' @param xlabel,ylabel x- and y-axis labels.
#' @return a \code{ggplot} object
#' @examples
#' data(copy_number)
#'
#' segments <- copy_number_segments(copy_number)
#'
#' genome_copy_number_plot(copy_number, segments, sample = "X17222", ylabel = "relative copy number")
#'
#' absolute_copy_numbers <- 0:8
#' relative_copy_numbers <- absolute_to_relative_copy_number(absolute_copy_numbers, ploidy = 4.01, cellularity = 0.81)
#' copy_number_steps <- data.frame(absolute_copy_number = absolute_copy_numbers, copy_number = relative_copy_numbers)
#'
#' genome_copy_number_plot(
#'   copy_number,
#'   segments,
#'   sample = "X17222",
#'   copy_number_steps = copy_number_steps,
#'   min_copy_number = 0.25,
#'   max_copy_number = 2.5,
#'   xlabel = NULL,
#'   ylabel = "relative copy number")
#'
#' # filter for specific sample and convert relative copy numbers to log2 ratios
#' log2_ratio <- copy_number[copy_number$sample == "X17222", ]
#' log2_ratio$copy_number <- log2(log2_ratio$copy_number)
#' log2_ratio$segmented <- log2(log2_ratio$segmented)
#'
#' log2_ratio_segments <- copy_number_segments(log2_ratio)
#'
#' log2_ratio_steps <- copy_number_steps
#' log2_ratio_steps$copy_number <- log2(log2_ratio_steps$copy_number)
#'
#' genome_copy_number_plot(
#'   log2_ratio,
#'   log2_ratio_segments,
#'   copy_number_steps = log2_ratio_steps,
#'   min_copy_number = -2,
#'   max_copy_number = 3,
#'   xlabel = NULL,
#'   ylabel = expression(log[2]~ratio))
#'
#' # filter for specific sample and convert relative copy numbers to absolute copy numbers
#' absolute_copy_number <- copy_number[copy_number$sample == "X17222", ]
#' absolute_copy_number$copy_number <- relative_to_absolute_copy_number(absolute_copy_number$copy_number, ploidy = 4.01, cellularity = 0.81)
#' absolute_copy_number$segmented <- relative_to_absolute_copy_number(absolute_copy_number$segmented, ploidy = 4.01, cellularity = 0.81)
#'
#' absolute_segments <- copy_number_segments(absolute_copy_number)
#'
#' genome_copy_number_plot(
#'   absolute_copy_number,
#'   absolute_segments,
#'   min_copy_number = 0,
#'   max_copy_number = 10,
#'   copy_number_breaks = 0:10,
#'   ylabel = "absolute copy number") +
#'   ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "grey60"))
#'
#' @import tidyr dplyr ggplot2
#' @export
genome_copy_number_plot <- function(copy_number,
                                    segments,
                                    sample = NULL,
                                    chromosome_lengths = NULL,
                                    copy_number_steps = NULL,
                                    max_points_to_display = Inf,
                                    min_copy_number = NULL, max_copy_number = NULL, copy_number_breaks = NULL,
                                    point_colour = "black", point_alpha = 0.15, point_size = 0,
                                    segment_colour = "red", segment_alpha = 1, segment_line_size = 0.75,
                                    copy_number_step_colour = "blue", copy_number_step_alpha = 0.35, copy_number_step_line_size = 0.75,
                                    xlabel = "chromosome", ylabel = "copy number") {

  stopifnot(is.data.frame(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("start" %in% names(copy_number), is.numeric(copy_number$start))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))
  stopifnot("copy_number" %in% names(copy_number), is.numeric(copy_number$copy_number))

  stopifnot(is.data.frame(segments))
  stopifnot("chromosome" %in% names(segments))
  stopifnot("start" %in% names(segments), is.numeric(segments$start))
  stopifnot("end" %in% names(segments), is.numeric(segments$end))
  stopifnot("copy_number" %in% names(segments), is.numeric(segments$copy_number))

  if (is.null(chromosome_lengths)) {
    chromosome_lengths <- chromosome_lengths(copy_number)
  } else {
    stopifnot(is.data.frame(chromosome_lengths))
    stopifnot("chromosome" %in% names(chromosome_lengths))
    stopifnot("length" %in% names(chromosome_lengths), is.numeric(chromosome_lengths$length))
  }

  # compute offsets and genome coordinates for each chromosome
  chromosomes <- chromosome_lengths %>%
    mutate(offset = lag(cumsum(length), default = 0)) %>%
    mutate(start = offset + 1, end = offset + length) %>%
    mutate(mid = offset + round(length / 2))

  offsets <- select(chromosomes, chromosome, offset)

  # filter copy number data for the specified sample
  if (!is.null(sample)) {
    stopifnot("sample" %in% names(copy_number))
    stopifnot("sample" %in% names(segments))
    selected_sample <- sample
    copy_number <- filter(copy_number, sample == selected_sample)
    segments <- filter(segments, sample == selected_sample)
  }

  # filter out missing and non-finite values
  copy_number <- filter(copy_number, is.finite(copy_number))
  segments <- filter(segments, is.finite(copy_number))

  # compute mid-point position for the copy number bins
  copy_number <- mutate(copy_number, position = (start + end) / 2)

  # convert to genome coordinates
  copy_number <- copy_number %>%
    left_join(offsets, by = "chromosome") %>%
    mutate(position = position + offset) %>%
    select(-offset)

  segments <- segments %>%
    left_join(offsets, by = "chromosome") %>%
    mutate_at(vars(start, end), ~ . + offset) %>%
    select(-offset)

  if (is.null(min_copy_number)) {
    min_copy_number <- min(copy_number$copy_number, segments$copy_number)
  } else {
    stopifnot(is.numeric(min_copy_number), length(min_copy_number) == 1, !is.na(min_copy_number))
    copy_number <- filter(copy_number, copy_number >= min_copy_number)
    segments <- filter(segments, copy_number >= min_copy_number)
  }

  if (is.null(max_copy_number)) {
    max_copy_number <- max(copy_number$copy_number, segments$copy_number)
  } else {
    stopifnot(is.numeric(max_copy_number), length(max_copy_number) == 1, !is.na(max_copy_number))
    copy_number <- filter(copy_number, copy_number <= max_copy_number)
    segments <- filter(segments, copy_number <= max_copy_number)
  }

  stopifnot(is.numeric(max_points_to_display), length(max_points_to_display) == 1, !is.na(max_points_to_display))
  if (max_points_to_display < nrow(copy_number))
    copy_number <- sample_n(copy_number, max_points_to_display)

  segment_lines <- segments %>%
    mutate(segment_number = row_number()) %>%
    select(segment_number, start, end, copy_number) %>%
    pivot_longer(c(start, end), names_to = "type", values_to = "position") %>%
    arrange(segment_number)

  xmin <- min(chromosomes$start)
  xmax <- max(chromosomes$end)

  plot <- ggplot(data = copy_number, mapping = aes(x = position, y = copy_number)) +
    geom_vline(xintercept = chromosomes$end, colour = "grey90")

  if (!is.null(copy_number_steps)) {
    xmax <- xmin + (xmax - xmin) * 1.04

    copy_number_steps <- copy_number_steps %>%
      filter(copy_number >= min_copy_number, copy_number <= max_copy_number) %>%
      arrange(desc(absolute_copy_number))

    if (nrow(copy_number_steps) > 0) {
      plot <- plot +
        geom_hline(yintercept = copy_number_steps$copy_number, colour = copy_number_step_colour, alpha = copy_number_step_alpha, size = copy_number_step_line_size) +
        geom_label(data = copy_number_steps, mapping = aes(x = xmin + 0.98 * (xmax - xmin), y = copy_number, label = absolute_copy_number)) +
        theme(panel.grid = element_blank())
    }
  }

  if (is.null(copy_number_breaks)) copy_number_breaks = waiver()

  plot <- plot +
    geom_point(colour = point_colour, alpha = point_alpha, size = point_size) +
    geom_line(data = segment_lines, mapping = aes(x = position, y = copy_number, group = segment_number), colour = segment_colour, alpha = segment_alpha, size = segment_line_size) +
    scale_x_continuous(limits = c(xmin, xmax), expand = expansion(mult = 0), breaks = chromosomes$mid, labels = chromosomes$chromosome) +
    scale_y_continuous(limits = c(min_copy_number, max_copy_number), breaks = copy_number_breaks, expand = expansion(mult = 0)) +
    labs(x = xlabel, y = ylabel) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 11),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(size = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )

  if (!is.null(copy_number_steps) && nrow(copy_number_steps) > 0) {
    plot <- plot +
      theme(panel.grid = element_blank())
  }

  return(plot)
}

#' Copy number plot for a single chromosome
#'
#' Chromosome copy number plot with gene annotation.
#'
#' \code{copy_number} values can be on the relative or absolute scales or can
#' be log2 ratios but the same scale should be used consistently in each of the
#' \code{copy_number}, \code{segments} and \code{copy_number_steps} data frames.
#'
#' @param copy_number a data frame containing \code{chromosome}, \code{start},
#' \code{end} and \code{copy_number} columns where there is a row for each copy
#' number bin; may optionally contain a \code{sample} column if the data frame
#' contains data for multiple samples.
#' @param segments a data frame containing \code{chromosome}, \code{start},
#' \code{end} and \code{copy_number} columns; may optionally contain a
#' \code{sample} column if the data frame contains data for multiple samples.
#' @param sample the sample (required if the \code{copy_number} and
#' \code{segments} contain data for multiple samples)
#' @param chromosome the chromosome to display.
#' @param start the start coordinate within the specified chromosome.
#' @param end the end coordinate within the specified chromosome.
#' @param copy_number_steps a data frame containing \code{absolute_copy_number}
#' and \code{copy_number} columns.
#' @param genes a data frame containing \code{name}, \code{chromosome},
#' \code{start} and \code{end} columns.
#' @param max_points_to_display maximum number of copy number points to display
#' (downsampling may be carried out if there are more copy number values than
#' this number).
#' @param min_copy_number the minimum \code{copy_number} to display.
#' @param max_copy_number the maximum \code{copy_number} to display.
#' @param copy_number_breaks breaks at which grid lines will be displayed.
#' @param point_colour the colour of the copy number points.
#' @param point_alpha the transparency of the copy number points.
#' @param point_size = the size of the copy number points.
#' @param segment_colour the colour of the copy number segments.
#' @param segment_alpha the transparency of the copy number segments.
#' @param segment_line_size the size of the lines for copy number segments.
#' @param copy_number_step_colour the colour of the copy number step lines.
#' @param copy_number_step_alpha the transparency of the copy number step lines.
#' @param copy_number_step_line_size the size of the lines for the copy number
#' steps.
#' @param gene_colour the colour of the gene bars.
#' @param gene_alpha the transparency of the gene bars.
#' @param position_scale the scaling factor for the position x-axis.
#' @param xlabel,ylabel x- and y-axis labels.
#' @return a \code{ggplot} object
#' @examples
#' data(copy_number)
#' data(genes)
#'
#' segments <- copy_number_segments(copy_number)
#'
#' absolute_copy_numbers <- 0:8
#' relative_copy_numbers <- absolute_to_relative_copy_number(absolute_copy_numbers, ploidy = 4.01, cellularity = 0.81)
#' copy_number_steps <- data.frame(absolute_copy_number = absolute_copy_numbers, copy_number = relative_copy_numbers)
#'
#' chromosome_copy_number_plot(
#'   copy_number,
#'   segments,
#'   sample = "X17222",
#'   chromosome = 3,
#'   copy_number_steps = copy_number_steps,
#'   genes = genes,
#'   min_copy_number = 0.25, max_copy_number = 2.5)
#'
#' # filter for specific sample and convert relative copy numbers to log2 ratios
#' log2_ratio <- copy_number[copy_number$sample == "X17222", ]
#' log2_ratio$copy_number <- log2(log2_ratio$copy_number)
#'
#' log2_ratio_segments <- copy_number_segments(log2_ratio)
#'
#' log2_ratio_steps <- copy_number_steps
#' log2_ratio_steps$copy_number <- log2(log2_ratio_steps$copy_number)
#'
#' chromosome_copy_number_plot(
#'   log2_ratio,
#'   log2_ratio_segments,
#'   chromosome = 17, start = 7250000, end = 7750000,
#'   copy_number_steps = log2_ratio_steps,
#'   genes = genes,
#'   min_copy_number = -2, max_copy_number = 2,
#'   position_scale = 1,
#'   xlabel = "position",
#'   ylabel = expression(log[2]~ratio))
#'
#' @import tidyr dplyr ggplot2 scales
#' @export
chromosome_copy_number_plot <- function(copy_number,
                                        segments,
                                        sample = NULL,
                                        chromosome, start = NULL, end = NULL,
                                        copy_number_steps = NULL,
                                        genes = NULL,
                                        max_points_to_display = Inf,
                                        min_copy_number = NULL, max_copy_number = NULL, copy_number_breaks = NULL,
                                        point_colour = "black", point_alpha = NULL, point_size = NULL,
                                        segment_colour = "red", segment_alpha = 1, segment_line_size = 0.75,
                                        copy_number_step_colour = "blue", copy_number_step_alpha = 0.35, copy_number_step_line_size = 0.75,
                                        gene_colour = "darkgreen", gene_alpha = 0.25,
                                        position_scale = 1e-6,
                                        xlabel = "position (Mbp)", ylabel = "copy number") {

  stopifnot(is.data.frame(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("start" %in% names(copy_number), is.numeric(copy_number$start))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))
  stopifnot("copy_number" %in% names(copy_number), is.numeric(copy_number$copy_number))

  stopifnot(is.data.frame(segments))
  stopifnot("chromosome" %in% names(segments))
  stopifnot("start" %in% names(segments), is.numeric(segments$start))
  stopifnot("end" %in% names(segments), is.numeric(segments$end))
  stopifnot("copy_number" %in% names(segments), is.numeric(segments$copy_number))

  # filter copy number data for the specified sample
  if (!is.null(sample)) {
    stopifnot("sample" %in% names(copy_number))
    stopifnot("sample" %in% names(segments))
    selected_sample <- sample
    copy_number <- filter(copy_number, sample == selected_sample)
    segments <- filter(segments, sample == selected_sample)
  }

  # filter copy number data for specified chromosome
  selected_chromosome <- chromosome
  copy_number <- filter(copy_number, chromosome == selected_chromosome)
  segments <- filter(segments, chromosome == selected_chromosome)

  # filter out missing and non-finite values
  copy_number <- filter(copy_number, is.finite(copy_number))
  segments <- filter(segments, is.finite(copy_number))

  # compute mid-point position for the copy number bins
  copy_number <- mutate(copy_number, position = (start + end) / 2)

  # filter copy number data for specified start and end
  if (!is.null(start)) {
    stopifnot(is.numeric(start), length(start) == 1, !is.na(start))
    selected_start <- start
    copy_number <- filter(copy_number, position >= selected_start)
    segments <- segments %>%
      filter(end > selected_start) %>%
      mutate(start = pmax(start, selected_start))
  }

  if (!is.null(end)) {
    stopifnot(is.numeric(end), length(end) == 1, !is.na(end))
    selected_end <- end
    copy_number <- filter(copy_number, position <= selected_end)
    segments <- segments %>%
      filter(start < selected_end) %>%
      mutate(end = pmin(end, selected_end))
  }

  if (is.null(min_copy_number)) {
    min_copy_number <- min(copy_number$copy_number, segments$copy_number)
  } else {
    stopifnot(is.numeric(min_copy_number), length(min_copy_number) == 1, !is.na(min_copy_number))
    copy_number <- filter(copy_number, copy_number >= min_copy_number)
    segments <- filter(segments, copy_number >= min_copy_number)
  }

  if (is.null(max_copy_number)) {
    max_copy_number <- max(copy_number$copy_number, segments$copy_number)
  } else {
    stopifnot(is.numeric(max_copy_number), length(max_copy_number) == 1, !is.na(max_copy_number))
    copy_number <- filter(copy_number, copy_number <= max_copy_number)
    segments <- filter(segments, copy_number <= max_copy_number)
  }

  if (max_points_to_display < nrow(copy_number))
    copy_number <- sample_n(copy_number, max_points_to_display)

  number_of_dots <- nrow(copy_number)
  if (is.null(point_size)) point_size <- 0 + max(2500 - number_of_dots, 0) * 1.5 / 2500
  if (is.null(point_alpha)) point_alpha <- 0.15 + max(1000 - number_of_dots, 0) * 0.35 / 1000

  segment_lines <- segments %>%
    mutate(segment_number = row_number()) %>%
    select(segment_number, start, end, copy_number) %>%
    pivot_longer(c(start, end), names_to = "type", values_to = "position") %>%
    arrange(segment_number)

  xmin <- start
  xmax <- end
  if (is.null(xmin)) xmin <- min(copy_number$position, segment_lines$position)
  if (is.null(xmax)) xmax <- max(copy_number$position, segment_lines$position)
  limits <- c(xmin, xmax)

  plot <- ggplot()

  if (!is.null(genes)) {
    genes <- filter(genes, chromosome == selected_chromosome, start <= xmax, end >= xmin)
    if (nrow(genes) > 0) {
      gene_boundaries <- genes %>%
        pivot_longer(c(start, end), names_to = "type", values_to = "position") %>%
        filter(position >= xmin, position <= xmax) %>%
        distinct(position)
      genes <- genes %>%
        mutate(start = pmax(start, xmin)) %>%
        mutate(end = pmin(end, xmax))
      plot <- plot +
        geom_rect(
          data = genes,
          mapping = aes(xmin = start, xmax = end, ymin = min_copy_number, ymax = max_copy_number),
          fill = gene_colour,
          alpha = gene_alpha
        ) +
        geom_vline(
          data = gene_boundaries,
          mapping = aes(xintercept = position),
          colour = gene_colour,
          alpha = gene_alpha
        ) +
        geom_text(
          data = genes,
          mapping = aes(x = (start + end) / 2, y = min_copy_number + 0.975 * (max_copy_number - min_copy_number), label = name),
          colour = gene_colour,
          size = 4
        )
    }
  }

  if (!is.null(copy_number_steps)) {
    limits <- c(xmin, xmin + (xmax - xmin) * 1.04)

    copy_number_steps <- copy_number_steps %>%
      filter(copy_number >= min_copy_number, copy_number <= max_copy_number) %>%
      arrange(desc(absolute_copy_number))

    if (nrow(copy_number_steps) > 0) {
      plot <- plot +
        geom_hline(yintercept = copy_number_steps$copy_number, colour = copy_number_step_colour, alpha = copy_number_step_alpha, size = copy_number_step_line_size) +
        geom_label(data = copy_number_steps, mapping = aes(x = limits[1] + 0.98 * (limits[2] - limits[1]), y = copy_number, label = absolute_copy_number))
    }
  }

  if (is.null(copy_number_breaks)) copy_number_breaks = waiver()

  plot <- plot +
    geom_point(data = copy_number, mapping = aes(x = position, y = copy_number), size = point_size, colour = point_colour, alpha = point_alpha) +
    geom_line(data = segment_lines, mapping = aes(x = position, y = copy_number, group = segment_number), colour = segment_colour, alpha = segment_alpha, size = segment_line_size) +
    scale_x_continuous(limits = limits, expand = expansion(mult = 0), labels = scales::unit_format(scale = position_scale, big.mark = ",", unit = "", sep = "")) +
    scale_y_continuous(limits = c(min_copy_number, max_copy_number), breaks = copy_number_breaks, expand = expansion(mult = 0)) +
    labs(x = xlabel, y = ylabel) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 11),
      axis.ticks.y = element_line(size = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank()
    )

  if (!is.null(copy_number_steps) && nrow(copy_number_steps) > 0) {
    plot <- plot +
      theme(panel.grid = element_blank())
  }

  return(plot)
}

#' Copy number density plot
#'
#' Creates a ggplot for the density of the given copy numbers values.
#'
#' @param copy_numbers a numeric vector containing copy number values.
#' @param copy_number_steps a data frame containing \code{copy_number}
#' and \code{absolute_copy_number} columns.
#' @param min_copy_number the minimum relative copy number to display.
#' @param max_copy_number the maximum relative copy number to display.
#' @param copy_number_step_colour the colour of the copy number step lines.
#' @param copy_number_step_alpha the transparency of the copy number step lines.
#' @param copy_number_step_line_size the size of the lines for the copy number
#' steps.
#' @param xlabel,ylabel x- and y-axis labels.
#' @return a \code{ggplot} object
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#'
#' copy_number_density_plot(copy_number$segmented)
#'
#' absolute_copy_numbers <- 0:8
#' relative_copy_numbers <- absolute_to_relative_copy_number(absolute_copy_numbers, ploidy = 4.01, cellularity = 0.81)
#' copy_number_steps <- data.frame(absolute_copy_number = absolute_copy_numbers, copy_number = relative_copy_numbers)
#'
#' copy_number_density_plot(copy_number$segmented, copy_number_steps, 0, 2.5)
#' @import dplyr ggplot2
#' @export
copy_number_density_plot <- function(copy_numbers,
                                     copy_number_steps = NULL,
                                     min_copy_number = NULL, max_copy_number = NULL,
                                     copy_number_step_colour = "blue", copy_number_step_alpha = 0.35, copy_number_step_line_size = 0.75,
                                     xlabel = "copy number", ylabel = "density") {

  copy_number_density <- copy_number_density(copy_numbers, min_copy_number, max_copy_number)

  if (is.null(min_copy_number)) min_copy_number <- min(copy_number_density$copy_number)
  if (is.null(max_copy_number)) max_copy_number <- max(copy_number_density$copy_number)

  plot <- ggplot(data = copy_number_density, mapping = aes(x = copy_number, y = density)) +
    geom_density(stat = "identity") +
    scale_x_continuous(limits = c(min_copy_number, max_copy_number), expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = xlabel, y = ylabel) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      panel.grid = element_blank()
    )

  if (!is.null(copy_number_steps)) {
    copy_number_steps <- copy_number_steps %>%
      filter(copy_number >= min_copy_number, copy_number <= max_copy_number) %>%
      arrange(desc(absolute_copy_number))

    if (nrow(copy_number_steps) > 0) {
      plot <- plot +
        geom_vline(xintercept = copy_number_steps$copy_number, colour = copy_number_step_colour, alpha = copy_number_step_alpha, size = copy_number_step_line_size) +
        geom_label(data = copy_number_steps, mapping = aes(x = copy_number, y = 1.1 * max(copy_number_density$density), label = absolute_copy_number))
    }
  }

  return(plot)
}

#' Heat map for given distances for grid of ploidies and cellularities
#'
#' Heat map of the given distances for a grid of ploidy and cellularity
#' values.
#'
#' @param distances a data frame containing \code{ploidy}, \code{cellularity}
#' and \code{distance} columns.
#' @param low_distance_colour the colour for low distances.
#' @param high_distance_colour the colour for high distances.
#' @param xlabel,ylabel x- and y-axis labels.
#' @return a \code{ggplot} object
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' segments <- copy_number_segments(copy_number)
#'
#' distances <- absolute_copy_number_distance_grid(segments$copy_number,
#'                                                 segments$weight,
#'                                                 distance_function = "MAD")
#' distance_heatmap(distances)
#' @import dplyr ggplot2
#' @export
distance_heatmap <- function(distances,
                             low_distance_colour = "red",
                             high_distance_colour = "blue",
                             xlabel = "cellularity", ylabel = "ploidy") {

  distances <- mutate(distances, goodness_of_fit = 1 / distance)

  ggplot(data = distances, mapping = aes(x = cellularity, y = ploidy)) +
    geom_raster(mapping = aes(fill = goodness_of_fit), interpolate = TRUE) +
    scale_fill_gradient(low = high_distance_colour, high = low_distance_colour) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    labs(x = xlabel, y = ylabel) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
}

#' Distance heat map for fitting a set of relative copy numbers
#'
#' Heat map of the distance of the given relative copy numbers to whole numbers
#' after scaling to absolute copy numbers for various ploidies and
#' cellularities.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @param low_distance_colour the colour for low distances.
#' @param high_distance_colour the colour for high distances.
#' @param min_ploidy the minimum ploidy to display.
#' @param max_ploidy the minimum ploidy to display.
#' @param xlabel,ylabel x- and y-axis labels.
#' @return a \code{ggplot} object
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' segments <- copy_number_segments(copy_number)
#'
#' absolute_copy_number_distance_heatmap(
#'   segments$copy_number,
#'   segments$weight,
#'   distance_function = "MAD"
#' )
#' @import dplyr ggplot2
#' @export
absolute_copy_number_distance_heatmap <- function(relative_copy_numbers, weights = NULL,
                                                  distance_function = "MAD",
                                                  min_ploidy = 1.5, max_ploidy = 5.5,
                                                  low_distance_colour = "red",
                                                  high_distance_colour = "blue",
                                                  xlabel = "cellularity", ylabel = "ploidy") {
  distances <- absolute_copy_number_distance_grid(
    relative_copy_numbers, weights,
    min_ploidy = min_ploidy, max_ploidy = max_ploidy,
    distance_function = distance_function)

  distance_heatmap(distances,
                   low_distance_colour, high_distance_colour,
                   xlabel, ylabel)
}
