# Compute the scaling factor for converting relative copy numbers to absolute
# copy numbers for the given ploidy and cellularity.
absolute_copy_number_scaling_factor <- function(ploidy, cellularity) {
  ploidy + (2 / cellularity) - 2
}

#' Convert relative copy numbers to absolute copy numbers
#'
#' Convert relative copy numbers to absolute copy numbers based on the given
#' ploidy and cellularity.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are
#' from the tumour.
#' @return a numeric vector of absolute copy numbers.
#' @examples
#' relative_to_absolute_copy_number(c(0.98, 1.6, 1.23), 4.01, 0.77)
#' @export
relative_to_absolute_copy_number <- function(relative_copy_numbers, ploidy, cellularity) {

  stopifnot(is.numeric(relative_copy_numbers))
  stopifnot(is.numeric(ploidy))
  stopifnot(is.numeric(cellularity))

  ploidy + absolute_copy_number_scaling_factor(ploidy, cellularity) * (relative_copy_numbers - 1)
}

#' Convert absolute copy numbers to relative copy numbers
#'
#' Convert absolute copy numbers to relative copy numbers based on the given
#' ploidy and cellularity.
#'
#' @param absolute_copy_numbers a numeric vector of absolute copy numbers.
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @return a numeric vector containing relative copy numbers, i.e. ratios of
#' copy numbers to the average copy number.
#' @examples
#' absolute_to_relative_copy_number(4, 4.01, 0.77)
#' absolute_to_relative_copy_number(1:10, 4.01, 0.77)
#' @export
absolute_to_relative_copy_number <- function(absolute_copy_numbers, ploidy, cellularity) {

  stopifnot(is.numeric(absolute_copy_numbers))
  stopifnot(is.numeric(ploidy))
  stopifnot(is.numeric(cellularity))

  ((absolute_copy_numbers - ploidy) / absolute_copy_number_scaling_factor(ploidy, cellularity)) + 1
}

#' Compute tumour DNA fraction for the given absolute copy number and
#' cellularity
#'
#' Compute the tumour DNA fraction for the given absolute copy number(s) in the
#' tumour and the cellularity, i.e. the fraction of all cells that are tumour
#' cells.
#'
#' @param absolute_copy_number the absolute copy number(s) (a numeric value).
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @return the fraction of DNA that originates from tumour cells.
#' @examples
#' tumour_fraction(3, 0.82)
#' @export
tumour_fraction <- function(absolute_copy_number, cellularity) {

  stopifnot(is.numeric(absolute_copy_number))
  stopifnot(is.numeric(cellularity))

  tumour <- absolute_copy_number * cellularity
  normal <- 2 * (1 - cellularity)
  tumour / (tumour + normal)
}

#' Copy number density estimation
#'
#' Obtain density estimates for the distribution of the given copy number
#' values.
#'
#' @param copy_numbers a numeric vector of relative copy number values.
#' @param min_copy_number the lower copy number value in the range of values to
#' estimate the density.
#' @param max_copy_number the upper copy number value in the range of values to
#' estimate the density.
#' @param n the number of equally-spaced points for which the density will be
#' estimated (a smaller number may be returned if \code{min_copy_number} and/or
#' \code{max_copy_number} are specified).
#' @return a data frame with copy number and density columns
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#'
#' density <- copy_number_density(copy_number$segmented)
#'
#' density <- copy_number_density(copy_number$segmented,
#'                                min_copy_number = 0,
#'                                max_copy_number = 2.5)
#' @import tibble dplyr
#' @export
copy_number_density <- function(copy_numbers, min_copy_number = NULL, max_copy_number = NULL, n = 512) {

  stopifnot(is.numeric(copy_numbers), length(copy_numbers) > 0)

  copy_number <- tibble(copy_number = copy_numbers) %>%
    filter(!is.na(copy_number))

  if (!is.null(min_copy_number)) {
    stopifnot(is.numeric(min_copy_number), length(min_copy_number) == 1)
    copy_number <- filter(copy_number, copy_number >= min_copy_number)
  }

  if (!is.null(max_copy_number)) {
    stopifnot(is.numeric(max_copy_number), length(max_copy_number) == 1)
    copy_number <- filter(copy_number, copy_number <= max_copy_number)
  }

  stopifnot(is.numeric(n))

  if (nrow(copy_number) < 2) return(tibble(copy_number = double(), density = double()))

  density <- density(copy_number$copy_number, n = n)
  density <- tibble(copy_number = density$x, density = density$y)

  if (!is.null(min_copy_number)) density <- filter(density, copy_number >= min_copy_number)
  if (!is.null(max_copy_number)) density <- filter(density, copy_number <= max_copy_number)

  density
}

#' Obtain maxima in a distribution of relative copy numbers.
#'
#' This function obtains peaks in the density for the given set of relative
#' copy numbers.
#'
#' @param copy_numbers a numeric vector of absolute copy number values.
#' @param min_copy_number the lower copy number value in the range of
#' values to estimate the density.
#' @param max_copy_number the upper copy number value in the range of
#' values to estimate the density.
#' @param lower_threshold the lower threshold for a maximum as a fraction of
#' the density of the largest maximum.
#' @return a data frame with copy number and density columns and a row for each
#' maximum in the distribution.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' maxima <- copy_number_maxima(copy_number$segmented, lower_threshold = 0.1)
#' @import dplyr
#' @export
copy_number_maxima <- function(copy_numbers, min_copy_number = 0, max_copy_number = 3, lower_threshold = 0) {

  stopifnot(is.numeric(copy_numbers), length(copy_numbers) > 0)
  stopifnot(is.numeric(min_copy_number), length(min_copy_number) == 1)
  stopifnot(is.numeric(max_copy_number), length(max_copy_number) == 1)
  stopifnot(is.numeric(lower_threshold), length(lower_threshold) == 1)

  copy_number_density(copy_numbers, min_copy_number, max_copy_number) %>%
    mutate(up = density > lag(density)) %>%
    filter(up & !lead(up)) %>%
    filter(density >= lower_threshold * max(density)) %>%
    select(-up)
}

#' Collapse copy number bins into segments
#'
#' Collapses copy number segments for the given copy number data frame.
#'
#' Adjacent bins with the same segmented copy number (\code{segmented} column)
#' for the same sample and on the same chromosome are treated as belonging to
#' the same segment.
#'
#' The returned data frame contains a row for each segment with columns for
#' the segment number (\code{segment}), the number of bins or rows in the
#' provided copy number data frame that were collapsed to form the segment
#' (\code{bin_count}), the sum of lengths of those bins
#' (\code{sum_of_bin_lengths}) and the weight for the segment \code{weight}).
#' The segment weight is the median-scaled sum of bin lengths.
#'
#' @param copy_number a data frame containing \code{sample}, \code{chromosome},
#'     \code{start}, \code{end} and \code{segmented} columns.
#' @return a data frame containing copy number segments with
#' \code{segment}, \code{sample}, \code{chromosome}, \code{start}, \code{end},
#' \code{copy_number}, \code{bin_count}, \code{sum_of_bin_lengths} and
#' \code{weight}.
#' @examples
#' data(copy_number)
#'
#' # segments for all samples
#' segments <- copy_number_segments(copy_number)
#'
#' # segments for single sample
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' segments <- copy_number_segments(copy_number)
#' @import dplyr
#' @export
copy_number_segments <- function(copy_number) {

  stopifnot(is.data.frame(copy_number))
  stopifnot("sample" %in% names(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("start" %in% names(copy_number), is.numeric(copy_number$start))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))
  stopifnot("segmented" %in% names(copy_number), is.numeric(copy_number$segmented))

  copy_number %>%
    filter(!is.na(segmented)) %>%
    mutate(length = end - start + 1) %>%
    arrange(sample, chromosome, start) %>%
    mutate(new_segment = row_number() == 1 | !(sample == lag(sample) & chromosome == lag(chromosome) & segmented == lag(segmented))) %>%
    mutate(segment = cumsum(new_segment)) %>%
    group_by(segment) %>%
    summarize(
      sample = first(sample),
      chromosome = first(chromosome),
      start = first(start),
      end = last(end),
      copy_number = first(segmented),
      bin_count = n(),
      sum_of_bin_lengths = sum(length),
      weight = sum(length) / median(length)
    )
}

#' Obtain chromosome lengths and offsets
#'
#' Obtains chromosome lengths and offset for the given copy number data.
#'
#' Offsets are used to convert chromosome coordinates into genomic coordinates.
#'
#' @param copy_number a data frame containing \code{chromosome} and \code{end}
#' columns.
#' @return a data frame containing \code{chromosome}, \code{length} and
#' \code{offset} columns.
#' @examples
#' data(copy_number)
#' chromosome_offsets(copy_number)
#' @import dplyr
#' @export
chromosome_offsets <- function(copy_number) {

  stopifnot(is.data.frame(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))

  copy_number %>%
    group_by(chromosome) %>%
    summarize(length = as.numeric(max(end))) %>%
    mutate(offset = lag(cumsum(length), default = 0))
}

#' Convert chromosome coordinates to genomic coordinates
#'
#' Converts the specified columns from chromosome coordinates to genomic
#' coordinates.
#'
#' Genomic coordinates are obtained by adding the chromosome coordinate to the
#' sum of the lengths of all preceding chromosomes.
#'
#' @param copy_number a data frame containing \code{chromosome} and one or more
#' coordinate columns, e.g. \code{position}, \code{start} or \code{end}.
#' @param column_names a vector of column names to convert to genomic
#' coordinates.
#' @param offsets a data frame containing \code{chromosome} and \code{offset}
#' columns (optional, if not given then offsets will be calculated for the
#' given copy number).
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' copy_number <- convert_to_genomic_coordinates(copy_number, c("start", "end"))
#' @import dplyr
#' @export
convert_to_genomic_coordinates <- function(copy_number,
                                           column_names,
                                           offsets = NULL) {
  stopifnot(is.data.frame(copy_number))
  stopifnot("chromosome" %in% names(copy_number))

  stopifnot(is.character(column_names), length(column_names) > 0)
  stopifnot(all(column_names %in% names(copy_number)))

  if (is.null(offsets)) {
    offsets <- chromosome_offsets(copy_number)
  } else {
    stopifnot(is.data.frame(offsets))
    stopifnot("chromosome" %in% names(offsets))
    stopifnot("offset" %in% names(offsets), is.numeric(offsets$offset))
  }

  copy_number %>%
    left_join(select(offsets, chromosome, offset), by = "chromosome") %>%
    mutate_at(vars(column_names), ~ . + offset) %>%
    select(-offset)
}
