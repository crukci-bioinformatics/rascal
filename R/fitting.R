#' Compute distance function for fitting relative copy numbers to absolute
#'
#' Computes a distance function on fitting relative copy numbers to whole number
#' absolute values with the given ploidy and cellularity.
#'
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' absolute_copy_number_distance(3, 0.67, copy_number$segmented)
#' @export
absolute_copy_number_distance <- function(ploidy, cellularity,
                                          relative_copy_numbers, weights = NULL,
                                          distance_function = c("MAD", "RMSD")) {

  distance <- 1e10
  if (cellularity < 0.0 || cellularity > 1.0) return(distance)

  absolute_copy_numbers <- relative_to_absolute_copy_number(relative_copy_numbers, ploidy, cellularity)

  absolute_copy_number_steps <- round(absolute_copy_numbers)

  differences <- abs(absolute_copy_numbers - absolute_copy_number_steps)
  # differences[which(absolute_copy_numbers > 0 & differences > 0.5)] <- 0.5

  if (is.null(weights)) weights <- rep(1, length(relative_copy_numbers))

  if (is.character(distance_function)) {
    if (distance_function[1] == "MAD") {
      distance <- weighted.mean(differences, weights, na.rm = TRUE)
    } else if (distance_function[1] == "RMSD") {
      distance <- sqrt(weighted.mean(differences ^ 2, weights, na.rm = TRUE))
    }
  }

  distance
}

#' Compute distance function for fitting relative copy numbers for a grid of
#' ploidies and cellularities
#'
#' Computes a distance function on fitting relative copy numbers to whole number
#' absolute values for a grid of ploidies and cellularities.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param min_ploidy,max_ploidy the range of ploidies.
#' @param ploidy_step the stepwise increment of ploidies along the grid.
#' @param min_cellularity,max_cellularity the range of cellularities.
#' @param cellularity_step the stepwise increment of cellularities along the
#' grid.
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' segments <- copy_number_segments(copy_number)
#' distances <- absolute_copy_number_distance_grid(segments$copy_number, segments$weight)
#' @import dplyr
#' @importFrom purrr pmap_dbl
#' @export
absolute_copy_number_distance_grid <- function(relative_copy_numbers, weights = NULL,
                                               min_ploidy = 1.5, max_ploidy = 5.5, ploidy_step = 0.01,
                                               min_cellularity = 0.01, max_cellularity = 1.0, cellularity_step = 0.01,
                                               distance_function = c("MAD", "RMSD")) {

  grid <- expand_grid(
    ploidy = seq(min_ploidy, max_ploidy, ploidy_step),
    cellularity = seq(min_cellularity, max_cellularity, cellularity_step)
  )

  distances <- pmap_dbl(grid,
                        absolute_copy_number_distance,
                        relative_copy_numbers = relative_copy_numbers,
                        weights = weights,
                        distance_function = distance_function)

  mutate(grid, distance = distances)
}

#' Find local minimum in distance function on fitting relative copy numbers
#'
#' Finds the local minimum in the distance function from scaling relative copy
#' numbers to absolute copy number starting from the given ploidy and
#' cellularity.
#'
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' find_minimum(3, 0.67, copy_number$segmented)
#' @import tibble
#' @export
find_minimum <- function(ploidy, cellularity,
                         relative_copy_numbers, weights = NULL,
                         distance_function = c("MAD", "RMSD")) {
  minimum <- nlm(
    function(params, relative_copy_numbers, weights, distance_function) { absolute_copy_number_distance(params[1], params[2], relative_copy_numbers, weights, distance_function) },
    c(ploidy, cellularity),
    relative_copy_numbers,
    weights,
    distance_function
  )
  tibble(ploidy = minimum$estimate[1], cellularity = minimum$estimate[2], distance = minimum$minimum)
}

#' Determine if a given ploidy and cellularity for fitting copy numbers is
#' acceptable
#'
#' Determines if the given ploidy and cellularity are compatible with the
#' provided relative copy numbers based on the proportion of copy number values
#' that correspond to or are below the zero copy number state and on the
#' proportion that are sufficiently close to a whole number copy number state.
#'
#' @param ploidy the tumour ploidy.
#' @param cellularity the cellularity, i.e. the fraction of cells that are from
#' the tumour.
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param max_proportion_zero the maximum proportion of fitted absolute copy
#' number values in the zero copy number state.
#' @param min_proportion_close_to_whole_number the minimum proportion of fitted
#' absolute copy number values sufficiently close to a whole number.
#' @param max_distance_from_whole_number the maximum distance from a whole
#' number that a fitted absolute copy number can be to be considered
#' sufficiently close.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#' is_acceptable_ploidy_and_cellularity(3, 0.67, copy_number$segmented)
#' @export
is_acceptable_ploidy_and_cellularity <- function(ploidy, cellularity,
                                                 relative_copy_numbers, weights = NULL,
                                                 max_proportion_zero = 0.05,
                                                 min_proportion_close_to_whole_number = 0.5,
                                                 max_distance_from_whole_number = 0.15) {

  absolute_copy_numbers <- relative_to_absolute_copy_number(relative_copy_numbers, ploidy, cellularity)
  absolute_copy_number_steps <- round(absolute_copy_numbers)

  if (is.null(weights)) weights <- rep(1, length(absolute_copy_numbers))
  sum_of_weights <- sum(weights)

  # filter based on proportion at or below copy number state 0
  prop_zero <- (sum(weights[which(absolute_copy_number_steps <= 0)]) / sum_of_weights)
  if (prop_zero > max_proportion_zero) return(FALSE)

  # filter based on proportion not sufficiently close to whole number copy number state
  prop_state <- (sum(weights[which(abs(absolute_copy_numbers - absolute_copy_number_steps) < max_distance_from_whole_number)]) / sum_of_weights)
  if (prop_state < min_proportion_close_to_whole_number) return(FALSE)

  return(TRUE)
}

#' Find best fit solutions for fitting copy numbers through grid-based search
#'
#' Find best fit solutions for a grid-based search through the given ranges of
#' ploidies and cellularities.
#'
#' @param relative_copy_numbers a numeric vector containing relative copy
#' numbers, i.e. ratios of copy numbers to the average copy number.
#' @param weights a numeric vector of weights to apply to each copy number value
#' (should be same length as relative_copy_numbers)
#' @param min_ploidy,max_ploidy the range of ploidies.
#' @param ploidy_step the stepwise increment of ploidies along the grid.
#' @param min_cellularity,max_cellularity the range of cellularities.
#' @param cellularity_step the stepwise increment of cellularities along the
#' grid.
#' @param distance_function the distance function to use, either "MAD" for the
#' mean absolute difference or "RMSD" for the root mean square difference, where
#' differences are between the fitted absolute copy number values and the
#' nearest whole number.
#' @param distance_filter_scale_factor the distance threshold above which
#' solutions will be discarded as a multiple of the solution with the smallest
#' distance.
#' @param max_proportion_zero the maximum proportion of fitted absolute copy
#' number values in the zero copy number state.
#' @param min_proportion_close_to_whole_number the minimum proportion of fitted
#' absolute copy number values sufficiently close to a whole number.
#' @param max_distance_from_whole_number the maximum distance from a whole
#' number that a fitted absolute copy number can be to be considered
#' sufficiently close.
#' @param solution_proximity_threshold how close two solutions can be before one
#' will be filtered; reduces the number of best fit solutions where there are
#' many minima in close proximity.
#' @param keep_all set to \code{TRUE} to return all solutions but with
#' additional \code{best_fit} column to indicate which are the local minima that
#' are acceptable solutions (may be useful to avoid computing the distance grid
#' twice)
#' @return the distance in the fitted absolute copy numbers to whole numbers.
#' @examples
#' data(copy_number)
#' copy_number <- copy_number[copy_number$sample == "X17222", ]
#'
#' segments <- copy_number_segments(copy_number)
#'
#' solutions <- find_best_fit_solutions(
#'   segments$copy_number, segments$weight,
#'   min_ploidy = 1.5, max_ploidy = 5.5,
#'   distance_function = "MAD")
#' @import dplyr
#' @export
find_best_fit_solutions <- function(relative_copy_numbers, weights = NULL,
                                    min_ploidy = 1.5, max_ploidy = 5.5, ploidy_step = 0.01,
                                    min_cellularity = 0.2, max_cellularity = 1.0, cellularity_step = 0.01,
                                    distance_function = c("MAD", "RMSD"),
                                    distance_filter_scale_factor = 1.25,
                                    max_proportion_zero = 0.05,
                                    min_proportion_close_to_whole_number = 0.5,
                                    max_distance_from_whole_number = 0.15,
                                    solution_proximity_threshold = 5,
                                    keep_all = FALSE) {

  # compute distance function for grid-based range of ploidies and cellularities
  distances <- absolute_copy_number_distance_grid(
    relative_copy_numbers, weights,
    min_ploidy, max_ploidy, ploidy_step,
    min_cellularity, max_cellularity, cellularity_step,
    distance_function)

  # add identifier for each solution to help with marking up the best fit
  # solutions if returning all solutions
  distances <- distances %>%
    mutate(id = row_number())

  # find minima, i.e. those grid points which have a smaller distance than all
  # adjacent points
  distances <- distances %>%
    mutate(x = dense_rank(cellularity)) %>%
    mutate(y = dense_rank(ploidy))

  solutions <- select(distances, id, x, y, distance)

  for (xdelta in -1:1) {
    for (ydelta in -1:1) {
      if (xdelta != 0 || ydelta != 0) {
        solutions <- solutions %>%
          mutate(xc = x + xdelta, yc = y + ydelta) %>%
          left_join(select(distances, xc = x, yc = y, dc = distance), by = c("xc", "yc")) %>%
          filter(is.na(dc) | distance <= dc) %>%
          select(id, x, y, distance)
      }
    }
  }

  distances <- select(distances, -x, -y)

  solutions <- solutions %>%
    select(id) %>%
    left_join(distances, by = "id")

  # only retain solutions with distances no more than
  # distance_filter_scale_factor times the the minimum value
  if (is.numeric(distance_filter_scale_factor) && nrow(solutions) > 1) {
    solutions <- solutions %>%
      filter(distance < distance_filter_scale_factor * min(distance))
  }

  # retain solutions with acceptable ploidies and cellularities
  solutions <- solutions %>%
    rowwise() %>%
    filter(
      is_acceptable_ploidy_and_cellularity(
        ploidy, cellularity,
        relative_copy_numbers, weights,
        max_proportion_zero = max_proportion_zero,
        min_proportion_close_to_whole_number = min_proportion_close_to_whole_number,
        max_distance_from_whole_number = max_distance_from_whole_number
      )
    ) %>%
    ungroup()

  # remove redundant solutions that have very similar ploidy and cellularity to
  # another solution with a smaller distance
  if (is.numeric(solution_proximity_threshold) && nrow(solutions) > 1) {

    solution_proximity_threshold2 <- solution_proximity_threshold ^ 2

    solutions <- solutions %>%
      arrange(distance) %>%
      mutate(rank = row_number()) %>%
      mutate(x = cellularity * 100) %>%
      mutate(y = (ploidy - min(ploidy)) * 100 / (max(ploidy) - min(ploidy)))

    pairs <- expand_grid(rank1 = 1:nrow(solutions), rank2 = 1:nrow(solutions)) %>%
      filter(rank2 > rank1) %>%
      left_join(select(solutions, rank1 = rank, x1 = x, y1 = y, distance1 = distance), by = "rank1") %>%
      left_join(select(solutions, rank2 = rank, x2 = x, y2 = y, distance2 = distance), by = "rank2") %>%
      mutate(xy2 = (x2 - x1) ^ 2 + (y2 - y1) ^ 2)

    to_remove <- tibble(rank = integer(0))

    while (nrow(pairs) > 0) {

      rank <- min(pairs$rank1)

      proximal <- pairs %>%
        filter(rank1 == rank, xy2 <= solution_proximity_threshold2) %>%
        select(rank = rank2)

      to_remove <- bind_rows(to_remove, proximal)

      pairs <- pairs %>%
        filter(rank1 != rank) %>%
        anti_join(proximal, by = c("rank1" = "rank"))
    }

    solutions <- solutions %>%
      anti_join(to_remove, by = "rank") %>%
      select(id, ploidy, cellularity, distance)
  }


  # return either the entire distance grid or just the selected solutions
  if (keep_all) {
    distances %>%
      mutate(best_fit = id %in% solutions$id) %>%
      select(ploidy, cellularity, distance, best_fit)
  } else {
    select(solutions, ploidy, cellularity, distance)
  }
}
