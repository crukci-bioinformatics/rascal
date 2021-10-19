#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--input-file"), dest = "input_file",
              help = "RDS, CSV or tab-delimited file containing copy number table with sample, chromosome, start, end copy_number and segmented columns or RDS file containing QDNAseqCopyNumbers object"),

  make_option(c("-o", "--output-file"), dest = "output_file",
              help = "Comma-separated value (CSV) output file"),

  make_option(c("-s", "--sample"), dest = "sample",
              help = "Name of sample to perform copy number fitting on (all samples if unset)"),

  make_option(c("--min-ploidy"), dest = "min_ploidy", default = 1.25,
              help = "The minimum ploidy to consider (default: %default)"),

  make_option(c("--max-ploidy"), dest = "max_ploidy", default = 5.25,
              help = "The maximum ploidy to consider (default: %default)"),

  make_option(c("--min-cellularity"), dest = "min_cellularity", default = 0.2,
              help = "The minimum cellularity to consider (default: %default)"),

  make_option(c("--max-cellularity"), dest = "max_cellularity", default = 1.0,
              help = "The maximum cellularity to consider (default: %default)")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

opt <- parse_args(option_parser, args)

input_file <- opt$input_file
output_file <- opt$output_file
sample <- opt$sample
min_ploidy <- opt$min_ploidy
max_ploidy <- opt$max_ploidy
min_cellularity <- opt$min_cellularity
max_cellularity <- opt$max_cellularity

if (is.null(input_file)) stop("Copy number RDS/CSV/TSV file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

library(tibble)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(rascal)

# function for extracting the copy number for a given sample
# can handle QDNAseqCopyNumbers object or a copy number data frame
copy_number_for_sample <- function(copy_number, sample) {
  if (any(class(copy_number) == "QDNAseqCopyNumbers")) {
    copy_number <- copy_number[,sample]
    copy_number_values <- Biobase::assayDataElement(copy_number, "copynumber")[,1]
    segmented_values <- Biobase::assayDataElement(copy_number, "segmented")[,1]
    Biobase::fData(copy_number) %>%
      rownames_to_column(var = "id") %>%
      as_tibble() %>%
      select(id, chromosome, start, end) %>%
      mutate(across(c(start, end), as.integer)) %>%
      mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
      mutate(sample = sample) %>%
      mutate(copy_number = copy_number_values) %>%
      mutate(segmented = segmented_values) %>%
      select(sample, chromosome, start, end, copy_number, segmented)
  } else {
    filter(copy_number, sample == !!sample)
  }
}

# read copy number data from input file
if (str_detect(input_file, "\\.rds$")) {
  copy_number <- readRDS(input_file)
} else if (str_detect(input_file, "\\.csv(\\.gz)?$")) {
  copy_number <- read_csv(input_file, col_types = cols(sample = "c", chromosome = "f", start = "i", end = "i", copy_number = "d", segmented = "d"))
} else {
  copy_number <- read_tsv(input_file, col_types = cols(sample = "c", chromosome = "f", start = "i", end = "i", copy_number = "d", segmented = "d"))
}

# check contents are correct and obtain sample names
if (any(class(copy_number) == "data.frame")) {
  required_columns <- c("sample", "chromosome", "start", "end", "segmented")
  missing_columns <- setdiff(required_columns, colnames(copy_number))
  if (length(missing_columns) > 0) stop("missing columns in", input_file, ": ", str_c(missing_columns, collapse = ", "))
  samples <- sort(unique(copy_number$sample))
} else if (any(class(copy_number) == "QDNAseqCopyNumbers")) {
  if (!requireNamespace("QDNAseq", quietly = TRUE)) stop("QDNAseq package needs to be installed")
  samples <- sort(Biobase::sampleNames(copy_number))
} else {
  stop(input_file, " should contain either a data frame or a QDNAseqCopyNumbers object")
}

# check specified sample is found in the copy number data
if (!is.null(sample)) {
  if (sample %in% samples) {
    samples <- sample
  } else {
    stop(sample, " not found in ", input_file)
  }
}

all_solutions <- NULL

number_of_samples <- length(samples)
append <- FALSE

for (sample_index in 1:number_of_samples) {
  sample <- samples[sample_index]

  sample_copy_number <- copy_number_for_sample(copy_number, sample)

  # copy number fitting requires relative copy numbers where values are relative
  # to the average copy number across the genome - using the median segmented
  # copy number
  relative_copy_number <- mutate(sample_copy_number, across(c(copy_number, segmented), ~ . / median(segmented, na.rm = TRUE)))

  segments <- copy_number_segments(relative_copy_number)

  solutions <- find_best_fit_solutions(
    segments$copy_number, segments$weight,
    min_ploidy = min_ploidy, max_ploidy = max_ploidy, ploidy_step = 0.01,
    min_cellularity = min_cellularity, max_cellularity = max_cellularity, cellularity_step = 0.01,
    distance_function = "MAD")

  message(sample_index, "/", number_of_samples, " ", sample, " ", nrow(solutions))

  if (nrow(solutions) == 0) next

  solutions <- solutions %>%
    transmute(sample = sample, ploidy, cellularity, distance)

  solutions %>%
    mutate(across(c(ploidy, cellularity), round, digits = 2)) %>%
    mutate(across(distance, round, digits = 3)) %>%
    write_csv(output_file, append = append)

  append <- TRUE

  all_solutions <- bind_rows(all_solutions, solutions)
}

message("Solutions found for ", length(unique(all_solutions$sample)), " of ", length(samples), " samples")

all_solutions %>%
  count(sample) %>%
  print(n = Inf)

