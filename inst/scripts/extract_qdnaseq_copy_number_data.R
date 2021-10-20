#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), dest = "input_file",
              help = "RDS file containing QDNAseqCopyNumbers object"),

  make_option(c("-o", "--output"), dest = "output_file",
              help = "Comma-separated value (CSV) output file"),

  make_option(c("-a", "--append"), dest = "append", action = "store_true", default = FALSE,
              help="Append to existing output file"),

  make_option(c("-s", "--sample"), dest = "sample",
              help = "Name of sample to extract copy number data for (all samples if unset)"),

  make_option(c("-m", "--median-scale"), dest = "median_scale", action = "store_true", default = FALSE,
              help = "Scale copy numbers by dividing by the median segmented value for each sample"),

  make_option(c("-d", "--digits"), dest = "digits", default = 3,
              help = "Number of digits to round copy number values to (default: %default)")
)

option_parser <- OptionParser(usage = "usage: %prog [options]", option_list = option_list, add_help_option = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- "--help"

opt <- parse_args(option_parser, args)

input_file <- opt$input_file
output_file <- opt$output_file
append <- opt$append
sample <- opt$sample
median_scale <- opt$median_scale
digits <- opt$digits

if (is.null(input_file)) stop("QDNAseqCopyNumbers RDS file must be specified")
if (is.null(output_file)) stop("Output file must be specified")

library(tibble)
library(readr)
library(dplyr)
library(Biobase)

copy_number <- readRDS(input_file)
if (!is(copy_number, "QDNAseqCopyNumbers")) stop(input_file, "does not contain a QDNAseqCopyNumbers object")

samples <- sort(sampleNames(copy_number))
if (!is.null(sample)) {
  if (!sample %in% samples) stop("Sample ", sample, " not found in copy number data")
  samples <- sample
}

locations <- fData(copy_number) %>%
  rownames_to_column(var = "id") %>%
  as_tibble() %>%
  select(id, chromosome, start, end) %>%
  mutate(across(c(start, end), as.integer)) %>%
  mutate(chromosome = factor(chromosome, levels = unique(chromosome)))

number_of_samples <- length(samples)
if (number_of_samples > 1) message("Samples: ", number_of_samples)

for (sample_index in 1:number_of_samples) {
  sample <- samples[sample_index]
  message(sample_index, "/", number_of_samples, " ", sample)

  copy_number_for_sample <- copy_number[, sample]

  copy_number_values <- assayDataElement(copy_number_for_sample, "copynumber")
  if (!all(rownames(copy_number_values) == locations$id)) stop("Error extracting copy number")

  segmented_values <- assayDataElement(copy_number_for_sample, "segmented")
  if (!all(rownames(segmented_values) == locations$id)) stop("Error extracting segmented copy number")

  copy_number_for_sample <- locations %>%
    mutate(sample = sample, copy_number = copy_number_values[,1], segmented = segmented_values[,1])

  if (median_scale) {
    copy_number_for_sample <- copy_number_for_sample %>%
      mutate(across(c(copy_number, segmented), ~ . / median(segmented, na.rm = TRUE)))
  }

  copy_number_for_sample <- copy_number_for_sample %>%
    mutate(across(c(copy_number, segmented), round, digits)) %>%
    mutate(copy_number = ifelse(copy_number == 0, 0, copy_number)) %>%
    mutate(segmented = ifelse(segmented == 0, 0, segmented)) %>%
    select(sample, chromosome, start, end, copy_number, segmented)

  print(summary(copy_number_for_sample$segmented))

  write_csv(copy_number_for_sample, output_file, append = append)

  append <- TRUE
}

