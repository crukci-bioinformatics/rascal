#' Sample copy number data set
#'
#' A dataset containing copy number data for 2 samples obtained through shallow
#' whole genome sequencing and processing using 30kb bins.
#'
#' Each row represents one bin, i.e. a window of 30 kilobases. The dataset
#' contains both pre- and post-segmented copy number values. These are relative
#' copy numbers, i.e. ratios of copy numbers to the average copy number, and
#' have not been log2-transformed.
#'
#' @format A data frame with 206398 rows and 6 variables:
#' \describe{
#'   \item{sample}{the name of the sample}
#'   \item{chromosome}{the chromosome}
#'   \item{start}{the start position of the bin within the chromosome}
#'   \item{end}{the end position of the bin within the chromosome}
#'   \item{copy_number}{the relative copy number for the bin}
#'   \item{segmented}{the segmented relative copy number for the bin}
#' }
"copy_number"

#' Sample set of genes
#'
#' A set of genes and their locations to illustrate display of genes in the
#' chromosome copy number plot.
#'
#' @format A data frame with 4 rows and 4 variables:
#' \describe{
#'   \item{name}{the gene name or symbol}
#'   \item{chromosome}{the chromosome}
#'   \item{start}{the start position of the gene within the chromosome}
#'   \item{end}{the end position of the gene within the chromosome}
#' }
"genes"
