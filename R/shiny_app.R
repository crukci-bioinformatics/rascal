#' Shiny app for absolute copy number scaling
#'
#' Starts the shiny app in a new tab in the the default web browser.
#' If a browser tab isn't opened it may be necessary to manually open the
#' URL shown in the console.
#'
#' @param load_sample_data load a sample data set on starting the Shiny app.
#' @param launch_browser launch the application in a web browser on starting the
#' Shiny app.
#' @param port the TCP port that the application should listen on.
#' @return opens the Shiny app in a web browser.
#' @import shiny
#' @export
start_shiny_app <- function(load_sample_data = TRUE, launch_browser = TRUE, port = getOption("shiny.port")) {
  dir <- system.file("shiny", package = "rascal")
  if (dir == "") stop("Could not find shiny app directory in the rascal package")
  shiny::shinyOptions(load_sample_data = load_sample_data)
  shiny::runApp(dir, port = port, launch.browser = launch_browser, display.mode = "normal")
}
