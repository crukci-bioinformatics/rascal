library(shiny)
library(shinyjs)
library(colourpicker)
library(tibble)
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(DT)
library(rascal)

options(shiny.maxRequestSize = 1024 * 1024 * 1024)
options(scipen = 999)

ui <- fluidPage(
  useShinyjs(),
  tags$style(type = "text/css", "body { padding-top: 60px; }"),
  # tags$style(type = "text/css", ".navbar-default { background-color: #c2d1f0; }"),
  navbarPage(
    title = a(href="https://www.cruk.cam.ac.uk", target = "_blank", img(style = "width: 180px; margin-right: 10px;", src = "CRUK_CI_logo.png")),
    windowTitle = "rascal - absolute copy number scaling",
    position = "fixed-top",
    tabPanel(
      title = strong(style = "color: #231F7F; font-size: 120%;", HTML("<em>rascal</em> &ndash; absolute copy number scaling")),
      div(style = "margin-top: 10px;"),
      fluidRow(
        column(
          width = 3,
          fileInput(
            'copy_number_file',
            'Copy number file upload',
            accept = c(".rds", ".txt", ".tsv", ".csv", ".gz"),
            width = "100%"
          )
        ),
        column(
          width = 2,
          selectInput("sample", label = "Sample", choices = list())
        ),
        column(
          width = 2,
          selectInput("gene", label = "Gene", choices = list())
        ),
        column(
          width = 2,
          numericInput("ploidy", label = "Ploidy", value = NULL, min = 1.25, max = 5.25, step = 0.01)
        ),
        column(
          width = 2,
          numericInput("cellularity", label = "Cellularity", value = NULL, min = 0.2, max = 1.0, step = 0.01)
        ),
        column(
          width = 1,
          actionButton(style = "margin-top: 25px;", "clear_all_button", label = HTML("Clear all"))
        )
      ),
      fluidRow(
        column(
          width = 5,
          h5(htmlOutput("genome_copy_number_plot_label"))
        ),
        column(
          width = 1,
          downloadButton("save_genome_copy_plot", label = "PDF")
        ),
        column(
          width = 5,
          h5(htmlOutput("location_label"))
        ),
        column(
          width = 1,
          downloadButton("save_chromosome_copy_plot", label = "PDF")
        )
      ),
      fluidRow(
        column(
          width = 6,
          plotOutput(
            "genome_copy_number_plot",
            hover = hoverOpts(id = "genome_copy_number_plot_hover", delay = 50, delayType = "throttle"),
            click = "genome_copy_number_plot_click",
            height = "350px"
          ) #%>% addSpinner(spin = "dots", color = "grey")
        ),
        column(
          width = 6,
          plotOutput(
            "chromosome_copy_number_plot",
            hover = hoverOpts(id = "chromosome_copy_number_plot_hover", delay = 50, delayType = "throttle"),
            brush = "chromosome_copy_number_plot_brush",
            dblclick = "chromosome_copy_number_plot_dblclick",
            height = "350px"
          )
        )
      ),
      htmlOutput(
        style = "margin-top: 10px;",
        "copy_number_plot_hover_over_label"
      ),
      fluidRow(
        style = "margin-top: 10px;",
        column(
          width = 2,
          selectInput("distance_function", label = "distance function", choices = c("MAD", "RMSD"))
        ),
        column(
          width = 2,
          selectInput("copy_numbers_to_be_fitted", label = "applied to", choices = c("segments", "maxima"))
        ),
        column(
          width = 2,
          conditionalPanel(
            condition = "input.copy_numbers_to_be_fitted == 'maxima'",
            sliderInput("number_of_maxima", "number of maxima", value = 3, min = 2, max = 8, step = 1, ticks = FALSE)
          )
        ),
        column(
          width = 1,
          actionButton(style = "margin-top: 0px;", "cache_current_ploidy_and_cellularity_button", label = HTML("&#8595; Store")),
          actionButton(style = "margin-top: 10px;", "restore_cached_ploidy_and_cellularity_button", label = HTML("&#8593; Restore")),
          p()
        ),
        column(
          width = 5,
          htmlOutput(style = "margin-top: 7px;", "current_ploidy_and_cellularity_label"),
          htmlOutput(style = "margin-top: 25px;", "cached_ploidy_and_cellularity_label")
        )
      ),
      fluidRow(
        style = "margin-top: 10px;",
        column(
          width = 6,
          h5("Segmented copy number density"),
          plotOutput("segmented_copy_number_density_plot", height = "350px")
        ),
        column(
          width = 6,
          h5("Distance heat map"),
          plotOutput(
            "distance_heat_map_plot",
            click = "distance_heat_map_plot_click",
            dblclick = "distance_heat_map_plot_dblclick",
            hover = hoverOpts(id = "distance_heat_map_plot_hover", delay = 50, delayType = "throttle"),
            height = "350px"
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          HTML("&nbsp;")
        ),
        column(
          width = 6,
          htmlOutput("distance_heat_map_plot_hover_over_label")
        )
      ),
      fluidRow(
        style = "margin-top: 10px;",
        # column(
        #   width = 6,
        #   h5("Segmented copy number maxima"),
        #   DT::dataTableOutput("segmented_copy_number_maxima_table", width = "80%")
        # ),
        column(
          width = 6,
          offset = 3,
          align = "center",
          h5("Best fit solutions"),
          DT::dataTableOutput("distance_best_fit_solution_table", width = "80%")
        )
      )
    ),
    tabPanel(
      title = "Ploidy/cellularity cache",
      div(style = "margin-top: 10px;"),
      fluidRow(
        column(
          width = 4,
          fileInput(
            'ploidy_and_cellularity_file',
            'Ploidy/cellularity file upload',
            accept = c(".txt", ".tsv", ".csv"),
            width = "100%"
          )
        ),
        column(
          width = 1,
          actionButton(style = "margin-top: 25px;", "clear_cache_button", label = HTML("Clear"))
        ),
        column(
          width = 1,
          style = "margin-top: 25px;",
          downloadButton("save_cached_ploidies_and_cellularities", label = "Save")
        )
      ),
      fluidRow(
        column(
          width = 5,
          DT::dataTableOutput("cached_ploidy_and_cellularity_table")
        )
      )
    ),
    tabPanel(
      title = "Genes",
      div(style = "margin-top: 10px;"),
      fluidRow(
        column(
          width = 4,
          fileInput(
            'genes_file',
            'Genes file upload',
            accept = c(".txt", ".tsv", ".csv"),
            width = "100%"
          )
        ),
        column(
          width = 1,
          actionButton(style = "margin-top: 25px;", "clear_genes_button", label = HTML("Clear"))
        )
      ),
      fluidRow(
        column(
          width = 5,
          DT::dataTableOutput("genes_table")
        )
      )
    ),
    navbarMenu("More",
    tabPanel(
      title = "About",
      fluidRow(
        column(
          width = 8,
          h4(HTML("<em>rascal</em> (<u><b>r</b></u>elative to <u><b>a</b></u>bsolute copy number <u><b>scal</b></u>ing)")),
          em("Shiny app for scaling relative copy number data from shallow whole genome sequencing of cancer samples to absolute values and estimating the tumour ploidy and cellularity of the samples."),
          p(),
          "Several research groups at CRUK CI are using shallow whole genome sequencing as a relatively inexpensive method for obtaining copy number profiles for tumour samples, particularly as libraries from several samples can be multiplexed in a single lane of sequencing.",
          p(),
          "We are principally using",
          tags$a(href = "https://bioconductor.org/packages/release/bioc/html/QDNAseq.html", target = "_blank", "QDNAseq"),
          "for summing reads that align within genomic windows or bins, typically 30kb in size, and correcting for GC-content and mappability.",
          "This results in values that are relative to the average copy number within the sample for the GC and mappability of each bin. These relative copy numbers are smoothed and segmented and provide useful insight into genomic abnormalities in cancers.",
          p(),
          "For some research projects it is desirable to obtain absolute copy numbers.",
          "Normally this would require deeper whole genome sequencing from which allele fractions of germline SNPs can help determine the clonal architecture of a tumour sample.",
          "In the absence of such information, and noting the significant increase in cost for deeper sequencing, we can attempt to fit the relative copy number profiles to absolute copy numbers by evaluating various estimates of ploidy and cellularity.",
          p(),
          "The approach used in this application is based on concepts introduced in the",
          tags$a(href = "https://bioconductor.org/packages/release/bioc/html/ACE.html", target = "_blank", "ACE"),
          "package developed by Bauke Ylstra's group at Amsterdam UMC.",
          "The mathematics underpinning this approach assume a single dominant clone; estimating ploidy and cellularity for heterogenous tumour samples may prove difficult with this method.",
          p(),
          "This application was created using the R Shiny web application framework. It was developed by",
          tags$a(href = "https://www.cruk.cam.ac.uk/author/matthew-eldridge", target = "_blank", "Matt Eldridge"),
          "in the",
          tags$a(href = "https://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core", target = "_blank", "Bioinformatics Core"),
          "in collaboration with the",
          tags$a(href = "https://www.cruk.cam.ac.uk/research-groups/brenton-group", target = "_blank", "James Brenton's laboratory"),
          "at the",
          tags$a(href = "https://www.cruk.cam.ac.uk", target = "_blank", "Cancer Research UK Cambridge Institute.")
        )
      )
    ),
    tabPanel(
      title = "Help",
      fluidRow(
        column(
          width = 8,
          h4("User guide"),
          hr(),
          h4("Main page"),
          p(),
          "Upload a tab-delimited, CSV or R data object file (.rds) containing a copy number table (or data frame in the case of an .rds file) by clicking the", strong("Browse"), "button on the main page.",
          "The following columns are expected:",
          p(),
          tags$ul(
            tags$li(strong(em("sample")), "(optional)"),
            tags$li(strong(em("chromosome"))),
            tags$li(strong(em("start"))),
            tags$li(strong(em("end"))),
            tags$li(strong(em("copy_number")), "(optional)"),
            tags$li(strong(em("segmented")))
          ),
          "A single unnamed sample will be assumed if there is no", em("sample"), "column.",
          p(),
          "Each row in the table should correspond to a bin (or window) or a wider continuous copy number segment following segmentation.",
          "Values should be relative copy numbers that have not been log2-transformed.",
          "Segmented copy number values are required as these are used in fitting to absolute copy numbers.",
          "Copy number values for individual bins are optional but can be helpful in assessing how well the segmentation performed and showing the level of noise in the data.",
          p(),
          "Alternatively an R data object file (.rds) containing a QDNAseqCopyNumbers object obtained from processing shallow whole genome sequencing data using QDNAseq can also be uploaded.",
          p(),
          "Select a sample to view from the drop-down list.",
          p(),
          "Click on a chromosome in the whole genome copy number plot (left-hand side) to display the copy number for that chromosome on the right-hand side.",
          "Zoom in to a specific region on a chromosome by clicking and dragging to select the region in the chromosome copy number plot; double-click to zoom out again and view the whole chromosome.",
          p(),
          "Hover over a location to display the copy number, log2 ratio, fitted absolute copy number and tumour DNA fraction at this locus.",
          p(),
          "The tumour fraction is the fraction of tumour DNA at that location given the cellularity and absolute copy number. For example, a sample with cellularity 0.5 (50% tumour and 50% normal) would have a tumour fraction of 0.5 if the absolute copy number at that position in both the tumour and the normal is 2, or a fraction of 0.6 if the absolute copy number in the tumour is 3.",
          p(),
          "Select a ploidy and cellularity using the selectors at the top of the main page or by clicking on a point within the distance heatmap.",
          "The distance heatmap shows how well different choices of ploidy and cellularity scale the relative copy number data to whole numbers on the absolute copy number scale.",
          p(),
          "Best fit solutions are displayed as points in the heatmap and listed in the table below the heatmap.",
          "Select a solution to update the currently selected ploidy and cellularity.",
          p(),
          "Specify the distance function (mean absolute difference or root mean square difference) from the drop-down list.",
          "This is applied to segmented copy number values with the following options for which values to use:",
          p(),
          tags$ul(
            tags$li(strong(em("segments")), HTML("&mdash;"), "relative copy number values for each segment weighted by the size of the segment"),
            tags$li(strong(em("maxima")), HTML("&mdash;"), "relative copy number values for each peak in the segmented copy number density plot, each given equal weight; the selected number of the most frequently observed relative copy number states (maxima) are used")
          ),
          p(),
          "The selected ploidy and cellularity can be stored in a cache by clicking on the", strong("Store"), "button.",
          "Click on the", strong("Restore"), "button to select the ploidy and cellularity currently stored in the cache.",
          p(),
          "The copy number plots can be saved as PDF image files using the", strong("PDF"), "buttons.",
          hr(),
          h4("Ploidy/cellularity cache page"),
          p(),
          "The cached ploidies and cellularities for each sample are displayed on the", strong(em("Ploidy/cellularity cache")), "page.",
          "Cached ploidies and cellularities can be saved as a CSV file by clicking on the", strong("Save"), "button.",
          "Previously saved (or otherwise determined) ploidies and cellularities can be loaded from a tab-delimited or CSV file by clicking the", strong("Browse"), "button.",
          hr(),
          h4("Genes page"),
          p(),
          "A set of genes and their locations can be loaded on the", strong(em("Genes")), "page.",
          p(),
          "Genes are displayed on the chromosome copy number plot as vertical bars.",
          p(),
          "Selecting a gene from the table on this page or in the drop-down on the main page will display the copy number plot for the chromosome on which the gene is located.",
          "The tumour fraction for the selected gene will also be displayed alongside each of the best fit solutions to help in deciding which solution is most consistent with other supporting data, e.g. allele fraction for a homozygous variant in that gene from digital PCR or amplicon sequencing. However, this is only the case where there is a single absolute copy number fitted across the entire length of the gene.",
          hr(),
          h4("Settings page"),
          p(),
          "Various display settings can be adjusted on the", strong(em("Settings")), "page."
        )
      )
    ),
    tabPanel(
      title = "Settings",
      fluidRow(
        h4("Settings"),
        hr(),
        column(
          width = 3,
          h4("Copy number plots"),
          div(style = "margin-top: 20px;"),
          checkboxInput("fix_log2ratio_range", label = "Fix log2 ratio range", value = TRUE),
          conditionalPanel(
            condition = "input.fix_log2ratio_range",
            div(style = "margin-top: 20px;"),
            sliderInput("log2ratio_range", "log2 ratio range", value = c(-2, 2), min = -5, max = 5, step = 0.25, ticks = FALSE)
          ),
          checkboxInput("limit_number_of_points_to_display", label = "Limit number of points to display", value = TRUE),
          conditionalPanel(
            condition = "input.limit_number_of_points_to_display",
            div(style = "margin-top: 20px;"),
            sliderInput("max_points_to_display", "Maximum number of points to display", value = 10000, min = 0, max = 50000, step = 5000, ticks = FALSE)
          ),
          div(style = "margin-top: 20px;"),
          colourInput("bin_colour", "Colour of points", value = "black", returnName = TRUE),
          div(style = "margin-top: 20px;"),
          colourInput("segment_colour", "Colour of segments", value = "red", returnName = TRUE),
          div(style = "margin-top: 20px;"),
          colourInput("gene_colour", "Colour of genes", value = "darkgreen", returnName = TRUE),
          div(style = "margin-top: 30px;"),
          hr(),
          div(style = "margin-top: 30px;"),
          h4("PDF export"),
          div(style = "margin-top: 20px;"),
          numericInput("pdf_width", "Width (inches)", value = 10, min = 6, max = 12, step = 0.5),
          div(style = "margin-top: 20px;"),
          numericInput("pdf_height", "Height (inches)", value = 6, min = 4, max = 12, step = 0.5)
        ),
        column(
          offset = 1,
          width = 3,
          h4("Absolute copy number steps"),
          div(style = "margin-top: 20px;"),
          checkboxInput("show_absolute_copy_number", label = "Show absolute copy numbers", value = TRUE),
          conditionalPanel(
            condition = "input.show_absolute_copy_number",
            div(style = "margin-top: 20px;"),
            sliderInput("max_absolute_copy_number_step", "Largest displayed copy number step", value = 8, min = 4, max = 12, step = 1, ticks = FALSE),
            div(style = "margin-top: 20px;"),
            colourInput("absolute_copy_number_step_colour", "Colour of copy number steps", value = "blue", returnName = TRUE)
          ),
          div(style = "margin-top: 30px;"),
          hr(),
          div(style = "margin-top: 30px;"),
          h4("Filtering options for best fit solutions"),
          div(style = "margin-top: 20px;"),
          sliderInput("distance_filter_scale_factor", "Distance threshold as multiple of lowest observed value", value = 1.25, min = 1, max = 2, step = 0.05, ticks = FALSE),
          div(style = "margin-top: 20px;"),
          sliderInput("max_proportion_zero", "Maximum proportion of fitted copy numbers in the zero copy number state", value = 0.1, min = 0, max = 0.2, step = 0.01, ticks = FALSE),
          div(style = "margin-top: 20px;"),
          sliderInput("min_proportion_close_to_whole_number", "Minimum proportion of fitted copy numbers sufficiently close to whole number", value = 0.5, min = 0, max = 1, step = 0.05, ticks = FALSE),
          div(style = "margin-top: 20px;"),
          sliderInput("max_distance_from_whole_number", "Distance from whole number for fitted value to be considered sufficiently close", value = 0.15, min = 0, max = 0.25, step = 0.025, ticks = FALSE),
          div(style = "margin-top: 20px;"),
          sliderInput("solution_proximity_threshold", "Proximity threshold for two solutions below which one will be removed", value = 5, min = 0, max = 25, step = 1, ticks = FALSE)
        ),
        column(
          offset = 1,
          width = 3,
          h4("Copy number density plot"),
          div(style = "margin-top: 20px;"),
          sliderInput("relative_copy_number_range", "Relative copy number range", value = c(0, 2.5), min = 0, max = 5, step = 0.25, ticks = FALSE),
          div(style = "margin-top: 30px;"),
          hr(),
          div(style = "margin-top: 30px;"),
          h4("Distance heat map"),
          div(style = "margin-top: 20px;"),
          sliderInput("ploidy_range", "Range of ploidies", value = c(1.25, 5.25), min = 1, max = 10, step = 0.25, ticks = FALSE),
          div(style = "margin-top: 20px;"),
          sliderInput("cellularity_range", "Range of cellularities", value = c(0.25, 1.0), min = 0.05, max = 1.0, step = 0.05, ticks = FALSE),
          div(style = "margin-top: 20px;"),
          colourInput("heatmap_best_fit_point_colour", "Colour of best fit solutions", value = "black", returnName = TRUE),
          div(style = "margin-top: 20px;"),
          colourInput("heatmap_current_point_colour", "Colour of current solution", value = "orange", returnName = TRUE),
          div(style = "margin-top: 20px;"),
          colourInput("heatmap_low_distance_colour", "Low distance colour", value = "red", returnName = TRUE),
          div(style = "margin-top: 20px;"),
          colourInput("heatmap_high_distance_colour", "High distance colour", value = "blue", returnName = TRUE)
        )
      )
    )
    )
  ),
  tags$div(
    style = "clear:both",
    tags$div(style="line-height:200%;", br()),
    HTML("&copy;"),
    tags$script(type = "text/javascript", "var d = new Date(); document.write(d.getFullYear())"),
    "University of Cambridge",
    tags$div(style = "float:right", tags$a(href = "https://www.cruk.cam.ac.uk/terms-and-conditions", target = "_blank", "Terms and Conditions"))
  ),
  br()
)

server <- function(input, output, session) {

  reactive_values <- reactiveValues(
    copy_number_data = NULL,
    samples = NULL,
    sample = NULL,
    location = list(chromosome = NULL, start = NULL, end = NULL),
    ploidy = NA,
    cellularity = NA,
    ploidy_and_cellularity_file = NULL,
    ploidy_and_cellularity_cache = tibble(sample = character(), ploidy = numeric(), cellularity = numeric()),
    genes = tibble(name = character(), chromosome = character(), start = integer(), end = integer()),
    selected_gene = NULL
  )

  # handle load button click event by loading copy number data from file
  observe({
    file <- input$copy_number_file
    if (is.null(file)) return(NULL)
    load_copy_number_data(file)
  })

  # load copy number data from a file
  load_copy_number_data <- function(file, initialization = FALSE) {
    progress <- shiny::Progress$new()
    on.exit(progress$close())

    progress$set(value = 0.1, message = "Reading copy number data")

    # read copy number data from input file
    if (str_detect(file$name, "\\.rds$")) {
      copy_number_data <- readRDS(file$datapath)
    } else if (str_detect(file$name, "\\.csv(\\.gz)?$")) {
      copy_number_data <- read_csv(file$datapath, col_types = cols(sample = "c", chromosome = "f", start = "i", end = "i", copy_number = "d", segmented = "d"))
    } else {
      copy_number_data <- read_tsv(file$datapath, col_types = cols(sample = "c", chromosome = "f", start = "i", end = "i", copy_number = "d", segmented = "d"))
    }

    progress$set(value = 0.4, message = "Checking copy number data")

    # check contents are as expected and obtain sample names
    # if (is.data.frame(copy_number_data)) {
    if (any(class(copy_number_data) == "data.frame")) {

      expected_columns <- c("sample", "chromosome", "start", "end", "copy_number", "segmented")
      required_columns <- c("chromosome", "start", "end", "segmented")

      missing_columns <- setdiff(required_columns, colnames(copy_number_data))
      if (length(missing_columns) > 0) {
        showModal(modalDialog(title = "Error", strong(file$name), "is missing the following columns:", str_c(missing_columns, collapse = ", ")))
        return(NULL)
      }

      missing_sample_column <- !"sample" %in% colnames(copy_number_data)
      missing_copy_number_column <- !"copy_number" %in% colnames(copy_number_data)

      # handle situation where there isn't a sample column by assuming there is
      # just a single unnamed sample
      if (missing_sample_column) {
        copy_number_data <- mutate(copy_number_data, sample = "unknown sample")
      }

      # add copy_number column if it doesn't exist
      if (missing_copy_number_column) {
        copy_number_data <- mutate(copy_number_data, copy_number = NA)
      }

      copy_number_data <- select(copy_number_data, one_of(expected_columns))

      # check for missing values in sample, chromosome, start and end columns
      n <- nrow(copy_number_data)
      copy_number_data <- filter_at(copy_number_data, vars(sample, chromosome, start, end), all_vars(!is.na(.)))
      if (n != nrow(copy_number_data)) {
        showModal(modalDialog(title = "Error", strong(file$name), "contains rows with missing values for sample, chromosome, start and/or end."))
        return(NULL)
      }

      # check that the data frame contains some segmented copy number data
      if (nrow(filter(copy_number_data, !is.na(segmented))) == 0) {
        showModal(modalDialog(title = "Error", strong(file$name), "contains no segmented copy number data.."))
        return(NULL)
      }

      # sort within each sample by chromosome and position
      progress$set(value = 0.55, message = "Sorting copy number data")
      copy_number_data <- arrange(copy_number_data, sample, chromosome, start, end)

      # check for overlapping bins
      progress$set(value = 0.75, message = "Checking for overlapping bins")
      overlapping_bins <- filter(copy_number_data, row_number() > 1 & sample == lag(sample) & chromosome == lag(chromosome) & start <= lag(end))
      if (nrow(overlapping_bins) > 0) {
        message <- "contains overlapping or duplicate bins."
        if (missing_sample_column) {
          message <- str_c(message, " This could be because the data set contains multiple samples (note that the sample column is missing).")
        }
        showModal(modalDialog(title = "Error", strong(file$name), message))
        return(NULL)
      }

      samples <- sort(unique(copy_number_data$sample))
      chromosomes <- levels(copy_number_data$chromosome)

    } else if (class(copy_number_data) == "QDNAseqCopyNumbers") {
      if (!requireNamespace(package = "QDNAseq", quietly = TRUE)) {
        showModal(modalDialog(title = "Error", "The QDNAseq package needs to be installed in order to load a QDNAseqCopyNumbers object"))
        return(NULL)
      }

      assays <- Biobase::assayDataElementNames(copy_number_data)
      required_assays <- c("copynumber", "segmented")
      missing_assays <- setdiff(required_assays, assays)
      if (length(missing_assays) > 0) {
        showModal(modalDialog(title = "Error", strong(file$name), "is missing the following assay data elements:", str_c(missing_assays, collapse = ", ")))
        return(NULL)
      }

      samples <- sort(Biobase::sampleNames(copy_number_data))
      chromosomes <- unique(Biobase::fData(copy_number_data)$chromosome)

    } else {
      showModal(modalDialog(title = "Error", strong(file$name), "should contain either a data frame or a QDNAseqCopyNumbers object."))
      return(NULL)
    }

    progress$set(value = 0.9, message = "Updating current data")

    # update the reactive values
    reactive_values$copy_number_data <- copy_number_data
    reactive_values$samples <- samples

    first_sample <- samples[1]
    if (initialization) {
      reactive_values$sample <- NULL
    } else {
      reactive_values$sample <- first_sample
    }

    # defaults to first chromosome
    chromosome <- chromosomes[1]
    # if gene list exists see if a gene has been selected
    # otherwise select the chromosome of the first gene in the list
    genes <- isolate(reactive_values$genes)
    if (nrow(genes) > 0) {
      chromosome <- genes$chromosome[1]
      selected_gene <- isolate(reactive_values$selected_gene)
      if (!is.null(selected_gene)) {
        selected_gene <- filter(genes, name == selected_gene)
        if (nrow(selected_gene) == 1) chromosome <- selected_gene$chromosome
      }
    }
    reactive_values$location <- list(chromosome = chromosome, start = NULL, end = NULL)

    updateSelectInput(session, "sample", label = "Sample", choices = samples, selected = first_sample)

    cached_ploidy_and_cellularity <- get_cached_ploidy_and_cellularity(first_sample)
    update_ploidy_and_cellularity(cached_ploidy_and_cellularity$ploidy, cached_ploidy_and_cellularity$cellularity)
  }

  # clear all data including copy number, cached ploidies and cellularities, and genes
  observe({
    if (input$clear_all_button > 0) {
      reactive_values$copy_number_data <- NULL
      reactive_values$samples <- NULL
      reactive_values$sample <- NULL
      reactive_values$location <- list(chromosome = NULL, start = NULL, end = NULL)
      reactive_values$ploidy <- NA
      reactive_values$cellularity <- NA
      reactive_values$ploidy_and_cellularity_file <- NULL
      reactive_values$ploidy_and_cellularity_cache = tibble(sample = character(), ploidy = numeric(), cellularity = numeric())
      reactive_values$genes <- tibble(name = character(), chromosome = character(), start = integer(), end = integer())
      reactive_values$selected_gene <- NULL
      updateSelectInput(session, "sample", label = "Sample", choices = list())
      updateSelectInput(session, "gene", label = "Gene", choices = list())
      update_ploidy_and_cellularity(NA, NA)
    }
  })

  # respond to sample selection from drop down list
  observe({
    sample <- input$sample
    current_sample <- isolate(reactive_values$sample)
    if (sample != "" && (is.null(current_sample) || sample != current_sample)) {
      cached_ploidy_and_cellularity <- get_cached_ploidy_and_cellularity(sample)
      update_ploidy_and_cellularity(cached_ploidy_and_cellularity$ploidy, cached_ploidy_and_cellularity$cellularity)
      reactive_values$sample <- sample
    }
  })

  # respond to gene selection from drop down list
  observe({
    selected_gene <- input$gene
    if (selected_gene != "") {
      gene <- filter(isolate(reactive_values$genes), name == selected_gene)
      if (nrow(gene) == 1) {
        reactive_values$selected_gene <- selected_gene
        chromosome <- gene$chromosome
        chromosomes <- isolate(chromosomes_for_selected_sample())
        if (chromosome %in% chromosomes$chromosome) {
          reactive_values$location <- list(chromosome = chromosome, start = NULL, end = NULL)
        }
      }
    }
  })

  # the range of ploidies (handles case where minimum and maximum ploidy are the same)
  ploidy_range <- reactive({
    ploidy_range <- input$ploidy_range
    min_ploidy <- ploidy_range[1]
    max_ploidy <- ploidy_range[2]
    if (min_ploidy == max_ploidy) {
      min_ploidy <- min_ploidy - 0.125
      max_ploidy <- max_ploidy + 0.125
    }
    c(min_ploidy, max_ploidy)
  })

  # the range of cellularities (handles case where minimum and maximum cellularity are the same)
  cellularity_range <- reactive({
    cellularity_range <- input$cellularity_range
    min_cellularity <- cellularity_range[1]
    max_cellularity <- cellularity_range[2]
    if (min_cellularity == max_cellularity) {
      min_cellularity <- min_cellularity - 0.025
      max_cellularity <- max_cellularity + 0.025
    }
    c(min_cellularity, max_cellularity)
  })

  # respond to ploidy selection in numeric input control
  observe({
    ploidy <- input$ploidy
    current_ploidy <- isolate(reactive_values$ploidy)
    if (!is.na(ploidy) && (is.na(current_ploidy) || ploidy != current_ploidy)) {
      ploidy_range <- ploidy_range()
      min_ploidy <- ploidy_range[1]
      max_ploidy <- ploidy_range[2]
      if (ploidy < min_ploidy || ploidy > max_ploidy) {
        showModal(modalDialog(title = "Error", str_c("Ploidy values must be in the range ", min_ploidy, " - ", max_ploidy, ". The range can be adjusted in the settings page.")))
      } else {
        reactive_values$ploidy <- ploidy
      }
    }
  })

  # respond to cellularity selection in numeric input control
  observe({
    cellularity <- input$cellularity
    current_cellularity <- isolate(reactive_values$cellularity)
    if (!is.na(cellularity) && (is.na(current_cellularity) || cellularity != current_cellularity)) {
      cellularity_range <- cellularity_range()
      min_cellularity <- cellularity_range[1]
      max_cellularity <- cellularity_range[2]
      if (cellularity < min_cellularity || cellularity > max_cellularity) {
        showModal(modalDialog(title = "Error", str_c("Cellularity values must be in the range ", min_cellularity, " - ", max_cellularity, ". The range can be adjusted in the settings page.")))
      } else {
        reactive_values$cellularity <- cellularity
      }
    }
  })

  # update range of ploidies that can be selected
  observe({
    ploidy_range <- ploidy_range()
    min_ploidy <- ploidy_range[1]
    max_ploidy <- ploidy_range[2]
    current_ploidy <- isolate(reactive_values$ploidy)
    updateNumericInput(session, "ploidy", label = "Ploidy", value = current_ploidy, min = min_ploidy, max = max_ploidy, step = 0.01)
  })

  # update range of cellularities that can be selected
  observe({
    cellularity_range <- cellularity_range()
    min_cellularity <- cellularity_range[1]
    max_cellularity <- cellularity_range[2]
    current_cellularity <- isolate(reactive_values$cellularity)
    updateNumericInput(session, "cellularity", label = "Cellularity", value = current_cellularity, min = min_cellularity, max = max_cellularity, step = 0.01)
  })

  # return a list containing the cached ploidy and cellularity for the given sample
  get_cached_ploidy_and_cellularity <- function(sample) {
    selected_sample <- sample
    sample_cache <- filter(isolate(reactive_values$ploidy_and_cellularity_cache), sample == selected_sample)
    if (nrow(sample_cache) == 1)
      as.list(select(sample_cache, ploidy, cellularity))
    else
      list(ploidy = NA, cellularity = NA)
  }

  # update the current ploidy and cellularity selection
  update_ploidy_and_cellularity <- function(ploidy, cellularity) {
    ploidy_range <- isolate(ploidy_range())
    min_ploidy <- ploidy_range[1]
    max_ploidy <- ploidy_range[2]
    cellularity_range <- isolate(cellularity_range())
    min_cellularity <- cellularity_range[1]
    max_cellularity <- cellularity_range[2]
    if (!is.na(ploidy) && !is.na(cellularity)) {
      if (ploidy < min_ploidy || ploidy > max_ploidy || cellularity <= min_cellularity || cellularity > max_cellularity) {
        ploidy <- NA
        cellularity <- NA
      }
    } else {
      ploidy <- NA
      cellularity <- NA
    }
    reactive_values$ploidy <- ploidy
    reactive_values$cellularity <- cellularity
    updateNumericInput(session, "ploidy", label = "Ploidy", value = ploidy, min = min_ploidy, max = max_ploidy, step = 0.01)
    updateNumericInput(session, "cellularity", label = "Cellularity", value = cellularity, min = min_cellularity, max = max_cellularity, step = 0.01)
  }

  # the copy number data for the selected sample with segment details added
  copy_number_for_selected_sample <- reactive({

    copy_number_data <- reactive_values$copy_number_data
    if (is.null(copy_number_data)) return(NULL)

    selected_sample <- reactive_values$sample
    if (is.null(selected_sample)) return(NULL)

    if (class(copy_number_data) == "QDNAseqCopyNumbers") {
      copy_number <- copy_number_data[,selected_sample]
      copy_number_values <- Biobase::assayDataElement(copy_number, "copynumber")[,1]
      segmented_values <- Biobase::assayDataElement(copy_number, "segmented")[,1]
      copy_number <- Biobase::fData(copy_number) %>%
        rownames_to_column(var = "id") %>%
        as_tibble() %>%
        select(id, chromosome, start, end) %>%
        mutate_at(vars(start, end), as.integer) %>%
        mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
        arrange(chromosome, start) %>%
        mutate(sample = selected_sample) %>%
        mutate(copy_number = copy_number_values) %>%
        mutate(segmented = segmented_values) %>%
        select(sample, chromosome, start, end, copy_number, segmented)
    } else {
      copy_number <- filter(copy_number_data, sample == selected_sample)
    }

    # copy number fitting requires relative copy numbers where values are relative
    # to the average copy number across the genome - using the median segmented
    # copy number
    copy_number <- copy_number %>%
      mutate(copy_number = pmax(copy_number, 0)) %>%
      mutate(segmented = pmax(segmented, 0)) %>%
      mutate_at(vars(copy_number, segmented), ~ . / median(segmented, na.rm = TRUE))

    copy_number %>%
      mutate(position = (start + end) / 2) %>%
      mutate(log2ratio = log2(copy_number))
  })

  # chromosomes for the selected sample
  chromosomes_for_selected_sample <- reactive({
    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(NULL)
    chromosome_offsets(copy_number) %>%
      mutate(start = offset + 1, end = offset + length)
  })

  # segments for the selected sample
  segments_for_selected_sample <- reactive({
    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(NULL)
    copy_number_segments(copy_number) %>%
      mutate(log2ratio = log2(copy_number))
  })

  # relative copy number for whole absolute copy number steps for the
  # current ploidy and cellularity selection
  copy_number_steps <- reactive({
    ploidy <- reactive_values$ploidy
    cellularity <- reactive_values$cellularity
    if (is.na(ploidy) || is.na(cellularity)) return(NULL)
    tibble(absolute_copy_number = 0:input$max_absolute_copy_number_step) %>%
      mutate(relative_copy_number = absolute_to_relative_copy_number(absolute_copy_number, ploidy, cellularity)) %>%
      mutate(log2ratio = log2(relative_copy_number))
  })

  # copy number plot label
  output$genome_copy_number_plot_label <- renderUI({
    label <- reactive_values$sample
    if (is.null(label)) {
      label <- "Genome-wide copy number"
      copy_number_data <- isolate(reactive_values$copy_number_data)
      if (!is.null(copy_number_data)) {
        label <- "Loading copy number data..."
      }
    }
    HTML(label)
  })

  # selected location label
  get_location_label <- reactive({
    label <- NULL
    location <- reactive_values$location
    chromosome <- location$chromosome
    if (!is.null(location$chromosome)) {
      label <- location$chromosome
      if (!is.null(location$start)) {
        label <- str_c(label, ": ", prettyNum(location$start, big.mark = ","))
        if (!is.null(location$end)) {
          label <- str_c(label, " - ", prettyNum(location$end, big.mark = ","))
        }
      }
    }
    label
  })

  # selected location display label
  output$location_label <-renderUI({
    label <- NULL
    sample <- reactive_values$sample
    if (is.null(sample)) {
      copy_number_data <- isolate(reactive_values$copy_number_data)
      if (!is.null(copy_number_data)) {
        label <- "Loading copy number data..."
      }
    } else {
      label <- get_location_label()
      if (!is.null(label)) {
        label <- str_c("Chromosome ", str_replace(label, "-", "&ndash;"))
      }
    }
    if (is.null(label)) label <- "Copy number for selected chromosome"
    HTML(label)
  })

  # range of log2 ratios to display
  log2ratio_range <- reactive({
    if (input$fix_log2ratio_range) {
      input$log2ratio_range
    } else {
      copy_number <- copy_number_for_selected_sample()
      1.1 * log2(quantile(copy_number$copy_number, c(0.001, 0.999), na.rm = TRUE))
    }
  })

  # maximum number of points (bins) to display in copy number plot
  max_points_to_display <- reactive({
    if (input$limit_number_of_points_to_display) {
      input$max_points_to_display
    }
    else {
      Inf
    }
  })

  # genome copy number plot for the selected sample
  create_genome_copy_number_plot <- reactive({

    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(NULL)

    chromosomes <- chromosomes_for_selected_sample()

    copy_number <- copy_number %>%
      convert_to_genomic_coordinates("position", chromosomes)

    copy_number <- copy_number %>%
      select(position, copy_number = log2ratio)

    segments <- segments_for_selected_sample() %>%
      convert_to_genomic_coordinates(c("start", "end"), chromosomes) %>%
      select(start, end, copy_number = log2ratio)

    log2ratio_range <- log2ratio_range()

    copy_number_steps <- NULL
    if (input$show_absolute_copy_number) {
      copy_number_steps <- copy_number_steps()
      if (!is.null(copy_number_steps)) {
        copy_number_steps <- select(copy_number_steps, absolute_copy_number, copy_number = log2ratio)
      }
    }

    # cat("creating genome copy number plot\n")
    genome_copy_number_plot(
      copy_number,
      segments,
      chromosomes,
      copy_number_steps,
      max_points_to_display = max_points_to_display(),
      min_copy_number = log2ratio_range[1], max_copy_number = log2ratio_range[2],
      point_colour = input$bin_colour,
      segment_colour = input$segment_colour,
      copy_number_step_colour = input$absolute_copy_number_step_colour,
      xlabel = "chromosome", ylabel = expression(log[2]~ratio))
  })

  output$genome_copy_number_plot <- renderPlot({
    plot <- create_genome_copy_number_plot()
    if (is.null(plot)) plot <- ggplot()
    plot
  })

  # disable the save button if no genome copy number plot displayed
  observe({
    copy_number <- copy_number_for_selected_sample()
    plot <- create_genome_copy_number_plot()
    toggleState("save_genome_copy_plot", !is.null(copy_number) && !is.null(plot))
  })

  # save genome copy number plot as PDF file
  output$save_genome_copy_plot <- downloadHandler(
    filename = function() { str_c(isolate(reactive_values$sample), '.copy_number.pdf') },
    content = function(file) {
      ggsave(
        file,
        plot = create_genome_copy_number_plot() +
          labs(title = isolate(reactive_values$sample)) +
          theme(
            title = element_text(size = 8),
            axis.title = element_text(size = 9),
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 8)
          ),
        device = "pdf",
        width = input$pdf_width,
        height = input$pdf_height,
        units = "in"
      )
    }
  )

  # copy number plot for the selected sample and chromosome
  create_chromosome_copy_number_plot <- reactive({

    copy_number <- copy_number_for_selected_sample()
    location <- reactive_values$location
    if (is.null(copy_number) || is.null(location$chromosome)) return(NULL)

    copy_number <- copy_number %>%
      select(chromosome, position, copy_number = log2ratio)

    segments <- segments_for_selected_sample() %>%
      select(chromosome, start, end, copy_number = log2ratio)

    chromosomes <- chromosomes_for_selected_sample()

    log2ratio_range <- log2ratio_range()

    copy_number_steps <- NULL
    if (input$show_absolute_copy_number) {
      copy_number_steps <- copy_number_steps()
      if (!is.null(copy_number_steps)) {
        copy_number_steps <- select(copy_number_steps, absolute_copy_number, copy_number = log2ratio)
      }
    }

    position_scale <- 1e-6
    xlabel <- "position (Mbp)"
    xmin <- location$start
    if (is.null(xmin)) xmin <- 1
    xmax <- location$end
    if (is.null(xmax)) {
      xmax <- chromosomes %>%
        filter(chromosome == location$chromosome) %>%
        pull(length)
    }
    if ((xmax - xmin) < 5000000) {
      position_scale <- 1
      xlabel <- "position"
    }

    # cat("creating chromosome copy number plot", location$chromosome, "\n")
    chromosome_copy_number_plot(
      copy_number,
      segments,
      chromosome = location$chromosome,
      start = location$start,
      end = location$end,
      copy_number_steps,
      genes = reactive_values$genes,
      max_points_to_display = max_points_to_display(),
      min_copy_number = log2ratio_range[1], max_copy_number = log2ratio_range[2],
      point_colour = input$bin_colour,
      segment_colour = input$segment_colour,
      copy_number_step_colour = input$absolute_copy_number_step_colour,
      gene_colour = input$gene_colour,
      position_scale = position_scale,
      xlabel = xlabel, ylabel = expression(log[2]~ratio))
  })

  output$chromosome_copy_number_plot <- renderPlot({
    plot <- create_chromosome_copy_number_plot()
    if (is.null(plot)) plot <- ggplot()
    plot
  })

  # disable the save button if no chromosome copy number plot displayed
  observe({
    plot <- create_chromosome_copy_number_plot()
    toggleState("save_chromosome_copy_plot", !is.null(plot))
  })

  # save chromosome copy number plot as PDF file
  output$save_chromosome_copy_plot <- downloadHandler(
    filename = function() { str_c(isolate(reactive_values$sample), '.copy_number.pdf') },
    content = function(file) {
      ggsave(
        file,
        plot = create_chromosome_copy_number_plot() +
          labs(title = str_c(isolate(reactive_values$sample), "  chromosome ", isolate(get_location_label()))) +
          theme(
            title = element_text(size = 8),
            axis.title = element_text(size = 9),
            axis.text.x = element_text(size = 7),
            axis.text.y = element_text(size = 8)
          ),
        device = "pdf",
        width = input$pdf_width,
        height = input$pdf_height,
        units = "in"
      )
    }
  )

  # get the copy number segment or bin corresponding to the given chromosome position
  get_copy_number_at_chromosome_position <- function(chromosome, position) {

    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(NULL)

    segments <- segments_for_selected_sample()

    selected_chromosome <- chromosome

    segment <- filter(segments, chromosome == selected_chromosome, start <= position, end >= position)
    if (nrow(segment) == 0) {
      segment <- filter(copy_number, chromosome == selected_chromosome, start <= position, end >= position)
    }

    if (nrow(segment) != 1) return(NULL)

    select(segment, chromosome, start, end, copy_number, log2ratio)
  }

  # get the copy number segments or bins corresponding to the given chromosome range
  get_copy_number_for_chromosome_range <- function(chromosome, start, end) {

    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(NULL)

    segments <- segments_for_selected_sample()

    selected_chromosome <- chromosome
    selected_start <- start
    selected_end <- end

    segments <- filter(segments, chromosome == selected_chromosome, start <= selected_end, end >= selected_start)
    if (nrow(segments) == 0) {
      segments <- filter(copy_number, chromosome == selected_chromosome, start <= selected_end, end >= selected_start)
    }

    select(segments, chromosome, start, end, copy_number, log2ratio)
  }

  # get the copy number bin corresponding to the given genomic coordinate
  get_copy_number_at_genomic_position <- function(position) {

    chromosomes <- chromosomes_for_selected_sample()
    if (is.null(chromosomes)) return(NULL)

    chromosome <- filter(chromosomes, position >= start, position <= end)
    if (nrow(chromosome) != 1) return(NULL)

    position <- position - chromosome$offset

    get_copy_number_at_chromosome_position(chromosome$chromosome, position)
  }

  # display for the bin or segment hovered over
  output$copy_number_plot_hover_over_label <- renderUI({

    event <- input$genome_copy_number_plot_hover
    event2 <- input$chromosome_copy_number_plot_hover

    copy_number <- NULL
    if (!is.null(event)) {
      copy_number <- get_copy_number_at_genomic_position(event$x)
    } else if (!is.null(event2)) {
      location <- isolate(reactive_values$location)
      if (!is.null(location$chromosome)) {
        copy_number <- get_copy_number_at_chromosome_position(location$chromosome, event2$x)
      }
    }

    if (is.null(copy_number)) return(HTML("&nbsp;"))

    label <- paste0(
      "Chromosome ",
      as.character(copy_number$chromosome),
      " ",
      prettyNum(copy_number$start, big.mark = ","),
      "&mdash;",
      prettyNum(copy_number$end, big.mark = ",")
    )

    relative_copy_number <- copy_number$copy_number

    if (!is.na(relative_copy_number)) {
      label <- paste0(
        label,
        "&nbsp;&nbsp;&nbsp;log2 ratio ",
        strong(round(log2(relative_copy_number), digits = 2)),
        "&nbsp;&nbsp;&nbsp;relative copy number ",
        strong(round(relative_copy_number, digits = 2))
      )

      show_absolute_copy_number <- isolate(input$show_absolute_copy_number)
      if (show_absolute_copy_number) {
        ploidy <- isolate(reactive_values$ploidy)
        cellularity <- isolate(reactive_values$cellularity)
        if (!is.na(ploidy) && !is.na(cellularity)) {
          absolute_copy_number <- relative_to_absolute_copy_number(relative_copy_number, ploidy, cellularity)
          tumour_fraction <- tumour_fraction(absolute_copy_number, cellularity)
          label <- paste0(
            label,
            "&nbsp;&nbsp;&nbsp;absolute copy number ",
            strong(round(absolute_copy_number, digits = 2)),
            "&nbsp;&nbsp;&nbsp;tumour fraction ",
            strong(round(tumour_fraction, digits = 2))
          )
        }
      }
    }

    HTML(label)
  })

  # chromosome selection in genome copy number plot
  observe({
    event <- input$genome_copy_number_plot_click
    if (!is.null(event))
    {
      copy_number <- get_copy_number_at_genomic_position(event$x)
      if (!is.null(copy_number)) {
        reactive_values$location <- list(chromosome = copy_number$chromosome, start = NULL, end = NULL)
      }
    }
  })

  # zoom in on selected region in the chromosome copy number plot
  observe({
    event <- input$chromosome_copy_number_plot_brush
    if (!is.null(event)) {
      chromosome <- isolate(reactive_values$location$chromosome)
      start <- isolate(reactive_values$location$start)
      end <- isolate(reactive_values$location$end)
      if (is.null(start) || is.null(end) || (event$xmin - start) >= 1 || (end - event$xmax) >= 1) {
        reactive_values$location <- list(chromosome = chromosome, start = event$xmin, end = event$xmax)
        session$resetBrush("chromosome_copy_number_plot_brush")
      }
    }
  })

  # reset view to whole of selected chromosome on double click
  observe({
    event <- input$chromosome_copy_number_plot_dblclick
    if (!is.null(event)) {
      chromosome <- isolate(reactive_values$location$chromosome)
      reactive_values$location <- list(chromosome = chromosome, start = NULL, end = NULL)
    }
  })

  # segmented copy number maxima for the selected sample
  segmented_copy_number_maxima_for_selected_sample <- reactive({
    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(NULL)
    relative_copy_number_range <- input$relative_copy_number_range
    copy_number_maxima(copy_number$segmented, min_copy_number = relative_copy_number_range[1], max_copy_number = relative_copy_number_range[2], lower_threshold = 0)
  })

  # segmented copy number maxima to be used for fitting, limited by the
  # number selected with priority given to those with the greatest density
  segmented_copy_number_maxima_for_fitting <- reactive({
    maxima <- segmented_copy_number_maxima_for_selected_sample()
    if (is.null(maxima)) return(NULL)
    top_n(maxima, input$number_of_maxima, density)
  })

  # segmented copy number density plot
  output$segmented_copy_number_density_plot <- renderPlot({

    copy_number <- copy_number_for_selected_sample()
    if (is.null(copy_number)) return(ggplot())

    # cat("creating copy number density plot\n")

    relative_copy_number_range <- input$relative_copy_number_range

    copy_number_steps <- NULL
    if (input$show_absolute_copy_number) {
      copy_number_steps <- copy_number_steps()
      if (!is.null(copy_number_steps)) {
        copy_number_steps <- select(copy_number_steps, absolute_copy_number, copy_number = relative_copy_number)
      }
    }

    plot <- copy_number_density_plot(
      copy_number$segmented,
      copy_number_steps = copy_number_steps,
      min_copy_number = relative_copy_number_range[1],
      max_copy_number = relative_copy_number_range[2],
      copy_number_step_colour = input$absolute_copy_number_step_colour,
      xlabel = "relative copy number")

    if (input$copy_numbers_to_be_fitted == "maxima") {
      maxima_for_fitting <- segmented_copy_number_maxima_for_fitting()
      if (!is.null(maxima_for_fitting)) {
        plot <- plot +
          geom_point(data = maxima_for_fitting, aes(x = copy_number, y = density))
      }
    }

    plot
  })

  # table of segmented copy number maxima
  output$segmented_copy_number_maxima_table <- DT::renderDataTable(
    {
      if (input$copy_numbers_to_be_fitted == "maxima")
        maxima <- segmented_copy_number_maxima_for_fitting()
      else
        maxima <- segmented_copy_number_maxima_for_selected_sample()

      if (is.null(maxima)) {
        maxima <- tibble(relative_copy_number = numeric(), absolute_copy_number = numeric(), density = numeric())
      } else {
        if (input$copy_numbers_to_be_fitted != "maxima") maxima <- top_n(maxima, 10, density)

        maxima <- maxima %>%
          transmute(relative_copy_number = copy_number, absolute_copy_number = "", density)

        ploidy <- reactive_values$ploidy
        cellularity <- reactive_values$cellularity
        if (!is.na(ploidy) && !is.na(cellularity)) {
          maxima <- maxima %>%
            mutate(absolute_copy_number = relative_to_absolute_copy_number(relative_copy_number, ploidy, cellularity))
        }
      }

      datatable(
        maxima,
        colnames = c("Relative copy number", "Absolute copy number", "Density"),
        rownames = FALSE,
        selection = "single",
        options = list(dom = "t")
      ) %>%
        formatRound(columns = 1:3, digits = 2)
    },
    server = FALSE
  )

  # copy numbers to use for fitting
  copy_number_for_fitting <- reactive({

    copy_number_to_be_fitted <- NULL

    selected_copy_number_for_fitting <- input$copy_numbers_to_be_fitted

    if (selected_copy_number_for_fitting == "maxima") {
      maxima <- segmented_copy_number_maxima_for_fitting()
      if (!is.null(maxima)) copy_number_to_be_fitted <- tibble(copy_number = maxima$copy_number, weight = 1)
    } else if (selected_copy_number_for_fitting == "segments") {
      segments <- segments_for_selected_sample()
      if (!is.null(segments)) copy_number_to_be_fitted <- tibble(copy_number = segments$copy_number, weight = segments$weight)
    }

    if (is.null(copy_number_to_be_fitted) || nrow(copy_number_to_be_fitted) < 2) return(NULL)

    copy_number_to_be_fitted
  })

  # distances for absolute copy number fit for grid of ploidies and cellularities
  ploidy_and_cellularity_distances <- reactive({

    distance_function <- input$distance_function

    copy_number_to_be_fitted <- copy_number_for_fitting()
    if (is.null(copy_number_to_be_fitted)) return(NULL)

    ploidy_range <- ploidy_range()
    min_ploidy <- ploidy_range[1]
    max_ploidy <- ploidy_range[2]

    cellularity_range <- cellularity_range()
    min_cellularity <- cellularity_range[1]
    max_cellularity <- cellularity_range[2]

    distances <- find_best_fit_solutions(
      copy_number_to_be_fitted$copy_number, copy_number_to_be_fitted$weight,
      min_ploidy = min_ploidy, max_ploidy = max_ploidy, ploidy_step = 0.01,
      min_cellularity = min_cellularity, max_cellularity = max_cellularity, cellularity_step = 0.01,
      distance_function = distance_function,
      distance_filter_scale_factor = input$distance_filter_scale_factor,
      max_proportion_zero = input$max_proportion_zero,
      min_proportion_close_to_whole_number = input$min_proportion_close_to_whole_number,
      max_distance_from_whole_number = input$max_distance_from_whole_number,
      solution_proximity_threshold = input$solution_proximity_threshold,
      keep_all = TRUE
    )

    # note that the seq function occasionally gives values that are slightly out
    # hence the rounding belo
    distances %>%
      mutate_at(vars(ploidy, cellularity, distance), round, digits = 3)
  })

  # best fit solutions from grid search over ploidies and cellularities
  ploidy_and_cellularity_best_fit_solutions <- reactive({
    distances <- ploidy_and_cellularity_distances()
    if (is.null(distances)) return(NULL)
    distances %>%
      filter(best_fit) %>%
      select(-best_fit) %>%
      arrange(distance, ploidy)
  })

  # heat map representation of the distance function
  output$distance_heat_map_plot <- renderPlot({

    distances <- ploidy_and_cellularity_distances()
    if (is.null(distances)) return(ggplot())

    # cat("creating distance heat map\n")

    plot <- distance_heatmap(
      distances,
      low_distance_colour = input$heatmap_low_distance_colour,
      high_distance_colour = input$heatmap_high_distance_colour
    )

    best_fit_solutions <- filter(distances, best_fit)
    if (nrow(best_fit_solutions) > 0) {
      plot <- plot +
        geom_point(data = best_fit_solutions, aes(x = cellularity, y = ploidy, size = distance), colour = input$heatmap_best_fit_point_colour) +
        scale_size_continuous(limits = c(0.0, max(best_fit_solutions$distance)), range = c(2, 1.25))
    }

    ploidy <- reactive_values$ploidy
    cellularity <- reactive_values$cellularity
    if (!is.na(ploidy) && !is.na(cellularity)) {
      ploidy_range <- ploidy_range()
      min_ploidy <- ploidy_range[1]
      max_ploidy <- ploidy_range[2]
      cellularity_range <- cellularity_range()
      min_cellularity <- cellularity_range[1]
      max_cellularity <- cellularity_range[2]
      if (ploidy >= min_ploidy && ploidy <= max_ploidy && cellularity >= min_cellularity && cellularity <= max_cellularity) {
        plot <- plot +
          geom_point(data = tibble(ploidy = ploidy, cellularity = cellularity), aes(x = cellularity, y = ploidy), size = 2.5, colour = input$heatmap_current_point_colour)
      }
    }

    plot
  })

  # ploidy/cellularity selection in heatmap
  observe({
    event <- input$distance_heat_map_plot_click
    if (is.null(event)) return(NULL)

    copy_number <- isolate(copy_number_for_selected_sample())
    if (is.null(copy_number)) return(NULL)

    ploidy <- round(event$y, digits = 2)
    cellularity <- round(event$x, digits = 2)

    update_ploidy_and_cellularity(ploidy, cellularity)
  })

  # double-click selection in ploidy/cellularity heatmap results
  # in finding the local minimum closest to the clicked point
  observe({
    event <- input$distance_heat_map_plot_dblclick
    if (is.null(event)) return(NULL)

    distance_function <- isolate(input$distance_function)

    copy_number_to_be_fitted <- isolate(copy_number_for_fitting())
    if (is.null(copy_number_to_be_fitted)) return(NULL)

    minimum <- find_minimum(event$y, event$x, copy_number_to_be_fitted$copy_number, copy_number_to_be_fitted$weight, distance_function)

    ploidy <- round(minimum$ploidy, digits = 3)
    cellularity <- round(minimum$cellularity, digits = 3)

    update_ploidy_and_cellularity(ploidy, cellularity)
  })

  # display for the ploidy and cellularity hovered over in the heat map
  output$distance_heat_map_plot_hover_over_label <- renderUI({
    event <- input$distance_heat_map_plot_hover
    distance_function <- isolate(input$distance_function)
    copy_number_to_be_fitted <- isolate(copy_number_for_fitting())
    label <- ""
    if (!is.null(event) && !is.null(copy_number_to_be_fitted)) {
      cellularity <- round(event$x, digits = 2)
      ploidy <- round(event$y, digits = 2)
      distance <- absolute_copy_number_distance(ploidy, cellularity, copy_number_to_be_fitted$copy_number, copy_number_to_be_fitted$weight, distance_function)
      distance <- round(distance, digits = 3)
      label <- paste0(
        "ploidy ", ploidy,
        "&nbsp;&nbsp;cellularity ", cellularity,
        "&nbsp;&nbsp;distance ", distance
      )
    }
    HTML(label)
  })

  # table of segmented copy number maxima
  output$distance_best_fit_solution_table <- DT::renderDataTable(
    {
      solutions <- ploidy_and_cellularity_best_fit_solutions()

      column_names <- c("Ploidy", "Cellularity", "Distance")

      if (is.null(solutions)) {
        solutions <- tibble(ploidy = numeric(), cellularity = numeric(), distance = numeric())
      }

      selected_gene <- reactive_values$selected_gene
      if (!is.null(selected_gene)) {
        gene <- filter(reactive_values$genes, name == selected_gene)
        if (nrow(gene) == 1) {
          copy_number <- get_copy_number_for_chromosome_range(gene$chromosome, gene$start, gene$end)
          relative_copy_number <- unique(copy_number$copy_number)
          if (length(relative_copy_number) == 1) {
            solutions <- solutions %>%
              rowwise() %>%
              mutate(absolute_copy_number = relative_to_absolute_copy_number(relative_copy_number, ploidy, cellularity)) %>%
              mutate(tumour_fraction = tumour_fraction(absolute_copy_number, cellularity)) %>%
              mutate(tumour_fraction = round(tumour_fraction, digits = 2)) %>%
              ungroup() %>%
              select(-absolute_copy_number)
            column_names <- c(column_names, str_c("Tumour fraction (", gene$name, ")"))
          }
        }
      }

      colnames(solutions) <- column_names

      datatable(
        solutions,
        rownames = FALSE,
        selection = "single",
        options = list(
          pageLength = 10,
          dom = ifelse(nrow(solutions) > 10, "tp", "t")
        )
      ) %>%
        formatRound(columns = setdiff(column_names, "Distance"), digits = 2) %>%
        formatRound(columns = "Distance", digits = 3)
    },
    server = FALSE
  )

  # handle selection event in best fit solution table
  observe({
    selected_row <- input$distance_best_fit_solution_table_rows_selected
    if (!is.null(selected_row))
    {
      solutions <- ploidy_and_cellularity_best_fit_solutions()
      if (!is.null(solutions)) {
        selected_solution <- slice(solutions, selected_row)
        if (nrow(selected_solution) == 1) {
          update_ploidy_and_cellularity(selected_solution$ploidy, selected_solution$cellularity)
        }
      }
    }
  })

  # label displaying the current ploidy and cellularity and the resulting distance
  # in scaling to absolute copy numbers
  output$current_ploidy_and_cellularity_label <- renderUI({
    label <- paste0(strong("Current:"))

    ploidy <- reactive_values$ploidy
    cellularity <- reactive_values$cellularity
    if (is.na(ploidy) && is.na(cellularity)) return(HTML(label))

    label <- paste0(label, "&nbsp;&nbsp;ploidy ", ploidy, "&nbsp;&nbsp;cellularity ", cellularity)
    if (is.na(ploidy) || is.na(cellularity)) return(HTML(label))

    distance_function <- input$distance_function

    copy_number_to_be_fitted <- copy_number_for_fitting()
    if (is.null(copy_number_to_be_fitted)) return(HTML(label))

    distance <- absolute_copy_number_distance(ploidy, cellularity, copy_number_to_be_fitted$copy_number, copy_number_to_be_fitted$weight, distance_function)
    distance <- round(distance, digits = 3)
    label <- paste0(label, "&nbsp;&nbsp;distance ", distance)

    HTML(label)
  })

  # label displaying the cached ploidy and cellularity and the resulting distance
  # in scaling to absolute copy numbers
  output$cached_ploidy_and_cellularity_label <- renderUI({
    label <- paste0(strong("Cached:"))

    selected_sample <- reactive_values$sample
    cached_values <- reactive_values$ploidy_and_cellularity_cache
    if (is.null(selected_sample) || nrow(cached_values) == 0) return(HTML(label))

    cached_values_for_selected_sample <- filter(cached_values, sample == selected_sample)
    if (nrow(cached_values_for_selected_sample) != 1) return(HTML(label))

    ploidy <- cached_values_for_selected_sample$ploidy
    cellularity <- cached_values_for_selected_sample$cellularity
    label <- paste0(label, "&nbsp;&nbsp;ploidy ", ploidy, "&nbsp;&nbsp;cellularity ", cellularity)
    if (is.na(ploidy) || is.na(cellularity)) return(HTML(label))

    distance_function <- input$distance_function
    copy_number_to_be_fitted <- copy_number_for_fitting()
    if (is.null(copy_number_to_be_fitted)) return(HTML(label))

    distance <- absolute_copy_number_distance(ploidy, cellularity, copy_number_to_be_fitted$copy_number, copy_number_to_be_fitted$weight, distance_function)
    distance <- round(distance, digits = 3)
    label <- paste0(label, "&nbsp;&nbsp;distance ", distance)

    HTML(label)
  })

  # enable/disable cache button depending on whether there is a valid ploidy
  # and cellularity currently selected
  observe({
    enable_cache <- FALSE
    selected_sample <- reactive_values$sample
    ploidy <- reactive_values$ploidy
    cellularity <- reactive_values$cellularity
    if (!is.null(selected_sample) && !is.na(ploidy) && !is.na(cellularity)) {
      ploidy_range <- ploidy_range()
      min_ploidy <- ploidy_range[1]
      max_ploidy <- ploidy_range[2]
      cellularity_range <- cellularity_range()
      min_cellularity <- cellularity_range[1]
      max_cellularity <- cellularity_range[2]
      if (ploidy >= min_ploidy && ploidy <= max_ploidy && cellularity >= min_cellularity && cellularity <= max_cellularity) {
        enable_cache <- TRUE
      }
    }
    toggleState("cache_current_ploidy_and_cellularity_button", enable_cache)
  })

  # enable/disable restore button depending on whether there is a valid ploidy
  # and cellularity stored in the cache for the current sample
  observe({
    enable_restore <- FALSE
    selected_sample <- reactive_values$sample
    cached_values <- reactive_values$ploidy_and_cellularity_cache
    if (!is.null(selected_sample) && nrow(cached_values) > 0) {
      cached_values_for_selected_sample <- cached_values %>%
        filter(sample == selected_sample)
      if (nrow(cached_values_for_selected_sample) == 1) {
        ploidy <- cached_values_for_selected_sample$ploidy
        cellularity <- cached_values_for_selected_sample$cellularity
        if (!is.na(ploidy) && !is.na(cellularity)) {
          ploidy_range <- ploidy_range()
          min_ploidy <- ploidy_range[1]
          max_ploidy <- ploidy_range[2]
          cellularity_range <- cellularity_range()
          min_cellularity <- cellularity_range[1]
          max_cellularity <- cellularity_range[2]
          if (ploidy >= min_ploidy && ploidy <= max_ploidy && cellularity >= min_cellularity && cellularity <= max_cellularity) {
            enable_restore <- TRUE
          }
        }
      }
    }
    toggleState("restore_cached_ploidy_and_cellularity_button", enable_restore)
  })

  # cache current ploidy and cellularity
  observe({
    if (input$cache_current_ploidy_and_cellularity_button > 0) {
      selected_sample <- isolate(reactive_values$sample)
      ploidy <- isolate(reactive_values$ploidy)
      cellularity <- isolate(reactive_values$cellularity)
      cached_values <- isolate(reactive_values$ploidy_and_cellularity_cache)
      if (!is.null(selected_sample) && !is.na(ploidy) && !is.na(cellularity)) {
        ploidy_range <- isolate(ploidy_range())
        min_ploidy <- ploidy_range[1]
        max_ploidy <- ploidy_range[2]
        cellularity_range <- isolate(cellularity_range())
        min_cellularity <- cellularity_range[1]
        max_cellularity <- cellularity_range[2]
        if (ploidy >= min_ploidy && ploidy <= max_ploidy && cellularity >= min_cellularity && cellularity <= max_cellularity) {
          reactive_values$ploidy_and_cellularity_cache <- cached_values %>%
            filter(sample != selected_sample) %>%
            bind_rows(tibble(sample = selected_sample, ploidy = ploidy, cellularity = cellularity)) %>%
            arrange(sample)
        }
      }
    }
  })

  # restore cached ploidy and cellularity
  observe({
    if (input$restore_cached_ploidy_and_cellularity_button > 0) {
      selected_sample <- isolate(reactive_values$sample)
      cached_ploidy_and_cellularity <- get_cached_ploidy_and_cellularity(selected_sample)
      update_ploidy_and_cellularity(cached_ploidy_and_cellularity$ploidy, cached_ploidy_and_cellularity$cellularity)
    }
  })

  # table of segmented copy number maxima
  output$cached_ploidy_and_cellularity_table <- DT::renderDataTable(
    {
      ploidy_and_cellularity_cache <- reactive_values$ploidy_and_cellularity_cache
      datatable(
        ploidy_and_cellularity_cache,
        rownames = FALSE,
        selection = "single",
        options = list(
          pageLength = 10,
          dom = ifelse(nrow(ploidy_and_cellularity_cache) > 10, "ftip", "t")
        )
      )
    },
    server = FALSE
  )

  # load ploidies and cellularities
  observe({
    file <- input$ploidy_and_cellularity_file
    if (is.null(file)) return(NULL)
    load_ploidies_and_cellularities(file)
  })

  # load sample ploidies and cellularities from a file into the cache
  load_ploidies_and_cellularities <- function(file) {

    if (str_detect(file$name, "\\.csv$")) {
      ploidies_and_cellularities <- read_csv(file$datapath, col_types = cols(sample = "c", ploidy = "d", cellularity = "d"))
    } else {
      ploidies_and_cellularities <- read_tsv(file$datapath, col_types = cols(sample = "c", ploidy = "d", cellularity = "d"))
    }

    expected_columns <- c("sample", "ploidy", "cellularity")
    missing_columns <- setdiff(expected_columns, colnames(ploidies_and_cellularities))
    if (length(missing_columns) > 0) {
      showModal(
        modalDialog(
          title = "Error",
          strong(file$name), "is missing the following columns:", str_c(missing_columns, collapse = ", ")
        )
      )
      return(NULL)
    }

    # filter rows with missing values
    n <- nrow(ploidies_and_cellularities)
    ploidies_and_cellularities <- filter_at(ploidies_and_cellularities, vars(sample, ploidy, cellularity), all_vars(!is.na(.)))
    if (n != nrow(ploidies_and_cellularities)) {
      showModal(
        modalDialog(
          title = "Warning",
          "Some entries have missing values for sample, ploidy and/or cellularity and have been discarded."
        )
      )
    }

    # ensure only one entry for each gene
    n <- nrow(ploidies_and_cellularities)
    ploidies_and_cellularities <- distinct(ploidies_and_cellularities, sample, .keep_all = TRUE)
    if (n != nrow(ploidies_and_cellularities)) {
      showModal(
        modalDialog(
          title = "Warning",
          "Sample names should be unique - entries with duplicated names have been discarded."
        )
      )
    }

    ploidies_and_cellularities <- arrange(ploidies_and_cellularities, sample)

    reactive_values$ploidy_and_cellularity_cache <- ploidies_and_cellularities
    reactive_values$ploidy_and_cellularity_file <- file$name
  }

  # clear the ploidy/cellularity cache
  observe({
    if (input$clear_cache_button > 0) {
      reactive_values$ploidy_and_cellularity_cache = tibble(sample = character(), ploidy = numeric(), cellularity = numeric())
    }
  })

  # disable the save button if nothing in the cache
  observe({
    toggleState("save_cached_ploidies_and_cellularities", nrow(reactive_values$ploidy_and_cellularity_cache) > 0)
  })

  # save cached ploidies and cellularities
  output$save_cached_ploidies_and_cellularities <- downloadHandler(
    filename = function() {
      file <- reactive_values$ploidy_and_cellularity_file
      if (is.null(file)) file <- "ploidies_and_cellularities.csv"
      file
    },
    content = function(file) {
      cache <- isolate(reactive_values$ploidy_and_cellularity_cache)
      if (str_detect(file, "\\.csv$")) {
        write_csv(cache, file)
      }
      else {
        write_tsv(cache, file)
      }
    }
  )

  # handle selection event in ploidy/cellularity cache table by
  # switching to the selected sample if we have copy number data loaded
  # for that sample
  observe({
    selected_row <- input$cached_ploidy_and_cellularity_table_rows_selected
    if (!is.null(selected_row))
    {
      cache <- isolate(reactive_values$ploidy_and_cellularity_cache)
      selected_values <- slice(cache, selected_row)
      if (nrow(selected_values) == 1) {
        samples <- isolate(reactive_values$samples)
        if (!is.null(samples)) {
          selected_sample <- selected_values$sample
          if (selected_sample %in% samples) {
            updateSelectInput(session, "sample", label = "Sample", choices = samples, selected = selected_sample)
            update_ploidy_and_cellularity(selected_values$ploidy, selected_values$cellularity)
          }
        }
      }
    }
  })

  # load genes
  observe({
    file <- input$genes_file
    if (is.null(file)) return(NULL)
    load_genes(file)
  })

  # load genes from a file
  load_genes <- function(file) {

    if (str_detect(file$name, "\\.csv$")) {
      genes <- read_csv(file$datapath, col_types = cols(name = "c", chromosome = "c", start = "i", end = "i"))
    } else {
      genes <- read_tsv(file$datapath, col_types = cols(name = "c", chromosome = "c", start = "i", end = "i"))
    }

    expected_columns <- c("name", "chromosome", "start", "end")
    missing_columns <- setdiff(expected_columns, colnames(genes))
    if (length(missing_columns) > 0) {
      showModal(
        modalDialog(
          title = "Error",
          strong(file$name), "is missing the following columns:", str_c(missing_columns, collapse = ", ")
        )
      )
      return(NULL)
    }

    # filter rows with missing values
    n <- nrow(genes)
    genes <- filter_at(genes, vars(name, chromosome, start, end), all_vars(!is.na(.)))
    if (n != nrow(genes)) {
      showModal(
        modalDialog(
          title = "Warning",
          "Some entries have missing values for name, chromosome, start and/or end and have been discarded."
        )
      )
    }

    # ensure only one entry for each gene
    n <- nrow(genes)
    genes <- distinct(genes, name, .keep_all = TRUE)
    if (n != nrow(genes)) {
      showModal(
        modalDialog(
          title = "Warning",
          "Gene names should be unique - entries with duplicated names have been discarded."
        )
      )
    }

    reactive_values$genes <- genes

    selected_gene <- NULL
    if (nrow(genes) > 0) selected_gene <- genes$name[1]

    updateSelectInput(session, "gene", label = "Gene", choices = genes$name, selected = selected_gene)
  }

  # clear the genes table
  observe({
    if (input$clear_genes_button > 0) {
      reactive_values$genes <- tibble(name = character(), chromosome = character(), start = integer(), end = integer())
      reactive_values$selected_gene <- NULL
      updateSelectInput(session, "gene", label = "Gene", choices = list())
    }
  })

  # genes table
  output$genes_table <- DT::renderDataTable(
    {
      genes <- reactive_values$genes
      datatable(
        genes,
        rownames = FALSE,
        selection = "single",
        options = list(
          pageLength = 10,
          dom = ifelse(nrow(genes) > 10, "ftip", "t")
        )
      ) %>%
        formatStyle("chromosome", textAlign = "right") %>%
        formatRound(c("start", "end"), digits = 0, interval = 3, mark = ",")
    },
    server = FALSE
  )

  # handle selection event in genes table by updating the gene drop-down
  observe({
    selected_row <- input$genes_table_rows_selected
    if (!is.null(selected_row))
    {
      genes <- isolate(reactive_values$genes)
      selected_gene <- slice(genes, selected_row)
      if (nrow(selected_gene) == 1) {
        updateSelectInput(session, "gene", label = "Gene", choices = genes$name, selected = selected_gene$name)
      }
    }
  })

  # initialization with sample copy number data, ploidy and cellularity fit, and genes
  if (getShinyOption("load_sample_data", default = TRUE)) {
    load_ploidies_and_cellularities(list(name = "ploidies_and_cellularities.csv", datapath = "ploidies_and_cellularities.csv"))
    load_genes(list(name = "genes.csv", datapath = "genes.csv"))
    load_copy_number_data(list(name = "copy_number_data.rds", datapath = "copy_number_data.rds"), initialization = TRUE)
  }
}

shinyApp(ui, server)

