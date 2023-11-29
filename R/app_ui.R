#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyWidgets prettyRadioButtons sliderTextInput
#' @noRd
app_ui <- function(request) {
  # Create the page content
  fluidPage(
    # custom css
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # use shinyjs (for element activation/deactivation things)
    shinyjs::useShinyjs(),
    # A title
    titlePanel("ProxiExplorer"),
    p("Nick Sunderland - nicholas.sunderland@bristol.ac.uk"),
    # Horizontal line
    tags$hr(),
    # sidebar layout
    sidebarLayout(position = "right",
                  sidebarPanel(p(strong("Controls")),
                               width = 3,
                               hr(),
                               fluidRow(
                                 column(12, selectInput(inputId  = "data_input",
                                                        label    = "Data input:",
                                                        choices  = c("GLP1R - HbA1c vs. HF", "Custom input"),
                                                        selected = "GLP1R - HbA1c vs. HF"))
                               ),
                               fluidRow(
                                 column(6, fileInput(inputId = "exposure_file",
                                                     label   = "Exposure:")),
                                 column(6, fileInput(inputId = "outcome_file",
                                                     label   = "Outcome:"))
                               ),
                               fluidRow(
                                 column(6, actionButton(inputId = "import",
                                                        label   = "Import data")),
                                 column(6, actionButton(inputId = "clump",
                                                        label   = "Clump data")),
                               ),
                               hr(),
                               fluidRow(
                                 column(6, p(strong("Gene region")))
                               ),
                               fluidRow(
                                 column(6, selectInput(inputId  = "gene_chr",
                                                       label    = "Chrom:",
                                                       choices  = c(as.character(1:22),"X"),
                                                       selected = "6")),
                                 column(6, numericInput(inputId = "gene_flanks_kb",
                                                        label   = "Flanks (kb):",
                                                        value   = 250,
                                                        min     = 0,
                                                        max     = 1000))
                               ),
                               fluidRow(
                                 column(6, numericInput(inputId = "gene_start",
                                                        label   = "Start:",
                                                        value   = 39016574,
                                                        min     = 1,
                                                        max     = 250000000,
                                                        step    = 1000)),
                                 column(6, numericInput(inputId = "gene_end",
                                                        label   = "End:",
                                                        value   = 39055519,
                                                        min     = 1,
                                                        max     = 250000000,
                                                        step    = 1000))
                               ),
                               hr(),
                               fluidRow(
                                 column(6, p(strong("Clumping")))
                               ),
                               fluidRow(
                                 column(6, sliderTextInput(inputId  = "clump_p1",
                                                           label    = "p1",
                                                           choices  = c(5e-8,5e-6,5e-4,0.05,0.01,0.1,0.5,1.0),
                                                           selected = 5e-8,
                                                           grid     = TRUE)),
                                 column(6, sliderTextInput(inputId  = "clump_p2",
                                                           label    = "p2",
                                                           choices  = c(5e-8,5e-6,5e-4,0.05,0.01,0.1,0.5,1.0),
                                                           selected = 1.0,
                                                           grid     = TRUE))
                               ),
                               fluidRow(
                                 column(6, sliderTextInput(inputId  = "clump_r2",
                                                           label    = "r2",
                                                           choices  = c(0.0001,0.001,0.01,seq(0.1,1,0.1)),
                                                           selected = 0.001,
                                                           grid     = TRUE)),
                                 column(6, sliderTextInput(inputId  = "clump_kb",
                                                           label    = "kb",
                                                           choices  = seq(0,1000,50),
                                                           selected = 250,
                                                           grid     = TRUE))
                               ),
                               hr(),
                               fluidRow(
                                 column(12, p(strong("Mendelian Randomisation")))
                               ),
                               fluidRow(
                                 column(12, prettyRadioButtons(inputId      = "mr_method",
                                                               label        = "Method:",
                                                               choiceNames  = c("IVW", "IVW-LDcorr"),
                                                               choiceValues = c("ivw", "ivw_ldcorr"),
                                                               selected     = "ivw",
                                                               inline       = TRUE,
                                                               animation    = "smooth")
                                 )
                               )
                    ), # sidebar panel end
                    mainPanel(p(strong("Locus plot:")),
                              width = 9,
                              plotOutput(outputId = "locus_plot",
                                         height   = "400px"),
                              p(strong("Clumps:")),
                              plotOutput(outputId = "clumps_plot",
                                         height   = "400px")
                    ) # main panel end
      ) # sidebar layout end
    ) # fluid page end
} # ui function end






#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "ProxiExplorer"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
