# options(shiny.trace=TRUE)

# silence R CMD checks
BETA_exposure <- BETA_outcome <- BP <- BP_exposure <- CHR <- P <- P_cat <- NULL
P_exposure <- RSID_exposure <- SE_exposure <- SE_outcome <- index <- NULL

#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny ggplot2 ggrepel data.table genepi.utils
#' @importFrom shinyjs disable enable
#' @importFrom stats runif
#' @importFrom forcats fct_rev
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  # reactive values
  exposure_filepath <- reactiveVal(NULL)
  outcome_filepath  <- reactiveVal(NULL)
  pfile             <- reactiveVal(NULL)
  exposure_dat      <- reactiveVal(NULL)
  outcome_dat       <- reactiveVal(NULL)
  harmonised_dat    <- reactiveVal(NULL)

  # observe the import data button
  observeEvent(input$import, {
    c <- input$gene_chr
    s <- input$gene_start
    e <- input$gene_end
    f <- input$gene_flanks_kb * 1000

    # file path checks
    if(!file.exists(exposure_filepath())) return(NULL)
    d <- genepi.utils::standardise_gwas(exposure_filepath(), input_format="default")
    print(c)
    print(s)
    print(e)
    print(f)
    print(d)
    # TODO: column and data checks here
    d <- d[CHR==c & BP>=s-f & BP<=e+f, ]
    # calculate needed things
    d[, P_cat := cut(P, breaks=c(1,5e-2,5e-4,5e-5,5e-6,5e-8,5e-100))]
    # set the data
    exposure_dat(d)

    # TABLE STRUCTURE:
    #
  })

  # observe clump button
  observeEvent(input$clump, {

    # check data
    if(is.null(exposure_dat())) return(NULL)
    if(is.null(outcome_dat())) {
      o <- genepi.utils::standardise_gwas(outcome_filepath(), input_format="default")
      outcome_dat(o)
    }

    p1 <- input$clump_p1
    p2 <- input$clump_p2
    r2 <- input$clump_r2
    kb <- input$clump_kb
    plink2 <- system.file("resources", "plink2", package="ProxiExplorer")
    plink_ref <- pfile()
    clumped_dat <- genepi.utils::clump(exposure_dat(), p1, p2, r2, kb, plink2, plink_ref)
    # harmonise with outcome
    h <- harmonise(clumped_dat, outcome_dat(), gwas1_trait="exposure", gwas2_trait="outcome", merge=c("RSID"="RSID"))
    # set the data
    harmonised_dat(h)

    # TABLE STRUCTURE:
    #
  })



  # the main display output
  output$locus_plot <- renderPlot({

    # check data
    if(is.null(exposure_dat())) return(NULL)

    # create the locus plot
    p <- ggplot(data    = exposure_dat(),
                mapping = aes(x=BP, y=-log10(P))) +
      geom_point(color="lightgray") +
      annotate(geom = "rect", xmin=input$gene_start, xmax=input$gene_end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
      theme_classic() +
      labs(x     = paste0("Chromosome ", input$gene_chr, " position"),
           y     = expression(paste("-log"[10], plain(P))),
           subtitle = input$data_input)

    # if there is clumped data, plot
    if(!is.null(harmonised_dat())) {
      p <- p +
        geom_point(data = harmonised_dat()[!is.na(clump), ], mapping = aes(x=BP_exposure, y=-log10(P_exposure), color=clump, fill=clump), shape=23) +
        geom_vline(data = harmonised_dat()[index==TRUE, ],   mapping = aes(xintercept = BP_exposure), linetype="dotted", color="darkred") +
        geom_point(data = harmonised_dat()[index==TRUE, ],   mapping = aes(x=BP_exposure, y=-log10(P_exposure)), size=3, fill="red", color="red", shape=24) +
        geom_label(data = harmonised_dat()[index==TRUE, ],   mapping = aes(label = clump, x = BP_exposure, y = -0.5)) +
        geom_label_repel(data = harmonised_dat()[index==TRUE, ], mapping = aes(label = RSID_exposure, x = BP_exposure, y = -log10(P_exposure)))
    }

    # plot
    return(p)
  })


  # plot the clump analysis
  output$clumps_plot <- renderPlot({

    # check data
    if(is.null(harmonised_dat())) return(NULL)

    # create the locus plot
    p <- ggplot(data    = harmonised_dat(),
                mapping = aes(x=BETA_exposure, y=BETA_outcome)) +
      geom_point(aes(color=fct_rev(P_cat), fill=fct_rev(P_cat), alpha=fct_rev(P_cat))) +
      geom_point(data = harmonised_dat()[index==TRUE, ], size=3, fill="red", color="red", shape=24) +
      geom_smooth(method="lm", mapping = aes(weight = (1/SE_exposure)*(1/SE_outcome)), color="red") +
      viridis::scale_fill_viridis(name="Exposure P-value", discrete=TRUE) +
      viridis::scale_color_viridis(name="Exposure P-value", discrete=TRUE) +
      scale_alpha_discrete(name="Exposure P-value") +
      facet_wrap(~clump, ncol=6, scales = "free") +
      labs(x = expression('\u03B2'[exposure]),
           y = expression('\u03B2'[outcome]),
           title = paste0("Gene clumps - ", input$data_input),
           subtitle= "Regression line weighted by 1/SE")

    # show
    return(p)


  })




  # observing the data input select box
  observeEvent(input$data_input, {
    if(input$data_input == "Custom input") {

      # set random data, activate file input and wait for input data
      exposure_dat(data.table(CHR = rep("6", 10000), BP = sample(1:50000, 10000, replace = FALSE), P = runif(10000)))
      outcome_dat(data.table(CHR = rep("6", 10000), BP = sample(1:50000, 10000, replace = FALSE), P = runif(10000)))
      shinyjs::enable(id="exposure_file")
      shinyjs::enable(id="outcome_file")

    } else {

      # use example dataset, disable file input
      shinyjs::disable(id="exposure_file")
      shinyjs::disable(id="outcome_file")

      # assess which option was chosen
      if(input$data_input == "GLP1R - HbA1c vs. HF") {

        exposure_filepath(system.file("extdata", "hba1c_jurgens_2022_gwas.tsv.gz", package="ProxiExplorer"))
        outcome_filepath(system.file("extdata", "heartfailure_shah_2020_gwas.tsv.gz", package="ProxiExplorer"))
        pfile(sub(".pgen", "", system.file("extdata", "glp1r.pgen", package="ProxiExplorer")))

      }
      # other options

    }
  })
}
