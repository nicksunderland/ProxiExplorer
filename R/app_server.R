# options(shiny.trace=TRUE)
# options(repos=c(BiocManager::repositories(version = "3.18")))

# silence R CMD checks
BETA_exposure <- BETA_outcome <- BP <- BP_exposure <- CHR <- P <- P_cat <- NULL
P_exposure <- RSID_exposure <- SE_exposure <- SE_outcome <- index <- RSID <- NULL

#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny shinyWidgets ggplot2 ggrepel data.table genepi.utils
#' @importFrom shinyjs disable enable
#' @importFrom stats runif
#' @importFrom forcats fct_rev
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  # reactive values
  # global DRUG_PROXIES object, with default genes and settings, created in class_proxy.R
  exposure_filepath <- reactiveVal(NULL)
  outcome_filepath  <- reactiveVal(NULL)
  pfile             <- reactiveVal(NULL)
  exposure_dat      <- reactiveVal(NULL)
  outcome_dat       <- reactiveVal(NULL)
  clumped_dat       <- reactiveVal(NULL)
  harmonised_dat    <- reactiveVal(NULL)
  qtl_dat           <- reactiveVal(NULL)

  # observe the import data button
  observeEvent(input$import, {
    c <- input$gene_chr
    s <- input$gene_start
    e <- input$gene_end
    f <- input$gene_flanks_kb * 1000

    # file path checks
    if(!file.exists(exposure_filepath())) return(NULL)

    #########
    cat(file=stderr(), exposure_filepath(), "\n")
    #########

    d <- genepi.utils::standardise_gwas(exposure_filepath(), input_format="default")

    #########
    cat(file=stderr(), utils::str(d), "\n")
    #########

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

    # clump parameters from the GUI
    p1 <- input$clump_p1
    p2 <- input$clump_p2
    r2 <- input$clump_r2
    kb <- input$clump_kb

    # the plink2 executable
    if(get_os() == "osx") {
      plink2 <- system.file("resources", "plink2", package="ProxiExplorer")
    } else if(get_os() == "linux") {
      plink2 <- system.file("resources", "plink2_linux", package="ProxiExplorer")
      #########
      cat(file=stderr(), plink2, "\n")
      #########
    } else {
      message("Operating system: ", get_os())
      stop("operating system error")
    }


    # give the plink2 executable permission to run
    Sys.chmod(plink2, "777", use_umask = FALSE)
    # get the reference file
    plink_ref <- pfile()
    # run clumping
    clumped_dat <- genepi.utils::clump(exposure_dat(), p1, p2, r2, kb, plink2, plink_ref)
    # set clumped data
    clumped_dat(clumped_dat)

    # harmonise with outcome
    h <- harmonise(clumped_dat, outcome_dat(), gwas1_trait="exposure", gwas2_trait="outcome", merge=c("RSID"="RSID"))
    # set the data
    harmonised_dat(h)

    # TABLE STRUCTURE:
    #
  })

  # observe the QTL combo box
  observeEvent(input$qtl_source, {

    # check if requested
    if(input$qtl_source == "None") return(NULL)

    # get QTL data
    # if(input$qtl_source == )
    # qtl_source <- input$qtl_source
    # qtl_tissue <- input$qtl_tissue
    # s <- input$gene_start
    # e <- input$gene_end
    # f <- input$gene_flanks_kb * 1000
    #
    # # file path checks
    # if(!file.exists(exposure_filepath())) return(NULL)
    #
    #

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
           y     = expression(paste("-log[10]", plain(P))),
           subtitle = input$data_input)

    # if there is clumped data, plot
    if(!is.null(clumped_dat())) {
      p <- p +
        geom_point(data = clumped_dat()[!is.na(clump), ], mapping = aes(x=BP, y=-log10(P), color=clump, fill=clump), shape=23) +
        geom_vline(data = clumped_dat()[index==TRUE, ],   mapping = aes(xintercept = BP), linetype="dotted", color="darkred") +
        geom_point(data = clumped_dat()[index==TRUE, ],   mapping = aes(x=BP, y=-log10(P)), size=3, fill="red", color="red", shape=24) +
        geom_label(data = clumped_dat()[index==TRUE, ],   mapping = aes(label = clump, x = BP, y = -0.5)) +
        geom_label_repel(data = clumped_dat()[index==TRUE, ], mapping = aes(label = RSID, x = BP, y = -log10(P)))
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
           title = paste0("Variant clumps - ", input$data_input),
           subtitle= "Regression line weighted by 1/SE")

    # show
    return(p)
  })

  # plot the QTL data
  output$qtl_plot <- renderPlot({

    # check data
    if(is.null(qtl_dat())) return(NULL)

    # create the locus plot
    p <- ggplot(data    = qtl_dat(),
                mapping = aes(x=BP, y=-log10(P))) +
      geom_point(color="lightgray") +
      annotate(geom = "rect", xmin=input$gene_start, xmax=input$gene_end, ymin=0, ymax=Inf, fill="blue", alpha = 0.05) +
      theme_classic() +
      labs(x     = paste0("Chromosome ", input$gene_chr, " position"),
           y     = expression(paste("-log[10]", plain(P))),
           subtitle = input$data_input)

    # if there is clumped data, plot
    if(!is.null(clumped_dat())) {
      p <- p +
        geom_point(data = clumped_dat()[!is.na(clump), ], mapping = aes(x=BP, y=-log10(P), color=clump, fill=clump), shape=23) +
        geom_vline(data = clumped_dat()[index==TRUE, ],   mapping = aes(xintercept = BP), linetype="dotted", color="darkred") +
        geom_point(data = clumped_dat()[index==TRUE, ],   mapping = aes(x=BP, y=-log10(P)), size=3, fill="red", color="red", shape=24) +
        geom_label(data = clumped_dat()[index==TRUE, ],   mapping = aes(label = clump, x = BP, y = -0.5)) +
        geom_label_repel(data = clumped_dat()[index==TRUE, ], mapping = aes(label = RSID, x = BP, y = -log10(P)))
    }


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

      # get & set the defaults
      exposure_filepath(     system.file("extdata", DRUG_PROXIES[[input$data_input]]$exposure_file, package="ProxiExplorer"))
      outcome_filepath(      system.file("extdata", DRUG_PROXIES[[input$data_input]]$outcome_file,  package="ProxiExplorer"))
      pfile(sub(".pgen", "", system.file("extdata", DRUG_PROXIES[[input$data_input]]$pfile,         package="ProxiExplorer")))
      updateSelectInput(session,     inputId = "gene_chr",       selected = DRUG_PROXIES[[input$data_input]]$gene_chr)
      updateNumericInput(session,    inputId = "gene_flanks_kb", value    = DRUG_PROXIES[[input$data_input]]$gene_flanks_kb)
      updateNumericInput(session,    inputId = "gene_start",     value    = DRUG_PROXIES[[input$data_input]]$gene_start)
      updateNumericInput(session,    inputId = "gene_end",       value    = DRUG_PROXIES[[input$data_input]]$gene_end)
      updateSliderTextInput(session, inputId = "clump_p1",       selected = DRUG_PROXIES[[input$data_input]]$clump_p1)
      updateSliderTextInput(session, inputId = "clump_p2",       selected = DRUG_PROXIES[[input$data_input]]$clump_p2)
      updateSliderTextInput(session, inputId = "clump_r2",       selected = DRUG_PROXIES[[input$data_input]]$clump_r2)
      updateSliderTextInput(session, inputId = "clump_kb",       selected = DRUG_PROXIES[[input$data_input]]$clump_kb)

    }
  })
}
