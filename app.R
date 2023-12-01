# Launch the ShinyApp (Do not remove this comment)
# To deploy, run: rsconnect::deployApp()
# Or use the blue button on top of this file

# use these the check consistency of renv
# ?renv::status()
# renv::status()

# debugging
# cat(file=stderr(), ..., "\n")
# rsconnect::deployApp()
# rsconnect::showLogs(streaming = TRUE)
# Navigate to the application (on ShinyApps.io) in your browser, and watch the R console for output.
# to show the errors causing a crash use: rsconnect::showLogs()

pkgload::load_all(export_all = FALSE,helpers = FALSE,attach_testthat = FALSE)
options( "golem.app.prod" = TRUE)
ProxiExplorer::run_app() # add parameters here (if any)




