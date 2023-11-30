#' @title Drug Proxy Class
#' @description
#' An S3 class holding default settings for a drug proxy instrument.
#' @param drug_name .
#' @param gene_name .
#' @param gene_chr .
#' @param gene_start .
#' @param gene_end .
#' @param gene_flanks_kb .
#' @param clump_p1 .
#' @param clump_p2 .
#' @param clump_r2 .
#' @param clump_kb .
#'
#' @return a Proxy S3 object
#' @export
#'
Proxy <- function(drug_name,
                  gene_name,
                  exposure_name,
                  exposure_file,
                  gene_chr,
                  gene_start,
                  gene_end,
                  gene_flanks_kb,
                  clump_p1,
                  clump_p2,
                  clump_r2,
                  clump_kb,
                  pfile,
                  outcome_name = NULL,
                  outcome_file = NULL,
                  qtl_name     = NULL,
                  qtl_file     = NULL,
                  proxy_id     = NULL) {

  s <- structure(
    list(
      proxy_id       = proxy_id,
      drug_name      = drug_name,
      gene_name      = gene_name,
      exposure_name  = exposure_name,
      exposure_file  = exposure_file,
      gene_chr       = gene_chr,
      gene_start     = gene_start,
      gene_end       = gene_end,
      gene_flanks_kb = gene_flanks_kb,
      clump_p1       = clump_p1,
      clump_p2       = clump_p2,
      clump_r2       = clump_r2,
      clump_kb       = clump_kb,
      pfile          = pfile,
      outcome_name   = outcome_name,
      outcome_file   = outcome_file,
      qtl_name       = qtl_name,
      qtl_file       = qtl_file),
    class = "Proxy"
  )

  return(s)
}

# curate a list of defaults
DRUG_PROXIES <- list(
  Proxy(proxy_id      = "GLP1R - Heart failure incidence",
        drug_name     = "incretins",
        gene_name     = "GLP1R",
        exposure_name = "HbA1c (Jurgens, 2022)",
        exposure_file = "hba1c_jurgens_2022_gwas.tsv.gz",
        gene_chr="6", gene_start=39016574, gene_end=39055519, gene_flanks_kb=250, clump_p1=5e-8, clump_p2=1, clump_r2=0.001, clump_kb=250,
        pfile         = "glp1r.pgen",
        outcome_name  = "Heart failure incidence (Shah, 2020)",
        outcome_file  = "heartfailure_shah_2020_gwas.tsv.gz")
)

# set the names of the list to the gene drug - name
names(DRUG_PROXIES) <- sapply(DRUG_PROXIES, function(x) paste0(x$proxy_id))
