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
                  qtls         = list(),
                  genes_file   = NULL,
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
      qtls           = qtls,
      genes_file     = genes_file),
    class = "Proxy"
  )

  return(s)
}

# curate a list of defaults here
DRUG_PROXIES <- list(
  Proxy(proxy_id      = "GLP1R - Heart failure incidence",
        drug_name     = "incretins",
        gene_name     = "GLP1R",
        exposure_name = "GLP1R/HbA1c (Jurgens, 2022)",
        exposure_file = "GLP1R/hba1c_jurgens_2022_gwas.tsv.gz",
        gene_chr      ="6",
        gene_start    = 39016574,
        gene_end      = 39055519,
        gene_flanks_kb= 250,
        clump_p1      = 5e-8,
        clump_p2      = 1,
        clump_r2      = 0.001,
        clump_kb      = 250,
        pfile         = "GLP1R/glp1r.pgen",
        outcome_name  = "GLP1R/Heart failure incidence (Shah, 2020)",
        outcome_file  = "GLP1R/heartfailure_shah_2020_gwas.tsv.gz",
        qtls          = list("GTEx_v8" = "GLP1R/glp1r_gtex.tsv.gz"),
        genes_file    = "GLP1R/glp1r_gtex_genes.tsv.gz"
        )
)

# set the names of the list to the gene drug - name
names(DRUG_PROXIES) <- sapply(DRUG_PROXIES, function(x) paste0(x$proxy_id))





