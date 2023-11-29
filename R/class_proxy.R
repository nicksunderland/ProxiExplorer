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
                  gene_chr,
                  gene_start,
                  gene_end,
                  gene_flanks_kb,
                  clump_p1,
                  clump_p2,
                  clump_r2,
                  clump_kb) {

  s <- structure(
    list(
      drug_name      = drug_name,
      gene_name      = gene_name,
      gene_chr       = gene_chr,
      gene_start     = gene_start,
      gene_end       = gene_end,
      gene_flanks_kb = gene_flanks_kb,
      clump_p1       = clump_p1,
      clump_p2       = clump_p2,
      clump_r2       = clump_r2,
      clump_kb       = clump_kb),
    class = "Proxy"
  )

  return(s)
}

# curate a list of defaults
DRUG_PROXIES <- list(
  Proxy(drug_name="incretin class", gene_name="GLP1R", gene_chr="6", gene_start=39016574, gene_end=39055519, gene_flanks_kb=250, clump_p1=5e-8, clump_p2=1, clump_r2=0.001, clump_kb=250)
)

# set the names of the list to the gene drug - name
names(DRUG_PROXIES) <- sapply(DRUG_PROXIES, function(x) paste0(x$gene_name, " - ", x$drug_name))
