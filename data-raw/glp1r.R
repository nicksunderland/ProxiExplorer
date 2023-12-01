# produce the GLP1R files

library(data.table)

# gene definitions
gene_chrom=6
gene_start_pos=39016574
gene_end_pos=39055519
gene_flanks_kb=5e5


# 1) HbA1c GWAS - Jurgens et al. (PMID: 35177841)
#    - download GWAS
# TODO - find source from paper
hba1c_data <- fread("/Users/xx20081/Documents/local_data/Jurgens_Pirruccello_2022_GWAS_Sumstats/GWAS_sumstats/GWAS_sumstats_EUR__invnorm_glycatedhaemoglobinhba1c__TOTALsample.tsv.gz")
# I similar to this, but I had already standardised the column names
# hba1c <- hba1c_data[BP_b37 >= gene_start_pos-gene_flanks_kb &
#                     BP_b37 <= gene_end_pos+gene_flanks_kb &
#                     CHR==gene_chrom, ]

# 2) GTEx data V8
chr6_gtex <- fread("/Users/xx20081/Documents/local_data/gtex_v8/gtex_v8_chr6.tsv.gz")
gene_win <- chr6_gtex[BP_b37 >= gene_start_pos-gene_flanks_kb &
                      BP_b37 <= gene_end_pos+gene_flanks_kb &
                      CHR==gene_chrom, ]
fwrite(gene_win, "/Users/xx20081/git/ProxiExplorer/inst/extdata/glp1r_gtex.tsv.gz", sep="\t")

# 3) GTEx gene data
chr6_gtex_gene <- fread("/Users/xx20081/Documents/local_data/gtex_v8/gtex_v8_genes.tsv.gz")
genes_win <- chr6_gtex_gene[gene_start >= gene_start_pos-gene_flanks_kb &
                            gene_end <= gene_end_pos+gene_flanks_kb &
                            gene_chr==as.character(gene_chrom), ]
fwrite(genes_win, "/Users/xx20081/git/ProxiExplorer/inst/extdata/glp1r_gtex_genes.tsv.gz", sep="\t")

# 4) The plink reference
# cut the plink reference using ?extract command - need to relook up what I did; but essentially just chop out
# the gene region of interest and then --make-pgen

# curate a list of defaults
DRUG_PROXIES <- list(
  Proxy(proxy_id      = "GLP1R - Heart failure incidence",
        drug_name     = "incretins",
        gene_name     = "GLP1R",
        exposure_name = "GLP1r/HbA1c (Jurgens, 2022)",
        exposure_file = "GLP1r/hba1c_jurgens_2022_gwas.tsv.gz",
        gene_chr      ="6",
        gene_start    = 39016574,
        gene_end      = 39055519,
        gene_flanks_kb= 250,
        clump_p1      = 5e-8,
        clump_p2      = 1,
        clump_r2      = 0.001,
        clump_kb      = 250,
        pfile         = "GLP1r/glp1r.pgen",
        outcome_name  = "GLP1r/Heart failure incidence (Shah, 2020)",
        outcome_file  = "GLP1r/heartfailure_shah_2020_gwas.tsv.gz",
        qtl_name      = "GTEx_v8",
        qtl_file      = "GLP1r/glp1r_gtex.tsv.gz",
        genes_file    = "GLP1r/glp1r_gtex_genes.tsv.gz"
  )
)

# set the names of the list to the gene drug - name
names(DRUG_PROXIES) <- sapply(DRUG_PROXIES, function(x) paste0(x$proxy_id))









