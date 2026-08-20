# Run dN / dS for different subgroups and pathways

# LIBRARIES
library(stringr)
library(readr)
library(data.table)
library(dndscv)

# FUNCTIONS
source("src/mutation_timer.R") # Borrow code from mutation_timer -> parse_indel_vcf and parse_vcf are useful
source("src/dnds.R")
source("src/utils.R")

# Paths
outdir <- "outputs/01_landscape/dnds"
dir.create(outdir, recursive = TRUE)
mtr_snvdir <- "data/processed/mutation_timer/snvs/"
mtr_indeldir <- "data/processed/mutation_timer/indels/"
cnadir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA"
ccfdir <- "data/processed/PPCG_SNV_clustering"
purdir <- "data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/purity_ploidy"
driverdir <- "data/meta/drivers/"
ncores <- 10

############ Run dN / dS cancer genes all samples ########################
# NB: Only samples that have both indels and snvs

mtr_ids_id = fs::dir_ls(mtr_indeldir, glob = "*vcf") %>% stringr::str_extract("PPCG[^_]+")
mtr_ids_snv = fs::dir_ls(mtr_snvdir, glob = "*vcf") %>% stringr::str_extract("PPCG[^_]+")
mtr_ids = dplyr::intersect(mtr_ids_id, mtr_ids_snv) %>%
    dplyr::intersect(str_remove(qc_pass(), "_DNA"))

ids = mtr_ids[mtr_ids %in% str_remove(get_prim_smps(), "_DNA")]

# STEP 1: Gather all mutation calls
clonal_muts = gather_calls(mtr_snvdir, mtr_indeldir, ids = ids, subclonal = F, ncores = ncores)
write_delim(clonal_muts, file = file.path(outdir, "all_clonal_muts.tsv"), delim = "\t")
subclonal_muts = gather_calls(mtr_snvdir, mtr_indeldir, ids = ids, subclonal = T, ncores = ncores)
write_delim(subclonal_muts, file = file.path(outdir, "all_subclonal_muts.tsv"), delim = "\t")
all_muts = gather_calls(mtr_snvdir, mtr_indeldir, ids = ids, subclonal = NA, ncores = ncores)
write_delim(all_muts, file = file.path(outdir, "all_muts.tsv"), delim = "\t")

# STEP 2: Build sample metadata for primary samples passing QC
prim_smps <- get_prim_smps()
s_df <- data.frame(
    sample_id = prim_smps,
    sample    = str_remove(prim_smps, "_DNA"),
    gs_group  = get_gs_group(prim_smps)
)
s_df <- s_df[!is.na(s_df$gs_group), ]

# STEP 3: Load driver gene lists
driver_lists = list.files(driverdir)
driver_genes = data.table::rbindlist(lapply(driver_lists, function(fp) {
    genes = read_delim(file.path(driverdir, fp), delim = "\t")$gene
    return(data.frame(pathway = str_remove(fp, ".tsv"), gene = genes))
}))

# Genes in dndscv come from RefCDS and might have diff names
data("refcds_hg19", package = "dndscv")
refcds19_genes = sapply(RefCDS, function(x) x$gene_name)
driver_genes = driver_genes[driver_genes$gene %in% refcds19_genes, ]

gs_groups = c("all", "1", "2", "3", "4+")
driver_lists_names = unique(driver_genes$pathway)

base_outdir = file.path(outdir, "dnds_clonal_and_subclonal")
print("dnds for clonal and subclonal mutations...")
for (gs_group in gs_groups) {
    foreach(driver_list = driver_lists_names) %dopar% {
        run_dnds(all_muts, s_df, driver_genes,
                 pathway = driver_list, base_outdir = base_outdir,
                 ordering = "all", gs_group = gs_group, outcome = "all")
    }
}

# For subclonal variants
base_outdir = file.path(outdir, "dnds_subclonal")
print("dnds for subclonal mutations...")
for (gs_group in gs_groups) {
    foreach(driver_list = driver_lists_names) %dopar% {
        run_dnds(subclonal_muts, s_df, driver_genes,
                 pathway = driver_list, base_outdir = base_outdir,
                 ordering = "all", gs_group = gs_group, outcome = "all")
    }
}

# For clonal variants
base_outdir = file.path(outdir, "dnds_clonal")
print("dnds for clonal mutations...")
for (gs_group in gs_groups) {
    foreach(driver_list = driver_lists_names) %dopar% {
        run_dnds(clonal_muts, s_df, driver_genes,
                 pathway = driver_list, base_outdir = base_outdir,
                 ordering = "all", gs_group = gs_group, outcome = "all")
    }
}
