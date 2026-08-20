# Estimate the correlation between clock-like signatures and age at diagnosis in the PPCG Evolution cohort

# LIBRARIES
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(VariantAnnotation)
library(doParallel)
library(ggrepel)
library(patchwork)
library(MASS)
library(lemon)

cores = 10
registerDoParallel(ncores = cores)

# FUNCTIONS
source("src/real_timing.R")
source("src/plot_functions.R")
source("src/plot_theme.R")
source("src/utils.R")

get_pval_cor <- function(x, y, test = "pearson"){
    crt <- cor.test(x,y, method = test)
    p = crt$p.value
    if (p < 0.001){
        return("p < 0.001")
    } else {
        return(as.character(round(p, 3)))
    }
}

get_estimate_cor <- function(x, y, test = "pearson"){
    crt <- cor.test(x, y, method = test)
    return(round(crt$estimate, 3))
}

plot_facet <- function(df, assignment = "soft", age_var){
    df_long <- df %>% 
        pivot_longer(cols = ends_with(assignment), 
                     names_to = "signature", values_to = "mutations") 
    
    df_long$mutations <- as.numeric(df_long$mutations)

    df_long$signature <- str_replace(df_long$signature, "_", " ") %>% str_replace("sbs", "SBS") %>% 
        str_replace("clock", "Clock")

    cor_data <- df_long %>%
                dplyr::group_by(signature) %>%
                dplyr::summarise(cor = get_estimate_cor(age_at_tumour_collection, mutations, "spearman"), 
                          p_value = get_pval_cor(age_at_tumour_collection, mutations, "spearman")) %>%
                mutate(cor_label = paste0("rho = ", round(cor, 2), "; ", p_value))

    p <- ggplot(df_long, aes(x = age_at_tumour_collection, y = mutations)) + 
        geom_point(alpha = 0.3) + 
        geom_smooth(method = "rlm") +  
        geom_text(data = cor_data, aes(x = Inf, y = Inf, label = cor_label),
            hjust = 1.1, vjust = 1.5, size = 3, inherit.aes = FALSE) + 
        facet_wrap(~signature, nrow = 1, scales = "free_y") + 
        theme(
            strip.text = element_text(size = 10) 
        ) + 
        labs(x = "Age at tumour collection", y = "Number Mutations")
    p <- axes2lemon(p, bt = "both", lt = "both")
    return(p)
}

plot_mut2age <- function(df, muts_var, age_var, col_var, fig_fp, rm_out = F, low_rm = 0.05, up_rm = 0.95){
    if (rm_out){
        bl = df[[muts_var]] < quantile(df[[muts_var]], up_rm, na.rm = T) & 
             df[[muts_var]] > quantile(df[[muts_var]], low_rm, na.rm = T)
        df = df[which(bl),]
    }
    ypos <- max(df[[muts_var]], na.rm = T) * 0.95
    crr <- cor.test(df[[muts_var]], df[[age_var]], method = "pearson")
    p <- ggplot(df, aes_string(x = age_var, y = muts_var)) + 
         geom_point(aes_string(col = col_var)) + geom_smooth(method = "rlm") + 
    annotate(geom = "text", x = 50, y = ypos, 
            label = paste0("p-value=", round(crr$p.value, 3), 
                            "; rho=", round(crr$estimate, 3)))
    ggsave(fig_fp)

    return(p)
}


# genome = load_genome("/camp/svc/reference/Genomics/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa")


## PATHS
outdir = "outputs/01_landscape/real_timing/"
figdir = "figures/01_landscape/real_timing/"
dir.create(outdir, recursive = T)
dir.create(figdir, recursive = T)
meta_fp = "data/meta/PPCG_donors_clin_20241217.csv"
mtr_dir = "data/processed/mutation_timer/snvs"
# sigtimer_dir = "/nemo/project/proj-vanloo-secure/working/nmensah/projects/ppcg/data/legacy/ppcg_signature_timing/210904_ppcg_sbs40-5"
sbs_assign_fp = "data/processed/SNV_signature_annotations__2024-06-21/SNVs_all_signature_20240621.tsv"

sbs_a = read_delim(sbs_assign_fp)
meta = read_delim(meta_fp, ",")
vcf_fs = list.files(mtr_dir, pattern = ".vcf", full = T)
# qs_fs = list.files(sigtimer_dir, pattern = ".*_timing_results.qs", full = T)
prim_smps = get_prim_smps()
prim_smps = prim_smps[!duplicated(extract_ppcg_pt(prim_smps))]

prim_smps = dplyr::intersect(prim_smps, qc_pass())

# Estimate burden of mutations

# ALL mutations
filter_clone = NULL

# SBS1+SBS5/40 hard assignment
target_sigs = c("SBS1", "SBS5", "SBS40a", "SBS40b")
na_df = data.frame(NA, NA, NA, NA)
colnames(na_df) = target_sigs
sigs_hard = data.table::rbindlist(foreach(smp = prim_smps) %dopar% {
    vcf_fp = vcf_fs[grepl(smp, vcf_fs)]
    if (length(vcf_fp) == 0){return(na_df)}
    return(get_sbs_sig(
        vcf_fp, genome = NULL, filter_clone = filter_clone, 
        sbs_a = sbs_a, target_sigs = target_sigs, simplify = F
    ))
})


# SBS1+SBS5/40 soft assignment
sigs_soft  <- data.table::rbindlist(foreach(smp = prim_smps) %dopar% {
    vcf_fp = vcf_fs[grepl(smp, vcf_fs)]
    if (length(vcf_fp) == 0){return(na_df)}
    return(get_sbs_sig_sp(
        vcf_fp, genome = NULL, filter_clone = filter_clone, 
        sbs_a = sbs_a, target_sigs = target_sigs, simplify = F
    ))
})

# Association with age at diagnosis
df = data.frame(
    smp_id = prim_smps, 
    # cpgs = cpgs, 
    sbs1_soft = sigs_soft$SBS1, 
    sbs5_soft = sigs_soft$SBS5, 
    sbs40_soft = sigs_soft$SBS40a + sigs_soft$SBS40b, 
    sbs1_hard = sigs_hard$SBS1, 
    sbs5_hard = sigs_hard$SBS5, 
    sbs40_hard = sigs_hard$SBS40a + sigs_hard$SBS40b
)

df$record_id = extract_ppcg_pt(df$smp_id)
df = left_join(df, meta, by = "record_id")
df$age_at_tumour_collection = as.numeric(df$age_at_tumour_collection)

df$clock_hard <- df$sbs1_hard + df$sbs5_hard + df$sbs40_hard
df$clock_soft <- df$sbs1_soft + df$sbs5_soft + df$sbs40_soft
write_delim(df, file = file.path(outdir, "clock2age.tsv"), delim = "\t")


df = read_delim(file = file.path(outdir, "clock2age.tsv"), delim = "\t")

# remove MMRd patients that might be outliers
mmrd_sigs = c("SBS15", "SBS44", "SBS21")
is_mmrd = extract_ppcg_pt(sbs_a$Sample_ID[rowSums(sbs_a[mmrd_sigs]) > 0])
is_ok = !df$record_id %in% is_mmrd
df <- df[which(is_ok),]

# Plot linear regression for the different potential clock signatures
p1 <- plot_facet(df, "soft", "age_at_tumour_collection")
p2 <- plot_facet(df, "hard", "age_at_tumour_collection")
p <- (p1/p2) + plot_annotation(tag_levels = "a")
write_tsv(
  df %>% dplyr::select(age_at_tumour_collection,
                        sbs1_soft, sbs5_soft, sbs40_soft, clock_soft,
                        sbs1_hard, sbs5_hard, sbs40_hard, clock_hard),
  file.path(outdir, "SF2_source_data.tsv")
)
save_ggplot(p, file.path(figdir, "SF2_clock_sbs2age"), w = 200, h = 100)
