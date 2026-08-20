# Perform a Random Clinical Trial emulation akin to Thompson et al, Nat Gen 2025

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(survminer)
library(survival)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

# Lookup vector for drug classification
drug_dictionary <- c(
    "Abiraterone" = "ARPI",
    "Enzalutamide" = "ARPI", 
    "Apalutamide" = "ARPI",
    "Darolutamide" = "ARPI",
    "Cabazitaxel" = "Taxane", 
    "Docetaxel" = "Taxane",
    "Paclitaxel" = "Taxane",
    "Carboplatin" = "Platinum",
    "Cisplatin" = "Platinum",
    "Olaparib" = "PARP inhibitor",
    "Goserelin" = "ADT",
    "Leuprorelin" = "ADT", 
    "Radium-223" = "Radiopharmaceutical"
)

plot_km <- function(data, km_fit, title, xlab, ylab, filename){
    ggsurv <- ggsurvplot(km_fit, data = data,
        risk.table = TRUE,
        pval = TRUE,
        conf.int = FALSE,
        palette = c("darkolivegreen4", "grey50"),
        legend.title = "Evolutionary Subtype",
        title = title, xlab = xlab, ylab = ylab, 
        # Legend management: Move it INSIDE the plot to save precious canvas space
        legend = c(0.7, 0.8),
        legend.labs = c("ARPI", "Taxane"),

        # Make the legend background perfectly transparent so it doesn't block the lines
        ggtheme = theme_classic() + theme(
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.margin = margin(0, 0, 0, 0)
        ),

        # --- SHRINK ALL ELEMENTS FOR 60x80mm ---
        size = 0.4,                  # Thinner survival lines (default is too thick for small plots)
        censor.size = 2,             # Smaller censor ticks
        pval.size = 2.5,             # P-value text size (in mm)
        fontsize = 2.5,              # Risk table inner text size (in mm)
        
        # Scale down all thematic fonts (size in points: Nature uses ~6-8pt)
        font.main = 7,
        font.x = 7,
        font.y = 7,
        font.tickslab = 6,
        font.legend = 6,
        
        # --- FIX 3: RISK TABLE MANAGEMENT ---
        risk.table.height = 0.2,    # Use 30% of the vertical space for the table
        risk.table.y.text = FALSE,   # CRITICAL: Replaces y-axis text with color bars to save horizontal space
        risk.table.title = "At risk",
        
        # Clean up the table theme
        tables.theme = theme_cleantable() + 
            theme(
                plot.title = element_text(size = 7), # Table title size
                axis.text.x = element_blank(),       # Remove duplicate x-axis labels from table
                axis.title.y = element_blank(), 
                # This explicitly anchors the left and right margins to align with the plot above it
                plot.margin = margin(0, 5.5, 0, 5.5)
            )
    )
    pdf(file.path(figdir, filename), width = 55/25.3, height = 60/25.3)
    print(ggsurv)
    dev.off()
    # Kaplan-Meier survival analysis
    return(km_fit)
}

# PATHS
figdir <- "figures/05_clinical_utility"
outdir <- "outputs/05_clinical_utility"
hmf_treatment_response_path <- "data/controlled/hmf/hmf_ttf_annotation.rds"
hmf_clinical_metadata_path <- "data/controlled/hmf/metadata.tsv"
imf_activities_path <- "data/controlled/HMF_cluster_maxact_per_biopsy.rds"
hmf_evotypes_path <- "outputs/04_progression/HMF_sample_assignments_G3.tsv"
dir.create(figdir, showWarnings = FALSE)
dir.create(outdir, showWarnings = FALSE)

# LOAD DATA
hmf_treatment_response = readRDS(hmf_treatment_response_path)
hmf_clinical_metadata = read_delim(hmf_clinical_metadata_path)
imf_activities = readRDS(imf_activities_path)
evolutionary_subtypes = read_delim(hmf_evotypes_path)
evolutionary_subtypes$patientIdentifier <- gsub('TI*V*$', '', evolutionary_subtypes$sample_id, perl=TRUE)

######### CLINICAL CORRELATES #########
hmf_clinical_metadata$patientIdentifier <- gsub('TI*V*$', '', hmf_clinical_metadata$sampleId, perl=TRUE)
hmf_clinical_metadata = merge(hmf_clinical_metadata, evolutionary_subtypes, by = "patientIdentifier")
hmf_clinical_metadata$treatment_line <- sapply(hmf_clinical_metadata$patientIdentifier, function(pid){
    nline = hmf_treatment_response %>% dplyr::filter(patientIdentifier == pid & biopsyState == "Post") %>% arrange(treatment_line) %>% dplyr::pull(treatment_line) %>% unique() %>% head(n = 1)
    if (length(nline) == 0) return(NA)
    return(nline - 1)
})
hmf_clinical_metadata$evolutionary_subtype <- ifelse(hmf_clinical_metadata$evolutionary_subtype == "3", "O-III", 
                                              ifelse(hmf_clinical_metadata$evolutionary_subtype == "2", "O-II", "O-I"))
hmf_clinical_metadata$evolutionary_subtype <- factor(hmf_clinical_metadata$evolutionary_subtype, levels = c("O-I", "O-II", "O-III"))

######### RCT EMULATION #########

# Merge evolutionary subtypes with clinical data
hmf_treatment_response = merge(hmf_treatment_response, evolutionary_subtypes, by = "patientIdentifier")
# add cluster information
idx <- match(hmf_treatment_response$sample_id, rownames(imf_activities))
all(rownames(imf_activities)[idx] == hmf_treatment_response$sample_id, na.rm = TRUE)
hmf_treatment_response$imf_cluster <- imf_activities[idx, "max.clust"]
hmf_treatment_response$imf6_activity <- imf_activities[idx, "CIN.Cluster-6"]
hmf_treatment_response$imf8_activity <- imf_activities[idx, "CIN.Cluster-8"]
hmf_treatment_response$imf5_activity <- imf_activities[idx, "CIN.Cluster-5"]

hmf_treatment_response$is_arpi <- ifelse(hmf_treatment_response$name %in% names(drug_dictionary[drug_dictionary == "ARPI"]), 1, 0)
hmf_treatment_response$is_taxane <- ifelse(hmf_treatment_response$name %in% names(drug_dictionary[drug_dictionary == "Taxane"]), 1, 0)
hmf_treatment_response$ttf_event <- ifelse(hmf_treatment_response$reason == "Length of treatment", 0, 1)
hmf_treatment_response$evolutionary_subtype <- ifelse(hmf_treatment_response$evolutionary_subtype == "3", "O-III", 
                                              ifelse(hmf_treatment_response$evolutionary_subtype == "2", "O-II", "O-I"))
hmf_treatment_response$evolutionary_subtype <- factor(hmf_treatment_response$evolutionary_subtype, levels = c("O-I", "O-II", "O-III"))

# Analysis
# TTF in Taxane vs ARPI in 2nd line
# set arm to Taxane for subsequent analysis
taxane_2nd_line <- hmf_treatment_response %>% dplyr::filter(treatment_line == 2 & is_taxane == 1 & is_arpi == 0) %>% dplyr::mutate(arm = "Taxane")
# set arm to ARPI for subsequent analysis
arpi_2nd_line <- hmf_treatment_response %>% dplyr::filter(treatment_line == 2 & is_arpi == 1 & is_taxane == 0) %>% dplyr::mutate(arm = "ARPI")
drugs_2nd_line <- rbind(taxane_2nd_line, arpi_2nd_line)

# Significant interaction between evolutionary subtype O-III and treatment arm, suggesting predictive biomarker
coxph(Surv(TTF, ttf_event) ~ (arm), data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-III",]) %>% summary()
coxph(Surv(TTF, ttf_event) ~ (arm), data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-II",]) %>% summary()
coxph(Surv(TTF, ttf_event) ~ (arm), data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-I",]) %>% summary()
coxph(Surv(TTF, ttf_event) ~ evolutionary_subtype + (arm) + evolutionary_subtype:arm, data = drugs_2nd_line) %>% summary() 


# Cox model predicts better TTF with ARPI even within patients with low IMF6 activity in O-III
coxph(Surv(TTF, ttf_event) ~ (arm), data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-III" & drugs_2nd_line$imf_cluster != "CIN.Cluster-6",])
coxph(Surv(TTF, ttf_event) ~ (arm), data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-II" & drugs_2nd_line$imf_cluster != "CIN.Cluster-6",])
coxph(Surv(TTF, ttf_event) ~ (arm), data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-I" & drugs_2nd_line$imf_cluster != "CIN.Cluster-6",])


write_tsv(
  drugs_2nd_line %>% dplyr::filter(evolutionary_subtype == "O-III") %>% dplyr::select(TTF, ttf_event, arm),
  file.path(outdir, "EXDF10d_source_data.tsv")
)
km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-III",])
plot_km(
    drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-III",], km_fit,
    "",
    "Days from treatment initiation", "Survival rate",
    "EXDF10d_TTF_in_Taxane_vs_ARPI_2nd_line_O-III.pdf"
)
write_tsv(
  drugs_2nd_line %>% dplyr::filter(evolutionary_subtype == "O-II") %>% dplyr::select(TTF, ttf_event, arm),
  file.path(outdir, "EXDF10c_source_data.tsv")
)
km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-II",])
plot_km(
    drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-II",], km_fit,
    "",
    "Days from treatment initiation", "Survival rate",
    "EXDF10c_TTF_in_Taxane_vs_ARPI_2nd_line_O-II.pdf"
)
write_tsv(
  drugs_2nd_line %>% dplyr::filter(evolutionary_subtype == "O-I") %>% dplyr::select(TTF, ttf_event, arm),
  file.path(outdir, "EXDF10b_source_data.tsv")
)
km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-I",])
plot_km(
    drugs_2nd_line[drugs_2nd_line$evolutionary_subtype == "O-I",], km_fit,
    "",
    "Days from treatment initiation", "Survival rate",
    "EXDF10b_TTF_in_Taxane_vs_ARPI_2nd_line_O-I.pdf"
)


write_tsv(
  drugs_2nd_line %>% dplyr::filter(imf_cluster != "CIN.Cluster-6", evolutionary_subtype == "O-III") %>% dplyr::select(TTF, ttf_event, arm),
  file.path(outdir, "EXDF10g_source_data.tsv")
)
km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = drugs_2nd_line[drugs_2nd_line$imf_cluster != "CIN.Cluster-6" & drugs_2nd_line$evolutionary_subtype == "O-III",])
plot_km(
    drugs_2nd_line[drugs_2nd_line$imf_cluster != "CIN.Cluster-6" & drugs_2nd_line$evolutionary_subtype == "O-III",], km_fit,
    "",
    "Days from treatment initiation", "Survival rate",
    "EXDF10g_TTF_in_Taxane_vs_ARPI_2nd_line_O-III_non_CIN.Cluster-6.pdf"
)

write_tsv(
  drugs_2nd_line %>% dplyr::filter(imf_cluster != "CIN.Cluster-6", evolutionary_subtype == "O-II") %>% dplyr::select(TTF, ttf_event, arm),
  file.path(outdir, "EXDF10f_source_data.tsv")
)
km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = drugs_2nd_line[drugs_2nd_line$imf_cluster != "CIN.Cluster-6" & drugs_2nd_line$evolutionary_subtype == "O-II",])
plot_km(
    drugs_2nd_line[drugs_2nd_line$imf_cluster != "CIN.Cluster-6" & drugs_2nd_line$evolutionary_subtype == "O-II",], km_fit,
    "",
    "Days from treatment initiation", "Survival rate",
    "EXDF10f_TTF_in_Taxane_vs_ARPI_2nd_line_O-II_non_CIN.Cluster-6.pdf"
)

write_tsv(
  drugs_2nd_line %>% dplyr::filter(imf_cluster != "CIN.Cluster-6", evolutionary_subtype == "O-I") %>% dplyr::select(TTF, ttf_event, arm),
  file.path(outdir, "EXDF10e_source_data.tsv")
)
km_fit <- survfit(Surv(TTF, ttf_event) ~ arm, data = drugs_2nd_line[drugs_2nd_line$imf_cluster != "CIN.Cluster-6" & drugs_2nd_line$evolutionary_subtype == "O-I",])
plot_km(
    drugs_2nd_line[drugs_2nd_line$imf_cluster != "CIN.Cluster-6" & drugs_2nd_line$evolutionary_subtype == "O-I",], km_fit,
    "",
    "Days from treatment initiation", "Survival rate",
    "EXDF10e_TTF_in_Taxane_vs_ARPI_2nd_line_O-I_non_CIN.Cluster-6.pdf"
)
