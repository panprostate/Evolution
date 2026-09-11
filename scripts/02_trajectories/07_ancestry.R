# Plot prevalence of ordering iii in african ancestry vs european ancestry patients

rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)
library(lemon)
library(janitor)
library(broom)

# FUNCTIONS
source("src/utils.R")
source("src/plot_functions.R")
source("src/plot_theme.R")

plot_trajectory_donut <- function(traj_table, cohort_label = "Cohort") {
  
  # 1. Convert the table object to a standard dataframe
  df <- as.data.frame(traj_table)
  # Standardize column names (table() usually outputs 'Var1' and 'Freq')
  colnames(df) <- c("trajectory", "count")
  
  # 2. Map names and calculate proportions
  df <- df %>%
    mutate(
      traj_name = case_when(
        trajectory == "Ordering 1" ~ "Canonical",
        trajectory == "Ordering 2" ~ "Alternative",
        trajectory == "Ordering 3" ~ "Gain-enriched",
        TRUE ~ as.character(trajectory) # Fallback just in case
      ),
      fraction = count / sum(count),
      # Only show percentage if it's > 2% to avoid text overlapping in tiny slices
      percent_label = ifelse(fraction > 0.02, paste0(round(fraction * 100, 1), "%"), "")
    ) %>%
    # Factor levels: drawing order so Canonical starts at the top
    mutate(traj_name = factor(traj_name, levels = c("Gain-enriched", "Alternative", "Canonical")))
  
  # Get total N for the center text
  total_n <- sum(df$count)
  
  # 3. Define the trajectory colors
  # Loaded from plot_theme.R, but re-assigning names here to ensure they match the biological names used in the figure
  names(trajectory_colors) <- c("Canonical", "Alternative", "Gain-enriched")

  # 4. Build the ggplot
  p <- ggplot(df, aes(x = 2, y = fraction, fill = traj_name)) +
    # The donut ring
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.25) +
    
    # The percentage text centered inside the slices
    geom_text(aes(label = percent_label), 
              position = position_stack(vjust = 0.5), 
              color = "white", size = 2.5, fontface = "bold") +
    
    # Wrap into polar coordinates
    coord_polar(theta = "y", start = 0) +
    
    # Leave the donut hole (x limits 0.5 to 2.5)
    xlim(0.5, 2.5) +
    
    # Center annotation text
    annotate("text", x = 0.5, y = 0, label = paste0(cohort_label, "\nn = ", total_n), 
             size = 3, fontface = "bold", color = "black", lineheight = 0.8) +
    
    # Apply colors
    scale_fill_manual(values = trajectory_colors) +

    # Minimal Nature-ready theme
    theme_void() + 
    theme(legend.position = "none")

  return(p)
}


## PATHS
outdir <- "outputs/02_trajectories"
figdir <- "figures/02_trajectories/ancestry"
clinfile <- "data/meta/PPCG_donors_clin_20241217.csv"
dir.create(outdir, recursive = TRUE); dir.create(figdir, recursive = TRUE) 
trajectories_fps <- list.files(outdir, pattern = "PPCG_Feb2026.*mergedseg_with_clonality.txt", full = T)
heroic_fp <- file.path(outdir, "heroic_assignments.tsv")

# LOAD DATA
tj_samples <- load_trajectories(trajectories_fps)
heroic_tj_samples <- read_delim(heroic_fp, delim = "\t")
clin_data <- read_delim(clinfile, delim = ",")
tj_samples$record_id <- str_remove(tj_samples$sample, "[a-z]_DNA")
tj_samples = merge(tj_samples, clin_data)
tj_samples %>% dplyr::select(sample, predicted_ancestry, trajectory) %>% 
  write_delim(file.path(outdir, "PPCG_trajectories_with_ancestry.txt"), delim = "\t")

# PLOT BY ANCESTRY
ppcg_eur_donuts = plot_trajectory_donut(table(tj_samples$trajectory[tj_samples$predicted_ancestry == "EUR"]), cohort_label = "PPCG European")
ppcg_afr_donuts = plot_trajectory_donut(table(tj_samples$trajectory[tj_samples$predicted_ancestry == "AFR"]), cohort_label = "PPCG African")
heroic_afr_donuts = plot_trajectory_donut(as.table(setNames(heroic_tj_samples$n, heroic_tj_samples$trajectory)), cohort_label = "HEROIC African")

write_tsv(tj_samples %>% dplyr::select(sample, trajectory, predicted_ancestry), file.path(outdir, "Fig3f_source_data.tsv"))
ggsave(file.path(figdir, "Fig3f_pie_trajectory_EUR.pdf"), ppcg_eur_donuts, height = 30/25.4, width = 30/25.4, dpi = 300)
ggsave(file.path(figdir, "Fig3f_pie_trajectory_AFR.pdf"), ppcg_afr_donuts, height = 30/25.4, width = 30/25.4, dpi = 300)
write_tsv(heroic_tj_samples, file.path(outdir, "Fig3g_source_data.tsv"))
ggsave(file.path(figdir, "Fig3g_pie_trajectory_HEROIC_AFR.pdf"), heroic_afr_donuts, height = 30/25.4, width = 30/25.4, dpi = 300)

# HEROIC VALIDATION
df = data.frame(
    n = c(36, 547, 6, 16, 13, 39), 
    cohort = c("PPCG", "PPCG", "PPCG", "PPCG", "HEROIC", "HEROIC"), 
    ancestry = c("EUR", "EUR", "AFR", "AFR", "AFR", "AFR"), 
    trajectory = c("Ordering 3", "Not Ordering 3", "Ordering 3", "Not Ordering 3", "Ordering 3", "Not Ordering 3")
)

# Robustness check 1: controlling for gleason grade group
tb = table(tj_samples$predicted_ancestry, tj_samples$gleason_grade_group)
tb / rowSums(tb)    

tj_samples$is_order3 <- ifelse(tj_samples$trajectory == "Ordering 3", TRUE, FALSE)
model <- glm(
  is_order3 ~ (predicted_ancestry == "AFR") + gleason_grade_group,
  data = tj_samples,
  family = binomial()
)
summary(model) 

# Robustness check 2: only explained by lower T2E
tcga_subtypes <- read_delim(file.path(outdir, "PPCG_trajectories_with_tcga_subtypes_and_you_subtypes.txt"), delim = "\t")
samples_t2e <- tcga_subtypes %>%
  # filter(tcga_subtype == "ERG fusion")
  filter(str_detect(tcga_subtype, "fusion"))

tj_samples$is_t2e <- ifelse(tj_samples$sample %in% samples_t2e$sample, TRUE, FALSE)

tj_samples_t2e <- tj_samples %>%
  filter(sample %in% samples_t2e$sample)

tj_samples_nt2e <- tj_samples %>%
  filter(!sample %in% samples_t2e$sample)

# T2E positive samples only
fisher.test(
    table(tj_samples_t2e$predicted_ancestry == "AFR", tj_samples_t2e$trajectory == "Ordering 3")
)

# T2E negative samples only
fisher.test(
    table(tj_samples_nt2e$predicted_ancestry == "AFR", tj_samples_nt2e$trajectory == "Ordering 3")
)


model <- glm(
  is_order3 ~ (predicted_ancestry == "AFR") + is_t2e,
  data = tj_samples,
  family = binomial()
)
summary(model) 

# Robustness check 3: consider CDK12 mutations
drivers_fp <- "data/raw/Somatic_drivers/DRIVER_ANNOTATIONS__2023-09-25_release4__SNV+indel__filteredVariants__driverGenes.txt"
drivers <- read_delim(drivers_fp)

sample_ids = colnames(drivers)[grep("^PPCG", colnames(drivers))]
cdk12_mutants <- sapply(drivers[drivers$symbol == "CDK12", sample_ids], function(x) x != ".")
tj_samples$is_cdk12 <- tj_samples$sample %in% names(cdk12_mutants[cdk12_mutants])

model <- glm(
  is_order3 ~ (predicted_ancestry == "AFR") + is_cdk12 + is_t2e + gleason_grade_group,
  data = tj_samples,
  family = binomial()
)
summary(model) 
