# Check overlap with Woodcock et al trajectories
rm(list = ls(all = TRUE))

# PACKAGES
library(tidyverse)

# FUNCTIONS
source("src/utils.R")

# PATHS
woodcock_subtypes = read_delim('data/meta/woodcock_subtypes.csv')

wood_class=read_csv('data/meta/ordering_classifications_woodcock2024.csv') %>% dplyr::mutate(
    Local_Patient=str_extract(Sample, "PD\\d+")
) %>% dplyr::select(-Sample)
wood_batch = read_csv('data/meta/ICGC_samples_w_batch_info.csv')
wood_pid = read_csv('data/meta/woodcock_pd_id.txt') %>% dplyr::mutate(
    Local_ID=Sample, Local_Patient=str_extract(Sample, "PD\\d+")
) %>% dplyr::select(-Sample)

# Create ordering classifications for woodcock samples
mtrack = read_csv('data/raw/PPCG_Data_Release/Metadata/PPCG_Master_Tracking_SV_Sheet_15_July_2020.csv') %>%
    dplyr::filter(Country=="UK") %>%
    dplyr::select(PPCG_Sample_ID, Local_ID) %>%
    dplyr::mutate(Local_Patient=str_extract(Local_ID, "PD\\d+")) %>%
    dplyr::mutate(ppevo_qcpass=PPCG_Sample_ID %in% qcp, in_woodcock=Local_Patient %in% wood_class$Local_Patient) %>%
    left_join(ptissue %>% dplyr::select(PPCG_Sample_ID, Tissue_Origin), by=c("PPCG_Sample_ID")) %>%
    left_join(wood_class, by=c("Local_Patient"))
ppevo_wc = mtrack %>% dplyr::filter(Tissue_Origin=="Primary" & ppevo_qcpass==T & in_woodcock==T)
write_csv(ppevo_wc, 'data/meta/240813_woodcock_orderings_ppcg.csv')

# Create mapping between PID and PPCG ID for Dan Woodcock
ppmeta = read_csv('data/raw/PPCG_Data_Release/Metadata/PPCG_Master_Tracking_SV_Sheet_15_July_2020.csv') %>%
    dplyr::filter(Country=="UK") %>%
    dplyr::select(PPCG_Sample_ID, Local_ID)
write_csv(ppmeta, 'data/export/240813_ppcg_ukpid.csv')
