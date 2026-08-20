# Generate final QC sample list: data/qc/qc_pass.csv
# This represents the final set of high-quality samples for subclonal reconstruction

# Load libraries
library(data.table)
library(dplyr)
library(readr)
library(stringr)

source("src/utils.R")

# Load datasets
qsamp = data.table::fread('data/qc/sample_level_qc.csv')
qblack = data.table::fread('data/qc/main_tracking_sheet_blacklist.csv', header = T)
qnrpcc = data.table::fread('data/qc/PPCG_Feb2026_NRPCC.txt')
qcmet = data.table::fread('data/qc/ppcg_mutational_signatures_sample_feature_chart.csv')
failedfit = data.table::fread('data/qc/PPCG_samples_no_acceptable_fit_found_Feb2026.txt', header = F)
qctissue = data.table::fread('data/meta/Sample_Donor_Tissue_Origin_2024.csv')

# Cleanup datasets as required
qsamp = qsamp %>% mutate(pass_sample_qc=ifelse(QC_Assessment!="Excluded", T,F)) %>% dplyr::select(PPCG_Sample_ID, pass_sample_qc)
qblack = qblack %>% filter(Country!='') %>% select(PPCG_Sample_ID=WGS_AssayID) %>% mutate(pass_blacklist=F) %>%
    mutate(PPCG_Sample_ID=ifelse(PPCG_Sample_ID=="PPCG0412_DNA", "PPCG0412a_DNA", PPCG_Sample_ID))
qnrpcc = qnrpcc %>% mutate(pass_nrpcc_qc=ifelse(NRPCC>=10, T, F), PPCG_Sample_ID = extract_ppcg_id(Barcode, full = FALSE)) %>% select(PPCG_Sample_ID, pass_nrpcc_qc, nrpcc = NRPCC)
qcmet = qcmet %>% filter(selected_one_sample_per_donor==T) %>% select(PPCG_Sample_ID=WGS_AssayID)
# Merge relevant fields into a new dataset `qc`
datasets = list(qsamp, qblack, qnrpcc %>% dplyr::select(-nrpcc), qcmet)
datasets = lapply(datasets, as.data.frame)
qc = Reduce(function(x, y) merge(x, y, by='PPCG_Sample_ID', all.x=T), datasets)
qc = data.table(qc) %>% mutate(pass_blacklist = ifelse(is.na(pass_blacklist), T, F))

# Replace all NAs with FALSE but dont use funs as it was depreciated
qc = qc %>% mutate_all( ~ ifelse(is.na(.), F, .) )

# Generate final qc selection accounting for all selection criteria
qc = qc %>% mutate(PASS=ifelse(pass_sample_qc & pass_blacklist & pass_nrpcc_qc, T, F))

# Select one primary and met per patient from QC passing samples
## Identify primary and met samples and add patient ID
qctissue_cln = qctissue %>% mutate(Patient_ID=str_extract(PPCG_Sample_ID, 'PPCG[0-9]+')) %>%
    mutate(Tissue_Origin=ifelse(Tissue_Origin=="Recurrence", "Metastasis", Tissue_Origin)) %>% dplyr::select(-Country)
qnrpcc_cln = qnrpcc %>% dplyr::select(PPCG_Sample_ID, nrpcc)
qc = merge(qc, qctissue_cln, by='PPCG_Sample_ID', all.x=T) %>% merge(., qnrpcc_cln, by='PPCG_Sample_ID', all.x=T)
## Make selection by group
qsel = qc %>% dplyr::filter(PASS==T) %>% group_by(Patient_ID, Tissue_Origin) %>% slice_max(nrpcc, n=1) %>% dplyr::pull(PPCG_Sample_ID)
qc = qc %>% dplyr::mutate(ppcg_evo_sel_one=ifelse(PPCG_Sample_ID %in% qsel, T, F)) 

## Count sample by tissue type that passed filter
qc %>% dplyr::filter(ppcg_evo_sel_one==T) %>% dplyr::count(Tissue_Origin) # 58 metastatic samples and 631 primary samples
qc[(pass_sample_qc & pass_blacklist & pass_nrpcc_qc) & ppcg_evo_sel_one==F]

## Cleanup
qc = qc %>% dplyr::select(-PPCG_Donor_ID, -Patient_ID)
qc$PASS = NULL
qc$PASS = qc$ppcg_evo_sel_one

# Write
data.table::fwrite(qc, 'data/qc/qc_pass.csv')


# Check numbers to create a CONSORT diagram as Supplementary Figure

# Numbers after sample QC and blacklist filtering
n_distinct(unique(extract_ppcg_pt(qc$PPCG_Sample_ID))) 
n_distinct(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist])  # 1197 samples
n_distinct(unique(extract_ppcg_pt(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist]))) # 989 patients
table(qc$Tissue_Origin[qc$pass_sample_qc & qc$pass_blacklist]) # 12 BPH, 37 Normal, 147 Metastasis, 1001 Primary

# Numbers after filtering to primaries
n_distinct(extract_ppcg_pt(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist & qc$Tissue_Origin == "Primary"])) # 950 patients
n_distinct(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist & qc$Tissue_Origin == "Primary"]) # 1,001 primary tumour samples

# Numbers after filtering by NRPCC > 10
n_distinct(extract_ppcg_pt(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist & qc$pass_nrpcc_qc & qc$Tissue_Origin == "Primary"])) # 638 patients
n_distinct(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist & qc$pass_nrpcc_qc & qc$Tissue_Origin == "Primary"]) # 658 samples

# Numbers after selecting only highest NRPCC sample per patient
n_distinct(extract_ppcg_pt(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist & qc$pass_nrpcc_qc & qc$Tissue_Origin == "Primary" & qc$ppcg_evo_sel_one])) # 638 patients
n_distinct(qc$PPCG_Sample_ID[qc$pass_sample_qc & qc$pass_blacklist & qc$pass_nrpcc_qc & qc$Tissue_Origin == "Primary" & qc$ppcg_evo_sel_one]) # 638 samples

qc$gs_grade <- get_gs_group(qc$PPCG_Sample_ID)

table(qc$gs_grade[qc$pass_sample_qc & qc$pass_blacklist & qc$pass_nrpcc_qc & qc$Tissue_Origin == "Primary" & qc$ppcg_evo_sel_one]) 


data.table::fwrite(qc %>% dplyr::filter(ppcg_evo_sel_one), 'data/qc/ST1.csv')