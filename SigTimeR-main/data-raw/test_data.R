## code to prepare test datasets:
## - test_clonal
## - test_subclonal
## - test_sigs

# Load signatures
active_sigs = read.table(
system.file("extdata", "COSMIC_v3.2_SBS_GRCh38.txt", package="SigTimeR"),
header=T)[c(2,10,12)]
active_sigs = as.matrix(active_sigs)
  
# Generate test data
sim = SigTimeR::simulate(
    active_sigs,
    1000, 0.7,
    list(
    c(0.1, 0.3, 0.6), c(0.5, 0.2, 0.3)
))

# Set variables
test_clonal = sim$clonal
test_subclonal = sim$subclonal
test_sigs = sim$params$sigs

# Add test data to package
usethis::use_data(test_clonal, overwrite = TRUE)
usethis::use_data(test_subclonal, overwrite = TRUE)
usethis::use_data(test_sigs, overwrite = TRUE)
