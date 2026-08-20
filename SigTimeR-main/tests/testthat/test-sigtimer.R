test_that("simulation completes", {
  
  active_sigs = read.table(
    system.file("extdata", "COSMIC_v3.2_SBS_GRCh38.txt", package="SigTimeR"),
    header=T)[c(2,10,12)]
  active_sigs = as.matrix(active_sigs)
  
  sim = SigTimeR::simulate(
    active_sigs,
    1000, 0.7,
    list(
    c(0.1, 0.3, 0.6), c(0.5, 0.2, 0.3)
  ))
  
  st_sim = SigTimeR::STData(sim$clonal, sim$subclonal, sim$params$sigs)
  
  hdp_result = SigTimeR::time(st_sim, method="hdp")
  expect_true(class(hdp_result)=="STOut")
  expect_true(stmethod(hdp_result)=="hdp")
  
  nnls_result = SigTimeR::time(st_sim, method="nnls")
  expect_true(class(nnls_result)=="STOut")
  expect_true(stmethod(nnls_result)=="nnls")
  
  dm_result = SigTimeR::time(st_sim, method="dm")
  expect_true(class(dm_result)=="STOut")
  expect_true(stmethod(dm_result)=="dm")
  
})


