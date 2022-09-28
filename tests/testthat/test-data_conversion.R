library(simutils)
data <- scater::mockSCE()
test_that("data class", {
  ## Seurat
  suppressWarnings(data1 <- data_conversion(data, return_format = "Seurat"))
  expect_s4_class(data1, "Seurat")

  ## SingleCellExperiment
  suppressWarnings(data2 <- data_conversion(data, return_format = "SingleCellExperiment"))
  expect_s4_class(data2, "SingleCellExperiment")

  ## list
  suppressWarnings(data3 <- data_conversion(data, return_format = "list"))
  expect_identical(class(data3), "list")
})
