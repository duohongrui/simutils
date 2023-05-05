library(simutils)
set.seed(111)
data <- matrix(rpois(100*100, 2),
               nrow = 100,
               dimnames = list(paste0('Gene', 1:100),
                               paste0('Cell', 1:100)))
result <- cell_properties(data, verbose = TRUE)

test_that("cell proporties length", {
  expect_length(result$library_size, 100)
  expect_length(result$zero_fraction_cell, 100)
  expect_length(result$TMM_factor, 100)
  expect_length(result$effective_library_size, 100)
  expect_length(result$prop_outliers_cell, 1)
  expect_equal(class(result$cell_cor), c("matrix", "array"))
})

result_gene <- gene_properties(data, verbose = TRUE)
test_that("gene proporties length", {
  expect_length(result_gene$mean_expression, 100)
  expect_length(result_gene$sd, 100)
  expect_length(result_gene$cv, 100)
  expect_length(result_gene$zero_fraction_gene, 100)
  expect_length(result_gene$dispersion, 100)
  expect_length(result_gene$prop_outliers_gene, 1)
})
