library(rinfino)
library(dplyr)
library(readr)
context("Prep expression matrix")

test_that("prep_expression_matrix can process kallisto output", {
  root_dir <- system.file("testdata", "kallisto", "output", package = "rinfino")
  output_file <- tempfile()
  df1 <- rinfino::prep_expression_matrix(root = root_dir, output_arg = output_file)
  testthat::expect_is(df1, 'data.frame')
  testthat::expect_true(file.exists(output_file))
  df2 <- readr::read_tsv(output_file)
  testthat::expect_equal(as.tbl(df1), as.tbl(df2))
})

