library(rinfino)
library(dplyr)
context("Load expdata")

test_that("load_expdata can read data from a matrix", {
  df <- rinfino::load_expdata(system.file("testdata",
                                          "test_expression_matrix.tsv.gz",
                                          package="rinfino"))
  testthat::expect_is(df, class = "tbl_df")
})

