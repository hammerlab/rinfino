library(rinfino)
context("run_pca")

test_that("run_pca works with default options", {
  data("rcctils_expression")
  pca <- rcctils_expression %>%
    rinfino::filter_expdata(fun = function(x) {max(x)>0}) %>%
    rinfino::run_pca()
})

