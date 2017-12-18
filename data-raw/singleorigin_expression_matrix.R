
library(rinfino)
library(dplyr)
library(devtools)

singleorigin_expression <- rinfino::load_expdata(system.file("data-raw",
                                                         "singleorigin_data.all_genes.txt.gz",
                                                         package="rinfino"),
                                                 gene_col = 'gene_name')

singleorigin_sampleinfo <- readr::read_tsv(system.file('data-raw', 'singleorigin_data.map_sample_to_celltype.txt', package = 'rinfino')) %>%
  dplyr::mutate(sample_name = as.character(sample_name))

singleorigin_expression <- rinfino::update_sampleinfo(expdata = singleorigin_expression,
                                                      df = singleorigin_sampleinfo,
                                                      sample_id = 'sample_name')

devtools::use_data(singleorigin_expression, compress = 'xz', overwrite = T)
devtools::use_data(singleorigin_sampleinfo, compress = 'xz', overwrite = T)







