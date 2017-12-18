library(rinfino)
library(devtools)

rcctils_expression <- rinfino::load_expdata(system.file("testdata",
                                        "test_expression_matrix.tsv.gz",
                                        package="rinfino"))

rcctils_sampleinfo <- readr::read_tsv(system.file('data-raw', 'E-MTAB-5640.sdrf.txt', package = 'rinfino')) %>%
  dplyr::mutate(sample_id = `Comment[ENA_RUN]`,
                filename = basename(`Comment[FASTQ_URI]`),
                sample_description = stringr::str_c(`Characteristics[phenotype]`,
                                                    `Characteristics[cell type]`,
                                                    sep = '-'))

rcctils_expression <- rinfino::update_sampleinfo(expdata = rcctils_expression,
                           df = rcctils_sampleinfo,
                           sample_id = 'sample_id')

devtools::use_data(rcctils_expression, compress = 'xz', overwrite = T)
