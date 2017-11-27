library(rinfino)
library(readr)
library(dplyr)
library(stringr)

# load & save data
rcctils_sampleinfo <- readr::read_tsv(system.file('data-raw', 'E-MTAB-5640.sdrf.txt', package = 'rinfino')) %>%
  dplyr::mutate(sample_id = `Comment[ENA_RUN]`,
                filename = basename(`Comment[FASTQ_URI]`),
                sample_description = stringr::str_c(`Characteristics[phenotype]`,
                                                    `Characteristics[cell type]`,
                                                    sep = '-'))
devtools::use_data(rcctils_sampleinfo, compress = 'xz')
