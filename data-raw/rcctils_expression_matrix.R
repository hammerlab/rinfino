library(rinfino)
library(devtools)

rccutils_expression <- rinfino::load_expdata(system.file("testdata",
                                        "test_expression_matrix.tsv.gz",
                                        package="rinfino"))
devtools::use_data(rcctils_expression, pkg = 'rinfino', compress = 'xz')
