
#' Run estimate on expression data
#' @import dplyr estimate readr glue tidyr
#' @importFrom magrittr %>%
#' @export
run_estimate <- function(expdf, tmpdir = tempdir(),
                         input_filename = file.path(tmpdir, 'estimate_input.gct'),
                         output_filename = file.path(tmpdir, 'estimate_output.gct')
                         ) {
	# get sample names
	sampleinfo <- attr(expdf, 'sampleinfo')
	samplenames <- sampleinfo %>% dplyr::select(`_SAMPLE_ID`) %>% unlist()
    # write out estimate input
	input_file <- file(input_filename, 'w')
	writeLines("#1.2", input_file)
	writeLines(glue::glue("{nrow(expdf)}\t{ncol(expdf)}"), input_file)
	exp_gct <- expdf %>%
	               dplyr::rename(`NAME` = `_GENE`) %>%
				   dplyr::mutate(`Description`=`NAME`) %>%
				   dplyr::select(`NAME`, `Description`, dplyr::one_of(samplenames))
	writeLines(paste(colnames(exp_gct), collapse="\t"), input_file)
	write.table(exp_gct, quote=FALSE, sep="\t", file=input_file, row.names=FALSE, col.names=FALSE)
	close(input_file)
    # run ESTIMATE and save the results into another GCT file
    # for some reason this needs to be in namespace (isn't found otherwise)
    SI_geneset <- estimate::SI_geneset
    estimate::estimateScore(input_filename, output.ds = output_filename, platform = "illumina")
    # read results back into dataframe
    res <- readr::read_tsv(output_filename, comment="#", skip=1) %>%
                    dplyr::select(-`NAME`) %>%
                    as.data.frame
    # reshape result
    res %>% 
        tidyr::gather(`_SAMPLE_ID`, value, -`Description`) %>% 
        tidyr::spread(Description, value)
}

