
#' Run estimate on expression data
#' @param expdf (tbl_df or data.frame) GxS expression data with sampleinfo attr
#' @param tmpdir (chr) optional path to working directory in which to execute estimate
#' @param input_filename (chr) optional full path to input file for estimate
#' @param output_filename (chr) optional full path to output file for estimate (raw, unprocessed)
#' @param trans (function) optional transform to apply to expression data prior to computing ESTIMATE
#' @param filter_fun (function) optional filter to apply to expression data, based on row-wise values
#' @import dplyr estimate readr glue tidyr
#' @importFrom magrittr %>%
#' @export
run_estimate <- function(expdf, tmpdir = tempdir(),
                         input_filename = file.path(tmpdir, 'estimate_input.gct'),
                         output_filename = file.path(tmpdir, 'estimate_output.gct'),
                         trans = NULL,
                         filter_fun = NULL
                         ) {
    # possibly transform expdata
    expdf <- filter_expdata(expdf, trans = trans, fun = filter_fun)
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
    data(SI_geneset, package='estimate') # has to be in globalenv to work, apparentlty
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

