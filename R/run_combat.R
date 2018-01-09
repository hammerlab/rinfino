#' loads expression data in GxS format, with one row per gene & one column per sample
#' returns a tbl_df also in GxS format, with sample-batch-assignment in attr `sampleinfo`
#'
#' @note Expects data in format prepared by prep_gene_expression_matrix.R
#' @note Genes are filtered to remove those with 0 expression across all samples
#'       these can be filled in with 0 if expressed in at least one batch on combination
#'
#' @param path (chr) path to expression matrix
#' @param gene_col (chr) name of column containing gene symbol
#' @param data_loader (function) function to use when loading data
#' @param batch (chr) optional description to give this batch
#' @param expression_filter (dbl) optional value by which to filter expression (if NULL, no filter)
#' @param ... extra args to data_loader
#' @return tbl_df containing expression data, with columns `_SAMPLE_ID` and `_BATCH`
#' @import dplyr readr stringr purrrlyr
#' @importFrom magrittr %>%
#' @export
load_expdata <- function(path, gene_col='Gene_symbol', data_loader=readr::read_tsv, batch=1,
                         filter=NULL, trans=NULL, ...) {
    df <- data_loader(path, ...) %>%
        dplyr::rename(`_GENE` = !!gene_col)

    # filter and/or transform data on load
    if (!is.null(filter) || !is.null(trans)) {
        df <- filter_expdata(df, fun = filter, trans = trans)
    }

    # annotate sampleinfo
    sampleinfo <- dplyr::tbl_df(list(`_SAMPLE_ID` = df %>% dplyr::select(-`_GENE`) %>% names())) %>%
        dplyr::mutate(sample_id = stringr::str_replace(`_SAMPLE_ID`, pattern = 'X([[:digit:]].*)', replacement='\\1'),
                      batch = !!batch
                      )
    structure(df, sampleinfo = sampleinfo)
}



#' Transform sample_id in both expdata & in attributes
#' @param expdata (tbl_df or data.frame) expression data in GxS format, with sampleinfo attribute
#' @param fun (function) function to apply to sample-id column names, and attributes
#' @return transformed expdata object
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
transform_sampleid <- function(expdata, fun) {
  names(expdata)[-1] <- fun(names(expdata)[-1]) # TODO filter by name, not position
  sampleinfo <- attr(expdata, 'sampleinfo') %>% 
        dplyr::mutate(`_SAMPLE_ID` = fun(`_SAMPLE_ID`))
  attr(expdata, 'sampleinfo') <- sampleinfo
  expdata
}


#' Update sampleinfo affiliated with expdata
#' @param expdata (tbl_df or data.frame) sample expression data in GxS format, with attribute sampleinfo
#' @param df (tbl_df or data.frame) sample-specific info to be merged in with sampleinfo
#' @param sample_id (str) name of column in new sampleinfo containing `_SAMPLE_ID`
#' @import dplyr
#' @importFrom magrittr %>%
#' @export 
update_sampleinfo <- function(expdata, df, sample_id = '_SAMPLE_ID') {
    sampleinfo <- attr(expdata, 'sampleinfo')
	sampleinfo2 <- sampleinfo %>% 
    	dplyr::left_join(df %>% dplyr::mutate_(.dots = list(`_SAMPLE_ID` = sample_id)),
                         suffix = c('.expdata', ''),
                         by = '_SAMPLE_ID')
    attr(expdata, 'sampleinfo') <- sampleinfo2
	expdata
}


#' load several expdata files all at once
#'
#' @note this function is a lightweight wrapper around \code{load_expdata}.
#'       Please see ?merge_expdata for more info on parameters
#'
#' @param path (vector or list) filepaths containing expression matrices
#' @param batch (vector or list) batch name corresponding to each filepath
#' @param ... other params to \code{load_expdata}, each of which can be a vector or singleton
#' @return list of imported dfs returned by load_expdata
#' @import purrr dplyr
#' @importFrom magrittr %>%
#' @export
load_all_expdata <- function(path, batch=NULL, ...) {
    if (is.null(batch))
        batch <- seq_len(length(path))
    dplyr::tbl_df(list(path=path, batch=batch, ...)) %>%
     purrr::pmap(purrr::lift_dl(load_expdata)) %>%
     merge_expdata(.)
}


#' merge expdata dfs together, combining sampleinfo
#' @param dflist (list) list of dataframes imported using load_expdata
#' @param fill (value or NULL) how to fill in NA values, if at all
#' @param filter (function) function to apply to rows, filtering genes/transcripts from result
#' @import purrr dplyr
#' @importFrom magrittr %>%
#' @export
merge_expdata <- function(dflist, fill=NULL, filter=function(x) {max(x) > 0}) {
    df <- dflist %>%
        purrr::reduce(dplyr::full_join, by = '_GENE')
    if ((!is.null(fill) && !is.na(fill)) || !is.null(filter))
        df <- filter_expdata(df, trans = function(x) {ifelse(is.na(x), fill, x)},
                             fun = filter)

    sampleinfo <- dflist %>% purrr::map_df(~ attr(., 'sampleinfo'))
    structure(df, sampleinfo = sampleinfo)
}


#' convert expdaat (as GxS matrix) to tbl_df format
#' @param expmat (matrix) expression data in GxS format
#' @param sampleinfo (tbL_df or data.frame) sampleinfo with `_SAMPLE_ID` column corresponding to column names in expmat
#' @return expdata (tbl_df) containing expression data, with sampleinfo attribute
#' @import tibble
#' @importFrom magrittr %>%
#' @export
expdata_as_df <- function(expmat, sampleinfo) {
    expdata <- as.data.frame(expmat) %>%
        tibble::rownames_to_column('_GENE')
    structure(expdata, sampleinfo = sampleinfo)
}

#' convert expdata (as GxS tbl_df) to matrix format
#' @param df (tbl_df or data.frame) a GxS data frame containing expression data
#' @param trans (function) optional transformation to apply to values in matrix
#' @import dplyr
#' @importFrom magrittr %>%
#' @return matrix of expression data, posssibly transformed values
#' @export
expdata_as_matrix <- function(df, trans=NULL) {
    expmat <- as.matrix(df %>% dplyr::select(-`_GENE`))
    colnames(expmat) <- df %>% dplyr::select(-`_GENE`) %>% names()
    rownames(expmat) <- df %>% dplyr::select(`_GENE`) %>% unlist()
    if (!is.null(trans))
        expmat <- trans(expmat)
    expmat
}

#' helper function to filter gxs expdata more efficiently
#' than can be done using data.frame tools
#'
#' @note Default behavior with both fun & trans NULL is to transform df to matrix & back
#'       This is retained for purposes of performance testing
#'
#' @param df (tbl_df or data.frame) GxS expression data with a column `_GENE`
#' @param fun (function) optional function to be applied to each row (GENE) of expression data
#'               returns a boolean value (TRUE: keep, FALSE: discard)
#' @param trans (function) optional transformation to apply to expression data AFTER filtering, for convenience
#' @import dplyr tibble
#' @importFrom magrittr %>%
#' @return filtered df (tbl_df) containing filtered and/or transformed GxS expression data
filter_gxs_expdata <- function(df, fun=function(x) {max(x)>0}, trans=NULL) {
    expmat <- as.matrix(df %>% dplyr::select(-`_GENE`))
    colnames(expmat) <- df %>% dplyr::select(-`_GENE`) %>% names()
    rownames(expmat) <- df %>% dplyr::select(`_GENE`) %>% unlist()
    sampleinfo <- attr(df, 'sampleinfo')
    if (!is.null(fun)) {
        # todo: improve to support hadley-style quosures (`~ max(.) > 1`)
        keep_exp <- apply(expmat, MARGIN = 1, FUN = fun)
        expmat <- expmat[keep_exp,]
    }
    # more efficient to do transformation here
    if (!is.null(trans))
        expmat <- trans(expmat)
    
    # this is most expensive line in this code
    expdata_as_df(expmat, sampleinfo = sampleinfo)
}


#' helper function to filter gxs expdata more efficiently
#' than can be done using data.frame tools
#'
#' @param df (tbl_df or data.frame) GxS expression data with a column `_GENE`
#' @param fun (function) optional function to be applied to each row (GENE) of expression data
#'               returns a boolean value (TRUE: keep, FALSE: discard)
#' @param trans (function) optional transformation to apply to expression data AFTER filtering, for convenience
#' @import dplyr tibble
#' @importFrom magrittr %>%
#' @export
#' @return filtered df (tbl_df) containing filtered and/or transformed GxS expression data
filter_expdata <- function(df, fun = function(x) {max(x)>0}, trans = NULL) {
    if (is.null(fun) && is.null(trans)) {
        df
    } else {
        filter_gxs_expdata(df, fun = fun, trans = trans)
    }
}

#' helper function to transpose a dataframe of expression data
#' The expectation is that that dataframe will have all-columns-but-one be of the same type
#' and so a matrix-based transpose will be more efficient than one that uses tidyr::gather/spread
#' (preliminary testing on a single expression matrix indicated a 5x improvement)
#'
#' Additionally, while data are in matrix form, one can apply a transformation function
#' more efficiently than while operating on the full dataframe. Thus, `trans` is supported.
#'
#' @param df (tbl_df or data.frame) containing expression data to be transposed
#' @param trans (function) optional function to apply to expression data during transpose
#' @param current_colname (string) name of field to create after transform, holding values of current colnames
#' @param current_rowname (string) name of field existing in df, holding current rownames
#' @import dplyr tibble
#' @importFrom magrittr %>%
#' @return transposed df, of type tbl_df
#' @export
transpose_expdata <- function(df, trans=NULL, current_colname, current_rowname) {
    expmat <- as.matrix(df %>% dplyr::select(-one_of(current_rowname)))
    colnames(expmat) <- df %>% dplyr::select(-one_of(current_rowname)) %>% names()
    rownames(expmat) <- df %>% dplyr::select(one_of(current_rowname)) %>% unlist()
    if (!is.null(trans))
        expmat <- trans(expmat)
    transmat <- t(expmat)
    df <- as.data.frame(transmat) %>%
        tibble::rownames_to_column(current_colname)
}


#' Function to transpose an SxG df of expression data into GxS
#' Optionally applies a numerical transformation (`trans`) to each expression value
#'
#' @param SxG df (tbl_df or data.frame) containing expression data
#' @param trans (function) optional numerical function to apply to expression data during transform
#' @return GxS df (tbl_df or data.frame) containing possibly-transformed expression data
#' @export
transpose2gxs <- function(df, trans=NULL) {
    transpose_expdata(df, trans = trans, current_rowname = '_SAMPLE_ID', current_colname = '_GENE')
}


#' Filter GxS expression data to specific list of genes
#' If `include_random` is true, will optionally select random genes
#' to round out the comparison
#'
#' @param df (tbl_df) SxG df with expression data, where sampleinfo stored in sampleinfo attr
#' @param n_genes (int) desired total number of genes (irrelevant if include_random == FALSE)
#' @param genelist (tbl_df) contains column `_GENE` containing list of genes to filter to
#' @param include_random (bool) whether to sample from remaining genes to reach `n_genes`
#' @import dplyr
#' @importFrom magrittr %>%
#' @return tbl_df containing filtered expression data, with sampleinfo attr
#' @export
filter_genes <- function(df,
                         n_genes = 900,
                         genelist = rinfino::genelist,
                         include_random = TRUE) {
    if (is.null(genelist))
        return(df)
    dff <- df %>% dplyr::semi_join(genelist, by = '_GENE')
    if (include_random) {
        # select random set of genes to get to desired `n_genes`
        n_random <- n_genes - nrow(dff)
        random_genes <- df %>%
            dplyr::select(`_GENE`) %>%
            dplyr::anti_join(genelist, by = '_GENE') %>%
            dplyr::sample_n(n_random)
        genelist <- dplyr::bind_rows(genelist, random_genes)
        dff <- df %>% dplyr::semi_join(genelist, by = '_GENE')
    }
    structure(dff, sampleinfo = attr(df, 'sampleinfo'))
}


#' Run combat to normalize for batch effects
#'
#' @param df (tbl_df) SxG expression data, with batch stored in sampleinfo attribute
#' @param mod (formula) adjustment for any linear covariates in sampleinfo
#' @param trans (formula) optional formula to apply to expression data
#' @import sva
#' @export
run_combat <- function(df, mod = ~ 1, trans=log1p, batch = 'batch') {
    expmat <- expdata_as_matrix(df, trans=trans)
    sampleinfo <- attr(df, 'sampleinfo')
    mm <- model.matrix(mod, data=sampleinfo)
    # re-order to match order of columns in expmat
    sampleinfo <- sampleinfo[match(colnames(expmat), sampleinfo$`_SAMPLE_ID`),]
    # result is a GxS matrix 
    res <- sva::ComBat(as.matrix(expmat), sampleinfo[[batch]], mm)
    expdata_as_df(res, sampleinfo = sampleinfo)
}


#' Test cibersort using one of the expression matrices provided
#' @export
test_cibersort <- function(tmpdir = tempdir(), output_file = file.path(tmpdir, 'cibersort_output.tsv'), test_case = 'ExampleMixtures-GEPs.txt',
                           cibersort_home = '~/cibersort', input_file = file.path(cibersort_home, test_case)) {
    execute_cibersort(tmpdir = tmpdir, cibersort_home = cibersort_home, input_file = input_file,
                      output_file = output_file)
}


#' Execute cibersort using Rserve
#' @param input_file path to input file on disk
#' @param output_file path where output should be stored (on disk)
#' @param tmpdir directory to use for intermediate results
#' @param cibersort_home directory into which CIBERSORT_package.zip has been unzipped
#' @param training_input path to gene signatures for cell types (defaults to LM22)
#' @param rserve_log path in which to store Rserve stderr (contents printed on error)
#' @param error_file path in which to store cibersort stderr (contents printed on error)
#' @import e1071 glue Rserve readr
#' @return cibersort results in a data.frame
#' @note used by \code\{run_cibersort\} to call out to cibersort.
#' @seealso \code\{run_cibersort\} which takes an expdata object rather than a filepath
#' @export
execute_cibersort <- function(input_file, 
                              output_file,
                              absolute_mode = FALSE,
                              tmpdir = tempdir(),
                              cibersort_home = '~/cibersort',
                              training_input = file.path(cibersort_home, 'LM22.txt'),
                              rserve_log = file.path(tmpdir, "Rserve_log.txt"),
                              error_file = file.path(tmpdir, 'cibersort_errors.txt')) {
    cat(glue::glue("Running cibersort in {tmpdir}"))
    cat("\n")
    shell(command = 'R', cmdargs = glue::glue('CMD Rserve --no-save > {rserve_log} 2>&1'))
    cmdargs <- glue::glue("-Xmx10g -Xms3g",
                          " -jar {cibersort_home}/CIBERSORT.jar",
                          " -M {input_file} -B {training_input}",
                          "{ifelse(absolute_mode == TRUE, ' -A', '')}",
                          " > {output_file}",
                          " 2> {error_file}")
    # uses e1071 package
    return_code <- try(shell(command = "java",
                             cmdargs = cmdargs))
    shell(command = 'pkill', cmdargs = "-9 Rserve")
    if (inherits(return_code, 'try-error') || return_code != 0) {
        cat("cibersort stderr:\n")
        try(writeLines(readLines(error_file)))
        cat("rserve log:\n")
        try(writeLines(readLines(rserve_log)))
        stop("Error executing cibersort")
    }
    cibersort <- suppressWarnings(readr::read_delim(output_file, delim = '\t', comment=">", trim_ws = TRUE,
                                                    col_types = cols(
                                                                    .default = col_double(),
                                                                    Column = col_integer(),
                                                                    X27 = col_character()
                                                                    ), col_names = T))
    if (nrow(cibersort) == 0) {
        cat("cibersort stderr:\n")
        writeLines(readLines(error_file))
        cat("rserve log:\n")
        writeLines(readLines(rserve_log))
        stop("No cibersort results.")
    }
    cat("Cibersort completed.\n")
    cibersort <- as.data.frame(cibersort) %>% dplyr::select(-`X27`, -`Column`)
    cibersort
}


#' Transform 'Absolute' mode of cibersort results to relative proportions
#' @param df (tbl_df or data.frame) cibersort results formatted by execute_cibersort
#' @return df containing relative proprtions (with absolute score) for each immune cell type
#' @import dplyr tidyr
#' @importFrom magrittr %>%
#' @export
transform_cibersort_to_relative <- function(df) {
    if (!("Absolute score" %in% names(df)))
		stop("Column 'Absolute score' not found - are you sure this was run using `absolute_mode = TRUE`?")

	cibersort_long <- 
	  df %>% 
	  tidyr::gather(cell_type, proportion, -`_SAMPLE_ID`, -`Absolute score`, -`P-value`, -`Pearson Correlation`)

	# Note that proportions no longer sum to 1; instead they sum to 'Absolute score'
	stopifnot(cibersort_long %>%
				  dplyr::group_by(`_SAMPLE_ID`, `Absolute score`) %>%
				  dplyr::summarise(total = sum(proportion)) %>% 
				  dplyr::ungroup() %>% 
				  dplyr::filter(total != `Absolute score`) %>% nrow()
			  == 0)

	# derive relative proportions for each cell type
	cibersort_long <- cibersort_long %>%
		dplyr::mutate(relative_proportion = proportion / `Absolute score`)

	cibersort <- cibersort_long %>%
		dplyr::select(-proportion) %>%
		tidyr::spread(key = cell_type, value = relative_proportion)
}



#' Run cibersort on expression data
#' @import Rserve glue readr dplyr 
#' @importFrom magrittr %>%
#' @export
run_cibersort <- function(df,
                          output_file = file.path(tmpdir, 'cibersort_output.tsv'),
                          error_file = file.path(tmpdir, 'cibersort_errors.txt'),
                          absolute_mode = FALSE,
                          cibersort_home = '~/cibersort',
                          training_data = NULL,
                          genelist = NULL,
                          trans = NULL, 
                          tmpdir = tempdir(),
                          filter_fun = NULL) {
    # default options if provided data are NULL
    if (!is.null(genelist)) {
        if (absolute_mode == FALSE) {
            genelist <- load_lm22_genes()
        }
    }
    if (is.null(training_data)) {
        training_input <- file.path(cibersort_home, 'LM22.txt')
    } else {
        stop('custom training data not yet supported')
        ## TODO generate signature matrix (call another function)
        ## save signature matrix as training_input
    }

    # prepare input matrix, filtered & possibly transformed
    # write matrix to file, in order to call cibersort
    input_file = file.path(tmpdir, 'cibersort_input.tsv')
    expmat <- df %>%
        filter_genes(genelist = genelist, include_random = FALSE) %>%
        filter_expdata(fun = filter_fun, trans = trans) %>%
        dplyr::rename(Gene_symbols = `_GENE`) %>%
        readr::write_tsv(path = input_file)
    # execute cibersort
    cibersort <- execute_cibersort(tmpdir=tmpdir, cibersort_home=cibersort_home, 
                      input_file=input_file, training_input=training_input,
                      output_file=output_file, error_file=error_file, absolute_mode=absolute_mode)
    # read & format results 
    cibersort <- cibersort %>% 
        dplyr::mutate(`_SAMPLE_ID` = colnames(expmat)[-1]) %>%
        dplyr::select(`_SAMPLE_ID`, one_of(names(cibersort)))
    if (absolute_mode == TRUE) {
        # reformat cibersort values to be proportions for consistency
        cibersort <- transform_cibersort_to_relative(cibersort)
    }
    cibersort
}

#' Helper function to load lm22 genes from cibersort_home
#' and prepare in format for `filter_genes`
#' @param cibersort_home (str) path to cibersort home dir
#' @importFrom magrittr %>%
#' @import dplyr readr
#' @export
load_lm22_genes <- function(cibersort_home = '~/cibersort') {
        lm22_path <- file.path(cibersort_home, 'LM22.txt')
        lm22_genes <- readr::read_tsv(lm22_path) %>% 
            dplyr::mutate(`_GENE` = `Gene symbol`)
}



#' Function to transpose an GxS df of expression data into SxG
#' Optionally applies a numerical transformation (`trans`) to each expression value
#'
#' @param GxS df (tbl_df or data.frame) containing expression data
#' @param trans (function) optional numerical function to apply to expression data during transform
#' @return SxG df (tbl_df or data.frame) containing possibly-transformed expression data
#' @export
transpose2sxg <- function(df, trans=NULL) {
    transpose_expdata(df, trans = trans, current_rowname = '_GENE', current_colname = '_SAMPLE_ID')
}

#' run PCA analysis on expression data
#'
#' @param df (tbl_df) GxS expression data, with columns indicating sample_ids
#' @param plot (boolean) if TRUE, also plot PCA results
#' @param trans (function) if provided, transform expression data prior to PCA
#' @return pca object for each sample run (columns in the provided df)
#' @import dplyr ggplot2 txtplot
#' @importFrom magrittr %>%
#' @export
run_pca <- function(df, use_ggplot=TRUE, trans=log1p, group = 'batch', colour = 'batch', ...) {
    # transform df to SxG orientation, possibly applying `trans` function
    tdf <- transpose2sxg(df, trans=trans)

    # merge in sampleinfo (sample_id & batch)
    sampleinfo <- attr(df, 'sampleinfo')

    # run pca
    pca <- prcomp(tdf %>% dplyr::select(-`_SAMPLE_ID`),
                  center = TRUE,
                  scale. = TRUE)

    if (use_ggplot) {
        pcadf <- tbl_df(pca$x) %>%
            bind_cols(tdf %>% dplyr::select(`_SAMPLE_ID`)) %>%
            dplyr::left_join(sampleinfo, by = '_SAMPLE_ID')
        pcaplot <- ggplot2::ggplot(pcadf, ggplot2::aes_string(x = 'PC1', y = 'PC2', group = group, colour = colour, ...)) +
            ggplot2::geom_point()
        print(pcaplot)
    } else {
        # produce ascii plot
        txtplot::txtplot(pca$x[,'PC1'], pca$x[,'PC2'], xlab = 'PC1', ylab = 'PC2')
    }
    pca
}

