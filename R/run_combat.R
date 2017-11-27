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
#' @export
load_expdata <- function(path, gene_col='Gene_symbol', data_loader=readr::read_tsv, batch=1,
                         filter=NULL, trans=NULL, ...) {
    df <- data_loader(path, ...) %>%
        dplyr::rename(`_GENE` = !!gene_col)

    # filter and/or transform data on load
    if (!is.null(filter) || !is.null(trans)) {
        df <- filter_gxs_expdata(df, fun = filter, trans = trans)
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
#' @param sampleinfo (tbl_df or data.frame) sample-specific info to be merged in with sampleinfo
#' @param sample_id (str) name of column in new sampleinfo containing `_SAMPLE_ID`
#' @import dplyr
#' @export 
update_sampleinfo <- function(expdata, sampleinfo, sample_id = '_SAMPLE_ID') {
	sampleinfo2 <- sampleinfo %>% 
    	dplyr::left_join(df %>% dplyr::mutate(`_SAMPLE_ID` = !!sample_id),
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


#' convert expdata (as GxS tbl_df) to matrix format
#' @param df (tbl_df or data.frame) a GxS data frame containing expression data
#' @param trans (function) optional transformation to apply to values in matrix
#' @import dplyr
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
    df <- as.data.frame(expmat) %>%
        tibble::rownames_to_column('_GENE')
    structure(df, sampleinfo = sampleinfo)
}
filter_expdata <- filter_gxs_expdata


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
#' @return tbl_df containing filtered expression data, with sampleinfo attr
#' @export
filter_genes <- function(df,
                         n_genes = 900,
                         genelist = rinfino::genelist,
                         include_random = TRUE) {
    dff <- df %>% dplyr::semi_join(genelist, by = '_GENE')
    if (include_random) {
        # select random set of genes to get to desired `n_genes`
        n_random <- n_genes - nrow(dff)
        random_genes <- df %>%
            dplyr::select(`_GENE`) %>%
            dplyr::anti_join(genelist, by = `_GENE`) %>%
            dplyr::sample_n(n_random)
        genelist <- dplyr::bind_rows(genelist, random_genes)
        dff <- df %>% dplyr::semi_join(genelist, by = `_GENE`)
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
run_combat <- function(df, mod = ~ 1, trans=log1p) {
    expmat <- expdata_as_matrix(df, trans=trans)
    sampleinfo <- attr(df, 'sampleinfo')
    mm <- model.matrix(mod, data=sampleinfo)
    # re-order to match order of columns in expmat
    sampleinfo <- sampleinfo[match(colnames(expmat), sampleinfo$`_SAMPLE_ID`),]
    sva::ComBat(as.matrix(expmat), sampleinfo$batch, mm)
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
#' @export
run_pca <- function(df, use_ggplot=TRUE, trans=log1p, group = 'batch', colour = 'batch', ...) {
    # transform df to SxG orientation, possibly applying `trans` function
    tdf <- transpose2sxg(df, trans=trans)

    # merge in sampleinfo (sample_id & batch)
    sampleinfo <- attr(df, 'sampleinfo')
    tdf <- tdf %>%
        dplyr::left_join(sampleinfo, by = '_SAMPLE_ID')

    # run pca
    pca <- prcomp(tdf %>% dplyr::select(-`_SAMPLE_ID`, -sample_id, -batch),
                  center = TRUE,
                  scale. = TRUE)

    if (use_ggplot) {
        pcadf <- tbl_df(pca$x) %>%
            bind_cols(sampleinfo)
        pcaplot <- ggplot2::ggplot(pcadf, ggplot2::aes_(x = 'PC1', y = 'PC2', group = group, colour = colour, ...)) +
            ggplot2::geom_point()
        print(pcaplot)
    } else {
        # produce ascii plot
        txtplot::txtplot(pca$x[,'PC1'], pca$x[,'PC2'], xlab = 'PC1', ylab = 'PC2')
    }
    pca
}

