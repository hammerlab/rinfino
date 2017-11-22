#' Prepare expression matrix of transcript- or gene-specific data
#' given a directory containing kallisto output (one folder per sample)
#'
#' @param root (chr) filepath to 'root' directory containing kallisto output
#' @param output_arg (chr) argument to `output_writer`, typically the filename to be written
#' @param skip_missing (bool) whether to report & ignore directories not containing output
#' @param metric (chr) which metric to capture from kallisto (see `?tximport` for details)
#' @param by (chr) description / name for column in output matrix containing gene or transcript
#' @param by_gene (bool) whether to organize results by gene (TRUE) or transcript (FALSE)
#' @param suffix (chr) file suffix to read in (either `tsv` or `h5` -- h5 gives errors on some platforms)
#' @param filename (chr) name of kallisto output file to ingest
#' @param output_writer (function) function to use when saving results.
#' @import dplyr readr glue stringr lubridate biomaRt tximport tibble purrr
#' @export
prep_expression_matrix <- function(root, output_arg,
                                   skip_missing=TRUE,
                                   cohort = basename(root),
                                   metric = 'scaledTPM',
                                   by = 'Gene_symbol',
                                   by_gene = TRUE,
                                   suffix = 'tsv',
                                   filename = glue::glue('abundance.{suffix}'),
                                   output_writer = readr::write_tsv) {

    # ---- Helper functions ----
    tmpdir <- tempdir()
    dir.create(tmpdir, recursive=T, showWarnings=F)
    shell <- function(cmdargs, command="bash") {
        # We need the full path to the main program
        command_path <- Sys.which(command)
        system(glue::glue('{command_path} {cmdargs}'))
    }

    # ---- index of kallisto output ----
    folders <- dir(root, full.names = T)
    df <- tbl_df(list(folder = folders, cohort = cohort)) %>%
        dplyr::mutate(run_id = basename(folder),
                      result_path = file.path(path.expand(folder), filename)
                      ) %>%
        dplyr::mutate(file_exists = purrr::map_chr(result_path, file.exists))

    # report how many files don't exist
    n_nonexist <- df %>% dplyr::filter(file_exists == FALSE) %>% nrow()
    if (n_nonexist>0)
        cat(glue::glue('Note: {n_nonexist} folders in {root} do not contain results. Skipping these.'))
    n_exist <- df %>% dplyr::filter(file_exists == TRUE) %>% nrow()
    cat(glue::glue('{n_exist} sample results identified in output folder {root}'))

    # optionally remove incomplete runs
    if (skip_missing)
        df <- df %>% dplyr::filter(file_exists == TRUE)

    # make sure all files exist
    if (df %>%
        dplyr::mutate(file_exists = purrr::map_chr(result_path, file.exists)) %>%
        dplyr::filter(file_exists != TRUE) %>%
        nrow() > 0)
        stop('Not all files exist - either resolve this or run with `skip_missing`')

    # summarise number of patients by cohort
    cohort_counts <- df %>%
                     dplyr::group_by(cohort) %>%
                     dplyr::summarize(`Number of Runs` = n_distinct(run_id))
    cohort_counts

    # copy kallisto files to tmpdir, so they are in one directory
    cmd_returns <- df %>%
                   dplyr::mutate(save_as = stringr::str_c(tmpdir, "/", run_id, '-', filename)) %>% # we want this to be easily parsed
                   dplyr::mutate(cp_args = stringr::str_c(result_path, " ", save_as)) %>%
                   dplyr::mutate(retcode = purrr::map_chr(cp_args, shell, command = 'cp'))

    stopifnot(all(cmd_returns$retcode == '0'))

    # get tx -> gene mapping
    mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host='ensembl.org')
    t2g <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "external_gene_name"), mart=mart)
    t2g <- dplyr::rename(t2g, target_id=ensembl_transcript_id, ext_gene=external_gene_name)
    head(t2g)

    # list of result files
    result_files <- dir(tmpdir, pattern=glue::glue("*\\.{suffix}"))
    result_files

    # Import and reduce
    txi <- tximport::tximport(paste(tmpdir, result_files, sep="/"), type="kallisto", countsFromAbundance = metric, txOut = TRUE)

    if (by_gene) {
        txi.sum <- tximport::summarizeToGene(txi, t2g)
        expression <- txi.sum$counts
    } else {
        expression <- txi$counts
    }

    # -abundance.tsv part of filename is redundant and can be omitted
    colnames(expression) <- result_files %>% stringr::str_replace(., pattern = glue::glue("-{filename}"), replacement = "")

    expdf <- expression %>%
        data.frame() %>%
        tibble::rownames_to_column(by)

    if (!is.null(output_arg) && !is.null(output_writer))
      expdf %>%
          output_writer(output_arg)

    cat(glue::glue('Output written to {output_arg}'), "\n")
    cat("Completed\n")

    expdf
}

