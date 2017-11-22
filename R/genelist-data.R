#' Selected genes for use with infino models
#' 
#' Used as a default list of genes when running `filter_genes`. Gene list is
#' comprised of known marker genes for immune deconvolution (e.g. from Bindea et al),
#' as well as novel genes determined to be differentially expressed among cell types
#' in the default _infino_ training set.
#'
#' @docType data
#' 
#' @usage data(genelist)
#' 
#' @format an object of class \code{"tbl_df"}
#' 
#' @keywords datasets
#'
#' @references Bindea et al (2013) Immunity Oct 17;39(4):782-95 
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24138885}{PubMed})
"genelist"
