#' Import Gene Expression Quantification Data for Single-Cell RNA-Seq
#'
#' Imports the gene x cell matrix output from either Alevin, Alevin-fry, Cellranger, or Kallisto and returns a SingleCellExperiment.
#'
#' @param quant_dir Path to directory where output files are located.
#' @param tool Type of tool used to create files (alevin, alevin-fry, cellranger, or kallisto).
#' @param include_unspliced Whether or not to include the unspliced reads in the counts matrix.
#'   If TRUE, the main "counts" assay will contain unspliced reads and spliced reads and an additional "spliced"
#'   assay will contain spliced reads only. If TRUE, requires that data has been aligned to a reference contianing
#'   spliced and unspliced reads.
#'   Default is TRUE.
#' @param usa_mode Logical indicating if Alevin-fry was used, if the USA mode was invoked.
#'   Default is FALSE.
#' @param filter Logical indicating whether or not to filter the counts matrix.
#'   Filtering is performed using DropletUtils::emptyDrops and cannot be performed with Cellranger.
#' @param fdr_cutoff FDR cutoff to use for DropletUtils::emptyDrops.
#'   Default is 0.01.
#' @param tech_version Technology or kit used to process library (i.e. 10Xv3, 10Xv3.1).
#' @param ... Any arguments to be passed into DropletUtils::emptyDrops.
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # read in single cell RNA seq data processed using cellranger
#' import_quant_data(quant_dir,
#'   tool = "cellranger"
#' )
#'
#' # read in single-cell RNA seq data processed using alevin-fry in USA mode
#' # with alignment to cDNA only and including counts for spliced cDNA only
#' import_quant_data(quant_dir,
#'   tool = "alevin-fry",
#'   include_unspliced = FALSE,
#'   usa_mode = TRUE
#' )
#'
#' # read in single-nuclei RNA-seq data processed using alevin-fry in
#' # USA mode with alignment to cDNA + introns and including counts for
#' # unspliced cDNA and perform filtering
#' import_quant_data(quant_dir,
#'   tool = "alevin-fry",
#'   usa_mode = TRUE,
#'   filter = TRUE
#' )
#'
#' # read in single-nuclei RNA-seq data processed using kallisto with
#' # alignment to cDNA + introns and including counts for unspliced cDNA
#' # and perform filtering
#' import_quant_data(quant_dir,
#'   tool = "kallisto",
#'   filter = TRUE
#' )
#' }
#'
import_quant_data <- function(quant_dir,
                              tool = c("cellranger", "alevin", "alevin-fry", "kallisto"),
                              include_unspliced = TRUE,
                              usa_mode = FALSE,
                              filter = FALSE,
                              fdr_cutoff = 0.01,
                              tech_version = NULL,
                              ...) {
  which_counts <- match.arg(which_counts)

  if (!(tool %in% c("cellranger", "alevin", "alevin-fry", "kallisto"))) {
    stop("Tool must be either cellranger, alevin, alevin-fry, or kallisto.")
  }

  # checks for intron_mode and usa_mode
  if (!is.logical(include_unspliced)) {
    stop("include_unspliced must be set as TRUE or FALSE")
  }
  if (!is.logical(usa_mode)) {
    stop("usa_mode must be set as TRUE or FALSE")
  }
  if (!is.logical(filter)) {
    stop("filter must be set as TRUE or FALSE")
  }

  # check that usa_mode and intron_mode are used with the proper tools
  if (usa_mode & tool %in% c("cellranger", "alevin", "kallisto")) {
    stop("USA mode only compatible with alevin-fry.")
  }
  if (include_unspliced & tool %in% c("cellranger")) {
    stop("Include unspliced not compatible with cellranger.")
  }

  # check that filter is not used with cellranger
  if (filter & tool == "cellranger") {
    stop("Cannot perform emptyDrops filtering on cellranger output.")
  }
  if (filter) {
    if (!(is.numeric(fdr_cutoff))) {
      stop("fdr_cutoff is not a number.")
    }
    if (fdr_cutoff < 0 | fdr_cutoff > 1) {
      stop("fdr_cutoff must be a number between 0 - 1.")
    }
  }

  if (tool %in% c("alevin-fry", "alevin")) {
    sce <- read_alevin(quant_dir, include_unspliced, usa_mode, tech_version)
  } else if (tool == "kallisto") {
    sce <- read_kallisto(quant_dir, include_unspliced)
  } else if (tool == "cellranger") {
    sce <- read_cellranger(quant_dir)
  }

  if (filter) {
    sce <- filter_counts(sce, fdr_cutoff, ...)
  }

  return(sce)
}
