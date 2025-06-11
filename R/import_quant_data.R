#' Import Gene Expression Quantification Data for Single-Cell RNA-Seq
#'
#' Imports the gene x cell matrix output from either Alevin, Alevin-fry, Cell Ranger, or Kallisto
#'   and returns a SingleCellExperiment.
#'
#' @param quant_dir Path to directory where output files are located.
#'   Used when tool is alevin, alevin-fry, or kallisto.
#' @param cellranger_h5_file Optional path to H5 file output by Cell Ranger.
#'   Used when tool is cellranger.
#' @param tool Type of tool used to create files (alevin, alevin-fry, cellranger, or kallisto).
#' @param include_unspliced Whether or not to include the unspliced reads in the counts matrix.
#'   If TRUE, the main "counts" assay will contain unspliced reads and spliced reads and an additional "spliced"
#'   assay will contain spliced reads only. If TRUE, requires that data has been aligned to a reference containing
#'   spliced and unspliced reads.
#'   Default is TRUE.
#' @param usa_mode Logical indicating if Alevin-fry was used, if the USA mode was invoked.
#'   Default is FALSE.
#' @param filter Logical indicating whether or not to filter the counts matrix.
#'   Filtering is performed using DropletUtils::emptyDrops and cannot be performed with Cell Ranger.
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
#' import_quant_data(cellranger_h5_file,
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
                              cellranger_h5_file = "",
                              tool = c("cellranger", "alevin", "alevin-fry", "kallisto"),
                              include_unspliced = TRUE,
                              usa_mode = FALSE,
                              filter = FALSE,
                              fdr_cutoff = 0.01,
                              tech_version = NULL,
                              ...) {
  which_counts <- match.arg(which_counts)

  stopifnot(
    "Tool must be one of cellranger, alevin, alevin-fry, or kallisto." =
      tool %in% c("cellranger", "alevin", "alevin-fry", "kallisto"),
    "h5_file must be provided when using cellranger as tool" = tool == "cellranger" & file.exists(cellranger_h5_file),
    "include_unspliced must be set as TRUE or FALSE" = is.logical(include_unspliced),
    "usa_mode must be set as TRUE or FALSE" = is.logical(usa_mode),
    "filter must be set as TRUE or FALSE" = is.logical(filter),
    "USA mode only compatible with alevin-fry." = !(usa_mode && tool %in% c("cellranger", "alevin", "kallisto")),
    "Include unspliced not compatible with cellranger." = !(include_unspliced && tool %in% c("cellranger")),
    "Cannot perform emptyDrops filtering on cellranger output." = !(filter && tool == "cellranger")
  )

  if (filter) {
    stopifnot(
      "fdr_cutoff must be a number." = is.numeric(fdr_cutoff),
      "fdr_cutoff must be a number between 0 - 1." = fdr_cutoff >= 0 && fdr_cutoff <= 1
    )
  }

  if (tool %in% c("alevin-fry", "alevin")) {
    sce <- read_alevin(quant_dir, include_unspliced, usa_mode, tech_version)
  } else if (tool == "kallisto") {
    sce <- read_kallisto(quant_dir, include_unspliced)
  } else if (tool == "cellranger") {
    sce <- read_cellranger(cellranger_h5_file)
  }

  if (filter) {
    sce <- filter_counts(sce, fdr_cutoff, ...)
  }

  return(sce)
}
