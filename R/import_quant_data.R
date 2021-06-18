#' Import Gene Expression Quantification Data for Single-Cell RNA-Seq
#'
#' Imports the gene x cell matrix output from either Alevin, Alevin-fry, Cellranger, or Kallisto and returns a SingleCellExperiment.
#'
#' @param quant_dir Full path to directory where output files are located.
#' @param tool Type of tool used to create files (Alevin, Alevin-fry, Cellranger, or Kallisto).
#' @param intron_mode Boolean indicating if the files included alignment to intronic regions. Default is FALSE.
#' @param usa_mode Boolean indicating if Alevin-fry was used, if the USA mode was invoked. Default is FALSE.
#' @param which_counts The type of counts (cDNA or intron) to include if alignment to intronic regions is TRUE. Default is FALSE.
#' @param intron_metadata Full path to a two column tsv file containing gene names for both spliced and intronic regions.
#'           Only required if intron_mode = TRUE and usa_mode = FALSE.
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' import_quant_data(quant_dir,
#' tool = "cellranger")
#'
#' import_quant_data(quant_dir,
#' tool = "alevin-fry",
#' intron_mode = FALSE,
#' usa_mode = TRUE)
#' }
#'
import_quant_data <- function(quant_dir, tool = c("cellranger", "alevin", "alevin-fry", "kallisto"),
                              intron_mode = FALSE, usa_mode = FALSE, which_counts = c("cDNA", "intron"), intron_metadata_path) {

  tool <- match.arg(tool)
  which_counts <- match.arg(which_counts)

  if (tool %in% c("alevin-fry", "alevin")){
    sce <- read_alevin(quant_dir, intron_mode, usa_mode, which_counts, intron_metadata_path)
  } else if (tool == "kallisto") {
    sce <- read_kallisto(quant_dir, intron_mode, which_counts, intron_metadata_path)
  } else if (tool == "cellranger") {
    sce <- read_cellranger(quant_dir)
  }

  return(sce)
}
