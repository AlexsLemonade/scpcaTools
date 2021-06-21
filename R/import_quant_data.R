#' Import Gene Expression Quantification Data for Single-Cell RNA-Seq
#'
#' Imports the gene x cell matrix output from either Alevin, Alevin-fry, Cellranger, or Kallisto and returns a SingleCellExperiment.
#'
#' @param quant_dir Full path to directory where output files are located.
#' @param tool Type of tool used to create files (Alevin, Alevin-fry, Cellranger, or Kallisto).
#' @param intron_mode Boolean indicating if the files included alignment to intronic regions. Default is FALSE.
#' @param usa_mode Boolean indicating if Alevin-fry was used, if the USA mode was invoked. Default is FALSE.
#' @param which_counts If intron_mode is TRUE, which type of counts should be included,
#'        only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'        Default is "spliced".
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
                              intron_mode = FALSE, usa_mode = FALSE, which_counts = c("spliced", "unspliced")) {

  tool <- match.arg(tool)
  which_counts <- match.arg(which_counts)

  # checks for intron_mode and usa_mode
  if(!is.boolean(intron_mode)){
    stop("intron_mode must be set as TRUE or FALSE")
  }
  if(!is.boolean(usa_mode)){
    stop("usa_mode must be set as TRUE or FALSE")
  }

  ## check that intron_metadata_path is provided
  # if(intron_mode == TRUE & usa_mode == FALSE){
  #   if(!file.exists(intron_metadata_path)){
  #     stop("Missing intron metadata file, check that file path is correct")
  #   } else {
  #     intron_metadata <- readr::read_tsv(intron_metadata_path)
  #     if(colnames(intron_metadata != c("spliced", "intron"))) {
  #       stop("Incorrect column names for intron metadata")
  #     }
  #   }
  # }

  if (tool %in% c("alevin-fry", "alevin")){
    sce <- read_alevin(quant_dir, intron_mode, usa_mode, which_counts, intron_metadata_path)
  } else if (tool == "kallisto") {
    sce <- read_kallisto(quant_dir, intron_mode, which_counts, intron_metadata_path)
  } else if (tool == "cellranger") {
    sce <- read_cellranger(quant_dir)
  }

  return(sce)
}
