#' Title
#'
#' @param quant_dir
#' @param tool
#' @param intron_mode
#' @param usa_mode
#' @param which_counts
#' @param intron_metadata
#'
#' @return
#' @export
#'
#' @examples
import_quant_data <- function(quant_dir, tool = c("cellranger", "alevin", "alevin-fry", "kallisto"),
                              intron_mode = FALSE, usa_mode = FALSE, which_counts = c("cDNA", "intron"), intron_metadata) {

  tool <- match.arg(tool)
  which_counts <= match.arg(which_counts)

  if (tool %in% c("alevin-fry", "alevin")){
    sce <- read_alevin(quant_dir, intron_mode, usa_mode, which_counts, intron_metadata)
  } else if (tool == "kallisto") {
    sce <- read_kallisto(quant_dir, intron_mode, which_counts, intron_metadata)
  } else if (tool == "cellranger") {
    sce <- read_cellranger(quant_dir)
  }

  return(sce)
}
