#' Generate a QC report from a SingleCellExperiment object
#'
#' @param sample_name The name of the sample for report headers
#' @param unfiltered_sce A SingleCellExperiment object that the report will describe
#' @param filtered_sce An optional filtered single cell experiment derived from first
#' @param output The output file path that will be created.
#'   If the file name does not include an extension, ".html" will be added automatically.
#'   Any directories in the path will be created as needed.
#'
#' @return The full path of the output file
#' @export
#'
#'
#' @examples
#' \dontrun{
#' generate_qc_report("Sample 1", my_sce, output = "reports/sample1_report.html")
#' }
#'
generate_qc_report <- function(sample_name,
                               unfiltered_sce,
                               filtered_sce = NULL,
                               output = NULL){
  if(!inherits(unfiltered_sce, "SingleCellExperiment")){
    stop("`unfiltered_sce` must be a SingleCellExperiment object.")
  }
  # check that the filtered sce is as expected
  if(!is.null(filtered_sce)){
    if (!inherits(filtered_sce, "SingleCellExperiment")){
      stop("`filtered_sce` must be a SingleCellExperiment object.")
    }
    if (ncol(unfiltered_sce) < ncol(filtered_sce)){
      stop("`filtered_sce` should have fewer cells than `unfiltered_sce`")
    }
    if (!all(colnames(filtered_sce) %in% colnames(unfiltered_sce))){
      "Some cells in `filtered_sce` are not present in `unfiltered_sce`, from which it should be derived."
    }
  }

  if(is.null(output)){
    output_file = glue::glue("{sample_name}_qc_report")
    output_dir = "."
  } else {
    output_file = basename(output)
    output_dir = dirname(output)
  }

  rmd <- system.file(file.path("rmd", "qc_report.rmd"), package = "scpcaTools")
  suppressPackageStartupMessages(rmarkdown::render(
    rmd,
    output_file = output_file,
    output_dir = output_dir,
    params = list(
      sample = sample_name,
      unfiltered_sce = unfiltered_sce,
      filtered_sce = filtered_sce
    ),
    envir = new.env(),
    quiet = TRUE
  ))
}
