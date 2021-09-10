#' Generate a QC report from a SingleCellExperiment object
#'
#' @param sce A SingleCellExperiment object that the report will describe
#' @param sample_name The name of the sample
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
#' generate_qc_report(my_sce, "Sample 1", output = "reports/sample1_report.html")
#' }
#'
generate_qc_report <- function(sce,
                               sample_name,
                               filtered_sce = NULL,
                               output = NULL){
  if(!inherits(sce, "SingleCellExperiment")){
    stop("`sce` must be a SingleCellExperiment object.")
  }
  # check that the filtered sce is as expected
  if(!is.null(filtered_sce)){
    if (!inherits(filtered_sce, "SingleCellExperiment")){
      stop("`filtered_sce` must be a SingleCellExperiment object.")
    }
    if (ncol(sce) < ncol(filtered_sce)){
      stop("`filtered_sce` should have fewer cells than `sce`")
    }
    if (!all(colnames(filtered_sce) %in% colnames(sce))){
      "Some cells in `filtered_sce` are not present in `sce`, from which it should be derived."
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
  rmarkdown::render(
    rmd,
    output_file = output_file,
    output_dir = output_dir,
    params = list(
      sce = sce,
      sample = sample_name,
      filtered_sce = filtered_sce
    ),
    envir = new.env()
  )
}
