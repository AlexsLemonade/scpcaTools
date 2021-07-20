#' Generate a QC report from a SingleCellExperiment object
#'
#' @param sce A SingleCellExperiment object that the report will describe
#' @param sample_name The name of the sample
#' @param outfile The output file path that will be created.
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
                               output = NULL){
  if(!inherits(sce,"SingleCellExperiment")){
    stop("`sce` must be a SingleCellExperiment object.")
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
      sample = sample_name
    ),
    envir = new.env()
  )
}
