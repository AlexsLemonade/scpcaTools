#' Generate a QC report from a SingleCellExperiment object
#'
#' @param sce
#' @param sample_name
#' @param outfile
#' @param ...
#'
#' @return The path of the output file
#' @export
#'
#'
#' @examples
#'
generate_qc_report <- function(sce,
                               sample_name,
                               output = NULL,
                               ...){
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
