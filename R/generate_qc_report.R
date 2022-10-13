#' Generate a QC report from a SingleCellExperiment object
#'
#' @param library_id The name of the library_id for report headers
#' @param unfiltered_sce A SingleCellExperiment object that the report will describe
#' @param filtered_sce An optional filtered single cell experiment derived from first
#' @param processed_sce An optional single cell experiment that has been normalized and
#'   contains PCA and UMAP embeddings
#' @param output The output file path that will be created.
#'   If the file name does not include an extension, ".html" will be added automatically.
#'   Any directories in the path will be created as needed.
#' @param ... Additional arguments to pass to rmarkdown::render()
#'
#' @return The full path of the output file
#' @export
#'
#'
#' @examples
#' \dontrun{
#' generate_qc_report("Library 1", my_sce, output = "reports/sample1_report.html")
#' }
#'
generate_qc_report <- function(library_id,
                               unfiltered_sce,
                               filtered_sce = NULL,
                               processed_sce = NULL,
                               output = NULL,
                               ...){
  ### Check dependencies in report template
  required_packages <- c("ggplot2", "kableExtra", "scater", "scran")

  installed_status <- sapply(required_packages,
                             requireNamespace,
                             quietly = TRUE)

  if(!all(installed_status)){
    missing_packages <- required_packages[!installed_status]|>
      paste(collapse = ", ")
    stop(paste("The following packages are required for generating a QC report:",
               missing_packages,
               sep = "\n       "))
  }

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
      stop("Some cells in `filtered_sce` are not present in `unfiltered_sce`, from which it should be derived.")
    }
  }

  if(is.null(output)){
    output_file = glue::glue("{library_id}_qc_report")
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
      library = library_id,
      unfiltered_sce = unfiltered_sce,
      filtered_sce = filtered_sce,
      processed_sce = processed_sce
    ),
    intermediates_dir = tempdir(),
    knit_root_dir = tempdir(),
    envir = new.env(),
    quiet = TRUE,
    ...
  ))
}
