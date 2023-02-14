#' Generate a QC report from a SingleCellExperiment object
#'
#' @param library_id The name of the library_id for report headers
#' @param unfiltered_sce A SingleCellExperiment object that the report will describe
#' @param filtered_sce An optional filtered single cell experiment derived from first
#' @param processed_sce An optional single cell experiment that has been normalized and
#'   contains PCA and UMAP embeddings
#' @param rmd_file An optional path to the rmd file to be rendered, if no file is provided,
#'   the default `qc_report.rmd` file present in the package will be used.
#' @param extra_params An optional named list of additional parameters to use when
#'   rendering the provided rmd file.
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
                               rmd_file = NULL,
                               extra_params = c(),
                               output = NULL,
                               ...){

  ### Check dependencies for generating report
  required_packages <- c("ggplot2", "kableExtra", "rmarkdown",
                         "scater", "scran")

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

  # create list of parameters
  params_list <- list(
    library = library_id,
    unfiltered_sce = unfiltered_sce,
    filtered_sce = filtered_sce,
    processed_sce = processed_sce
  )

  # define rmd file if none provided
  if(is.null(rmd_file)){
    rmd_file <- system.file(file.path("rmd", "qc_report.rmd"), package = "scpcaTools")
  } else {
    # rmd file is provided, check that it exists
    if(!file.exists(rmd_file)){
      stop("provided rmd file does not exist.")
    }

    # if additional rmd file is provided, add extra params
    params_list <- append(params_list,
                          extra_params)
  }

  suppressPackageStartupMessages(rmarkdown::render(
    rmd_file,
    output_file = output_file,
    output_dir = output_dir,
    params = params_list,
    intermediates_dir = tempdir(),
    knit_root_dir = tempdir(),
    envir = new.env(),
    quiet = TRUE,
    ...
  ))
}
