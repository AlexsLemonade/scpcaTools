% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_qc_report.R
\name{generate_qc_report}
\alias{generate_qc_report}
\title{Generate a QC report from a SingleCellExperiment object}
\usage{
generate_qc_report(
  library_id,
  unfiltered_sce,
  filtered_sce = NULL,
  processed_sce = NULL,
  report_template = NULL,
  extra_params = NULL,
  output = NULL,
  ...
)
}
\arguments{
\item{library_id}{The name of the library_id for report headers}

\item{unfiltered_sce}{A SingleCellExperiment object that the report will describe}

\item{filtered_sce}{An optional filtered single cell experiment derived from first}

\item{processed_sce}{An optional single cell experiment that has been normalized and
contains PCA and UMAP embeddings}

\item{report_template}{An optional path to a template rmd file to be rendered,
if no file is provided, the default `qc_report.rmd` file present in the package will be used.}

\item{extra_params}{An optional named list of additional parameters to use when
rendering the provided rmd file.}

\item{output}{The output file path that will be created.
If the file name does not include an extension, ".html" will be added automatically.
Any directories in the path will be created as needed.}

\item{...}{Additional arguments to pass to rmarkdown::render()}
}
\value{
The full path of the output file
}
\description{
Generate a QC report from a SingleCellExperiment object
}
\examples{
\dontrun{
generate_qc_report("Library 1", my_sce, output = "reports/sample1_report.html")
}

}
