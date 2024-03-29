% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrate_sces.R
\name{integrate_fastmnn}
\alias{integrate_fastmnn}
\title{Integrate a list of SingleCellExperiment (SCE) objects using the fastMNN
approach from the `batchelor` package.
Source: http://www.bioconductor.org/packages/release/bioc/html/batchelor.html}
\usage{
integrate_fastmnn(merged_sce, batch_column, ...)
}
\arguments{
\item{merged_sce}{A merged SCE object as prepared by `scpcaTools::merge_sce_list()`.}

\item{batch_column}{The variable in `merged_sce` indicating batches.}

\item{...}{Additional arguments to pass into `batchelor::fastMNN()`}
}
\value{
An integrated SCE object produced by `fastMNN`
}
\description{
Integrate a list of SingleCellExperiment (SCE) objects using the fastMNN
approach from the `batchelor` package.
Source: http://www.bioconductor.org/packages/release/bioc/html/batchelor.html
}
