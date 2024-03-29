% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_ilisi.R
\name{calculate_ilisi}
\alias{calculate_ilisi}
\title{Calculate iLISI (integration Local Inverse Simpson's Index) scores}
\usage{
calculate_ilisi(merged_sce, pc_name, batch_column = "library_id")
}
\arguments{
\item{merged_sce}{The merged SCE object containing data from multiple batches}

\item{pc_name}{The name that allows access to the PCs. Example: fastMNN_PCA}

\item{batch_column}{The variable in `merged_sce` indicating the grouping of interest.
Generally this is either batches or cell types. Default is "library_id".}
}
\value{
Data frame with three columns, one row per cell. Columns are `ilisi_score`,
 `ilisi_score_norm`, `cell_barcode`, and `batch_id`. The column `ilisi_score_norm`
 is a normalized version if `ilisi_score` where values range from [0,1], inclusive
}
\description{
Calculate iLISI (integration Local Inverse Simpson's Index) scores
}
