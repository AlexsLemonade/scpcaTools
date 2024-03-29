% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sce_to_seurat.R
\name{sce_to_seurat}
\alias{sce_to_seurat}
\title{Convert SingleCellExperiment object to Seurat object}
\usage{
sce_to_seurat(sce, assay_name = "counts")
}
\arguments{
\item{sce}{SingleCellExperiment object}

\item{assay_name}{The assay name (default "counts") to include in
the Seurat object. This name will be applied as the assay name in
the Seurat object. If the default "counts" assay is used, then the
assay name will instead be "RNA," consistent with Seurat defaults.}
}
\value{
Seurat object
}
\description{
Convert SingleCellExperiment object to Seurat object
}
\examples{
\dontrun{
sce_to_seurat(sce = sce_object)
}
}
