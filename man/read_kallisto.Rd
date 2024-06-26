% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_kallisto.R
\name{read_kallisto}
\alias{read_kallisto}
\title{Read in counts data processed with Kallisto}
\usage{
read_kallisto(quant_dir, include_unspliced = TRUE, round_counts = TRUE)
}
\arguments{
\item{quant_dir}{Path to directory where output files are located.}

\item{include_unspliced}{Whether or not to include the unspliced reads in the counts matrix.
If TRUE, the main "counts" assay will contain unspliced reads and spliced reads and an additional "spliced"
assay will contain spliced reads only. If TRUE, requires that data has been aligned to a reference containing
spliced and unspliced reads.
Default is TRUE.}

\item{round_counts}{Logical indicating in the count matrix should be rounded to integers on import.
Default is TRUE.}
}
\value{
SingleCellExperiment of unfiltered gene x cell counts matrix
}
\description{
Read in counts data processed with Kallisto
}
\examples{
\dontrun{

# import output files processed with kallisto with alignment to cDNA + introns,
# including all unspliced cDNA counts in final counts matrix
read_kallisto(quant_dir,
  include_unspliced = TRUE,
  which_counts = "unspliced"
)
}
}
