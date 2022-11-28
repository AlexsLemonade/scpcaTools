% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_sce_list.R
\name{merge_sce_list}
\alias{merge_sce_list}
\title{Merge a list of SCEs as preparation for formal integration}
\usage{
merge_sce_list(
  sce_list = list(),
  batch_column = "library_id",
  barcode_column = "barcode",
  retain_coldata_cols = c("sum", "detected", "total", "subsets_mito_sum",
    "subsets_mito_detected", "subsets_mito_percent", "miQC_pass", "prob_compromised"),
  preserve_rowdata_cols = NULL
)
}
\arguments{
\item{sce_list}{A list of SingleCellExperiment objects. The list may optionally
be named with batch information. If no names are provided, names will be generated
based on the SCE's index.}

\item{batch_column}{A character value giving the resulting colData column name
to differentiate originating SingleCellExperiment objects. Often these values
are library IDs. Default value is `library_id`.}

\item{barcode_column}{A character value giving the colData column name that
stores cell barcodes. If this column does not yet exist it will be created.
Default value is  `barcode`.}

\item{retain_coldata_cols}{A vector of colData columns which should be retained
in the the final merged SCE object.}

\item{preserve_rowdata_cols}{A vector of column names that appear in originating
SCE objects' rowData slots which should not be renamed with
the given SCE object name or index is name is not given. These are generally
columns which are not specific to the given library's preparation or statistics.
For example, such a vector might contain items like "Gene", "ensembl_ids", etc.}
}
\value{
A SingleCellExperiment object containing all SingleCellExperiment objects
  present in the inputted list
}
\description{
This function takes an optionally-named (if named, ideally by a form of
 library ID) list of SingleCellExperiment (SCE) objects and merges them into
 one SCE object. The resulting colData slot will contain information in the
 given batch column to differentiate originating SCE objects based on the SCE
 names, or index if no names are provided. If not already present, cell barcodes
 will also be created as a colData column. Only features (genes) that are
 present in all provided SCE objects will be retained. To preserve original
 feature-level information, all final rowData slot column names will be
 appended with the given SCE's name, as `<column_name>-<sce_name>` sxcept for
 columns indicated to preserve. Currently any present altExps are not retained.
}