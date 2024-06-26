% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_sce_list.R
\name{get_altexp_metadata}
\alias{get_altexp_metadata}
\title{Helper function to get altExp metadata from an SCE that may not have the altexp
Returns main experiment metadata if the altExp is not present in the SCE}
\usage{
get_altexp_metadata(sce, altexp_name)
}
\arguments{
\item{sce}{SCE object to extract metadata from}

\item{altexp_name}{Name of the altExp to extract metadata for}
}
\value{
List with fields `library_id` and `sample_id`
}
\description{
Helper function to get altExp metadata from an SCE that may not have the altexp
Returns main experiment metadata if the altExp is not present in the SCE
}
