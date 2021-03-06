% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_counts_matrix_from_annotation.R
\name{rename_counts_matrix_from_annotation}
\alias{rename_counts_matrix_from_annotation}
\title{rename_counts_matrix_from_annotation}
\usage{
rename_counts_matrix_from_annotation(
  counts_matrix,
  annotation_table,
  new_columnname,
  old_columnname = ""
)
}
\arguments{
\item{counts_matrix}{a matrix of gene x cells counts}

\item{annotation_table}{annotation table containing at least one column with the counts matrix column names, and another column we want to use to rename the counts matrix.}

\item{new_columnname}{the columnname in the annotation_table, we want to use as the new counts matrix column names}

\item{old_columnname}{(optional) the columnname in the annotation_table containing the current counts matrix columnnames}
}
\value{
a renamed counts matrix
}
\description{
Renames the counts matrix column names from a specified column in a provided annotation table
}
