% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer_signif_annot_cols.R
\name{infer_signif_annot_cols}
\alias{infer_signif_annot_cols}
\title{infer_signif_annot_cols}
\usage{
infer_signif_annot_cols(normalised_counts_matrix, annotation_table)
}
\arguments{
\item{normalised_counts_matrix}{a filtered, normalized, and log transformed counts matrix}

\item{annotation_table}{a table of cell annotations, with  rows corresponding to colnames in normalised_counts_matrix}
}
\value{
a vector of annotation_table column names
}
\description{
infers what Principal Components should be analysed. Determines what columns in
annotation table represent the highest variance in the normalised counts matrix
}
