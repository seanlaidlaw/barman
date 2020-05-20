#' @title rename_counts_matrix_from_annotation
#'
#' @description Renames the counts matrix column names from a provided annotation table
#'
#' @param counts_matrix the featurecounts matrix
#' @param annotation_table annotation table containing at least one column with the featurecounts column names, and another column we want to use to rename the counts matrix.
#' @param dna_cell_id the columnname in the annotation_table, we want to use as the new featurecounts column names
#' @param rna_cell_id (optional) the columnname in the annotation_table containing the current featurecounts columnnames
#'
#' @return a renamed counts matrix
#'
#' @export

rename_counts_matrix_from_annotation <- function(counts_matrix, annotation_table, dna_cell_id, rna_cell_id="") {
	# we need to match the header of the featurecounts matrix with a Cell_ID of some sort, from the annotation table so that we can have the same ID for our DNA and RNA data.
	# this renames the counts matrix headers to the specified Cell_ID column from the annotation.
	sample_colnames = colnames(counts_matrix)[!colnames(counts_matrix) %in% c("Geneid","Chr","Start","End","Strand","Length")]

	# if rna_cell_id not given then determine from annotation_table
	if (rna_cell_id == "") {
		rna_cell_id = ""
		for (annotation_colname in colnames(annotation_table)) {
			if (all(sample_colnames %in% annotation_table[[annotation_colname]])) {
				rna_cell_id = annotation_colname
			}
		}
		if (rna_cell_id == "") {
			stop("Couldn't automatically determine which column of annotation corresponded to colnames of counts matrix. please manually specify with the rna_cell_id argument")
		}
	}

	new_sample_colnames = annotation_table[[dna_cell_id]][match(sample_colnames, annotation_table[[rna_cell_id]])]
	colnames(counts_matrix) = append(c("Geneid","Chr","Start","End","Strand","Length"), new_sample_colnames)
	return(counts_matrix)
}
