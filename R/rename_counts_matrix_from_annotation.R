#' @title rename_counts_matrix_from_annotation
#'
#' @description Renames the counts matrix column names from a specified column in a provided annotation table
#'
#' @param counts_matrix a matrix of gene x cells counts
#' @param annotation_table annotation table containing at least one column with the counts matrix column names, and another column we want to use to rename the counts matrix.
#' @param new_columnname the columnname in the annotation_table, we want to use as the new counts matrix column names
#' @param old_columnname (optional) the columnname in the annotation_table containing the current counts matrix columnnames
#'
#' @return a renamed counts matrix
#'
#' @export

rename_counts_matrix_from_annotation <- function(counts_matrix, annotation_table, new_columnname, old_columnname="") {
	# get current non-metadata column names
	sample_colnames = colnames(counts_matrix)[!colnames(counts_matrix) %in% c("Geneid","Chr","Start","End","Strand","Length")]

	# if old_columnname not given then determine from annotation_table
	if (old_columnname == "") {
		old_columnname = ""
		for (annotation_colname in colnames(annotation_table)) {
			if (all(sample_colnames %in% annotation_table[[annotation_colname]])) {
				old_columnname = annotation_colname
			}
		}
		if (old_columnname == "") {
			stop("Couldn't automatically determine which column of annotation corresponded to colnames of counts matrix. please manually specify with the old_columnname argument")
		}
	}

	new_sample_colnames = annotation_table[[new_columnname]][match(sample_colnames, annotation_table[[old_columnname]])]

	# make sure provided new column names are unique and equal in length to old column names
	if (length(annotation_table[[new_columnname]]) != length(annotation_table[[old_columnname]])) {
		stop("Length of new_columnname not equal to length of old_columnname")
	}
	if (length(unique(annotation_table[[new_columnname]])) != length((annotation_table[[new_columnname]]))) {
		stop("new_columnname column contains non-unique values")
	}


	# type conversion to character important to not have factors become digits when adding to this vector
	new_sample_colnames = c("Geneid","Chr","Start","End","Strand","Length", as.character(new_sample_colnames))

	# rename columns of counts_matrix dataframe to new_sample_colnames
	counts_matrix = setNames(counts_matrix, new_sample_colnames)

	return(counts_matrix)
}
