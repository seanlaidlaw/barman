#' @title plot_rna_clusters
#'
#' @description plots UMAP of a filtered and normalized SingleCellExperiement object
#'
#' @param normalised_counts_matrix a filtered, normalized, and log transformed counts matrix
#' @param annotation_table a table of cell annotations, with  rows corresponding to colnames in normalised_counts_matrix
#' @param umap_vars a vector of 2 colnames from annotation table, to manually specifiy the plot's
#' colour and shape values respectively
#'
#' @return ggplot2 object of UMAP
#'
#' @export


plot_rna_clusters <- function(normalised_counts_matrix, annotation_table, umap_vars) {
	annotation_table_cols = colnames(annotation_table)

	# trim metadata columns from counts matrix
	trimmed_counts_matrix = normalised_counts_matrix
	rownames(trimmed_counts_matrix) = trimmed_counts_matrix$Geneid
	trimmed_counts_matrix$Geneid = NULL
	trimmed_counts_matrix$Chr = NULL
	trimmed_counts_matrix$Start = NULL
	trimmed_counts_matrix$End = NULL
	trimmed_counts_matrix$Strand = NULL
	trimmed_counts_matrix$Length = NULL

	# trim annotation table to be same as counts matrix
	annotation_table = annotation_table[rownames(annotation_table) %in% colnames(trimmed_counts_matrix),]
	annotation_table = as.data.frame(annotation_table)
	colnames(annotation_table) = annotation_table_cols


	# construct SingleCellExperiment object from normalized tables
	normalised_SCE_object <- SingleCellExperiment::SingleCellExperiment(
		assays = list(logcounts = as.matrix(trimmed_counts_matrix)),
		colData = annotation_table
	)


	# if only one var given
	if (length(umap_vars) == 1) {
			pc1_major_var = umap_vars

			plt = scater::plotUMAP(
				normalised_SCE_object,
				colour_by = pc1_major_var)
		} else {
			pc1_major_var = umap_vars[1]
			pc2_major_var = umap_vars[2]

			plt = scater::plotUMAP(
				normalised_SCE_object,
				colour_by = pc1_major_var,
				shape_by = pc2_major_var)
	}

	return(plt)
}
