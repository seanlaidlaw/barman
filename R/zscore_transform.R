#' @title zscore_transform
#'
#' @description apply a zscore transformation to a list of columns based on average and SD of a column names provided by all_cells_to_normalise arg. Zscore should only be applied on counts matrix containing values that follow a normal distribution.
#'
#' @param all_cells_to_normalise list of column names that we want to zscore transform
#' @param counts_matrix counts matrix with normally distrubted data
#'
#' @return zscore transformed counts matrix
#'
#' @export
zscore_transform <- function(all_cells_to_normalise, counts_matrix) {

	# sanity check given arguments
	# if (!any(control_groups %in% all_cells_to_normalise)) {
	# 	print("Error: none of the control group cells are in the all_cells_to_normalise list.")
	# 	return(NA)
	# }
	if (!any(all_cells_to_normalise %in% colnames(counts_matrix))) {
		print("Error: none of the cells in all_cells_to_normalise list are counts_matrix column names.")
		return(NA)
	}

	zscore_counts <- counts_matrix
	rownames(zscore_counts) <- zscore_counts$Geneid
	zscore_counts[zscore_counts == 0] <- NA

	# filter list to only those in matrix
	all_cells_to_normalise <- all_cells_to_normalise[all_cells_to_normalise %in% colnames(zscore_counts)]

	# calculate Z scores
	zscore_counts$mean <- rowMeans(as.matrix(zscore_counts[, all_cells_to_normalise]))
	zscore_counts$stdD <- matrixStats::rowSds(as.matrix(zscore_counts[, all_cells_to_normalise]))
	zscore_counts <- as.data.frame(t(apply(zscore_counts[, c("mean", "stdD", all_cells_to_normalise)], 1, function(x) ((x - x[["mean"]]) / x[["stdD"]]))))
	zscore_counts$mean <- NULL
	zscore_counts$stdD <- NULL

	# readd old metadata columns
	zscore_counts$Geneid <- rownames(zscore_counts)
	zscore_counts <- merge(zscore_counts, counts_matrix[, c("Chr", "Geneid", "Start", "End", "Strand", "Length")], by = "Geneid")
	zscore_counts <- zscore_counts[zscore_counts$Geneid != "", ] # remove empty geneid rows

	# rearrange column order
	zscore_counts <- zscore_counts[, c("Geneid", "Chr", "Start", "End", "Strand", "Length", all_cells_to_normalise)]

	return(zscore_counts)
}
