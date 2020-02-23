#' @title barman
#'
#' @description apply a zscore transformation to a list of columns based on average and SD of a control list of column names. Zscore should only be applied on counts matrix containing values that follow a normal distribution.
#'
#' @param all_cells_to_normalise list of column names that we want to zscore transform
#' @param counts_matrix counts matrix with normally distrubted data
#' @param control_groups (optional) list of column names corresponding to samples we want to use for mean and SD calculations, in case we only want to calculate average expression from a control group.
#'
#' @return Nothing
#'
#' @export
zscore_transform <- function(all_cells_to_normalise, control_groups, counts_matrix) {

	# sanity check given arguments
	if (!any(control_groups %in% all_cells_to_normalise)) {
		print("Error: none of the control group cells are in the all_cells_to_normalise list.")
		return(NA)
	}
	if (!any(all_cells_to_normalise %in% counts_matrix)) {
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
	zscore_counts <- as.data.frame(t(apply(zscore_counts[, c("mean", "stdD", control_groups)], 1, function(x) ((x - x[["mean"]]) / x[["stdD"]]))))
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
