#' @title logR_by_ref_group
#'
#' @description converts raw counts matrix to log fold change values, based on fold change
#' between mean of specified experiemental group and specified reference group.
#'
#' @param raw_counts_matrix a counts matrix with no normalisation
#' @param reference_cells a list of strings corresponding to column names in counts matrix. These are
#' the names of cells we will use as a reference group
#' @param manual_filter optional list of values to pass to the internal filter_and_normalise_scRNA
#' function to filter out poor quality samples before calculating log fold change
#'
#' @return a counts matrix containing log fold change values between specified reference cells and
#' all others.
#'
logR_by_ref_group <- function(raw_counts_matrix, reference_cells, manual_filter = FALSE)
 {

	# convert raw counts to raw FPKM (use same filters as for other analysis)
	fpkm_counts <- barman::filter_and_normalise_scRNA(raw_counts_matrix, manual_filter = manual_filter)

	rownames(fpkm_counts) = fpkm_counts$Geneid

	# remove 0 counts
	fpkm_counts[fpkm_counts == 0] <- NA

	# add column with median expression of reference cells
	fpkm_counts$norm_factor <- rowMeans(as.matrix(fpkm_counts[, !(colnames(fpkm_counts) %in% c("Geneid","Chr","Start","End","Strand","Length"))]), na.rm = T)

	# calculate logR of expression between reference cells, and all other cells
	log_r_dataframe <- t(apply(fpkm_counts, 1, function(x) (log2((x + 1) / (x[["norm_factor"]] + 1)))))
	log_r_dataframe <- subset(log_r_dataframe, select = c(-norm_factor))
	log_r_dataframe <- as.data.frame(log_r_dataframe)
	cell_ids = colnames(log_r_dataframe)

	# readd columns back
	log_r_dataframe$Geneid = rownames(log_r_dataframe)
	merged_logr = merge(log_r_dataframe, raw_counts_matrix[,c("Geneid","Chr","Start","End","Strand","Length")]
	, by = "Geneid")
	merged_logr = merged_logr[,c("Geneid","Chr","Start","End","Strand","Length",cell_ids)]

	return(merged_logr)
}
