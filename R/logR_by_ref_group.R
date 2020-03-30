#' @title logR_by_ref_group
#'
#' @description converts raw counts matrix to log fold change values, based on fold change
#' between mean of specified experiemental group and specified reference group.
#'
#' @param fpkm_counts a normalized counts matrix
#' @param reference_cells a list of strings corresponding to column names in counts matrix. These are
#' the names of cells we will use as a reference group
#'
#' @return a counts matrix containing log fold change values between specified reference cells and
#' all others.
#'
#' @export

logR_by_ref_group <- function(fpkm_counts, reference_cells) {

	rownames(fpkm_counts) = fpkm_counts$Geneid

	trimmed_fpkm_counts <- as.matrix(fpkm_counts[, !(colnames(fpkm_counts) %in% c("Geneid","Chr","Start","End","Strand","Length"))])

	# remove 0 counts
	trimmed_fpkm_counts[trimmed_fpkm_counts == 0] <- NA

	# add column with median expression of reference cells
	trimmed_fpkm_counts = as.data.frame(trimmed_fpkm_counts)
	trimmed_fpkm_counts$norm_factor = rowMeans(trimmed_fpkm_counts[,colnames(trimmed_fpkm_counts) %in% reference_cells], na.rm = T)


	# calculate logR of expression between reference cells, and all other cells
	log_r_dataframe <- t(apply(trimmed_fpkm_counts, 1, function(x) (log2((x + 1) / (as.numeric(x[["norm_factor"]]) + 1)))))
	log_r_dataframe = as.data.frame(log_r_dataframe)

	log_r_dataframe <- subset(log_r_dataframe, select = c(-norm_factor))
	log_r_dataframe <- as.data.frame(log_r_dataframe)

	# read columns back
	cell_ids = colnames(log_r_dataframe)
	log_r_dataframe$Geneid = rownames(log_r_dataframe)
	merged_logr = merge(log_r_dataframe, raw_counts_matrix[,c("Geneid","Chr","Start","End","Strand","Length")]
	, by = "Geneid")
	merged_logr = merged_logr[,c("Geneid","Chr","Start","End","Strand","Length",cell_ids)]

	return(merged_logr)
}
