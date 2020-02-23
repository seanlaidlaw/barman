#' @title barman
#'
#' @description create a segment table of RNA expression, using the same segments as a provided CNV segment file
#'
#' @param cell_id the column name in counts_matrix of the cell we want to get avg expression for each segment
#' @param counts_matrix a featureCounts matrix
#' @param segment_file the CNV segment file from scCNV or converted from ASCAT with ascat2sccnvseg()
#'
#' @return a segmentation matrix of gene expression
#'
#' @export

get_expression_by_segment <- function(cell_id, counts_matrix, segment_file) {

	if (length(counts_matrix[[cell_id]]) == 0) {
		print("Error: counts matrix either doesnt have column:'", cell_id, "' or it doesnt contain any values.")
		return(NA)
	}

	counts_matrix_subset = counts_matrix[,c("Chr", "Start", "End", cell_id)]

	for (numeric_column in c("Start", "End")) {
		if (typeof(counts_matrix[[numeric_column]]) != "integer") {
			print("Error: ",numeric_column," column is not an integer, can not proceed with segmentation. Please verify column only contains integers.")
			return(NA)
		}
	}


	segmented_expression_table = data.frame(matrix(rep(NA,nrow(segment_file)*4),nrow(segment_file),4))
	colnames(segmented_expression_table) = c("Chr", "Start", "SegmentExp", "End")

	for (i in 1:nrow(segment_file)) {
		sgmt_chrom = (segment_file[i,1])
		sgmt_start = (segment_file[i,2])
		sgmt_end = (segment_file[i,5])

		sgmt_avg = mean(counts_matrix_subset[[cell_id]][(counts_matrix_subset$Chr == sgmt_chrom & counts_matrix_subset$Start >= sgmt_start & counts_matrix_subset$Start <= sgmt_end)])

		segmented_expression_table[i,1] = sgmt_chrom
		segmented_expression_table[i,2] = sgmt_start
		segmented_expression_table[i,3] = sgmt_avg
		segmented_expression_table[i,4] = sgmt_end
	}

	# remove NA values from table
	segmented_expression_table = segmented_expression_table[!is.na(segmented_expression_table$SegmentExp),]


	PCFexp = copynumber::pcf(subset(segmented_expression_table, select = c("Chr", "Start", "SegmentExp")), gamma = 15)

	# Set segmented value for each range that pcf() found
	segmented_expression_table$PCFexp = NA

	for (i in 1:nrow(PCFexp)) {
		sgmt_chrom = (PCFexp[i,2])
		sgmt_start = (PCFexp[i,4])
		sgmt_end = (PCFexp[i,5])
		sgmt_cn = (PCFexp[i,7])

		segmented_expression_table$PCFexp[(segmented_expression_table$Chr == sgmt_chrom & segmented_expression_table$Start >= sgmt_start & segmented_expression_table$Start <= sgmt_end)] <- sgmt_cn
	}

	segmented_expression_table = segmented_expression_table[,c("Chr", "Start",  "SegmentExp", "PCFexp", "End")]
	return(segmented_expression_table)
}
