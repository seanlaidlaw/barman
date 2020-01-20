#' @title barman
#'
#' @description Converts ASCAT summary segmentation format to scCNV format, for plotting ASCAT data with G_T_chr_plots
#'
#' @param ascat_summary a segment matrix loaded in from an ASCAT.summary.csv file
#'
#' @return a segment matrix in scCNV format
#'
#' @export
ascat2sccnvseg <- function(ascat_summary) {

	if (length(colnames(tolly_dna_segment_file)) != 8 ) {
		print(paste0("ERROR: expected 8 columns but only got: ", length(colnames(tolly_dna_segment_file)), ". Are you sure you provided an ASCAT summary table?"))
		return(NA)
	}

	# from ascat confluence (https://confluence.sanger.ac.uk/display/IT/ASCAT+NGS) we want col 7 as CN
	colnames(ascat_summary) = c("Segment number (just an incremental counter)","Chromosome", "Segment start (0 based)", "Segment end (1 based)", "Normal sample total copy number", "Normal sample minor allele copy number", "Tumour sample total copy number", "Tumour sample minor allele copy number")
	ascat_summary$`Segment number (just an incremental counter)` = NULL
	ascat_summary$`Normal sample total copy number` = NULL
	ascat_summary$`Normal sample minor allele copy number` = NULL
	ascat_summary$`Tumour sample minor allele copy number` = NULL

	# make colnames the same
	colnames(ascat_summary) = c("Chr", "Pos","end","rawCN")

	# add mean segmentation
	sgmts <- copynumber::pcf(subset(ascat_summary, select = c("Chr", "Pos","rawCN")), gamma = 15)
	sgmts <- subset(sgmts, select = c(chrom, start.pos, end.pos, mean))
	colnames(sgmts) = c("Chr", "Pos","end","segmentedCN")


	# Set segmented value for each range that pcf() found
	ascat_summary$segmentedCN = NA

	for (i in 1:nrow(sgmts)) {
		sgmt_chrom = (sgmts[i,1])
		sgmt_start = (sgmts[i,2])
		sgmt_end = (sgmts[i,3])
		sgmt_cn = (sgmts[i,4])

		ascat_summary$segmentedCN[(ascat_summary$Chr == sgmt_chrom & ascat_summary$Pos >= sgmt_start & ascat_summary$Pos <= sgmt_end)] <- sgmt_cn
	}

	ascat_summary = ascat_summary[,c("Chr", "Pos",  "rawCN", "segmentedCN", "end")]
	return(ascat_summary)
}