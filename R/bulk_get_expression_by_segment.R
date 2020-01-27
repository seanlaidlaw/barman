#' @title barman
#'
#' @description runs the get_expression_by_segment function for a whole directory of scCNV segment files
#'
#' @param counts_matrix a featureCounts matrix
#' @param dna_segments_dir folder of CNV segment files from scCNV or converted from ASCAT with ascat2sccnvseg()
#' @param output_dir an output directory to store RNAsegment files
#'
#' @return Nothing
#'
#' @export
bulk_get_expression_by_segment = function(counts_matrix, dna_segments_dir, output_dir) {

	# get list of samplenames from counts_matrix header
	cell_sample_names = colnames(counts_matrix)[!(colnames(counts_matrix) %in% c("Geneid","Chr","Start","End","Strand","Length"))]

	# make sure dna_segments_dir isnt empty
	if (length(list.files(dna_segments_dir)) < 1) {
		print("Error: dna_segments_dir argument is an empty folder")
		return(NA)
	}

	# make sure some files in dna_segments_dir match colnames of counts matrix
	anyfile = F
	for (i in length(cell_sample_names)) {
		if (any(grepl(paste0(".*?",cell_sample_names[i], ".*"), list.files(dna_segments_dir)))) {
			anyfile = T
		}
	}
	if (!anyfile) {
		print("Error: No match was found between filenames in dna_segments_dir and column names in counts_matrix")
		return(NA)
	}


	# create output directory
	dir.create(file.path(output_dir), showWarnings = FALSE)

	# for each sample_name find a dna_segment file that includes sample_name in filename
	for (sample_name in cell_sample_names) {
		if (any(grepl(paste0(".*?",sample_name, ".*"), list.files(dna_segments_dir)))) {

			dna_seg_filename = list.files(dna_segments_dir)[grepl(paste0(".*?",sample_name, ".*"), list.files(dna_segments_dir))]
			dna_seg = read.table(paste0(dna_segments_dir,"/",dna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)

			# get RNA expression for each segment of the segmented DNA
			rna_seg = get_expression_by_segment(cell_id = sample_name, counts_matrix = counts_matrix, segment_file = dna_seg)

			# write to output folder for further analysis later
			write.table(rna_seg, paste0(output_dir, "/RNAsegments_", sample_name, "_from_DNA_segments.tsv"), sep = "\t", row.names = F, col.names = T)
		}
	}
}
