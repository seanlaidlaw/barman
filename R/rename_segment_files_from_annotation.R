#' @title rename_segment_files_from_annotation
#'
#' @description Renames a folder of segmentation files to remove the bloat from filenames, leaving only the cell id and a ".txt" suffix
#'
#' @param segment_files_dir the folder of files to rename
#' @param annotation_table annotation table containing a column of unique cell IDs that can be used to differentiate each segment file
#' @param dna_cell_id the columnname in the annotation_table of the column containg the unique cell IDs
#'
#' @return Nothing
#'
#' @export

rename_segment_files_from_annotation <- function(segment_files_dir, annotation_table, dna_cell_id) {
	# do a first pass to check regex is specifc to just one file per dna_id and that those files exist
	nb_of_detected_files = 0
	for (dna_id in annotation_table[[dna_cell_id]]) {
		# look for any file in folder containing dna_cell_id and containing ".[a-z]{3}" (e.g. .txt or .tsv or .csv)
		cell_id_file_finder_regex = paste0(".*?",dna_id,".*?\\....")
		segment_files = list.files(segment_files_dir, pattern=cell_id_file_finder_regex)

		if (length(segment_files) > 1) {
			stop(paste0("regex of '", cell_id_file_finder_regex, "' is not specific enough to only find 1 file for this dna_id: '", dna_id, "'."))
		} else if (length(segment_files) == 1) {
			nb_of_detected_files = nb_of_detected_files + 1
		}
	}

	if (nb_of_detected_files == 0) {
		stop(paste0("No files with '", dna_id, "' regex were detected in folder :'", segment_files_dir, "'."))
	}

	# if didnt fail due to non-specific regex or absent files, then rename files to dna_id.txt
	for (dna_id in annotation_table[[dna_cell_id]]) {
		# look for any file in folder containing dna_cell_id and containing ".[a-z]{3}" (e.g. .txt or .tsv or .csv)
		cell_id_file_finder_regex = paste0(".*?",dna_id,".*?\\....")
		segment_files = list.files(segment_files_dir, pattern=cell_id_file_finder_regex)

		if (length(segment_files) == 1) {
			file.rename(from = paste0(segment_files_dir,"/",segment_files), to = paste0(segment_files_dir,"/",dna_id, ".txt"))
		}
	}

}
