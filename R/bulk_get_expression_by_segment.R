#' @title bulk_get_expression_by_segment
#'
#' @description runs the get_expression_by_segment function for a whole directory of scCNV segment files
#'
#' @param counts_matrix a featureCounts matrix
#' @param dna_segments_dir folder of CNV segment files from scCNV or converted from ASCAT with ascat2sccnvseg()
#' @param output_dir an output directory to store RNAsegment files
#' @param threads number of threads to run concurrently (defaults to number of avaliable cores minus 1)
#'
#' @return Nothing
#'
#' @export
bulk_get_expression_by_segment = function(counts_matrix, dna_segments_dir, output_dir, threads) {

	# get list of samplenames from counts_matrix header
	cell_sample_names = colnames(counts_matrix)[!(colnames(counts_matrix) %in% c("Geneid","Chr","Start","End","Strand","Length"))]

	# make sure dna_segments_dir isnt empty
	if (length(list.files(dna_segments_dir)) < 1) {
		print("Error: dna_segments_dir argument is an empty folder")
		return(NA)
	}

	# make sure some files in dna_segments_dir match colnames of counts matrix
	anyfile = F
	for (i in 1:length(cell_sample_names)) {
		if (any(grepl(paste0(".*?",cell_sample_names[i], ".*"), list.files(dna_segments_dir)))) {
			anyfile = T
		}
	}
	if (!anyfile) {
		print("Error: No match was found between filenames in dna_segments_dir and column names in counts_matrix, verify that cell_ids are the same in column names of counts matrix and in DNA segment files.")
		return(NA)
	}


	# create output directory
	dir.create(file.path(output_dir), showWarnings = FALSE)

	# tell user how many RNA to DNA matches were found
	present_in_dna = c()
	for (sample_name in cell_sample_names) {
		if (any(grepl(paste0(".*?",sample_name, ".*"), list.files(dna_segments_dir)))) {
			present_in_dna = append(present_in_dna, values = sample_name, after = length(present_in_dna))
		}
	}
	if (length(present_in_dna) == length(cell_sample_names)) {
		print("Found all RNA names in DNA segments folder")
	} else {
		print(paste0("Found ", length(present_in_dna), " DNA samples corresponding to ", length(cell_sample_names), " RNA names"))
	}


	# Define a function where the heavy lifting happens, that can be easily inserted into a foreach loop to parallelise
	parallel_seg_counter = function(sample_name) {
		print(sample_name)

		# read in DNA segment file
		dna_seg_filename = list.files(dna_segments_dir)[grepl(paste0(".*?",sample_name, ".*"), list.files(dna_segments_dir))]
		dna_seg = read.table(paste0(dna_segments_dir,"/",dna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)

		# get RNA expression for each segment of the segmented DNA
		rna_seg = get_expression_by_segment(cell_id = sample_name, counts_matrix = counts_matrix, segment_file = dna_seg)

		# write to output folder for further analysis later
		write.table(rna_seg, paste0(output_dir, "/", "RNAsegments_", sample_name, "_from_DNA_segments.tsv"), sep = "\t", row.names = F, col.names = T)
	}

	#setup parallel backend to use many processors
	total_cores = parallel::detectCores()

	# if threads option not given then run all detected but 1
	if (missing(threads)) {
		threads = total_cores[1]-1
	}
	cores_cluster <- parallel::makeCluster(threads) #not to overload your computer
	doParallel::registerDoParallel(cores_cluster)

	library(doParallel)
	foreach::foreach(i=1:(length(present_in_dna))) %dopar% {
		library(barman)

		parallel_seg_counter(present_in_dna[i])
	}

	parallel::stopCluster(cores_cluster)

}
