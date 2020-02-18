#' @title barman
#'
#' @description runs the G_T_chr_plot function for a whole directory of scCNV segment files and RNA segment
#'
#' @param dna_segments_dir folder of CNV segment files from scCNV or converted from ASCAT with ascat2sccnvseg()
#' @param rna_segments_dir folder of expression segment files, from get_expression_by_segment()
#' @param threads number of threads to run concurrently (defaults to number of avaliable cores minus 1)
#'
#' @return Nothing
#'
#' @export
bulk_G_T_chr_plots = function(dna_segments_dir, rna_segments_dir, threads) {

	# make sure dna_segments_dir isnt empty
	if (length(list.files(dna_segments_dir)) < 1) {
		print("Error: dna_segments_dir argument is an empty folder")
		return(NA)
	}
	# make sure rna_segments_dir isnt empty
	if (length(list.files(rna_segments_dir)) < 1) {
		print("Error: dna_segments_dir argument is an empty folder")
		return(NA)
	}



	# tell user how many RNA to DNA matches were found
	rna_sample_ids = c()
	for (file in list.files(rna_segments_dir)) {
		rna_sample_ids = append(rna_sample_ids, gsub("_from_DNA_segments.tsv", "", gsub("RNAsegments_","",file)), length(rna_sample_ids))
	}
	print(paste0("Generating plots for samples : ", unlist(rna_sample_ids)))

	# Define a function where the heavy lifting happens, that can be easily inserted into a foreach loop to parallelise
	parallel_seg_counter = function(sample_name) {

		# read in DNA segment file
		dna_seg_filename = list.files(dna_segments_dir)[grepl(paste0(".*?",sample_name, ".*"), list.files(dna_segments_dir))]
		dna_seg = read.table(paste0(dna_segments_dir,"/",dna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)

		# read in RNA segment file
		rna_seg_filename = list.files(rna_segments_dir)[grepl(paste0(".*?",sample_name, ".*"), list.files(rna_segments_dir))]
		rna_seg = read.table(paste0(rna_segments_dir,"/",rna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)


		plot_title = paste0(sample_name)


		plt = G_T_chr_plot(cnv_data = dna_seg, exp_data = rna_seg, title = plot_title, save = T)
		return(plt)
	}


	#setup parallel backend to use many processors
	total_cores = parallel::detectCores()

	# if threads option not given then run all detected but 1
	if (missing(threads)) {
		threads = total_cores[1]-1
	}
	cores_cluster <- parallel::makeCluster(threads) #not to overload your computer
	doParallel::registerDoParallel(cores_cluster)

	foreach::foreach(i=1:(length(rna_sample_ids))) %dopar% {
		library(barman)

		plt = parallel_seg_counter(rna_sample_ids[i])
		print(plt)
	}

	parallel::stopCluster(cores_cluster)

}
