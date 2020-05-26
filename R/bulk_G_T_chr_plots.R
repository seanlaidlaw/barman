#' @title bulk_G_T_chr_plots
#'
#' @description runs the G_T_chr_plot function for a whole directory of scCNV segment files and RNA segment
#'
#' @param dna_segments_dir folder of CNV segment files from scCNV or converted from ASCAT with ascat2sccnvseg()
#' @param rna_segments_dir folder of expression segment files, from get_expression_by_segment()
#' @param output_dir directory where the plots will be saved
#' @param threads number of threads to run concurrently (defaults to number of avaliable cores minus 1)
#'
#' @return Nothing
#'
#' @export
bulk_G_T_chr_plots <- function(dna_segments_dir, rna_segments_dir, output_dir="./", chromosomes=NA, threads) {
  if (!missing(dna_segments_dir)) {
    # make sure dna_segments_dir isnt empty
    if (length(list.files(dna_segments_dir)) < 1) {
      stop("dna_segments_dir argument is an empty folder")
    }
  }

  if (!missing(rna_segments_dir)) {
    # make sure rna_segments_dir isnt empty
    if (length(list.files(rna_segments_dir)) < 1) {
      stop("rna_segments_dir argument is an empty folder")
    }
  }



  # get sample_ids from folders provided
	if ((!missing(dna_segments_dir)) & (!missing(rna_segments_dir))) {
		# if we have both DNA and RNA segments, the RNA are more aggressivly filtered so use that list
		segments_dir = rna_segments_dir
	} else if (!missing(rna_segments_dir)) {
		segments_dir = rna_segments_dir
	} else if (!missing(dna_segments_dir)) {
		segments_dir = dna_segments_dir
	}

  sample_ids <- c()
  for (file in list.files(segments_dir)) {
  	sample_ids <- append(sample_ids, gsub("(_from_DNA_segments.tsv)|(-PCF-kmin.*copynumber.refLOCUS.NA.txt)", "", gsub("(SEGMENTSlogR.GCcorrected-M30-)|(RNAsegments_)", "", file)), length(sample_ids))
  }
  message(paste0("Generating plots for samples : ", paste(names(sample_ids),sample_ids,sep=" ",collapse=";")))




	# Define a function where the heavy lifting happens, that can be easily inserted into a foreach loop to parallelise
	parallel_seg_counter <- function(sample_name, run_DNA=TRUE, run_RNA=TRUE) {

		plot_title <- paste0(sample_name)

		# with both we can plot both
		if (run_DNA & run_RNA) {
			dna_seg_filename <- list.files(dna_segments_dir)[grepl(paste0(".*?", sample_name, ".*"), list.files(dna_segments_dir))]
			dna_seg <- read.table(paste0(dna_segments_dir, "/", dna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)
			rna_seg_filename <- list.files(rna_segments_dir)[grepl(paste0(".*?", sample_name, ".*"), list.files(rna_segments_dir))]
			rna_seg <- read.table(paste0(rna_segments_dir, "/", rna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)

			plt <- G_T_chr_plot(cnv_data = dna_seg, exp_data = rna_seg,output_dir=output_dir, title = plot_title, save = T)

		} else if (run_DNA) {
			# read in DNA segment file
			dna_seg_filename <- list.files(dna_segments_dir)[grepl(paste0(".*?", sample_name, ".*"), list.files(dna_segments_dir))]
			dna_seg <- read.table(paste0(dna_segments_dir, "/", dna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)
			plt <- G_T_chr_plot(cnv_data = dna_seg,output_dir=output_dir, title = plot_title, save = T)
		} else if (run_RNA) {
			# read in RNA segment file
			rna_seg_filename <- list.files(rna_segments_dir)[grepl(paste0(".*?", sample_name, ".*"), list.files(rna_segments_dir))]
			rna_seg <- read.table(paste0(rna_segments_dir, "/", rna_seg_filename), header = T, sep = "\t", stringsAsFactors = F)
			plt <- G_T_chr_plot(exp_data = rna_seg,output_dir=output_dir, title = plot_title, save = T)
		}

	 return(plt)
  }


  # setup parallel backend to use many processors
  total_cores <- parallel::detectCores()

  # if threads option not given then run all detected but 1
  if (missing(threads)) {
    threads <- total_cores[1] - 1
  }
  cores_cluster <- parallel::makeCluster(threads) # not to overload your computer
  doParallel::registerDoParallel(cores_cluster)


  # if RNA or DNA segment folder not provided then generate plots from the only provided segment folder
  run_DNA = TRUE
  run_RNA = TRUE

  if (missing(rna_segments_dir)) {
  	run_RNA = FALSE
  }
  if (missing(dna_segments_dir)) {
  	run_DNA = FALSE
  }


  suppressPackageStartupMessages(library(doParallel))
  foreach::foreach(i = 1:(length(sample_ids))) %dopar% {
	suppressPackageStartupMessages(library(barman))

	parallel_seg_counter(sample_ids[i],run_DNA=run_DNA, run_RNA=run_RNA)
  }

  parallel::stopCluster(cores_cluster)
}
