#' @title filter_and_normalise_scRNA__scran
#'
#' @description filters out poor quality cells and features, and converts counts matrix to scran normalized matrix
#' Expects ERCC features to be named 'ERCC-*' and mitochondrial genes to be named 'MT-*'
#'
#' @param counts_matrix a raw counts matrix from featureCounts
#' @param output_dir output directory to save plots
#' @param manual_filter vector of 4 values, manually specifying upper cutoffs to apply for total counts. e.g. c(20000,6000000,20,50) for upper limits of total counts, total_features_by_counts, pct mt, and pct ercc respectively. Also accepts lists instead of ints, where the first element of list is lower cutoff and second is upper cutoff.
#' @param filter_only boolean to return only filtered (i.e. non-fpkm normalized) matrix
#'
#' @return scran normalized counts matrix (returns 'logcounts(sce)')
#'
#' @export
filter_and_normalise_scRNA__scran <- function(counts_matrix, output_dir="./", manual_filter=FALSE, filter_only=FALSE) {

	# remove metadata (length, strand,etc.) columns from counts matrix to obtain a gene x cell matrix
	trimmed_counts_matrix = counts_matrix
	rownames(trimmed_counts_matrix) = trimmed_counts_matrix$Geneid
	trimmed_counts_matrix$Geneid = NULL
	trimmed_counts_matrix$Chr = NULL
	trimmed_counts_matrix$Start = NULL
	trimmed_counts_matrix$End = NULL
	trimmed_counts_matrix$Strand = NULL
	length_list = trimmed_counts_matrix
	length_list = subset(length_list, select = Length)
	trimmed_counts_matrix$Length = NULL

	# create a SingleCellExperiment object from matrix
	fc_sce <- SingleCellExperiment::SingleCellExperiment(
		assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix)+1)),
	)

	# define non-endogenous genes, as they need separate normalisations
	SingleCellExperiment::isSpike(fc_sce, "ERCC") = grepl("ERCC-", rownames(fc_sce))
	SingleCellExperiment::isSpike(fc_sce, "MT") =  grepl("MT-", rownames(fc_sce))

	# calculate QC statistics
	fc_sce <- scater::calculateQCMetrics(fc_sce,
		feature_controls = list(
		MT = SingleCellExperiment::isSpike(fc_sce, "MT"),
		ERCC = SingleCellExperiment::isSpike(fc_sce, "ERCC")
		)
	)



	# remove genes not expressed in any cell
	keep_feature <- rowSums(SingleCellExperiment::counts(fc_sce) > 0) > 0
	fc_sce <- fc_sce[keep_feature, ]


	if (typeof(manual_filter) != "list") {
		fc_sce <- scater::runPCA(
			fc_sce,
			use_coldata = TRUE,
			detect_outliers = TRUE
		)

		fc_sce$use <- !(fc_sce$outlier)

		qc_pca_plot <- scater::plotPCA(
			fc_sce,
			colour_by = "use")

		ggplot2::ggsave(paste0(output_dir,"/PCA_filter.png"),
				plot = qc_pca_plot,
				width = 7,
				height = 7,
				device="png",
				dpi="retina",
				units="in")

	} else {
		# apply manual cutoffs

		# if two numbers given then use the first as a lower and second as a higher cutoff
		if (length(manual_filter[[1]]) > 1) {
			total_features_by_counts = as.numeric(fc_sce$total_features_by_counts) >= as.numeric(manual_filter[[1]][1]) & fc_sce$total_features_by_counts <= as.numeric(manual_filter[[1]][2])
		} else {
			total_features_by_counts = as.numeric(fc_sce$total_features_by_counts) <= as.numeric(manual_filter[[1]])
		}

		if (length(manual_filter[[2]]) > 1) {
			total_counts = as.numeric(fc_sce$total_counts) >= as.numeric(manual_filter[[2]][1]) & fc_sce$total_counts <= as.numeric(manual_filter[[2]][2])
		} else {
			total_counts = as.numeric(fc_sce$total_counts) <= as.numeric(manual_filter[[2]])
		}

		if (length(manual_filter[[3]]) > 1) {
			pct_counts_MT = as.numeric(fc_sce$pct_counts_MT) >= as.numeric(manual_filter[[3]][1]) & fc_sce$pct_counts_MT <= as.numeric(manual_filter[[3]][2])
		} else {
			pct_counts_MT = as.numeric(fc_sce$pct_counts_MT) <= as.numeric(manual_filter[[3]])
		}

		if (length(manual_filter[[4]]) > 1) {
			pct_counts_ERCC = as.numeric(fc_sce$pct_counts_ERCC) >= as.numeric(manual_filter[[4]][1]) & fc_sce$pct_counts_ERCC <= as.numeric(manual_filter[[4]][2])
		} else {
			pct_counts_ERCC = as.numeric(fc_sce$pct_counts_ERCC) <= as.numeric(manual_filter[[4]])
		}

		fc_sce$use <-(pct_counts_ERCC & pct_counts_MT & total_counts & total_features_by_counts)

		fc_sce <- scater::runPCA(
			fc_sce,
			use_coldata = TRUE,
		)

		qc_pca_plot <- scater::plotPCA(
			fc_sce,
			colour_by = "use")

		ggplot2::ggsave(paste0(output_dir,"/PCA_filter.png"),
				plot = qc_pca_plot,
				width = 7,
				height = 7,
				device="png",
				dpi="retina",
				units="in")
	}

	# remove features that have less than 1 transcript in 2 or less cells
	keep_feature <- scater::nexprs(
		fc_sce[,SingleCellExperiment::colData(fc_sce)$use],
		byrow = TRUE,
		detection_limit = 1
	) >= 2

	# subset object to apply the filters we defined
	fc_sce <- fc_sce[keep_feature,]
	fc_sce <- fc_sce[,fc_sce$use]


	# calculate length for remaining features (required if we change normalisation to FPKM)
	length_list = merge(SingleCellExperiment::counts(fc_sce), length_list, by = "row.names")
	length_list = subset(length_list, select = Length)
	length_list = as.numeric(length_list$Length)

	if (filter_only) {
		# remove non-endogenous genes
		fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "ERCC")]
		fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "MT")]

		rownames(counts_matrix) = counts_matrix$Geneid
		counts_matrix$Geneid = NULL
		counts_matrix = merge(counts_matrix[,c("Chr","Start","End","Strand","Length")], SingleCellExperiment::counts(fc_sce), by = "row.names")
		counts_matrix$Geneid = counts_matrix$Row.names
		counts_matrix$Row.names = NULL
		counts_matrix = counts_matrix[,append("Geneid", colnames(counts_matrix)[1:length(colnames(counts_matrix))-1])]

		return(counts_matrix)
	}


	# scran Normalization
	suppressPackageStartupMessages(library(scran))

	# computes size factors that are used to scale the counts in each cell.
	sce= fc_sce
	rm(fc_sce)

	sce <- computeSumFactors(sce)
	summary(sizeFactors(sce))

	# we want to use the deconvolution size factors for the endogenous genes, and the
	# spike-in-based size factors for the spike-in transcripts so we use 'general.use=F'
	sce <- computeSpikeFactors(sce, general.use=FALSE)

	# apply normalisation
	sce <- normalize(sce)
	fc_sce = sce
	rm(sce)

	# remove non-endogenous genes
	fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "ERCC")]
	fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "MT")]
	norm_matrix = SingleCellExperiment::logcounts(fc_sce)


	fpkm_counts_matrix = counts_matrix
	rownames(fpkm_counts_matrix) = fpkm_counts_matrix$Geneid
	fpkm_counts_matrix$Geneid = NULL

	fpkm_counts_matrix = merge(fpkm_counts_matrix[,c("Chr","Start","End","Strand","Length")], norm_matrix, by = "row.names")
	fpkm_counts_matrix$Geneid = fpkm_counts_matrix$Row.names
	fpkm_counts_matrix$Row.names = NULL
	fpkm_counts_matrix = fpkm_counts_matrix[,append("Geneid", colnames(fpkm_counts_matrix)[1:length(colnames(fpkm_counts_matrix))-1])]

	return(fpkm_counts_matrix)
}
