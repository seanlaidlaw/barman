#' @title barman
#'
#' @description filters out poor quality cells and features, and converts counts matrix to FPKM
#'
#' @param counts_matrix a raw counts matrix from featureCounts
#' @param annotation_table an annotation table
#'
#' @return Nothing
#'
#' @export

filter_and_normalise_scRNA <- function(counts_matrix, annotation_table=FALSE, return_sce=FALSE) {

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


	if (annotation_table) {
		fc_sce <- SingleCellExperiment::SingleCellExperiment(
			assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix))),
			colData = annotation_table
		)
	} else {
		fc_sce <- SingleCellExperiment::SingleCellExperiment(
			assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix))),
		)
	}



	SingleCellExperiment::isSpike(fc_sce, "ERCC") <- grepl("ERCC-", rownames(fc_sce))
	SingleCellExperiment::isSpike(fc_sce, "MT") <- rownames(fc_sce) %in% c("MT-TF","MT-RNR1","MT-TV","MT-RNR2","MT-TL1","MT-ND1","MT-TI","MT-TQ","MT-TM","MT-ND2","MT-TW","MT-TA","MT-TN","MT-TC","MT-TY","MT-CO1","MT-TS1","MT-TD","MT-CO2","MT-TK","MT-ATP8","MT-ATP6","MT-CO3","MT-TG","MT-ND3","MT-TR","MT-ND4L","MT-ND4","MT-TH","MT-TS2","MT-TL2","MT-ND5","MT-ND6","MT-TE","MT-CYB","MT-TT","MT-TP")

	fc_sce <- scater::calculateQCMetrics(fc_sce,
										 feature_controls = list(
										 	MT = SingleCellExperiment::isSpike(fc_sce, "MT"),
										 	ERCC = SingleCellExperiment::isSpike(fc_sce, "ERCC")
										 )
	)

	# remove non-endogenous genes
	fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "ERCC")]
	fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "MT")]


	fc_sce <- scater::runPCA(
		fc_sce,
		use_coldata = TRUE,
		detect_outliers = TRUE
	)

	fc_sce$use <- !(fc_sce$outlier)

	scater::plotReducedDim(
		fc_sce,
		use_dimred = "PCA_coldata",
		size_by = "total_features_by_counts",
		shape_by = "use",
		colour_by = "outlier"
	)


	keep_feature <- scater::nexprs(
		fc_sce[,SingleCellExperiment::colData(fc_sce)$use],
		byrow = TRUE,
		detection_limit = 1
	) >= 4

	fc_sce <- fc_sce[keep_feature,]
	fc_sce <- fc_sce[,fc_sce$use]


	length_list = merge(SingleCellExperiment::counts(fc_sce), length_list, by = "row.names")
	length_list = subset(length_list, select = Length)
	length_list = as.numeric(length_list$Length)

	scater::fpkm(fc_sce) = scater::calculateFPKM(fc_sce, effective_length = length_list, exprs_values = "counts")

	if (return_sce) {
		return(fc_sce)
	}

	fpkm_counts_matrix = counts_matrix
	rownames(fpkm_counts_matrix) = fpkm_counts_matrix$Geneid
	fpkm_counts_matrix$Geneid = NULL
	fpkm_counts_matrix = merge(fpkm_counts_matrix[,c("Chr","Start","End","Strand","Length")], scater::fpkm(fc_sce), by = "row.names")
	fpkm_counts_matrix$Geneid = fpkm_counts_matrix$Row.names
	fpkm_counts_matrix$Row.names = NULL
	fpkm_counts_matrix = fpkm_counts_matrix[,append("Geneid", colnames(fpkm_counts_matrix)[1:length(colnames(fpkm_counts_matrix))-1])]

	return(fpkm_counts_matrix)
}
