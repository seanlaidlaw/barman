#' @title barman
#'
#' @description filters out poor quality cells and features, and converts counts matrix to FPKM
#'
#' @param counts_matrix a raw counts matrix from featureCounts
#' @param annotation_table an annotation table
#' @param return_sce boolean to return SCE object instead of matrix
#' @param manual_filter vector of 4 values, manually specifying upper cutoffs to apply for total counts. e.g. c(20000,6000000,20,50) for upper limits of total counts, total_features_by_counts, pct mt, and pct ercc respectively. Also accepts lists instead of ints, where the first element of list is lower cutoff and second is upper cutoff.
#' @param filter_only boolean to return only filtered (i.e. non-fpkm normalized) matrix
#'
#' @return normalized counts matrix
#'
#' @export
filter_and_normalise_scRNA <- function(counts_matrix, annotation_table=FALSE, return_sce=FALSE, manual_filter=FALSE, filter_only=FALSE) {

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


	if (missing(annotation_table)) {
		fc_sce <- SingleCellExperiment::SingleCellExperiment(
			assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix) + 1)),
		)
	} else {
		fc_sce <- SingleCellExperiment::SingleCellExperiment(
			assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix) + 1)),
			colData = annotation_table
		)
	}

	mt_genes = c("MT-TF","MT-RNR1","MT-TV","MT-RNR2","MT-TL1","MT-ND1","MT-TI","MT-TQ","MT-TM","MT-ND2","MT-TW",
				 "MT-TA","MT-TN","MT-TC","MT-TY","MT-CO1","MT-TS1","MT-TD","MT-CO2","MT-TK","MT-ATP8","MT-ATP6",
				 "MT-CO3","MT-TG","MT-ND3","MT-TR","MT-ND4L","MT-ND4","MT-TH","MT-TS2","MT-TL2","MT-ND5","MT-ND6",
				 "MT-TE","MT-CYB","MT-TT","MT-TP", "ENSG00000210049","ENSG00000211459","ENSG00000210077","ENSG00000210082",
				 "ENSG00000209082","ENSG00000198888","ENSG00000210100","ENSG00000210107","ENSG00000210112","ENSG00000198763",
				 "ENSG00000210117","ENSG00000210127","ENSG00000210135","ENSG00000210140","ENSG00000210144","ENSG00000198804",
				 "ENSG00000210151","ENSG00000210154","ENSG00000198712","ENSG00000210156","ENSG00000228253","ENSG00000198899",
				 "ENSG00000198938","ENSG00000210164","ENSG00000198840","ENSG00000210174","ENSG00000212907","ENSG00000198886",
				 "ENSG00000210176","ENSG00000210184","ENSG00000210191","ENSG00000198786","ENSG00000198695","ENSG00000210194",
				 "ENSG00000198727","ENSG00000210195","ENSG00000210196")


	SingleCellExperiment::isSpike(fc_sce, "MT") = grepl("ERCC-", rownames(fc_sce))
	SingleCellExperiment::isSpike(fc_sce, "ERCC") = rownames(fc_sce) %in% mt_genes

	fc_sce <- scater::calculateQCMetrics(fc_sce,
		feature_controls = list(
		MT = SingleCellExperiment::isSpike(fc_sce, "MT"),
		ERCC = SingleCellExperiment::isSpike(fc_sce, "ERCC")
		)
	)


	# remove non-endogenous genes
	fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "ERCC")]
	fc_sce = fc_sce[!SingleCellExperiment::isSpike(fc_sce, "MT")]

	if (typeof(manual_filter) != "list") {
		fc_sce <- scater::runPCA(
			fc_sce,
			use_coldata = TRUE,
			detect_outliers = TRUE
		)

		fc_sce$use <- !(fc_sce$outlier)

		qcplt4 <- scater::plotReducedDim(
			fc_sce,
			use_dimred = "PCA_coldata",
			size_by = "total_features_by_counts",
			colour_by = "use")
		print(qcplt4)

	} else {
		# apply manual cutoffs


		if (length(manual_filter[[1]]) > 1) {
			total_features_by_counts = fc_sce$total_features_by_counts >= manual_filter[[1]][1] & fc_sce$total_features_by_counts <= manual_filter[[1]][2]
		} else {
			total_features_by_counts = fc_sce$total_features_by_counts <= manual_filter[[1]]
		}
		if (length(manual_filter[[2]]) > 1) {
			total_counts = fc_sce$total_counts >= manual_filter[[2]][1] & fc_sce$total_counts <= manual_filter[[2]][2]
		} else {
			total_counts = fc_sce$total_counts <= manual_filter[[2]]
		}
		pct_counts_MT = fc_sce$pct_counts_MT <= manual_filter[[3]]
		pct_counts_ERCC = fc_sce$pct_counts_ERCC <= manual_filter[[4]]

		fc_sce$use <-(pct_counts_ERCC & pct_counts_MT & total_counts & total_features_by_counts )

		fc_sce <- scater::runPCA(
			fc_sce,
			use_coldata = TRUE,
		)

		qcplt4 <- scater::plotReducedDim(
			fc_sce,
			use_dimred = "PCA_coldata",
			size_by = "total_features_by_counts",
			colour_by = "use")
		print(qcplt4)
	}


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

	if (filter_only) {
		if (return_sce) {
			return(fc_sce)
		} else {
			rownames(counts_matrix) = counts_matrix$Geneid
			counts_matrix$Geneid = NULL
			counts_matrix = merge(counts_matrix[,c("Chr","Start","End","Strand","Length")], SingleCellExperiment::counts(fc_sce), by = "row.names")
			counts_matrix$Geneid = counts_matrix$Row.names
			counts_matrix$Row.names = NULL
			counts_matrix = counts_matrix[,append("Geneid", colnames(counts_matrix)[1:length(colnames(counts_matrix))-1])]

			return(counts_matrix)
		}
	}


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
