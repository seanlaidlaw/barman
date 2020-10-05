#' @title qc_plots
#'
#' @description plots total_counts total_features_by_counts, mt% and ercc% for given counts matrix
#'
#' @param counts_matrix a raw counts matrix from featureCounts, where feature names correspond to
#' hgnc symbols
#'
#' @return ggplot2 object of violin plot
#'
#' @export


qc_plots <- function(counts_matrix) {

	# if counts matrix is a featureCounts matrix with metadata columns then remove those
	if (any(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% colnames(counts_matrix))) {
		trimmed_counts_matrix = counts_matrix
		trimmed_counts_matrix = trimmed_counts_matrix[!duplicated(trimmed_counts_matrix[,c('Geneid')]),]
		trimmed_counts_matrix = trimmed_counts_matrix[trimmed_counts_matrix$Geneid != "",]
		rownames(trimmed_counts_matrix) = trimmed_counts_matrix$Geneid
		trimmed_counts_matrix$Geneid = NULL
		trimmed_counts_matrix$Chr = NULL
		trimmed_counts_matrix$Start = NULL
		trimmed_counts_matrix$End = NULL
		trimmed_counts_matrix$Strand = NULL
		trimmed_counts_matrix$Length = NULL
	} else {
		# if counts matrix is already gene x cell matrix then use it directly
		trimmed_counts_matrix = counts_matrix
	}


	fc_sce <- SingleCellExperiment::SingleCellExperiment(
		assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix) + 1)),
	)


	seurat_HCC38_BL <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(fc_sce))

	seurat_HCC38_BL[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_HCC38_BL, pattern = "MT-")
	seurat_HCC38_BL[["percent.ercc"]] <- Seurat::PercentageFeatureSet(seurat_HCC38_BL, pattern = "ERCC-")

	plt = Seurat::VlnPlot(seurat_HCC38_BL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ercc"), ncol = 4)
	print(plt)
}
