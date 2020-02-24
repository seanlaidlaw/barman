#' @title qc_plots
#'
#' @description plots total_counts total_features_by_counts, mt% and ercc% for given counts matrix
#'
#' @param counts_matrix a raw counts matrix from featureCounts
#'
#' @return ggplot2 object of violin plot
#'
#' @export


qc_plots <- function(counts_matrix) {

	if (length(grep("ENSG[0-9]*", counts_matrix$Geneid)) > 0) {
		mart72.hs <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
		new_gene_ids = biomaRt::getBM(
			attributes=c("ensembl_gene_id","hgnc_symbol"),
			filters="ensembl_gene_id",
			values=counts_matrix$Geneid,
			mart=mart72.hs)

		new_gene_ids$Geneid = new_gene_ids$ensembl_gene_id
		new_gene_ids$ensembl_gene_id = NULL
		new_gene_ids = new_gene_ids[!duplicated(new_gene_ids[,c('Geneid')]),]
		new_gene_ids = new_gene_ids[new_gene_ids$hgnc_symbol != "",]

		new_table = merge(new_gene_ids, counts_matrix, by = "Geneid")
		new_table$Geneid = NULL
		new_table$Geneid = new_table$hgnc_symbol
		new_table$hgnc_symbol = NULL
		counts_matrix = new_table
	}


	trimmed_counts_matrix = counts_matrix
	trimmed_counts_matrix = trimmed_counts_matrix[!duplicated(trimmed_counts_matrix[,c('Geneid')]),]
	trimmed_counts_matrix = trimmed_counts_matrix[trimmed_counts_matrix$Geneid != "",]
	rownames(trimmed_counts_matrix) = trimmed_counts_matrix$Geneid
	trimmed_counts_matrix$Geneid = NULL
	trimmed_counts_matrix$Chr = NULL
	trimmed_counts_matrix$Start = NULL
	trimmed_counts_matrix$End = NULL
	trimmed_counts_matrix$Strand = NULL
	length_list = trimmed_counts_matrix
	length_list = subset(length_list, select = Length)
	trimmed_counts_matrix$Length = NULL


	fc_sce <- SingleCellExperiment::SingleCellExperiment(
		assays = list(counts = as.matrix(trimmed_counts_matrix), logcounts = log2(as.matrix(trimmed_counts_matrix) + 1)),
	)


	seurat_HCC38_BL <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(fc_sce))

	seurat_HCC38_BL[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_HCC38_BL, pattern = "MT-")
	seurat_HCC38_BL[["percent.ercc"]] <- Seurat::PercentageFeatureSet(seurat_HCC38_BL, pattern = "ERCC-")

	plt = Seurat::VlnPlot(seurat_HCC38_BL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ercc"), ncol = 4)
	print(plt)
}
