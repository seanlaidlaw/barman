#' @title barman
#'
#' @description Convert FPKM normalized object SCE to Seurat Object
#'
#' @param counts_matrix a raw counts matrix from featureCounts
#' @param annotation_table an annotation table with at least one of the columns containing the colnames of the counts matrix
#' @param SCE_obj a SingleCellExperiment object created by running filter_and_normalise_scRNA() with the return_sce option set to TRUE
#'
#' @return a normalized Seurat object
#'
#' @export
filtered_SCE_to_Seurat <- function(SCE_obj, counts_matrix, annotation_table){


	if (!("fpkm" %in% names(SCE_obj@assays$data))) {
		print("Error: passed SingleCellExperiment object does not have normalized data, please run filter_and_normalise_scRNA() first")
		return(NA)
	}

	#make annotation as Seurat expects
	annotation_table_seurat = annotation_table
	sample_colnames = colnames(counts_matrix)[!colnames(counts_matrix) %in% c("Geneid","Chr","Start","End","Strand","Length")]
	for (annotation_colname in colnames(annotation_table_seurat)) {
	    if (all(sample_colnames %in% annotation_table_seurat[[annotation_colname]])) {
	        rna_cell_id = annotation_colname
	    }
	}
	rownames(annotation_table_seurat) = annotation_table_seurat[[rna_cell_id]]
	annotation_table_seurat[[rna_cell_id]] = NULL

	Seurat_obj <- CreateSeuratObject(counts = fpkm(SCE_obj), meta.data =  annotation_table_seurat)
	Seurat_obj <- FindVariableFeatures(object = Seurat_obj,dispersion.function = LogVMR)
	Seurat_obj <- ScaleData(object = Seurat_obj, display.progress=F)
	Seurat_obj = RunPCA(object = Seurat_obj,  npcs = 30, verbose = FALSE)
	Seurat_obj = RunUMAP(object = Seurat_obj,  npcs = 30, verbose = FALSE, dims = c(1,2))

	return(Seurat_obj)
}
