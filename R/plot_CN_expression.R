#' @title plot_CN_expression
#'
#' @description Plots a boxplot of average expression for each DNA copy number for each gene between
#' specified groups of cells
#'
#' @param CN_expression_dir directory of tables correlating expression to copy number for each
#' @param annotation_table an annotation table where rownames are the same cell ids as the prefix in filenames in
#' CN_expression_dir, and where at least one column of matrix determines cell group / type
#' @param annotation_column name of column in annotation_table corresponding to cell group / type
#'
#' @return a ggplot object of boxplot of average expression per copy number (optionally per cell
#' group)
#'
#' @export
plot_CN_expression = function(CN_expression_dir, annotation_table=NA, annotation_column=NA) {
	cumulative_cn_expr =  data.frame(Geneid=character(), Cell=character(), Expression=numeric(), rawCN=numeric(), pcfCN=numeric())

	for (cn_expr_file in list.files(CN_expression_dir)) {
			cn_expr_name = sub("_expr_CN.tsv$", "", cn_expr_file)
			cn_expr_path = paste0(CN_expression_dir, "/", cn_expr_file)
			cn_expr <- read.table(cn_expr_path, header = T, sep = "\t", stringsAsFactors = F)

			# add cell name as a column
			cn_expr$Cell = cn_expr_name

			# make table same column order each time so rbind works
			cn_expr = cn_expr[,c("Geneid", "Cell", "Expression", "rawCN", "pcfCN")]

			# bind to previously read in tables
			cumulative_cn_expr = rbind(cumulative_cn_expr, cn_expr)
	}


	cumulative_cn_expr = cumulative_cn_expr[!is.na(cumulative_cn_expr$Expression),]
	cumulative_cn_expr$Cell = as.factor(cumulative_cn_expr$Cell)
	cumulative_cn_expr$rawCN = as.factor(cumulative_cn_expr$rawCN)

	boxplot_colors = c("#354f6a","#e28140","#b13f46")

	if ((!is.na(annotation_table)) & (!is.na(annotation_column))) {
		annotation_table$Cell = rownames(annotation_table)
		cumulative_cn_expr = merge(cumulative_cn_expr, annotation_table[,c("Cell", annotation_column)], by = "Cell")
		colnames(cumulative_cn_expr)[length(colnames(cumulative_cn_expr))] = "Group"
		cumulative_cn_expr$Group = as.factor(cumulative_cn_expr$Group)

		plt = ggplot(cumulative_cn_expr[cumulative_cn_expr$rawCN %in% c(1,2,3),], aes(x = Group, y = Expression, fill = rawCN)) + geom_boxplot() + scale_fill_manual(values=boxplot_colors)
	} else {
		plt = ggplot(cumulative_cn_expr[cumulative_cn_expr$rawCN %in% c(1,2,3),], aes(x = rawCN, y = Expression, fill = rawCN)) + geom_boxplot() + scale_fill_manual(values=boxplot_colors)
	}


	return(plt)
}
