#' @title expression_boxplots
#'
#' @description Generates per chromosome boxplots for the mean expression for each gene, grouped by chromosome
#'
#' @param experimental_group list of column names corresponding to samples we want to use to show the mean expression of cells of interest
#' @param control_group  list of column names corresponding to samples we want to use to show the mean expression over control cells
#' @param counts_matrix counts matrix with expression we want to plot, needs to already be normalised
#'
#' @return ggplot2 object
#'
#' @export
expression_boxplots <- function(experimental_group, control_group, counts_matrix) {

  experimental_group <- experimental_group[experimental_group %in% colnames(counts_matrix)]
  control_group <- control_group[control_group %in% colnames(counts_matrix)]


  if (length(experimental_group) > 1) {
  	counts_matrix$Experiemental_mean <- rowMeans(as.matrix(counts_matrix[, experimental_group]), na.rm = T)
  } else {
  	counts_matrix$Experiemental_mean = counts_matrix[, experimental_group]
  }

  counts_matrix$Control_mean <- rowMeans(as.matrix(counts_matrix[, control_group]), na.rm = T)

  counts_matrix <- counts_matrix[, c("Chr", "Experiemental_mean", "Control_mean")]
  counts_matrix <- counts_matrix[counts_matrix$Chr %in% c(1:22, "X", "Y"), ]
  counts_matrix = as.data.frame(counts_matrix)


  melted_counts_matrix <- suppressMessages(reshape2::melt(counts_matrix))

  colnames(melted_counts_matrix)[2] <- "Group"

  # coerce value column into numeric only
  melted_counts_matrix$value <- as.numeric(melted_counts_matrix$value)
  melted_counts_matrix <- melted_counts_matrix[!is.na(melted_counts_matrix$value), ]

  # mark outliers
  is.outlier <- function(x) {
    x < quantile(x, .25) - 1.5 * IQR(x) |
      x > quantile(x, .75) + 1.5 * IQR(x)
  }

  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  melted_counts_matrix %>%
    group_by(Group) %>%
    dplyr::mutate(outlier = is.outlier(value)) %>%
    ungroup() -> melted_counts_matrix


  boxplot_color_palete <- c("#4252A1", "#006A75", "#E1433D") # blue, green, red

  # plot boxes + outlier points with statistical signficance brackets
  p <- ggplot2::ggplot(melted_counts_matrix, ggplot2::aes(x = Group, y = value, fill = Group)) +
  	ggplot2::facet_wrap(~Chr, nrow = 1) +
  	ggplot2::geom_boxplot(outlier.color = NA) +
  	ggplot2::scale_fill_manual(values = boxplot_color_palete) +
  	ggplot2::scale_colour_manual(values = boxplot_color_palete) +
  	ggplot2::xlab("Group of Cells by Chromosome") +
  	ggplot2::ylab("Log Fold Change in Expression") +
  	ggplot2::geom_point(data = melted_counts_matrix[melted_counts_matrix$outlier, ], ggplot2::aes(col = Group), size = 3) +
    ggpubr::geom_signif(
      test = "wilcox.test", comparisons = list(c("Experiemental_mean", "Control_mean")),
      map_signif_level = T,
      vjust = 0,
      textsize = 3,
      size = 0.5,
      step_increase = -0.05
    )

  return(p)
}
