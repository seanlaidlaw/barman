#' @title G_T_chr_plot
#'
#' @description Generates Genome CNV and Transcriptome Expression plots by chromsome from previously generated segmentation files.
#'
#' @param cnv_data a segmentation file produced by scCNV. Can also use an ASCAT.summary.csv if converted with the ascat2sccnvseg() function.
#' @param exp_data a segmentation file generated by calling the get_expression_by_segment() function on a featurecounts matrix and on a DNA segmentation file
#' @param title Plot title
#' @param output_dir output directory to store saved plots
#' @param save boolean to specify if the plot should be generated or saved as PNG. If set to TRUE, the plot will be saved as the specified title.
#' @param chromosomes optional list of chromosomes to plot, defaults to all chromosomes
#'
#' @return a plot object or nothing if save==TRUE
#'
#' @export
G_T_chr_plot <- function(cnv_data, exp_data, chromosomes=NA, output_dir="./", title, save) {

	if (missing(save)) {
		save = FALSE
	}

	if (missing(title)) {
		title = ""
	}

	if (is.na(chromosomes)) {
		# if no chromosome argument provided, use all
		chromosomes = paste0("chr",unique(cnv_data$Chr))
	} else {
		chromosomes = paste0("chr",unique(chromosomes))
	}


	if (!dir.exists(output_dir)) {
		dir.create(output_dir)
	}


	cnv_data$Chr <- paste0("chr", cnv_data$Chr)
	cnv_data_upper_cutoff <- round(mean(cnv_data$segmentedCN) + sd(cnv_data$segmentedCN) * 3)
	cnv_data_upper_cutoff = 4
	cnv_data <- cnv_data[cnv_data$rawCN <= cnv_data_upper_cutoff, ]

	if (!missing(exp_data)) {
	exp_data$Chr <- paste0("chr", exp_data$Chr)
	exp_data_upper_cutoff <- mean(exp_data$PCFexp) + signif(sd(exp_data$SegmentExp) * 3, digits = 2)
	exp_data <- exp_data[exp_data$SegmentExp <= exp_data_upper_cutoff, ]
	}

	dna_r0 <- 0
	dna_r1 <- 0.9
	rna_r0 <- 0
	rna_r1 <- 0.9

	if (save) {
		png(paste0(output_dir,"/", title, ".png"), width = 5120, height = 1600, units = "px", res = 300)
	}

	# by default use plot type 3 (2 data panels), but if we only have DNA or only RNA then only use 1 data panel (plot type 4)
	karyoplotr_plottype = 3
	if (missing(exp_data) | missing(cnv_data)) {
		karyoplotr_plottype = 4
	}

	# TODO: add kpAddLabels() for xlab and ylabs and add kpAddMainTitle() for main title
	kp <- karyoploteR::plotKaryotype(genome = "hg19", chromosomes = chromosomes, plot.type = karyoplotr_plottype, cex = 0.85, main = gsub(".*__", "", title))
	karyoploteR::kpDataBackground(kp, data.panel = 2, r0 = dna_r1, r1 = dna_r0)
	karyoploteR::kpAxis(karyoplot = kp, col = "gray50", data.panel = 2, ymin = 0, ymax = cnv_data_upper_cutoff, r1 = dna_r0, r0 = dna_r1, cex = 0.5, tick.pos = c(0,1,2,3,4), side = 2)
	karyoploteR::kpPoints(kp, chr = cnv_data$Chr, x = cnv_data$Pos, y = cnv_data$rawCN, col = "#00A6EDAA", data.panel = 2, ymin = 0, ymax = cnv_data_upper_cutoff, r0 = dna_r1, r1 = dna_r0, clipping = T)
	karyoploteR::kpSegments(karyoplot = kp, data.panel = 2, chr = cnv_data$Chr, x0 = cnv_data$Pos, x1 = cnv_data$end, y0 = cnv_data$segmentedCN, y1 = cnv_data$segmentedCN, ymin = 0, ymax = cnv_data_upper_cutoff, r0 = dna_r1, r1 = dna_r0, clipping = T, col = "black", lwd = 3)
#	karyoploteR::kpText(karyoplot = kp, data.panel = 2, chr = "chr1", y = 0.283, x = 0, labels = "Copy Number\n\n", srt = 90, cex = 0.8)

	if (!missing(exp_data)) {
		karyoploteR::kpDataBackground(kp, data.panel = 1, r0 = rna_r0, r1 = rna_r1)
		karyoploteR::kpAxis(karyoplot = kp, col = "gray50", data.panel = 1, ymin = round(min(exp_data$SegmentExp)), ymax = exp_data_upper_cutoff, r1 = rna_r1, r0 = rna_r0, cex = 0.5, numticks = 10, side = 2)
		karyoploteR::kpPoints(kp, chr = exp_data$Chr, x = exp_data$Start, y = exp_data$SegmentExp, col = "#FFBD07AA", data.panel = 1, ymin = round(min(exp_data$SegmentExp)), ymax = exp_data_upper_cutoff, r0 = rna_r0, r1 = rna_r1, clipping = T)
		karyoploteR::kpSegments(karyoplot = kp, data.panel = 1, chr = exp_data$Chr, x0 = exp_data$Start, x1 = exp_data$End, y0 = exp_data$PCFexp, y1 = exp_data$PCFexp, ymin = round(min(exp_data$SegmentExp)), ymax = exp_data_upper_cutoff, r0 = rna_r0, r1 = rna_r1, clipping = T, col = "black", lwd = 3)
	#	karyoploteR::kpText(karyoplot = kp, data.panel = 1, chr = "chr1", y = 0.3, x = 0, labels = "Segmented Expression\n\n", srt = 90, cex = 0.8)
	}

	if (save) {
		dev.off()
	}
}

