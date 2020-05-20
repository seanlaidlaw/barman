#' @title get_copynumber_for_gene
#'
#' @description Groups all segments that overlap each gene boundry and returns that to the
#' G_T_chr_plot which calls this internal function
#'
#' @return a list of values corresponding to one segment along with the PCF of the segments
#'
get_copynumber_for_gene = function(dna_segments, chr, start, end, gene_name) {

	if (!missing(gene_name)) {
	mart72.hs <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl", host = "grch37.ensembl.org")

		if (grepl("ENSG[0-9]*", gene_name)) {
			new_gene_ids = biomaRt::getBM(
				attributes=c("ensembl_gene_id","chromosome_name","start_position","end_position"),
				filters="ensembl_gene_id",
				values=gene_name,
				mart=mart72.hs)
		} else {
			new_gene_ids = biomaRt::getBM(
				attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"),
				filters="hgnc_symbol",
				values=gene_name,
				mart=mart72.hs)

		}
		if (nrow(new_gene_ids) < 1) {
			stop(paste0("No matches found on biomaRt for gene_name:", gene_name))
		}

		chr = new_gene_ids$chromosome_name
		start = new_gene_ids$start_position
		end = new_gene_ids$end_position
	} else {
		# if no gene+name given make sure that chr start and end are given
		if (missing(chr) | missing(start) | missing(end)) {
			stop("Please either provide gene_name arg or chr,start,end arguments")
		}
	}


	# look for segments that entirely contain gene
	n_entire_containing_segments = dna_segments[(dna_segments$Chr == as.character(chr) & dna_segments$Pos <= as.integer(start) & dna_segments$end >=as.integer(end)),]

	if (nrow(n_entire_containing_segments) > 0) {
		segStart = n_entire_containing_segments$Pos
		segEnd = n_entire_containing_segments$end
		segraw_1 = n_entire_containing_segments$rawCN
		segraw_2 = n_entire_containing_segments$rawCN
		segsegCN_1 = n_entire_containing_segments$segmentedCN
		segsegCN_2 = n_entire_containing_segments$segmentedCN
	} else {
		# if gene is spread over multiple segments look for the closest near start
		n_start_containing_segments = dna_segments[(dna_segments$Chr == as.character(chr)),]
		n_start_containing_segments = n_start_containing_segments[order(n_start_containing_segments$Pos),]
		n_start_containing_segments = tail(n_start_containing_segments[(n_start_containing_segments$Pos <= as.integer(start)),], 1)
		segStart = n_start_containing_segments$Pos
		segraw_1 = n_start_containing_segments$rawCN
		segsegCN_1 = n_start_containing_segments$segmentedCN

		n_end_containing_segments = dna_segments[(dna_segments$Chr == as.character(chr)),]
		n_end_containing_segments = n_end_containing_segments[order(n_end_containing_segments$end),]
		n_end_containing_segments = head(n_end_containing_segments[(n_end_containing_segments$end >= as.integer(end)),], 1)
		segEnd = n_end_containing_segments$end
		segraw_2 = n_end_containing_segments$rawCN
		segsegCN_2 = n_end_containing_segments$segmentedCN
	}

	return_values = list(chr, segStart, segEnd, segraw_1, segraw_2, segsegCN_1, segsegCN_2)
	return(return_values)
}


