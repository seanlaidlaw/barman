#' @title order_chromosome_list
#'
#' @description sorts a vector of strings/ints numerically then alphabetically to show chromosmes in
#' proper order
#'
#' @param random_order_chr_list a vector of unordered values either strings or integers
#'
#' @return a vector of sorted values first by numeric order then by alphabetical order
#'
#' @export
order_chromosome_list <- function(random_order_chr_list) {
	# determine if chromosome is numeric or alpha
	numeric_chrs = suppressWarnings(rand_chr_list[!is.na(as.numeric(random_order_chr_list))])
	alpha_chrs = suppressWarnings(random_order_chr_list[is.na(as.numeric(random_order_chr_list))])

	# sort order of chromosomes for numeric then alpha
	numeric_chrs = numeric_chrs[order(as.integer(numeric_chrs), decreasing = F)]
	alpha_chrs = alpha_chrs[order(alpha_chrs, decreasing = F)]

	# create a new list of numeric ordered chromosmes followed by alphabetical ordered ones
	ordered_chr_list = numeric_chrs
	ordered_chr_list = append(x = ordered_chr_list, values = alpha_chrs)

	# throw error if for some reason output doesnt contain all the input chromosomes
	if (!all(random_order_chr_list %in% ordered_chr_list)) {
		warning("not all input items passed to order_chromosome_list() are in output")
	}

	return(ordered_chr_list)
}
