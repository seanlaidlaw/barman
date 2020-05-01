#' @title remove_non_specific_annotation_columns
#'
#' @description removes columns from table where all values are the same or where all values are distinct strings
#'
#' @param annotation_table a table of annotations we want to prune
#' @param exclude a list of column names we want to keep even if they are non-specific
#'
#' @return a pruned annotation table
#'
#' @export

remove_non_specific_annotation_columns <- function(annotation_table, exclude) {
	# exclude argument can be given as c("filename") or as "filename". make sure provided argument is either list or char format
	if (!missing(exclude)) {
		if (typeof(exclude) == "character") {
			exclude = as.list(exclude)
		}

		if (typeof(exclude) != "list") {
			print("Error: function was passed an exclude option but it wasnt in list format: e.g. c('Col1','Col2')")
			return(NA)
		}
	} else {
		exclude = NA
	}


	# if a column in annotation table only contains one unique value then remove it
	for (i in colnames(annotation_table)) {
		if (!(i %in% exclude)) {
			if (length(unique(annotation_table[[i]])) == 1) {
					annotation_table[[i]] = NULL
			}
		}
	}


	# if a column in annotation table only contains unique values then remove it
	for (i in colnames(annotation_table)) {
		if (!(i %in% exclude)) {
			if (all(is.numeric(annotation_table[[i]])) == FALSE) {
				if (length(unique(annotation_table[[i]])) == length(annotation_table[[i]])) {
						annotation_table[[i]] = NULL
				}
			}
		}
	}

	return(annotation_table)
}
