# Barman
![](https://github.com/seanlaidlaw/barman/blob/master/img/barman_logo.png?raw=true)

Barman is an R library for easy data preprocessing and visualization of Genome and Transcriptome single-cell data.


## Installation
Barman is still under active development, as such installation is done via the actively maintained [git repository](https://github.com/seanlaidlaw/barman)

```
devtools::install_github("seanlaidlaw/barman")
```

## Standard workflow

Starting from a gene x cell counts matrix barman functions can be used in the following order to
achieve a standard workflow:

```
                +------------------------------+
                |       RNA Counts Matrix      |
                +--------------+---------------+
                               |
                               |
                               |                           +-----+
                               v                                 |
                +--------------+---------------+                 |
                | filter_and_normalise_scrna() |                 |
                +--------------+---------------+                 +
                               |                              counts
                               v                           normalisation
                +--------------+---------------+                 +
                |      logR_by_ref_group()     |                 |
                +---------+---------+----------+                 |
                          |         |                            |
                          |         |                      +-----+
                          |   OR    |
                          |         |                      +-----+
                          v         v                            |
+-------------------------+-+    +--+------------------------+   |
|   expression_boxplots()   |    |get_expression_by_segment()|   |
+---------------------------+    +-------------+-------------+   |
                                               |                 |
                                               |           visualisation
                                               v                 |
                                 +-------------+-------------+   |
                                 |       G_T_chr_plot()      |   |
                                 +---------------------------+   |
                                                                 |
                                                           +-----+
```

## Quick Start

Barman is designed for multiomic data, but provides functions for preprocessing DNA and RNA
independently.

Preprocessing the RNA requires a counts matrix, which can be generated from single-cell RNA bams
using a tool such as [FeatureCounts](dx.doi.org/10.1093/bioinformatics/btt656).

```
# Load the counts
counts_matrix = read.table("my_featurecounts_matrix.tsv", sep = "\t", header = TRUE)
```

Filtering can be done automatically using PCA outlier detection, and the resulting filtered matrix
is noramlized for gene length using FPKM.

```
# Filtering and Normalizing
counts_matrix = filter_and_normalise_scrna(counts_matrix) # automatic filtering
```

Filtering can also be manually specified using the `manual_filter` option and providing a list of
2 length vectors, specifying the lower and upper limits for `total_counts`, `total_features_by_counts`,
%MT genes, and %ERCC respectively.

```
# Filtering and Normalizing
fpkm_matrix = filter_and_normalise_scrna(counts_matrix, manual_filter = list(c(100000,100000000),
c(1000, 5000), 20, 50)) # filters out cells according to manual filters
```

We can see the difference in the QC by running the `qc_plots` function

```
qc_plots(fpkm_matrix)
```


We can also normalize based on a given reference group by using the `logR_by_ref_group` function to
calculate the log fold change between a group of specified reference cells and the rest:

```
reference_group_cells = rownames(fpkm_matrix)[fpkm_matrix$cell_type = "control",]
logFC_matrix = logR_by_ref_group(fpkm_matrix, reference_group_cells)
```


We can plot expression per chromosome of different cells by using the `expression_boxplots` function
which accepts two lists of cell ids, one as an experimental group, and one as a control. From this
it calcualtes the average expression per gene, and groupes by chromosome before plotting.
```
expression_boxplots(experimental_group = list("Cell_A01", "Cell_A02", "Cell_B03"), control_group =
list("Cell_D03", "Cell_D04", "Cell_D05", fpkm_matrix)
```



If scCNV has been run on the DNA bams then we can use the produced segmentation files, and a counts
matrix to produce segmentation files for the RNA expression data, thus allowing us to compare
segment to segment between DNA and RNA.

```
scCNV_segment_for_A01 = read.table("SEGMENT_scCNV_segment_for_A01.copynumber.refLocus.txt"
A01_expression_segments = get_expression_by_segment("Cell_A01", fpkm_matrix, scCNV_segment_for_A01)
```

This can also be done in bulk, without having to read in the scCNV produced segment files by using
the `bulk_get_expression_by_segment` function. this function uses all available cores except one to
process the files in parallel.

```
bulk_get_expression_by_segment(fpkm_matrix, "./Path_to_scCNV_outputs/", "./RNA_seg_output_folder")
```


Finally, we can see the Genome and Transcriptome segment plot, showing the DNA to RNA comparison by
running the `G_T_chr_plot` function.
```
G_T_chr_plot(cnv_data = scCNV_segment_for_A01, exp_data = A01_expression_segments, title = "Cell A01 Genome and Transcriptome")
```

