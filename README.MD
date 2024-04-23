# SpatialTouchstone

Welcome to `SpatialTouchstone`, an innovative tool designed to facilitate quality control for spatialomics technology data.

## Overview

Spatial omic technologies are revolutionizing how we understand biological tissues, enabling the mapping of molecular compositions at unprecedented spatial resolutions. However, with great power comes great complexity. The interpretation of spatial omics data requires robust quality control measures to ensure the reliability and reproducibility of the findings. `SpatialTouchstone` is developed to meet this need, providing a comprehensive suite of quality control metrics specifically tailored for spatial omics data.

## Functions

### readSpatial

Reads and preprocesses spatial transcriptomics data from a specified path for a given imaging based platform. The function supports 'Xenium' and 'CosMx' platforms, performing platform-specific loading and preprocessing steps. This includes loading the data, annotating it with sample metadata and processing cell metadata. For 'Xenium', it adds additional cell metadata and tissue coordinates as an embedding for custom plotting. For 'CosMx', it separates targeting and non-targeting probes, adds additional cell metadata and includes tissue coordinates as an embedding. 'Merscope' platform support is under development. Regardless of platform, data is stored in an assay named "RNA" for convenient and for each platform, this table will be used by subsequent functions.

### readTxMeta

Reads transcriptome metadata from a specified path, depending on the platform specified. Currently supports 'Xenium' and 'CosMx' platforms. For 'Xenium', it reads 'transcripts.csv.gz' and renames the 'feature_name' column to 'target' for consistency. For 'CosMx', it reads the appropriate transcriptome file matched by '\*tx_file.csv.gz'. 'Merscope' platform support is under development. For each platform, this table will be used by subsequent functions.

### getGlobalFDR

It calculate 'specificity' as Global False Discovery Rate (FDR) for specified features (or all features by default) in a Seurat object. It leverages transcript localization data to calculate FDR based on the proportion of negative control or blank barcode expressions compared to the expression of actual genes.

### getTxPerCell

This function calculates the mean number of transcripts per cell for a specified set of features (genes) in a Seurat object. If no features are specified, the function defaults to using all targets available within the RNA assay of the provided Seurat object. It's a useful metric for assessing the overall transcriptional activity within the sampled cells.

### getTxPerArea

Calculates the mean number of transcripts per unit area for specified features (genes) in a Seurat object. This calculation provides insight into the density of transcriptional activity relative to cell area, offering a normalized measure of gene expression that accounts for differences in cell size.

### getTxPerNuc

This function calculates the mean number of transcripts localized within the nucleus for a specified set of features (genes) across all cells in a given dataset. This measurement is specific to spatial transcriptomics data where the distinction between nuclear and cytoplasmic transcript localization can be made.

### getMeanExpression

Computes the mean expression levels for a specified list of features (genes) and control probes within a Seurat object. Computes the mean expression levels for a specified list of features (genes) and control probes within a Seurat object. This function separately calculates the mean expression for both target genes and control probes, then combines these results into a single data frame.

### getMeanSignalRatio

Computes the log-ratio of the mean expression levels of specified genes (or all genes if none are specified) to the mean expression levels of negative control probes within a Seurat object. This metric can provide insights into the overall signal strength relative to background noise levels.

### getCellTxFraction

Calculates the fraction of transcripts in cells for a given set of features or all features. The function works by reading transcript metadata for the given sample and platform, then filtering for the specified features (if any). For each platform, it identifies transcripts that are not assigned to any cell ('UNASSIGNED' for Xenium, 'None' for CosMx, etc.) and calculates the fraction of transcripts that are assigned to cells.

### getMaxRatio

Identifies the maximum mean expression value among the specified features (or all features if none are specified) and calculates its log-ratio to the mean expression value of negative control probes. This log-ratio reflects the dynamic range of the dataset, indicating the spread between the

### getMaxDetection

Calculates the distribution of maximal values by identifying the maximal expression values across the specified set of features (or all features if none are specified) within the dataset. This function illustrates the upper bounds of gene expression, providing crucial information for assessing the dataset's dynamic range and the sensitivity of detection methods used in the experiment. The function aggregates these maximal values and presents them in a data frame, facilitating further analysis of the expression distribution and detection efficiency across different genes.

### getMECR

Calculates the Mutually Exclusive Co-expression Rate (MECR) Implementation. This function is based on a predefined set of markers and their associated cell types, computing the rate at which pairs of markers are expressed in mutually exclusive patterns within cells. The calculation is limited to markers found in the dataset and can be adjusted to focus on a subset by specifying it. The result is a single MECR value that quantifies the overall mutually exclusive co-expression pattern, rounded to three decimal places. This is particularly useful for understanding cellular heterogeneity and the specificity of marker gene expression within a spatial transcriptomics dataset, offering insights into the cellular composition and the distinctiveness of cell-type-specific gene expression patterns.

### getMorans

Calculates spatial autocorrelation using Moran's I. This function processes gene-targeting probes by converting the input data into a `SingleCellExperiment` object, then into a `SpatialFeatureExperiment` object. Spatial coordinates are specified and used to further process the data, including normalization and nearest neighbor graph construction. Subsequently, Moran's I is calculated for each feature using the `Voyager` package. This process is repeated for control probes. The results for both gene-targeting and control probes are combined into a single data frame. The spatial autocorrelation analysis provides Moran's I values for each feature, indicating the degree of spatial clustering, and is performed separately for gene-targeting probes and control probes, allowing for a comparison of spatial patterns in gene expression and background noise.


### getSilhouetteWidth

Calculates silhouette width as cluster evaluation, providing a metric for assessing the quality of the clustering. By clustering the data and calculating the silhouette width, this function evaluates how similar an object is to its own cluster compared to other clusters. The silhouette width ranges from -1 to 1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. This measure is useful for determining the appropriate number of clusters in the data, guiding the interpretation of spatial patterns and cellular heterogeneity within the dataset.

### getSparsity

Shows the sparsity (as a count or proportion) of a matrix, which is a crucial aspect of high-dimensional data analysis, especially in single-cell and spatial transcriptomics. A high sparsity value (close to 1) means that the majority of the matrix elements are zero, indicating a large number of non-detected genes or features in the dataset. Understanding sparsity is important for choosing the appropriate data normalization and transformation techniques, which can significantly impact downstream analysis and interpretations. This function calculates the sparsity of the RNA counts matrix in a Seurat object, providing insights into the data's density and informing decisions on data preprocessing steps.

### getEntropy

Calculates the entropy of the expression matrix in a Seurat object, offering a measure of the distribution's uniformity across the detected genes. Entropy is a concept from information theory that measures the uncertainty or disorder within a dataset. In the context of transcriptomics, a higher entropy value can indicate a more heterogeneous expression profile across cells, reflecting the complexity and diversity of transcriptional states within the sample. By quantifying entropy, this function provides a summary statistic of the dataset's overall expression diversity, aiding in the understanding of cellular complexity and the identification of transcriptionally distinct cell populations.

### Plot function in utils.R

-   plotSampleLabel
-   plotPanelSize
-   plotCellCount
-   plotTxPerCell
-   plotTxPerArea
-   plotTxPerNuc
-   plotTxPerCellNorm
-   plotFractionTxInCell
-   plotSignalRatio
-   plotMeanExpression
-   plotMaxExpression
-   plotMECR
-   plotMorans
-   plotSilhouette
-   plotSparsity
-   plotEntropy
