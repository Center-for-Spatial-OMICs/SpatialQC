
### Loading Packages

# Loading packages
packages <- c("Seurat", "dplyr", "ggplot2", "shadowtext", "scales", "cowplot", "data.table", "Matrix", "matrixStats", "SingleCellExperiment", "SpatialExperiment", "SpatialFeatureExperiment", "bluster", "BiocParallel", "BioQC", "coop", "fs")

# Check and install packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Install the package from CRAN
    install.packages(pkg)

    # Check if it's a Bioconductor package
    if (pkg %in% c("SingleCellExperiment", "SpatialExperiment", "SpatialFeatureExperiment", "bluster", "BiocParallel", "BioQC")) {
      # Install the Bioconductor package
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    }

    # Load the package
    library(pkg, character.only = TRUE)
  } else {
    # Load the package
    library(pkg, character.only = TRUE)
  }
}

globalVariables(c("target", "celltype_pred", "cell_type", "FDR", "nCount_RNA", "cell_area", "cell_ID", "Area.um2", "metric","platform","value","cell","cor","sample_id","value_ref","group","prop","squish","type","cell_id","overlaps_nucleus","feature_name"))

#' @importFrom magrittr %>%
#' @description
#' Compute various metrics for a given set of samples
#' @details
#' Computes various metrics for spatial transcriptomic data for a given set of samples. This function processes a data.frame containing sample information, including paths to necessary files, to compute a variety of metrics for each sample.
#' The metrics outputted include the number of cells, specificity FDR, number of transcripts per cell, transcripts per area of segmented cell, transcripts per nucleus, signal to noise ratio, fraction of transcripts in segmented cells, Mutually Exclusive Co-Expression Rate (MECR), sparsity, and Shannon entropy.
#' To initiate vignette, the data frame must include the columns: sample_id, platform, expMat, tx_file, and cell_meta.
#' @title getAllMetrics
#' @param df_samples A data.frame that must contain the following columns:
#' \itemize{
#'   \item{sample_id}{Unique identifier for each sample.}
#'   \item{platform}{The platform used, either "Xenium" or "CosMx".}
#'   \item{expMat}{Path to the expression matrix file. For Xenium platform, this should be a path to the cell_feature_matrix folder.}
#'   \item{tx_file}{Path to the transcript file.}
#'   \item{cell_meta}{Path to the cell metadata file.}
#' }
#' @return A modified version of df_samples with the computed metrics added as new columns.
#' @export
getAllMetrics <- function(df_samples, features = NULL) {
  # Ordering columns to match to mapply()
  df_samples <- df_samples[, c('dataset','platform', 'expMat', 'tx_file', 'cellSegMeta')]

  gStartTime <- Sys.time()

  tryCatch({
    start.time <- Sys.time()
    print("Calculating number of cells ")
    df_samples$nCells <- NA
    df_samples$nCells <- mapply(getNcells, expMat = df_samples[, 3], platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Specificity
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Specificity ")
    df_samples$specificityFDR <- NA
    df_samples$specificityFDR <- mapply(getGlobalFDR, tx_file = df_samples[, 4],
                                        platform = df_samples[, 2], cellSegMeta =  df_samples[, 5],
                                        MoreArgs = list(features = features))
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Tx per Cell
  tryCatch({
    start.time <- Sys.time()
    print("Calculating number of transcripts per cell ")
    df_samples$TxPerCell <- NA
    df_samples$TxPerCell <- mapply(getTxPerCell, expMat = df_samples[, 3] , platform = df_samples[, 2],
                                   MoreArgs = list(features = features))
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Tx per Area
  tryCatch({
    start.time <- Sys.time()
    print("Calculating number transcripts per area of segmented cell ")
    df_samples$TxPerArea <- NA
    df_samples$TxPerArea <- mapply(getTxPerArea, platform = df_samples[, 2] ,
                                   cellSegMeta =  df_samples[, 5], MoreArgs = list(features = features) )
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # Tx per Nucleus
  tryCatch({
    start.time <- Sys.time()
    print("Calculating number transcripts per nuclei ")
    df_samples$TxPerNuc <- NA
    df_samples$TxPerNuc <- mapply(getTxPerNuc, tx_file = df_samples[, 4] , platform = df_samples[, 2],
                                  MoreArgs = list(features = features))
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Signal to noise ratio
  tryCatch({
    start.time <- Sys.time()
    print("Calculating signal to noise ratio ")
    df_samples$SigNoiseRatio <- NA
    df_samples$SigNoiseRatio <- mapply(getMeanSignalRatio, expMat = df_samples[, 3] , platform = df_samples[, 2] )
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Fraction of txs in cells
  tryCatch({
    start.time <- Sys.time()
    print("Calculating fraction of transcripts in segmented cells ")
    df_samples$CellTxFraction <- NA
    df_samples$CellTxFraction <- mapply(getCellTxFraction, tx_file = df_samples[, 4] ,
                                        platform = df_samples[, 2], MoreArgs = list(features = features) )
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  #Dynamic range
  tryCatch({
    start.time <- Sys.time()
    print("Calculating dynamic range ")
    df_samples$DynamicRange <- NA
    df_samples$DynamicRange <- mapply(getMaxRatio, expMat = df_samples[, 3] ,
                                      platform = df_samples[, 2], MoreArgs = list(features = features) )
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # MECR
  tryCatch({
    start.time <- Sys.time()
    print("Calculating MECR (Mutually Exclusive Co-Expression Rate) ")
    df_samples$MECR <- NA
    df_samples$MECR <- mapply(getMECR, expMat = df_samples[, 3] , platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # Sparsity
  tryCatch({
    start.time <- Sys.time()
    print("Calculating sparsity ")
    df_samples$sparsity <- NA
    df_samples$sparsity <- mapply(getSparsity, expMat = df_samples[, 3] , platform = df_samples[, 2],
                                  MoreArgs = list(features = features))
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Entropy
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Shannon entropy ")
    df_samples$entropy <- NA
    df_samples$entropy <- mapply(getEntropy, expMat = df_samples[, 3] , platform = df_samples[, 2],
                                 MoreArgs = list(features = features))
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # Complexity
  tryCatch({
    start.time <- Sys.time()
    print("Calculating complexity ")
    df_samples$Complexity <- NA
    df_samples$Complexity <- mapply(getComplexity, expMat = df_samples[, 3] , platform = df_samples[, 2],
                                    MoreArgs = list(features = features))
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # GlobalFDR
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Global FDR ")
    df_samples$GlobalFDR <- NA
    df_samples$GlobalFDR <- mapply(getGlobalFDR, tx_file = df_samples[, 4] , platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})

  # MeanExpression
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Mean Expression ")
    df_samples$MeanExpression <- NA
    df_samples$MeanExpression <- mapply(getMeanExpression,  expMat = df_samples[, 3] , platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # MaxDetection
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Max Expression ")
    df_samples$MaxExpression <- NA
    df_samples$MaxExpression <- mapply(getMaxDetection,  expMat = df_samples[, 3] , platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # Sparsity
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Sparsity ")
    df_samples$Sparsity <- NA
    df_samples$Sparsity <- mapply(getSparsity,  expMat = df_samples[, 3] , platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})


  # PanelSize
  tryCatch({
    start.time <- Sys.time()
    print("Calculating Panel Size ")
    df_samples$PanelSize <- NA
    df_samples$PanelSize <- mapply(getSparsity,  expMat = df_samples[, 3] , platform = df_samples[, 2])
    end.time <- Sys.time()
    print(end.time - start.time)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "")})
  gEndTime <- Sys.time()
  print(paste0("total time: ", round(gEndTime - gStartTime , digits = 2) ))
  return(df_samples)

}


#' @title readSpatial
#' @description
#' Reads data from either imaging based spatial transcriptomic (CosMx, MERSCOPE, Xenium)
#' and outputs a seurat object(seu_obj) with some common metadata for downstream comparison.
#' @details
#' Reads and preprocesses spatial transcriptomics data from a specified path for a given platform.
#' The function supports 'Xenium' and 'CosMx' platforms, performing platform-specific loading and
#' preprocessing steps. This includes loading the data, annotating it with sample metadata, and
#' processing cell metadata. For 'Xenium', it adds additional cell metadata and tissue coordinates
#' as an embedding for custom plotting. For 'CosMx', it fixes assays to separate targeting and
#' non-targeting probes, adds additional cell metadata, and includes tissue coordinates as an embedding.
#' 'Merscope' platform support is under development.
#' @param sample_id Identifier for the sample being processed.
#' @param path The file path from which to load the spatial transcriptomics data.
#' @param seurat Return a seurat object or not.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
##' @return Returns a Seurat object containing the loaded and processed spatial transcriptomics data,
#'         including sample metadata, cell metadata, and embeddings for tissue coordinates. The
#'         structure and content of the returned object vary depending on the platform.
#' @importFrom Seurat CreateSeuratObject
#' @importFrom data.table fread
#' @export

readSpatial <- function(sample_id, path, platform=NULL, seurat=FALSE){
  print(paste0("Reading: ", sample_id))

  # If platform not specified, try to guess which tech it is (from folder name)
  if(is.null(platform)) {
    xenium <- grep('*[X-x]enium*|*10x*', sample_id)
    cosmx <- grep('*[C-c]os[M-m]x*|*[N-n]anostring*', sample_id)
    mersc <- grep('*[M-m]erscope*|*MERSCOPE*|*[V-v]izgen*',sample_id)

    if(length(xenium) == 0 & length(cosmx) == 0 & length(mersc) == 0) {
      platform = 'Unknown'
    }

    if(length(xenium) > 0) { platform = 'Xenium'}
    if(length(cosmx) > 0) { platform = 'CosMx'}
    if(length(mersc) > 0) { platform = 'Merscope'}
  }

  ## Read and store tables as list outside of seurat
  if(seurat==FALSE) {
    # create empty list
    obj_list <- list()

    if(platform == 'Xenium') {
      obj_list[[sample_id]] <- list()
      obj_list[[sample_id]][['expMatrix']] <- Matrix::readMM(file.path(path, 'cell_feature_matrix/matrix.mtx.gz'))
      cols <- data.table::fread(file.path(path, 'cell_feature_matrix/barcodes.tsv.gz'), header = F)
      rows <- data.table::fread(file.path(path, 'cell_feature_matrix/features.tsv.gz'), header = F)
      rownames(obj_list[[sample_id]][['expMatrix']]) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(obj_list[[sample_id]][['expMatrix']]) <- cols$V1 ## this is the barcodes of cells

      obj_list[[sample_id]][['TxMatrix']] <- data.table::fread(file.path(path, 'transcripts.csv.gz'))
    }

    if(platform == 'CosMx') {
      obj_list[[sample_id]] <- list()
      pattern <- "exprMat_file.csv.gz$|exprMat*"
      file_name <- list.files(path = path, pattern = pattern, full.names = TRUE)
      obj_list[[sample_id]][['expMatrix']] <- data.table::fread(file_name)
      pattern <- "*tx_file.csv.gz$|*tx.csv$|tx_file_unique.csv.gz$"
      file_name <- list.files(path = path, pattern = pattern, full.names = TRUE)
      obj_list[[sample_id]][['TxMatrix']] <- data.table::fread(file_name[1])
    }

    # if(platform == 'Merscope') {
    #   obj_list[[sample_id]] <- list()
    #   pattern <- "exprMat_file.csv.gz$"
    #   file_name <- list.files(path = path, pattern = pattern, full.names = TRUE)
    #   obj_list[[sample_id]][['expMatrix']] <- data.table::fread(file_name)
    #   obj_list[[sample_id]][['TxMatrix']] <- data.table::fread(file.path(path, 'transcripts.csv.gz'))
    # }
    #

    return(obj_list)
  }


  if(platform == "Xenium"){
    print("Loading Xenium data")
    seu_obj <- LoadXenium(path, assay = "RNA")
    seu_obj@meta.data$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data

    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, "cells.csv.gz")
    )
    #Set up a few defined metadata columns
    seu_obj@meta.data$cell_area <- cell_meta$cell_area
    seu_obj@meta.data$nucleus_area <- cell_meta$nucleus_area
    seu_obj@meta.data$transcript_counts <- cell_meta$transcript_counts
    seu_obj@meta.data$negprobe_counts <- cell_meta$control_probe_counts

    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- GetTissueCoordinates(seu_obj)
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(coords, key="Tissue_", assay="RNA")

  } else if(platform == "CosMx"){
    print("Loading CosMx data")
    seu_obj <- LoadNanostring(path, fov="fov")
    seu_obj$sample_id <- sample_id
    seu_obj@meta.data$platform <- platform
    seu_obj@meta.data$path <- path #used in some functions to pull tx data

    #Fix assays to separate targeting and non-targeting probes
    ## Add negative control probe assay
    sys_probes <- grep("SystemControl", rownames(seu_obj), value=T)
    neg_probes <- grep("Negative", rownames(seu_obj), value=T)
    seu_obj[["ControlProbe"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[neg_probes,]
    )
    ## Make "Nanostring" assay
    tx_probes <- rownames(seu_obj)[!rownames(seu_obj) %in% c(sys_probes, neg_probes)]
    seu_obj[["RNA"]] <- CreateAssayObject(
      counts = seu_obj[["Nanostring"]]$counts[tx_probes,]
    )
    DefaultAssay(seu_obj) <- "RNA"
    seu_obj[["Nanostring"]] <- NULL

    ##### Cell Metadata
    print("Getting additional cell metadata")
    cell_meta <- data.table::fread(
      file.path(path, list.files(path, pattern="*metadata_file.csv.gz"))
    )
    #It's excessive, but we'll add all metadata into the object
    seu_obj@meta.data <- cbind(seu_obj@meta.data, cell_meta)
    seu_obj@meta.data$fov <- factor(paste0("FOV", seu_obj@meta.data$fov))
    seu_obj@meta.data$cell_area <- seu_obj$Area.um2
    seu_obj@meta.data$transcript_counts <- seu_obj$nCount_RNA
    seu_obj@meta.data$negprobe_counts <- seu_obj$nCount_ControlProbe

    ### Add tissue coordinates as an embedding for custom plotting
    print("Adding tissue coordinates as embedding")
    coords <- data.frame(
      Tissue_1 = cell_meta$CenterY_global_px,
      Tissue_2 = cell_meta$CenterX_global_px
    )
    coords <- as.matrix(coords[,1:2])
    colnames(coords) <- c("Tissue_1", "Tissue_2")
    rownames(coords) <- colnames(seu_obj)
    seu_obj[["tissue"]] <- CreateDimReducObject(m , key="Tissue_", assay="RNA")


  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()

  } else{
    print("Not a supported platform")
    stop()

  }

  return(seu_obj)
}

#' @title readTxMeta
#' @description
#' Reads the transcript localization and metadata table.
#' This table will be used by subsequent function in this package.
#' @details
#' #' This function reads transcriptome metadata from a specified path, depending on the platform specified.
#' Currently supports 'Xenium' and 'CosMx' platforms. For 'Xenium', it reads 'transcripts.csv.gz' and renames
#' the 'feature_name' column to 'target' for consistency. For 'CosMx', it reads the appropriate transcriptome
#' file matched by '*tx_file.csv.gz'. 'Merscope' platform support is under development.
#'
#' @param path The file path from which to read the transcriptome metadata.
#' @param platform  The platform for which the transcriptome metadata is being read. Valid options are
#'        'Xenium', 'CosMx', and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return Returns a data table with the transcriptome metadata. For 'Xenium', the 'feature_name' column is
#'         renamed to 'target'. For 'CosMx', the transcriptome file matching '*tx_file.csv.gz' is read.
#'         No return value for 'Merscope' as it stops execution with an error message.
#' @export
#' @importFrom data.table fread setnames
readTxMeta <- function(path, platform){
  if(platform == "Xenium"){
    df <- data.table::fread(file.path(path, "transcripts.csv.gz"))
  } else if(platform == "CosMx"){
    df <- data.table::fread(file.path(path,
                                      list.files(path, pattern = "*tx_file.csv.gz")))
    df <- unique(df)
  } else if(platform == "Merscope"){
    print("Working on support!")
    stop()
  } else{
    print("Platform not supported")
    stop()
  }
}


######## General Utilities ########

#' @title getPseudobulk
#' @description
#' Creates pseudobulk gene expression matrices from single-cell RNA-seq data.
#' @details
#' This function generates pseudobulk expression matrices by averaging the expression levels of cells within the same cell type.
#' It takes a Seurat object as input, along with a specification for the metadata column that defines cell types.
#' Each column in the output matrix represents a cell type, with row names corresponding to genes and cell counts averaged across cells of the same type.
#'
#' @param seu_obj A Seurat object containing single-cell RNA-seq data.
#' The object must have a counts matrix in the "RNA" assay and metadata that includes cell type annotations.
#' @param celltype_meta A string specifying the name of the metadata column in the Seurat object that contains cell type annotations.
#' Defaults to "cell_type".
#'
#' @return A matrix where each column represents a pseudobulk sample (averaged expression levels) for a specific cell type.
#' Row names correspond to genes.

getPseudobulk <- function(seu_obj, celltype_meta="cell_type") {


  celltype <- factor(seu_obj@meta.data[,celltype_meta])
  names(celltype) <- colnames(seu_obj)
  mat <- seu_obj[["RNA"]]$counts

  mat.summary <- do.call(cbind, lapply(levels(celltype), function(s) {
    cells <- names(celltype)[celltype==s]
    pseudobulk <- rowMeans(mat[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}



#' @title annotateData
#' @param seu_obj A Seurat object containing spatial transcriptomics data that you wish to annotate.
#' @param ref A Seurat object serving as the reference, containing single-cell RNA-seq data with cell type annotations.
#' @param celltype_meta A string specifying the name of the metadata column in the reference Seurat object that contains cell type annotations.
#' Defaults to "cell_type".
#'
#' @return The input Seurat object (`seu_obj`) with an additional metadata column (`celltype_pred`) containing the cell type predictions.
#' @description
#' Automates the annotation of spatial transcriptomics data using cell annotations from reference single-cell RNA-seq data.
#' It filters cells based on transcript counts and employs a machine learning model for cell type prediction.
#'
#' @details
#' `annotateData` leverages the power of single-cell reference data to annotate spatial transcriptomics data.
#' This approach is particularly useful in studies aiming to understand tissue composition and cellular localization
#' without the need for extensive manual annotation. The function requires the input spatial data (`seu_obj`) to be in
#' Seurat format and a similarly formatted single-cell RNA-seq dataset (`ref`) as a reference. The reference must contain
#' a cell type annotation column, which can be specified using the `celltype_meta` parameter. The function emphasizes the
#' importance of data pre-processing, as cells with low transcript counts are removed to improve the accuracy of the cell
#' type predictions made by the `insitutypeML` algorithm.
#' @importFrom InSituType insitutypeML
#' @export


annotateData <- function(seu_obj, ref, celltype_meta="cell_type"){
  print("Getting pseudobulk for reference")
  ref_mat <- getPseudobulk(ref)

  cells_keep <- colnames(seu_obj)[colSums(seu_obj[["RNA"]]$counts) > 10]
  seu_obj <- subset(seu_obj, cells = cells_keep)

  query_mat <- seu_obj[["RNA"]]$counts

  print("Annotated spatial data")
  insitutype_res <- InSituType :: insitutypeML(x = t(query_mat),
                                 neg = colMeans(seu_obj[["ControlProbe"]]$counts),
                                 reference_profiles = ref_mat)

  seu_obj$celltype_pred <- insitutype_res$clust
  return(seu_obj)
}

######## QC Metrics ########


#' @title getNcells
#' @description
#' Calculate the number of cells present in a given dataset, supporting both Seurat object input(seu_obj)
#' and direct matrix file paths for imaging based spatial transcriptomic platforms such as Xenium and CosMx.
#' @details
#' The function `getNcells` offers flexibility in determining cell counts, catering to different types of input. For `Xenium` platform datasets, it expects the path to a 'matrix.mtx.gz'
#' file within the specified `expMat` directory and calculates the number of cells by reading the column count of this sparse matrix. For `CosMx` datasets, it reads a CSV file specified by `expMat`,
#' calculating cell count based on the number of rows. This approach allows for integration into workflows that may start with raw data files or pre-processed Seurat objects, facilitating versatile
#' data handling practices in transcriptomic analyses.
#' @param seu_obj Optional; a Seurat object from which to calculate the number of cells. If not provided, `expMat` must be specified.
#' @param expMat A character string specifying the path to the expression matrix if `seu_obj` is not provided. This parameter is used in conjunction with `platform` to determine the method for reading the data.
#' @param platform Optional; specifies the platform ('Xenium' or 'CosMx') if `expMat` is used. Required if `expMat` is provided without `seu_obj`.
#' @return An integer indicating the total number of cells in the dataset.
#' @importFrom Matrix readMM
#' @importFrom data.table fread
#' @export


getNcells <- function(seu_obj = NULL, expMat = 'path_to_expMat', platform = NULL) {
  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {
        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        ncell <- ncol(Matrix::readMM(file.path(mtx_path)))
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
        ncell <- ncol(exp)

      }
      res <- ncell
    }
    if(platform == 'CosMx') {
      ncell <- nrow(data.table::fread(expMat))
    }

    res <- ncell

  }


  if(!is.null(seu_obj)) {
    ncell <- ncol(seu_obj)
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=ncell
    )
  }


  return(res)

}



#' @details
#' Computes the global FDR for specified features (or all features by default) in a Seurat object.
#' It leverages transcript localization data to calculate FDR based on the proportion of negative control
#' or blank barcode expressions compared to the expression of actual genes.
#' @title getGlobalFDR
#' @description
#' Calculate specificity as Global False Discovery Rate (FDR).
#' @param seu_obj A Seurat object containing spatial data, including a path and platform attribute.
#' @param features An optional vector of feature names (gene names) for which to calculate the FDR.
#'        If NULL (default), FDR is calculated for all features in the seu_obj.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param path Path to generate tx file.
#' @param tx_file Path to transcription file.
#' @export
#' @importFrom data.table fread .N
#' @importFrom Seurat CreateSeuratObject
#' @return A data frame with columns for sample_id, platform, and the
#' calculated mean FDR across the specified features.
getGlobalFDR <- function(seu_obj = NULL, features = NULL, tx_file ='path_to_txFile', cellSegMeta = 'path_to_cellMeta', platform = NULL, path) {
  # Initialize variable
  tx_df <- NULL

  # Load data based on the presence of the Seurat object
  if(is.null(seu_obj)) {
    tx_df <- data.table::fread(tx_file)

    if(platform == 'Xenium') {
      setnames(tx_df, "feature_name", "target")  # changing the colname to target to keep it consistent for all techs
    }
  } else {
    # If Seurat object is provided, extract the necessary information from it
    path <- unique(seu_obj$path)
    platform <- unique(seu_obj$platform)
    sample_id <- unique(seu_obj$sample_id)

    # Read Tx localization data
    tx_df <- readTxMeta(path, platform)
  }

  # Filter and process data
  negProbes <- tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*', tx_df$target)]
  allGenes <- unique(tx_df[!target %in% negProbes, target])  # list of unique genes (non control or blank probes) in panel

  # Create table with expression per each gene in panel
  expTableAll <- tx_df[, .(Count = .N), by = target]
  expTable <- expTableAll[expTableAll$target %in% allGenes, ]
  expNeg <- sum(expTableAll[!expTableAll$target %in% allGenes, ]$Count)  # sum of all negative control or blank or unassigned barcodes (i.e non specific)

  numGenes <- length(expTable$target)
  numNeg <- length(expTableAll[!expTableAll$target %in% allGenes, ]$target)

  expTable$FDR <- 1
  for(i in allGenes) {
    fdr = (expNeg / (expTable[expTable$target %in% i, ]$Count + expNeg)) * (numGenes / numNeg) * 1/100
    expTable[target == i, FDR := fdr]
  }

  # Decide what value to return
  if(is.null(features)) {
    mean_fdr <- mean(expTable$FDR)
  } else {
    mean_fdr <- mean(expTable$FDR[expTable$target %in% features])
  }

  # Return results as a dataframe
  res <- data.frame(
    sample_id = sample_id,
    platform = platform,
    value = mean_fdr
  )

  return(res)
}


#' This function calculates the mean number of transcripts per unit area for specified features (genes) in a Seurat object.
#' If no features are specified, the function defaults to using all targets within the RNA assay of the provided Seurat object.
#' This calculation provides insight into the density of transcriptional activity relative to cell area, offering a normalized
#' measure of gene expression that accounts for differences in cell size.
#' @title getTxPerArea
#' @description
#' Calculate the Mean Transcripts per Unit Area.
#' @param seu_obj A Seurat object.
#'        The object must have 'sample_id' and 'platform' metadata for identification and reporting. It is expected that
#'        `seu_obj$cell_area` contains the area information for each cell.
#' @param features An optional vector of feature names (e.g., gene symbols) for which to calculate transcripts per unit area.
#'        Defaults to NULL, which means the calculation uses all available features in the RNA assay of the Seurat object.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param tx_file Path to transcription file.
#' @param cellSegMeta Path to CellsegMeta file.
#' @return Returns a data frame with three columns: `sample_id`, `platform`, and `value`. The `value` column contains
#'         the computed mean number of transcripts per unit area for the selected features across all cells in the dataset.
#'         This data frame provides a concise summary of transcriptional activity normalized by cell size, which can be
#'         critical for downstream analyses, especially in studies where cell morphology and size are variable.
#' @import Seurat
#' @export
getTxPerArea <- function(seu_obj = NULL, features=NULL,
                         platform, cellSegMeta = 'path_to_cellMeta', tx_file = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      cell_meta <- data.table::fread(file.path(cellSegMeta))

      mean_tx_norm <- mean(cell_meta$total_counts / cell_meta$cell_area)

      # If features are specified - have to calculate differently
      if(!is.null(features)) {
        ## Have to read the transcripts file
        tx_df <- data.table::fread(tx_file)
        # subset the tx file with only existing assigned cells and features
        tx_df <- tx_df[tx_df$cell_id %in% cell_meta$cell_id,]
        tx_df <- tx_df[tx_df$feature_name %in% features,]
        # count number of features per cell - have to merge cell_meta and tx_df - to add the area information
        tx_df <- merge(tx_df, cell_meta, by = 'cell_id')
        tx_df <- tx_df %>% group_by(cell_id,cell_area) %>% tally()

        mean_tx_norm <- mean(tx_df$n / tx_df$cell_area)
      }

    }

    if(platform == 'CosMx') {
      cell_meta <- data.table::fread(file.path(cellSegMeta))
      cell_meta <- data.frame(cell_meta)

      if("nCount_RNA" %in% names(cell_meta) == FALSE) {
        # Creating nCount_RNA column
        nCount_RNA <- cell_meta %>%
          group_by(cell_ID) %>%
          tally() %>%
          data.frame()
        names(nCount_RNA)[2] <- "nCount_RNA"

        cell_meta <- merge(cell_meta, nCount_RNA, by="cell_ID")
      }

      if("Area.um2" %in% names(cell_meta) == TRUE) {
        names(cell_meta)[names(cell_meta$Area) %in% "Area.um2"] <- "Area"
      }

      #mean_tx_norm <- mean(cell_meta$nCount_RNA / cell_meta$Area.um2)
      mean_tx_norm <- mean(cell_meta$nCount_RNA / cell_meta$Area)

      if(!is.null(features)) {
        ## Have to read the transcripts file
        tx_df <- data.table::fread(tx_file)
        # subset the tx file with only existing assigned cells and features
        tx_df <- tx_df[cell_ID != 0 & target %in% features]
        # count number of features per cell - have to merge cell_meta and tx_df - to add the area information
        tx_df <- merge(tx_df, cell_meta, by = 'cell')
        tx_df <- tx_df %>% group_by(cell,Area.um2) %>% tally()

        mean_tx_norm <- mean(tx_df$n / tx_df$Area.um2)
      }

    }

    if(platform == 'Merscope') {

    }

    return(mean_tx_norm)

  }

  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    tx_count <- colSums(seu_obj[["RNA"]]$counts[features,])
    mean_tx_norm <- mean(tx_count / seu_obj$cell_area)
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=mean_tx_norm
    )
    return(res)
  }


}


#' @title getTxPerCell
#' @description
#' Calculate the average number of transcripts per cell for the given features.
#' @details
#' This function calculates the mean number of transcripts per cell for a specified set of features (genes)
#' in a Seurat object. If no features are specified, the function defaults to using all targets available
#' within the RNA assay of the provided Seurat object. It's a useful metric for assessing the overall
#' transcriptional activity within the sampled cells.
#' @param seu_obj A Seurat object with RNA assays. The object must have 'sample_id' and 'platform' metadata attributes for identification
#'        and reporting purposes.
#' @param features An optional vector of feature names (e.g., gene symbols) to include in the calculation.
#'        Defaults to NULL, in which case the calculation uses all available features in the RNA assay
#'        of the Seurat object.
#' @param expMat Path to exprMatrix file.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return A data frame with columns 'sample_id', 'platform', and 'value', where 'value' represents the
#'         mean number of transcripts per cell calculated across the specified features. This output can
#'         be useful for comparative analysis across samples or experimental conditions.
#' @export
#' @import Seurat
getTxPerCell <- function(seu_obj = NULL, features=NULL, expMat = 'path_to_exprMatrix',
                         platform){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {

        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]

        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      } else{
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells

    }

    if(platform == 'Merscope') {

    }

    if(is.null(features)) {
      features <- rownames(exp)
      # remove non specific probes
      # features <- features[-grep('Unassigned*|NegControl*|BLANK*|SystemControl*', features)]
      features <- features[-grep('Unassigned*|NegControl*|BLANK*|SystemControl*|NegPrb*', features)]



    }

    # Calculate average N of Txs per cell
    mean_tx <- mean(colSums(exp[features, ]))
    res <- mean_tx

  }

  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    mean_tx <- mean(colSums(seu_obj[["RNA"]]$counts[features,]))
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=mean_tx
    )

  }

  return(res)

}


#' @title getTxPerNuc
#' @description
#' Calculate the number of transcripts per nucleus.
#' @details
#' This function calculates the mean number of transcripts localized within the nucleus
#' for a specified set of features (genes) across all cells in a given dataset. This measurement
#' is specific to spatial transcriptomics data where the distinction between nuclear and cytoplasmic
#' transcript localization can be made. The function is designed to work with different platforms by
#' reading transcript localization data using a helper function `readTxMeta`. If no features are
#' specified, the function defaults to using all available genes.
#' @param seu_obj A Seurat object containing spatial transcriptomics data with metadata attributes
#'        including 'path' and 'platform' to identify the location and type of data. The object is
#'        expected to be annotated with sample identifiers and must have been preprocessed to include
#'        information on cellular localization of transcripts.
#' @param features An optional vector of feature names (genes) for which to calculate the mean transcripts
#'        per nucleus. If NULL (default), the function calculates this metric for all genes.
#' @param tx_file Path to tx file.
#' @param platform Platform.
#' @return A data frame with the columns: `sample_id`, `platform`, and `value`, where `value` represents
#'         the mean number of transcripts per nucleus across the specified features for the dataset.
#'         This provides a focused view of nuclear gene expression.
#' @export
#' @import Seurat
#' @importFrom dplyr filter group_by summarize

getTxPerNuc <- function(seu_obj=NULL, features=NULL, tx_file = NULL, platform = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      tx_df <- data.table::fread(tx_file)

      if(is.null(features)) {
        # subset the tx file with only existing assigned cells and features
        # remove neg control probes
        negProbes <- unique(tx_df$feature_name[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$feature_name)])
        # number of txs in nucleus
        nTx_nuc <- dim(tx_df[cell_id != 'UNASSIGNED' & overlaps_nucleus == 1 & !feature_name %in% negProbes])[1]
        # number of cells
        nCells <- length(unique(tx_df$cell_id))
      }
      if(!is.null(features)) {
        nTx_nuc <- dim(tx_df[cell_id != 'UNASSIGNED' & overlaps_nucleus == 1 & !feature_name %in% negProbes & feature_name %in% features])[1]

      }


    }

    if(platform == 'CosMx') {
      tx_df <- data.table::fread(tx_file)

      if(is.null(features)) {
        # subset the tx file with only existing assigned cells and features
        # remove neg control probes
        negProbes <- unique(tx_df$target[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tx_df$target)])
        # number of txs in nucleus
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & CellComp == 'Nuclear' & !target %in% negProbes])
        # number of cells
        nCells <- length(unique(tx_df$cell[tx_df$cell_ID != 0]))
      }
      if(!is.null(features)) {
        nTx_nuc <- nrow(tx_df[cell_ID != 0 & CellComp == 'Nuclear' & !target %in% negProbes & target %in% features])

      }

    }

    if(platform == 'Merscope') {

    }

    return(nTx_nuc / nCells)

  }

  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  path <- unique(seu_obj$path)
  platform <- unique(seu_obj$platform)

  # Read Tx localization data
  tx_df <- readTxMeta(path, platform)

  if(platform == "Xenium"){
    tx_df <- filter(tx_df, cell_id %in% colnames(seu_obj) &
                      overlaps_nucleus == 1 &
                      features %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())


  } else if(platform == "CosMx"){
    tx_df$cell_id <- paste(tx_df$cell_ID, tx_df$fov, sep="_")
    tx_df <- tx_df %>%
      filter(cell_id %in% colnames(seu_obj) &
               CellComp == "Nuclear" &
               target %in% features) %>%
      group_by(cell_id) %>%
      summarize(nuc_counts = n())

  } else if(platform == "Merscope"){
    print("Working on support")

  } else{
    print("Platform not supported")
  }

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=mean(tx_df$nuc_counts)
  )

  return(res)
}


#' @title getMeanExpression
#' @description
#' It calculated mean expression per probe.
#' @details
#' Computes the mean expression levels for a specified list of features (genes)
#' and control probes within a Seurat object. This function separately calculates the mean
#' expression for both target genes and control probes, then combines these results into a
#' single data frame. If no features are specified, it defaults to all genes in the RNA assay.
#' @param seu_obj A Seurat object.
#' This object must have an RNA assay for target genes and a ControlProbe assay for control probes.
#' @param features An optional vector of gene names for which to calculate mean expression.
#' If NULL, mean expression is calculated for all genes in the RNA assay.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param expMat Path to exprMatrix file.
#' @return A data frame with columns `target` (gene or control probe name), `value` (mean expression level),
#' `type` (indicating whether the row represents a Gene or Control), `platform`, and `sample_id`.
#' @export
#' @import Seurat
getMeanExpression <- function(seu_obj = NULL, features=NULL, expMat = 'path_to_expMatrix',
                              platform=NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {

        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
        # remove neg control
        exp <- exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),]
      } else{
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
        exp <- exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),]

      }
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp)

    }

    if(platform == 'Merscope') {

    }

    avg_exp_df <- as.data.frame(rowMeans(exp))
    colnames(avg_exp_df) <- 'MeanExpression'
    avg_exp <- mean(avg_exp_df$MeanExpression)

    if(is.null(features)) {
      res <- avg_exp
    } else {
      avg_exp_df_sub <- subset(avg_exp_df, rownames(avg_exp_df) %in% features)
      res <- mean(avg_exp_df_sub$MeanExpression)
    }

  }



  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    target_df <- data.frame(
      target = features,
      value = rowMeans(seu_obj[["RNA"]]$counts[features,]),
      type = "Gene"
    )

    control_df <- data.frame(
      target = rownames(seu_obj[["ControlProbe"]]$counts),
      value = rowMeans(seu_obj[["ControlProbe"]]$counts),
      type = "Control"
    )

    res <- rbind(target_df, control_df)
    res$platform <- unique(seu_obj$platform)
    res$sample_id <- unique(seu_obj$sample_id)


  }
  return(res)
}
### log-ratio of mean gene counts to mean neg probe counts
#' @title getMeanSignalRatio
#' @description
#' Calcuates the log-ratio of mean gene expression counts to mean negative probe counts.
#' @details
#' Computes the log-ratio of the mean expression levels of specified genes (or all genes if none are specified)
#' to the mean expression levels of negative control probes within a Seurat object. This metric can provide insights into
#' the overall signal strength relative to background noise levels.
#' @param seu_obj A Seurat object, expected to have both 'RNA' and 'ControlProbe'
#' assays for calculating mean expression levels of genes and negative control probes, respectively.
#' @param features An optional vector of gene names for which to calculate the mean log-ratio. If NULL, the calculation
#' is performed for all genes present in the 'RNA' assay of the seu_obj.
#' @return A data frame with three columns: `sample_id`, `platform`, and `value`, where `value` represents the calculated
#' mean log-ratio of gene expression to negative control probe expression for the selected genes. This summary can be
#' used to assess the signal-to-noise ratio in the dataset.
#' @export

getMeanSignalRatio <- function(seu_obj=NULL, features=NULL, platform=NULL, expMat = 'path_to_expMat'){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {

        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    noise <- exp[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ]
    exp <- exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),]
    #ratio <- mean( log10(rowMeans(exp + .1)) - log10(rowMeans(noise + .1)) )

    if(is.null(features)) {
      #return( suppressWarnings(mean( log10(rowMeans(exp + .1)) - log10(rowMeans(noise + .1)) )))
      res <- suppressWarnings(mean( log10(rowMeans(exp + .1)) - log10(rowMeans(noise + .1)) ))
    } else {
      #return(suppressWarnings(mean( log10(rowMeans(exp[rownames(exp) %in% features,] + .1)) - log10(rowMeans(noise + .1)) )))
      res <- suppressWarnings(mean( log10(rowMeans(exp[rownames(exp) %in% features,] + .1)) - log10(rowMeans(noise + .1)) ))
    }


  }

  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
    neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)

    ratio <- log10(tx_means) - log10(mean(neg_probe_means))
    ratio <- mean(ratio)

    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=ratio
    )


  }
  return(res)

}

#' @title getCellTxFraction
#' @description
#' Calculate fraction of transcripts in cells from spatial transcriptomic data.
#' @details
#' The function works by reading transcript metadata for the given sample and platform, then filtering for the specified features (if any).
#' For each platform, it identifies transcripts that are not assigned to any cell ('UNASSIGNED' for Xenium, 'None' for CosMx, etc.) and calculates the fraction of transcripts that are assigned to cells.
#' This metric is an important indicator of the quality of spatial transcriptomics data, reflecting the efficiency of transcript capture and cell assignment.
#' @param seu_obj A Seurat object containing spatial transcriptomics data, expected to include metadata fields for `path` and `platform` to facilitate reading transcript metadata with `readTxMeta`.
#' @param features An optional vector of gene identifiers for which to calculate the transcript fraction. If NULL, the calculation is performed using all available features in the dataset.
#' @param tx_file Path to transcription file.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param path path
#' @return A data frame containing the calculated fraction of transcripts in cells, along with `sample_id` and `platform` for context. This output provides insight into the efficiency of transcript capture and assignment within the given dataset.
#' @return A data frame with probe type (Gene or Control), the mean expression values, and the associated sample and platform identifiers.
#' @export
#' @import Seurat
getCellTxFraction <- function(seu_obj=NULL, features=NULL, tx_file = 'path_to_tx',
                              platform = NULL, path){

  if(is.null(seu_obj)) {

    tx_df <- data.table::fread(tx_file)
    total_tx_count <- nrow(tx_df)

    if(platform == 'Xenium') {
      if(is.null(features)) {
        unassigned_tx_count <- sum(tx_df$cell_id == 'UNASSIGNED')
      } else {
        tx_df <- tx_df[tx_df$feature_name %in% features,]
        total_tx_count <- nrow(tx_df)
        unassigned_tx_count <- sum(tx_df$cell_id == 'UNASSIGNED')

      }

    }

    if(platform == 'CosMx') {
      if(is.null(features)) {
        # unassigned_tx_count <- sum(tx_df$CellComp == '')
        unassigned_tx_count <- sum(tx_df$CellComp == 'None')
        total_tx_count <- nrow(tx_df)

      } else {
        tx_df <- tx_df[tx_df$target %in% features,]
        total_tx_count <- nrow(tx_df)
        # unassigned_tx_count <- sum(tx_df$CellComp == '')
        unassigned_tx_count <- sum(tx_df$CellComp == 'None')


      }

    }

    return( (total_tx_count - unassigned_tx_count) / total_tx_count )

  }


  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    path <- unique(seu_obj$path)
    platform <- unique(seu_obj$platform)

    tx_df <- readTxMeta(path, platform)

    if(platform == "Xenium"){
      #tx_df <- filter(tx_df, features %in% feature_name)
      tx_df <- tx_df[tx_df$feature_name %in% features, ]
      total_tx_count <- nrow(tx_df)
      unassigned_tx_count <- sum(tx_df$cell_id == "UNASSIGNED")

      cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count

    } else if(platform == "CosMx"){
      tx_df <- filter(tx_df, target %in% features)
      total_tx_count <- nrow(tx_df)
      unassigned_tx_count <- sum(tx_df$CellComp == "None")

      cell_tx_fraction <- (total_tx_count - unassigned_tx_count) / total_tx_count

    } else if(platform == "Merscope"){
      print("Working on support")

    } else{
      print("Platform not supported")
    }

    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=cell_tx_fraction
    )

    return(res)
  }



}




##### Dynamic Range
# Log-ratio of highest mean exp vs. mean noise
#' @title getMaxRatio
#' @description
#' Calculate the Log-ratio of highest mean expression vs. mean noise.
#' @details
#' The function identifies the maximum mean expression value among the specified features
#' (or all features if none are specified) and calculates its log-ratio to the mean expression value
#' of negative control probes. This log-ratio reflects the dynamic range of the dataset, indicating
#' the spread between the highest signal and background noise levels, which is critical for assessing
#' data quality and sensitivity of detection in spatial transcriptomics.
#' @param seu_obj A Seurat object including assays for target genes ('RNA')
#' and control probes ('ControlProbe').
#' @param features An optional vector of gene identifiers for which to perform the calculation. If NULL, the calculation
#' includes all genes in the RNA assay.
#' @param expMat Path to exprMatrix file.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return A data frame with `sample_id`, `platform`, and the calculated `value` of the log-ratio.
#' @export
#' @import Seurat
getMaxRatio <- function(seu_obj = NULL, features=NULL, expMat ='path_to_expMat', platform = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {


        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }
    tx_means <- rowMeans(exp[-grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)),])
    neg_probe_means <- rowMeans(exp[grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)), ])

    if(is.null(features)) {
      #return( log10(max(tx_means)) - log10(mean(neg_probe_means)) )
      res <- log10(max(tx_means)) - log10(mean(neg_probe_means))
    } else {
      #return( log10(max(tx_means[features])) - log10(mean(neg_probe_means)) )
      res <- log10(max(tx_means[features])) - log10(mean(neg_probe_means))
    }

  }


  if(!is.null(seu_obj)) {

    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features,])
    neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)

    ratio <- log10(max(tx_means)) - log10(mean(neg_probe_means))

    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=ratio
    )


  }
  return(res)
}

# Distribution of maximal values
#' @title getMaxDetection
#' @description
#' Calculate the distribution of maximal values.
#' @details This function identifies the maximal expression values across the specified set of features
#' (or all features if none are specified) within the dataset, illustrating the upper bounds of gene
#' expression. Such information is crucial for assessing the dataset's dynamic range and the sensitivity
#' of detection methods used in the experiment. The function aggregates these maximal values and presents
#' them in a data frame, facilitating further analysis of the expression distribution and detection
#' efficiency across different genes.
#'
#' @param seu_obj A Seurat object.
#' @param expMat Path to exprMatrix file.
#' @param features An optional vector of gene identifiers for which to analyze the distribution of maximal
#' expression values. If NULL, the calculation encompasses all genes in the dataset.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return A data frame summarizing the maximal expression values across the specified features or the
#' entire dataset, which can be used to analyze the upper limits of detection and expression within the
#' sample.
#' @return A data frame.
#' @export
#' @import Seurat

getMaxDetection <- function(seu_obj = NULL, features=NULL, expMat ='path_to_expMat', platform = NULL){

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {

        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }


    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp)
    }

    max_vals <- data.frame(matrixStats::rowMaxs(as.matrix(exp)))
    colnames(max_vals) <- 'MaxValue'
    max_val <- max(max_vals$MaxValue)

    if(is.null(features)) {
      res <- max_val
    } else {
      max_vals <- subset(max_vals, rownames(max_vals) %in% features)
      res <- max(max_vals$MaxValue)
    }

  }



  if(!is.null(seu_obj)) {
    if(is.null(features)){
      features <- rownames(seu_obj)
    } else{
      features <- features
    }

    max_vals <- matrixStats::rowMaxs(as.matrix(seu_obj[["RNA"]]$counts[features,]))

    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=max_vals,
      gene = features
    )

  }
  return(res)

}


##### Mutually Exclusive Co-expression Rate (MECR) Implementation
#' @title getMECR
#' @description
#' Calculate Mutually Exclusive Co-expression Rate (MECR) which is a metric used for determining specificity.
#' @details
#' Based on a predefined set of markers and their associated cell types, this function computes the rate
#' at which pairs of markers are expressed in mutually exclusive patterns within cells. The calculation is limited
#' to markers found in the dataset and can be adjusted to focus on a subset by specifying it. The result is a
#' single MECR value that quantifies the overall mutually exclusive co-expression pattern, rounded to three decimal places.
#' @param  seu_obj A seurat object.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param expMat Path to exprMatrix file.
#' @return A data frame containing the `sample_id`, `platform`, and the computed
#' MECR value. The MECR value is rounded to three decimal places and represents
#' the average mutually exclusive co-expression rate across the selected markers.
#' @export
#' @import Seurat
getMECR <- function(seu_obj=NULL, expMat = 'path_to_expMat', platform = NULL) {
  #This function comes from Hartman & Satija, bioRxiv, 2024
  #We are using a custom marker table. The original publication bases it on
  #scRNA-seq from matched tissue.
  marker_df <- data.frame(
    gene = c("EPCAM", "KRT19", "KRT8",
             "CD3E", "CD3D", "CD8A", "NKG7",
             "MS4A1", "CD79A",
             "PECAM1", "CLDN5", "VWF",
             "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
             "PDGFRA", "DPT", "COL1A1",
             "MYH11", "ACTG2"),
    cell_type = c("Epithelial", "Epithelial", "Epithelial",
                  "T", "T", "T", "T",
                  "B", "B",
                  "Endo", "Endo", "Endo",
                  "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                  "Fibro", "Fibro", "Fibro",
                  "Muscle", "Muscle")
  )
  rownames(marker_df) <- marker_df$gene

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {
        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }


    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp)
    }

    genes <- intersect(rownames(exp), rownames(marker_df))
    mtx <- as.matrix(exp[genes,])

  }

  if(!is.null(seu_obj)) {
    genes <- intersect(rownames(seu_obj), rownames(marker_df))
    mtx <- as.matrix(seu_obj[['RNA']]$counts[genes, ])
  }

  coexp.rates <- c()
  #print(paste0("Marker count: ", length(genes)))
  if (length(genes) > 25) { genes <- sample(genes, 25) }
  for (g1 in genes) {
    for (g2 in genes) {
      if ((g1 != g2) && (g1 > g2) && (marker_df[g1, "cell_type"] != marker_df[g2, "cell_type"])) {
        c1 <- mtx[g1, ]
        c2 <- mtx[g2, ]
        coexp.rates <- c(
          coexp.rates,
          sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0)) # >0 too liberal of an expression threshold?
      }
    }
  }



  if(is.null(seu_obj)) {
    return(round(mean(coexp.rates), digits=3))
  } else {
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=round(mean(coexp.rates), digits=3)
    )
    return(res)
  }

}
##### Distribution spatial autocorrelation

#' @title getMorans
#' @description
#' Calculate the Moran's I statistic for spatial autocorrelation of a dataset. Spatial autocorrelation is multi-directional and multi-dimensional, it has a value from -1 to 1.
#' @param  seu_obj A seurat object.
#' @param features A character vector specifying the features (genes or probes) to include in the analysis. If `NULL` (the default), all features in `seu_obj` are used.
#' @return A data frame.
#' @details
#' The function first processes gene-targeting probes by converting the input data into a `SingleCellExperiment` object, then into a `SpatialFeatureExperiment` object;
#' Spatial coordinates are specified and used to further process the data, including normalization and nearest neighbor graph construction. Subsequently, Moran's I is calculated for each feature using the `Voyager` package. The same process is repeated for control probes. The results for both gene-targeting and control probes are combined into a single data frame.
#' The spatial autocorrelation analysis provides Moran's I values for each feature, indicating the degree of spatial clustering. This is performed separately for gene-targeting probes and control probes, allowing for a comparison of spatial patterns in gene expression and background noise.
#' @return A data frame with Moran's I values for each feature (gene or control probe), along with sample ID, platform, gene/probe name, and type (Gene or Control). Each row corresponds to a feature, and columns include the sample identifier, platform used, Moran's I value, gene/probe name, and the type indicating whether it's a Gene-targeting probe or a Control probe.
#' @export
#' @import Seurat
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment counts
#' @importFrom SpatialExperiment toSpatialExperiment
#' @importFrom SpatialFeatureExperiment toSpatialFeatureExperiment
#' @importFrom Voyager runMoransI
#' @importFrom scater logNormCounts
#' @importFrom BiocParallel MulticoreParam

getMorans <- function(seu_obj,
                      features=NULL){
  #Requires SingleCellExperiment, SpatialFeatureExperiment, Voyager, scater

  if(is.null(features)){
    features <- rownames(seu_obj)
  } else{
    features <- features
  }

  #First run for gene-targeting probes
  print("Getting Moran's I for gene-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["RNA"]]$counts[features,]),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)

  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh",
                                                 dist_type = "idw", k = 20,
                                                 style = "W")

  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))

  spatial_cor <- as.data.frame(rowData(sfe))

  targeting <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Gene"
  )

  #Now run for control probes
  print("Getting Moran's I for non-targeting probes")
  sce <- SingleCellExperiment(list(counts=seu_obj[["ControlProbe"]]$counts),
                              colData = seu_obj@meta.data)
  colData(sce) <- cbind(colData(sce), Embeddings(seu_obj, 'tissue'))
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("Tissue_1",
                                                         "Tissue_2"))
  sfe <- toSpatialFeatureExperiment(spe)
  sfe <- sfe[, colSums(counts(sfe)) > 0]
  rowData(sfe)$means <- rowMeans(counts(sfe))
  rowData(sfe)$vars <- rowVars(counts(sfe))
  sfe <- scater::logNormCounts(sfe)

  #Nearest neighbor
  colGraph(sfe, "knn20") <- findSpatialNeighbors(sfe, method = "knearneigh",
                                                 dist_type = "idw", k = 20,
                                                 style = "W")
  #Moran's I
  sfe <- Voyager::runMoransI(sfe, colGraphName = "knn20", BPPARAM = MulticoreParam(8))

  spatial_cor <- as.data.frame(rowData(sfe))

  control <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=spatial_cor[,3], #morans I
    gene = rownames(spatial_cor),
    type = "Control"
  )

  res <- rbind(targeting, control)

  return(res)
}

##### Cluster evaluation: silhouette width
#' @title getSilhouetteWidth
#' @description
#' Calculate the silhouette width which provides a metric for evaluating consistency within data clusters.
#' @details The silhouette width calculation involves preprocessing steps including normalization and scaling of the data, PCA for dimensionality reduction, and clustering.
#' The silhouette width is then calculated for a downsampled subset of cells to ensure computational efficiency.
#' This metric helps in assessing the cohesion and separation of the identified clusters, with higher values indicating better defined clusters.
#' @param seu_obj A Seurat object.
#' @return A data frame containing the `sample_id`, `platform`,
#' and the average silhouette width across all clusters, rounded to three decimal places.
#' @export
#' @import Seurat
#' @importFrom bluster approxSilhouette
getSilhouetteWidth <- function(seu_obj){
  print("Clustering data")
  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData()
  VariableFeatures(seu_obj) <- rownames(seu_obj)
  seu_obj <- seu_obj %>%
    RunPCA(verbose=F) %>%
    FindNeighbors(dims=1:10) %>%
    FindClusters(resolution=0.2)

  #Downsample to 100 cells per cluster for silhouette calculation
  seu_obj <- subset(seu_obj, downsample = 10000)

  silhouette <- bluster::approxSilhouette(
    Embeddings(seu_obj, 'pca')[,1:10],
    clusters = seu_obj$seurat_clusters
  )

  silhouette <- as.data.frame(silhouette)

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(mean(silhouette$width), digits=3)
  )

}

##### Sparsity calculation #####
#Show the sparsity (as a count or proportion) of a matrix.
#For example, .99 sparsity means 99% of the values are zero. Similarly, a sparsity of 0 means the matrix is fully dense.
#' @title  getSparsity
#' @description
#' Calculate the sparsity (as a count or proportion) of a gene expression count matrix.
#' @details
#' For example, .99 sparsity means 99% of the values are zero. Similarly, a sparsity of 0 means the matrix is fully dense.
#' @param seu_obj A Seurat object.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param expMat Path to exprMatrix file.
#' @param features An optional vector of feature names (e.g., gene symbols) to include in the calculation.
#'        Defaults to NULL, in which case the calculation uses all available features in the RNA assay
#'        of the Seurat object.
#' @return A data frame.
#' @export
#' @importFrom BioQC entropy
#' @import Seurat
#' @importFrom coop sparsity
getSparsity <- function(seu_obj=NULL, features = NULL, expMat = 'path_to_expMat', platform = NULL) {

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {

        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      }  else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
    }

    if(is.null(features)) {
      #return(coop::sparsity(as.matrix(exp)))
      res <- coop::sparsity(as.matrix(exp))
    } else {
      #return(coop::sparsity(as.matrix(exp[rownames(exp) %in% features, ])))
      res <- coop::sparsity(as.matrix(exp[rownames(exp) %in% features, ]))
    }
  }

  if(!is.null(seu_obj)) {
    value = coop::sparsity(as.matrix(seu_obj@assays$RNA$counts))
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=round(value, digits=3)
    )

  }
  return(res)

}


#' @title getEntropy
#' @description
#' Calculate the entropy of the RNA count matrix within a Seurat object (seu_obj).
#' @details
#' Entropy is used to quantify the diversity or uniformity of gene expression across the dataset. This function calculates the entropy of the RNA count matrix, reflecting the distribution of gene expression levels.
#' A higher entropy value suggests a more uniform distribution across genes, while a lower value indicates concentration in a smaller number of genes.
#' This measure can provide insights into the complexity of the cellular composition and the heterogeneity within the sample.
#' @param seu_obj A Seurat object containing RNA count data.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @param expMat Path to exprMatrix file.
#' @param features An optional vector of feature names (e.g., gene symbols) to include in the calculation.
#'        Defaults to NULL, in which case the calculation uses all available features in the RNA assay
#'        of the Seurat object.
#' @return A data frame with the `sample_id`, `platform`, and the entropy value of the RNA count matrix, rounded to three decimal places.
#' @export
#' @importFrom BioQC entropy
getEntropy <- function(seu_obj=NULL, features = NULL, expMat = 'path_to_expMat', platform = NULL) {

  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {

      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {

        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
        bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]


        exp <- Matrix::readMM(file.path(mtx_path))
        cols <- data.table::fread(file.path(bar_path), header = F)
        rows <- data.table::fread(file.path(feat_path), header = F)
        rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
        colnames(exp) <- cols$V1 ## this is the barcodes of cells
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        exp <- do.call(rbind, exp) #bind all the lists
      }
    }

    if(platform == 'CosMx') {
      exp <- data.table::fread(file.path(expMat))
      remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
      exp <- data.frame(exp)
      exp <- exp[, !colnames(exp) %in% remove_cols]
      ## remove first 2 columns - usually FOV and Cell_ID information
      #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
      exp <- t(exp)
    }

    if(is.null(features)) {

      #return(BioQC::entropy(as.matrix(exp)))
      res <- BioQC::entropy(as.matrix(exp))
    } else {
      #return(BioQC::entropy(as.matrix(exp[rownames(exp) %in% features, ])))
      res <- BioQC::entropy(as.matrix(exp[rownames(exp) %in% features, ]))
    }
  }



  if(!is.null(seu_obj)) {
    value = BioQC::entropy(as.matrix(seu_obj@assays$RNA$counts))
    res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      value=round(value, digits=3)
    )

  }
  return(res)
}


### ref
#' @title getRCTD
#' @description
#' Performs cell type deconvolution on spatial transcriptomics data. It requires a Seurat object with spatial transcriptomics data (seu_obj) and a reference Seurat object (ref).
#' It filters out cells with fewer than 10 transcripts to improve accuracy and uses the Robust Cell Type Deconvolution (RCTD) metric for cell type deconvolution.
#' @details
#' The function first prepares the reference data by extracting counts, cell type annotations, and UMI counts. The reference dataset is then used to construct
#' a `Reference` object. For the query (spatial) data, cells with fewer than 10 transcript counts are filtered out, and tissue coordinates are prepared.
#' The function then runs RCTD, which uses spatial and reference data to deconvolve cell types present in the spatial dataset. The output is a modified RCTD object
#' with cell type weights normalized across each spatial location. Users should be aware that low transcript counts may affect the accuracy of deconvolution;
#' hence, filtering is a crucial step before analysis.
#' @param seu_obj A Seurat object containing spatial transcriptomics data.
#' @param ref A Seurat object containing single-cell RNA-seq data used as a reference for deconvolution.
#' @return A RCTD object containing the results of the deconvolution, including cell type weights for each spatial location.
#'
#' @importFrom dplyr filter
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix colSums
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom spacexr SpatialRNA
#' @importFrom spacexr Reference
#' @importFrom spacexr create.RCTD
#' @importFrom spacexr run.RCTD
#' @export
getRCTD <- function(seu_obj, ref){
  # Prep reference for RCTD
  counts <- ref[["RNA"]]$counts
  cluster <- ref$cell_type
  cluster <- factor(cluster)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)

  # Prep obj
  # NOTE: If tx counts are low, it will fail.
  # I've found that subseting to >10 counts tends to work.
  # Should decide if that's done here or in a separate function
  print("Filtering any cells with <10 tx counts")
  seu_obj <- subset(seu_obj, nCount_RNA > 10)
  counts <- seu_obj[["RNA"]]$counts
  coords <- GetTissueCoordinates(seu_obj)
  coords <- filter(coords, cell %in% colnames(seu_obj)) %>%
    column_to_rownames(var = "cell")
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))

  # Run RCTD in full mode
  print("Running RCTD")
  RCTD <- create.RCTD(query, reference, max_cores = 6,
                      UMI_min = 0, UMI_max = Inf, counts_MIN = 0,
                      UMI_min_sigma = 50)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  RCTD@results$weights <- RCTD@results$weights / rowSums(RCTD@results$weights)
  return(RCTD)
  #seu_obj$max_weight <- rowMaxs(RCTD@results$weights)

}

#' @title getMaxRCTD
#' @description
#' This function is designed to run Robust Cell Type Deconvolution (RCTD) on spatial transcriptomics data using a single-cell RNA-seq reference.
#' It filters cells based on transcript counts, prepares the data, runs RCTD.
#' From the RCTD metric, it then calculates the maximum cell type weight for each spatial location, indicating the most prevalent cell type. The results include sample identifiers, platforms, and the corresponding maximum weight values, rounded to three decimal places.
#' @details
#' Initially, the function prepares the reference data from a Seurat object by extracting RNA counts, cell types, and UMI counts to create a reference object for RCTD.
#' Spatial data is also preprocessed by filtering cells with low transcript counts and preparing coordinates. RCTD is then executed to deconvolve cell types.
#' Post RCTD, the function extracts the maximum weight for each cell type across spatial locations, which can be used to identify dominant cell types in specific areas.
#' The approach ensures a focus on significant cell type contributions, enhancing the interpretation of spatial transcriptomics data.
#' @param seu_obj A Seurat object containing spatial transcriptomics data.
#' @param ref A Seurat object containing single-cell RNA-seq data used as a reference for deconvolution.
#' @importFrom dplyr filter
#' @importFrom tibble column_to_rownames
#' @importFrom Matrix colSums
#' @import Seurat
#' @importFrom spacexr SpatialRNA Reference create.RCTD run.RCTD
getMaxRCTD <- function(seu_obj, ref){
  # Prep reference for RCTD
  counts <- ref[["RNA"]]$counts
  cluster <- ref$cell_type
  cluster <- factor(cluster)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)

  # Prep obj
  # NOTE: If tx counts are low, it will fail.
  # I've found that subseting to >10 counts tends to work.
  # Should decide if that's done here or in a separate function
  print("Filtering any cells with <10 tx counts")
  seu_obj <- subset(seu_obj, nCount_RNA > 10)
  counts <- seu_obj[["RNA"]]$counts
  coords <- GetTissueCoordinates(seu_obj)
  coords <- filter(coords, cell %in% colnames(seu_obj)) %>%
    column_to_rownames(var = "cell")
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))

  # Run RCTD in full mode
  print("Running RCTD")
  RCTD <- create.RCTD(query, reference, max_cores = 6,
                      UMI_min = 0, UMI_max = Inf, counts_MIN = 0,
                      UMI_min_sigma = 50)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  RCTD@results$weights <- RCTD@results$weights / rowSums(RCTD@results$weights)

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=round(rowMaxs(RCTD@results$weights), digits=3)
  )

  return(res)
}

#' @title getCorrelationExp
#' @description
#' Identifies common genes between a spatial transcriptomics dataset (seu_obj) and a single-cell RNA-seq reference dataset (ref).
#' It computes the mean expression levels of these genes in both datasets while including the mean signal of control probes from seu_obj.
#' This function is useful for preparing data for correlation analyses or visual comparison between datasets.
#' @details
#' The function begins by identifying genes present in both the spatial transcriptomics dataset and the single-cell reference.
#' For these common genes, it calculates the average expression levels in both datasets. Additionally, it computes the average signal
#' of control probes in the spatial dataset, which can serve as a background signal reference. The resulting data frame includes sample identifiers,
#' platforms, gene names, mean expression values in the sample and reference, and the negative control probe signal, facilitating easy comparison and correlation analysis.
#' @param seu_obj A Seurat object containing spatial transcriptomics data.
#' @param ref A Seurat object containing single-cell RNA-seq data used as a reference.
#' @return A data frame with columns for `sample_id`, `platform`, `gene`, `value_sample`, `value_ref`, and `neg_probe_signal`.
#' @export
getCorrelationExp <- function(seu_obj, ref){
  common_genes <- intersect(rownames(seu_obj), rownames(ref))
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    gene = common_genes,
    value_sample = rowMeans(seu_obj[["RNA"]]$counts[common_genes,]),
    value_ref = rowMeans(ref[["RNA"]]$counts[common_genes,]),
    neg_probe_signal = mean(seu_obj[["ControlProbe"]]$counts)
  )
  return(res)
}

#' @title getCorrelation
#' @description
#' Identifies common genes between a spatial transcriptomics dataset (seu_obj) and a single-cell RNA-seq reference dataset (ref).
#' From this objec  the Spearman correlation coefficient is computed based on the mean expression levels of these genes across the datasets.
#' The result is a single correlation coefficient value that quantifies the similarity between the two datasets at the expression level of shared genes.
#' @details
#' This function is particularly useful for assessing the overall concordance between spatial transcriptomics data and a reference single-cell RNA-seq dataset.
#' By focusing on common genes, it provides a metric of similarity that can help validate spatial data against established single-cell profiles or compare spatial datasets against different references.
#' @param seu_obj A Seurat object containing spatial transcriptomics data.
#' @param ref A Seurat object containing single-cell RNA-seq data used as a reference.
#' @return A data frame with columns for `sample_id`, `platform`, and the Spearman correlation coefficient `value`.
#' @export
getCorrelation <- function(seu_obj, ref){
  common_genes <- intersect(rownames(seu_obj), rownames(ref))
  cor_res <- cor(
    rowMeans(seu_obj[["RNA"]]$counts[common_genes,]),
    rowMeans(ref[["RNA"]]$counts[common_genes,]),
    method='spearman'
  )
  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=cor_res
  )
  return(res)
}

#' @title getCellTypeCor
#' @description
#' This function subsets both the spatial transcriptomics dataset (seu_obj) and the single-cell RNA-seq reference (ref) into cell type-specific groups.
#' From these cell-type specific group,
#' it computes the Spearman correlation coefficient for each cell type, providing a cell type-specific similarity metric between the datasets.
#' @details
#' Useful for detailed analysis of cell type-specific expression patterns, this function helps in understanding the concordance between spatial transcriptomics
#' data and single-cell RNA-seq references at the level of individual cell types. It requires that `seu_obj` includes cell type predictions and that `ref` contains
#' annotated cell types. The output is a list of correlation coefficients for each cell type, offering insights into which cell types show higher or lower similarity
#' between the two datasets.
#' @param seu_obj A Seurat object containing spatial transcriptomics data with cell type predictions (`celltype_pred`).
#' @param ref A Seurat object containing single-cell RNA-seq data with annotated cell types (`cell_type`).
#' @return A data frame with `cell_type` and the corresponding Spearman correlation coefficient for each cell type between `seu_obj` and
#' @export
getCellTypeCor <- function(seu_obj, ref){
  cell_types <- levels(ref$cell_type)

  cor_list <- list()
  for(i in 1:length(cell_types)){
    print(paste0("Correlating: ", cell_types[i]))
    cor_list[[i]] <- getCorrelation(
      subset(seu_obj, celltype_pred == cell_types[i]),
      subset(ref, cell_type == cell_types[i])
    )
  }
  cor_list <- do.call(rbind, cor_list)
  cor_list$cell_type <- cell_types
  return(cor_list)
}

#' @title getClusterMetrics
#' @description
#' Applies the clustering analysis to a Seurat object (seu_obj) to calculate the Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI)
#' between the resulting clusters and a specified cell type prediction.
#' ARI and NMI are metrics used to assess the quality of clustering and the agreement between predicted cell types and clustering outcomes.
#' @details
#' Clustering is performed on the `seu_obj` data, which is first normalized and scaled. PCA is then run, followed by neighbor finding
#' and cluster identification. The ARI and NMI are calculated to evaluate the clustering quality, with ARI measuring the similarity
#' between two data clusterings and NMI providing a normalized measure of the mutual dependence between the predicted cell types and
#' the clusters. These metrics offer insights into the coherence of cell type predictions and the overall quality of the clustering process.
#' @param seu_obj A Seurat object for which clustering metrics are to be calculated.
#' @param metadata_col The metadata column in `seu_obj` that contains cell type predictions for comparison against clustering results.
#'
#' @return A data frame with sample identifiers, platforms, and the calculated ARI and NMI metrics.
#' @importFrom Seurat NormalizeData ScaleData RunPCA FindNeighbors FindClusters
#' @importFrom dplyr data_frame
#' @importFrom bluster pairwiseRand
#' @importFrom NMI NMI
#' @export
getClusterMetrics <- function(seu_obj, metadata_col){
  print("Clustering data")
  seu_obj <- seu_obj %>%
    NormalizeData() %>%
    ScaleData()

  VariableFeatures(seu_obj) <- rownames(seu_obj)

  seu_obj <- seu_obj %>%
    RunPCA(verbose=F) %>%
    FindNeighbors(dims=1:20) %>%
    FindClusters(resolution=0.5)

  ari <- bluster::pairwiseRand(seu_obj$seurat_clusters,
                               seu_obj$celltype_pred,
                               mode=('index'),
                               adjusted=TRUE)

  nmi <- NMI(
    data.frame(cellid = colnames(seu_obj), cluster = seu_obj$seurat_clusters),
    data.frame(cellid = colnames(seu_obj), cluster = seu_obj$celltype_pred)
  )

  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=c(ari, nmi$value),
    metric = c("ARI", "NMI")
  )

}

#' @title getPanelSize
#' Get Panel Size from Expression Matrix or Seurat Object.
#' @description
#' Calculate the size of the panel by counting the number of features
#' in an expression matrix, excluding negative or control features. It can handle
#' different types of expression data input, including direct paths to expression
#' matrices or Seurat object with embedded expression data.
#' @param seu_obj An optional Seurat object containing the expression data. If provided,
#'        the function will calculate the panel size based on the number of features
#'        present in the object.
#' @param expMat The path to the expression matrix file or directory containing
#'        expression data files. This parameter is used if no Seurat object is provided.
#'        The function supports both flat files and Hierarchical Data Format (HDF5) files.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return Returns an integer indicating the number of valid features (excluding controls
#'         and system features) present in the expression data.
getPanelSize <- function(seu_obj = NULL, expMat = 'path_to_expMat', platform = NULL) {
  if(is.null(seu_obj)) {
    if(platform == 'Xenium') {
      expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      n <- length(expMat_checker)
      h5_file_path <- expMat_checker[1]
      for (i in 2:n) {
        one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
      }

      if(grepl('.h5', one_string_only) == FALSE) {
        mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]
        tmp <- fread(file.path(feat_path), header = F)
        nPanel <- nrow(tmp) - length(grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', tmp$V1))
      } else {
        mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
        mtx_path <- mtx_path[grepl('.h5', mtx_path)]
        exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
        tmp <- do.call(rbind, exp) #bind all the lists
        nPanel <-  nrow(tmp) - length(grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(tmp)))
      }
    }
    if(platform == 'CosMx') {
      tmp <- data.table::fread(expMat)
      nPanel <- ncol(tmp) - length(grep('Neg*|SystemControl*|Blank*|BLANK*', colnames(tmp))) - 2 ## - 2 because the first 2 columns are FOV and cell id
    }
  }

  if(!is.null(seu_obj)) {
    nPanel <- nrow(seu_obj)
  }

  return(nPanel)

}


#'@title getComplexity
#'@description
#'Calculate the complexity of expression data, determining the number
#' of features (genes) required to reach half of the total expression sum in a dataset.
#' @param seu_obj An optional Seurat object containing the expression data. If provided,
#'        the function will calculate complexity based on the data in the object.
#' @param features Optional character vector of features to focus the complexity calculation
#'        on a subset of all features. If NULL, all features are considered.
#' @param expMat The path to the expression matrix file or directory containing
#'        expression data files. This parameter is used if no Seurat object is provided.
#'        This function supports processing of '.mtx.gz' for sparse matrices, or '.h5'
#'        files for HDF5 format expression matrices.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @details
#' Calculate the total expression and determines the minimal set of features accounting for
#' at least half of this total, either across all features or a specified subset. If a Seurat
#' object is used, it derives the complexity from the RNA assay counts.
#' @export

getComplexity <- function(seu_obj=NULL, features = NULL, expMat = 'path_to_expMat', platform = NULL) {
if(is.null(seu_obj)) {
  if(platform == 'Xenium') {

    expMat_checker <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
    n <- length(expMat_checker)
    h5_file_path <- expMat_checker[1]
    for (i in 2:n) {
      one_string_only <- paste0(h5_file_path, expMat_checker[i], sep="")
    }

    if(grepl('.h5', one_string_only) == FALSE) {

      mtx_bar_feat_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      mtx_path <- mtx_bar_feat_path[grepl('matrix.mtx.gz', mtx_bar_feat_path)]
      bar_path <- mtx_bar_feat_path[grepl('barcodes.tsv.gz', mtx_bar_feat_path)]
      feat_path <- mtx_bar_feat_path[grepl('features.tsv.gz', mtx_bar_feat_path)]

      exp <- Matrix::readMM(file.path(mtx_path))
      cols <- data.table::fread(file.path(bar_path), header = F)
      rows <- data.table::fread(file.path(feat_path), header = F)
      rownames(exp) <- rows$V2 ## this is the gene symbol column of the dataframe rows
      colnames(exp) <- cols$V1 ## this is the barcodes of cells

    } else {
      mtx_path <- fs::dir_ls(expMat, recurse = TRUE, type = "file")
      mtx_path <- mtx_path[grepl('.h5', mtx_path)]
      exp <- Read10X_h5(filename = file.path(mtx_path), use.names = TRUE, unique.features = TRUE)
      exp <- do.call(rbind, exp) #bind all the lists
    }
  }

  if(platform == 'CosMx') {
    exp <- data.table::fread(file.path(expMat))
    remove_cols <- as.vector(c("V1","sampleID","slide","case", "fov", "cell_ID"))
    exp <- data.frame(exp)
    exp <- exp[, !colnames(exp) %in% remove_cols]
    ## remove first 2 columns - usually FOV and Cell_ID information
    #exp <- exp[, -c(1:2)] # "fov", "cell_ID" have been removed above
    exp <- t(exp) ## transposing for consistency  - row = genes, column= cells
  }

  total_sum <- sum(exp)


  if(is.null(features)) {

    row_sums <- rowSums(exp)
    cumulative_sums <- cumsum(row_sums)

    return(as.integer( which(cumulative_sums >= total_sum / 2)[1]))
    #cumulative_sums <- as.integer( which(cumulative_sums >= total_sum / 2)[1])


  }

  if(!is.null(features)) {
    total_sum_feat <- sum(exp[rownames(exp) %in% features, ])
    row_sums <- rowSums(as.data.frame(exp[rownames(exp) %in% features, ]))
    cumulative_sums <- cumsum(row_sums)
  }

  result <- which(cumulative_sums >= total_sum / 2)[1]
  weigh <- 31

  res <- as.integer(result) / (weigh * ( nrow(exp) - length(grep('Neg*|SystemControl*|Blank*|BLANK*|Unassigned*', rownames(exp)))) )
}

if(!is.null(seu_obj)) {

  total_sum <- sum(seu_obj@assays$RNA$counts)
  row_sums <- rowSums(seu_obj@assays$RNA$counts)
  cumulative_sums <- cumsum(row_sums)


  res <- data.frame(
    sample_id = unique(seu_obj$sample_id),
    platform = unique(seu_obj$platform),
    value=which(cumulative_sums >= total_sum / 2)[1]
  )

}

return(res)

}

######## Plotting ########

# All plots assume input is a tidy data frame with the following columns:
# 1) sample_id
# 2) platform
# 3) value (based on what is being plotted--from functions above)

#' @title plotSampleLabel
#' @description Performs Spearman correlation analysis for each cell type between spatial and reference datasets.
#' @param sample_meta sample_meta
#' @return Data frame with cell type-specific correlation values.
#' @import ggplot2
#' @export

plotSampleLabel <- function(sample_meta){
  df <- data.frame(
    samples = sample_meta$sample_id,
    platform = sample_meta$platform
  )

  p <- ggplot(df, aes(x="", y=samples)) +
    geom_text(aes(label=samples, color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) +
    theme_void() + theme(legend.position='none')
  return(p)
}

# plotPanelSize <- function(df){
#   p <- ggplot(df, aes(x="", y=sample_id)) +
#     geom_point(shape=21, color='black', alpha=0.8, stroke=1,
#                aes(size=value, fill=value)) +
#     geom_shadowtext(color = "black", size = 4, #fontface = "bold",
#                     bg.colour = "white", bg.r = .2,
#                     aes(label=scales::comma(value))) +
#     scale_fill_gradientn(colours=viridis::mako(100)) +
#     xlab("Panel size") + ylab("") +
#     scale_size(range = c(6,12)) +
#     scale_x_discrete(position='top',
#                      labels = c("")) +
#     theme_classic() +
#     theme(
#       legend.position='none',
#       axis.text.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.title.x = element_text(size=12),
#       axis.line.y = element_blank(),
#       axis.ticks.y = element_blank()
#     )
#
#   return(p)

#}

#' @title plotCellCount
#' @description Computes and compares cell type proportions between spatial and reference datasets.
#' @param df Dataframe.
#' @return Plot.
#' @import ggplot2
#' @export

plotCellCount <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_text(aes(label=scales::comma(value), color=platform), size=5, hjust=0.5) +
    scale_color_manual(values = c("#59C134", "#14B3E6")) +
    scale_x_discrete(position='top') +
    xlab("Cell count") + ylab("") +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(p)
}

#' @title plotTxPerCell
#' @description Creates a scatter plot visualizing the number of transcripts per cell across samples.
#' @param df Data frame with `sample_id` and `value` columns indicating sample IDs and the number of transcripts per cell, respectively.
#' @return A ggplot object visualizing the number of transcripts per cell for each sample.
#' @export
#' @import ggplot2

plotTxPerCell <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("Per cell") + ylab("") +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}
#' @title plotTxPerArea
#' @description Generates a scatter plot to show the number of transcripts per unit area across samples.
#' @param df Data frame with `sample_id` and `value` columns indicating sample IDs and the number of transcripts per unit area, respectively.
#' @return A ggplot object depicting the distribution of transcripts per unit area for each sample.
#' @export
#' @import ggplot2
plotTxPerArea <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\num^2") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title Plot Transcripts Per Nucleus
#' @description Produces a scatter plot that displays the number of transcripts per nucleus for each sample.
#' @param df Data frame with `sample_id` and `value` columns for sample IDs and the number of transcripts per nucleus, respectively.
#' @return A ggplot object illustrating the transcripts count per nucleus across samples.
#' @export
#' @import ggplot2
plotTxPerNuc <- function(df){
  df$column <- ""
  p <- ggplot(df, aes(x=column, y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per\nnucleus") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotTxPerCellNorm
#' @description Creates a plot visualizing normalized transcripts per cell across samples.
#' @param df Data frame with `sample_id` and normalized `value` columns.
#' @return A ggplot object visualizing normalized transcripts per cell for each sample.
#' @export
#' @import ggplot2
plotTxPerCellNorm <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    xlab("Per cell\nper target") + ylab("") +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    scale_size(range = c(7,12)) +
    scale_x_discrete(position='top') +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotFractionTxInCell
#' @description Generates a bar plot showing the fraction of transcripts located within cells for each sample.
#' @param df Data frame with `sample_id` and `value` columns indicating the fraction of transcripts in cells.
#' @return A ggplot object depicting the fraction of transcripts within cells across samples.
#' @export
#' @import ggplot2
plotFractionTxInCell <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90', stroke=1) +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Fraction Tx in Cells") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0),
                       breaks=c(0, 0.5, 1),
                       labels = c(0, 0.5, 1),
                       limits=c(0,1)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#plotTxPerCell_Intersect <- function()
#' @title plotSignalRatio
#' @description Creates a bar plot for the signal-to-noise ratio of gene expression across samples.
#' @param df Data frame with `sample_id` and `value` columns indicating the mean log10-ratio of expression over noise.
#' @return A ggplot object visualizing the signal-to-noise ratio for each sample.
#' @export
#' @import ggplot2
plotSignalRatio <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -.05) +
    xlab("Mean log10-ratio\nexpression over noise") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}
#' @import ggplot2
#' @title plotMeanExpression
#' @description Generates a jitter and boxplot visualizing the mean gene detection rate per cell across samples and cell types.
#' @param df Data frame with `sample_id`, `value`, and `type` columns.
#' @return A ggplot object depicting mean gene detection rates per cell, differentiated by cell type.
#' @export

plotMeanExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black",
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n detection per cell") + ylab("") +
    scale_x_log10(position='top', expand = c(0,0),
                  labels = label_log(digits = 2)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotMaxExpression
#' @description Creates a jitter and boxplot showing the maximal gene detection rate per cell across samples and platforms.
#' @param df Data frame with `sample_id`, `value`, and `platform` columns.
#' @return A ggplot object illustrating maximal gene detection rates per cell, differentiated by platform.
#' @export
#' @import ggplot2
plotMaxExpression <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=platform)) +
    geom_boxplot(color="black", fill="lightgrey",
                 alpha=0.5, outlier.size=0, outlier.colour = NA) +
    scale_colour_manual(values = c("#59C134", "#14B3E6")) +
    xlab("Maximal gene\ndetection per cell") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0),
                       limits = c(0,50), oob=squish) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

#' @title plotMECR
#' @description
#' Visualizes the Mutually Exclusive Co-expression Rate (MECR) values across samples with a scatter plot.
#' @param df Data frame with `sample_id` and `value` columns for MECR values.
#' @return A ggplot object showing MECR values for each sample.
#' @export
#' @import ggplot2
plotMECR <- function(df){
  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=scales::comma(value))) +
    scale_fill_gradientn(colours=RColorBrewer::brewer.pal(9, "YlOrRd"),
                         limits=c(0,0.1)) +
    xlab("MECR") + ylab("") +
    scale_size(range = c(7,12), limits = c(0, 0.1)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=10),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotMorans
#' @description Creates a plot for Moran's I spatial autocorrelation values across samples.
#' @param df Data frame with `sample_id`, `value`, and `type` columns for Moran's I values and sample types.
#' @return A ggplot object depicting Moran's I spatial autocorrelation values, differentiated by sample type.
#' @export
#' @import ggplot2
plotMorans <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_jitter(size=0.15, shape=16, aes(color=type),
                position = position_jitterdodge(jitter.width=0.25)) +
    geom_boxplot(color="black",
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 aes(fill=type)) +
    scale_fill_manual(values=c("lightgrey", "firebrick")) +
    scale_colour_manual(values=c("grey20", "firebrick")) +
    xlab("Mean gene\n spatial autocorrelation\n(Moran's I)") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotSilhouette
#' @description Generates a bar plot showing mean silhouette width values for clustering quality across samples.
#' @param df Data frame with `sample_id` and `value` columns for mean silhouette width values.
#' @return A ggplot object visualizing mean silhouette width for each sample.
#' @export
#' @import ggplot2
plotSilhouette <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_col(color='black', fill='grey90') +
    geom_text(aes(label=scales::comma(value)),
              hjust = 1, nudge_x = -0.001) +
    xlab("Mean\nsilhouette width\n(Louvain res=0.2)") + ylab("") +
    scale_x_continuous(position='top',
                       expand = c(0,0)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotCorrelation
#' @description Creates scatter plots to visualize correlation between spatial and reference dataset expression values across samples.
#' @param df Data frame with `sample_id`, `value_sample`, and `value_ref` columns indicating sample IDs and expression values in the spatial and reference datasets, respectively.
#' @return A ggplot object showing correlation scatter plots for each sample, with log-transformed axes.
#' @export
#' @import ggplot2
plotCorrelation <- function(df){
  #facet_wrap order is opposite to axis text orders, so we'll flip the levels
  df$sample_id <- factor(df$sample_id)
  df$sample_id <- factor(df$sample_id, levels = rev(levels(df$sample_id)))
  p <- ggplot(df, aes(x=value_ref, y=value_sample)) +
    geom_point(size=0.5, shape=16, stroke=0, alpha=0.75) +
    geom_abline(intercept = 0, slope = 1, linetype=2) +
    scale_x_log10(labels = label_log(digits = 2), limits=c(1e-4, 10),
                  oob=squish, position='top') +
    scale_y_log10(labels = label_log(digits = 2), limits=c(1e-4, 10), oob=squish) +
    xlab("snPATHO-seq\ncorrelation") + ylab("") +
    facet_wrap(~sample_id, ncol=1) +
    theme_bw() +
    theme(axis.text = element_text(size=10, color="black"),
          legend.position="none",
          axis.title.x = element_text(size=12),
          strip.text = element_blank(),
          strip.background = element_blank())
  return(p)
}

#' @title plotCellTypeCor
#' @description Visualizes correlation for each cell type across samples using jitter plots.
#' @param df Data frame with `sample_id`, `value`, and `cell_type` columns.
#' @return A ggplot object showing cell type-specific correlation for each sample.
#' @export
#' @import ggplot2
plotCellTypeCor <- function(df){
  p <- ggplot(df, aes(x=value, y = sample_id)) +
    geom_jitter(shape=21, color='black', aes(fill=cell_type),
                alpha=0.5, size=3, height=0.1) +
    stat_summary(fun = mean, geom = "crossbar", width=0.5) +
    scale_x_continuous(position='top', limits=c(0, 1),
                       expand=c(0,0), breaks=c(0, 0.5, 1)) +
    xlab("Cell type\ncorrelation") + ylab("") +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 6.5, 5.5, 5.5, unit='pt')
    )
  return(p)
}

#' @title plotCellTypeProportion
#' @description Creates bar plots to show the proportion of each cell type within samples.
#' @param df Data frame with `group`, `prop`, `cell_type`, and `sample_id` columns.
#' @return A ggplot object visualizing the proportion of each cell type within and across samples.
#' @export
#' @import ggplot2
plotCellTypeProportion <- function(df){
  p <- ggplot(df, aes(x=prop, y=group)) +
    geom_bar(position="stack", stat="identity",
             color='black', alpha=0.75,
             aes(fill=cell_type), width=0.7) +
    scale_x_continuous(position='top',
                       limits=c(0,1), expand = c(0,0),
                       breaks=c(0, 0.5, 1)) +
    ylab("") + xlab("Cell type\nproportion") +
    facet_wrap(~sample_id, ncol=1) +
    theme_classic() +
    theme(
      panel.spacing = grid::unit(0, "lines"),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position='none',
      axis.text.y = element_text(size=10, color="black"),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 6.5, 5.5, 5.5,
                           unit = "pt")
    )
  return(p)
}

#' @title plotARI
#' @description
#' Generates a scatter plot visualizing Adjusted Rand Index (ARI) values across samples to assess clustering performance.
#' @param cluster_metrics Data frame with `sample_id`, `value`, and `metric` columns, filtered for "ARI".
#' @return A ggplot object depicting ARI values for each sample.
#' @export
#' @import ggplot2
plotARI <- function(cluster_metrics){
  df <- cluster_metrics %>% filter(metric == "ARI")

  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=round(value, digits=2))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("ARI") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotNMI
#' @description
#' Creates a scatter plot to visualize Normalized Mutual Information (NMI) values across samples, indicating the quality of clustering.
#' @param cluster_metrics Data frame with `sample_id`, `value`, and `metric` columns, filtered for "NMI".
#' @return A ggplot object showing NMI values for each sample.
#' @export
#' @import ggplot2
plotNMI <- function(cluster_metrics){
  df <- cluster_metrics %>% filter(metric == "NMI")

  p <- ggplot(df, aes(x="", y=sample_id)) +
    geom_point(shape=21, color='black', alpha=0.8, stroke=1,
               aes(size=value, fill=value)) +
    geom_shadowtext(color = "black", size = 4, #fontface = "bold",
                    bg.colour = "white", bg.r = .2,
                    aes(label=round(value, digits=2))) +
    scale_fill_gradientn(colours=viridis::mako(100)) +
    xlab("NMI") + ylab("") +
    scale_size(range = c(6,12)) +
    scale_x_discrete(position='top',
                     labels = c("")) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(p)
}

#' @title plotRCTD
#' @description Generates boxplots for the maximum decomposition weights obtained from Robust Cell Type Deconvolution (RCTD) across samples.
#' @param df Data frame with `sample_id` and `value` columns for maximum decomposition weights.
#' @return A ggplot object showing boxplots of max decomposition weights for each sample.
#' @export
#' @import ggplot2
plotRCTD <- function(df){
  p <- ggplot(df, aes(x=value, y=sample_id)) +
    geom_boxplot(color="black",
                 alpha=0.5, outlier.size=0, outlier.colour = NA,
                 fill="lightgrey", width=0.5) +
    xlab("Max decomposition\nweight") + ylab("") +
    scale_x_continuous(position='top', expand = c(0,0),
                       limits=c(0, 1), breaks=c(0, 0.5, 1),
                       oob=squish) +
    #scale_x_log10(position='top', expand = c(0,0),
    #              labels = label_log(digits = 2)) +
    theme_classic() +
    theme(
      legend.position='none',
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10, color="black"),
      axis.title.x = element_text(size=12),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(5.5, 6.5, 5.5, 5.5, unit='pt')
    )
  return(p)
}
#' @import ggplot2
#' @description
#' Generates plots from technical metrics in spatial touchstone.
#' @title plotMetrics
#' @param metrics_list A list containing the metrics results.
#' @param PlotAutocorr Logical; if `TRUE`, plots the autocorrelation metric. Default is `TRUE`.
#' @param PlotSparsity Logical; if `TRUE`, plots the sparsity metric;Default is `TRUE`.
#' @param PlotEntropy Logical; if `TRUE`, plots the entropy metric; Default is `TRUE`.
#' @param PlotClusterSilhouette Logical; if `TRUE`, plots the cluster silhouette metric; Default is `TRUE`.
#' @param ncol Integer; the number of columns in the plot grid; Default is 14.
#' @param nrow Integer; the number of rows in the plot grid; Default is 1.
#' @param rel_widths Numeric vector; relative widths of columns in the plot grid.
#' @return A plot object created with cowplot::plot_grid.
#' @importFrom cowplot plot_grid
#' @export
#' @import ggplot2
# metrics_list <- list(
#   sample_meta = sample_meta,
#   cell_count = getNcells_test_1,
#   tx_per_cell = getTxPerCell_test_1,
#   tx_per_um2 = getTxPerArea_test_1,
#   tx_per_nuc = getTxPerNuc_test_1,
#   tx_per_cell_norm = getTxPerCell_test_1, #
#   tx_fraction_in_cell = getCellTxFraction_test_1, #
#   signal_ratio = getMeanSignalRatio_test_1,
#   mean_expression = getMeanExpression_test_1,
#   sparsity = getSparsity_test_1,
#   entropy = getEntropy_test_1,
#   silhouette = getSilhouetteWidth_test,
#   morans = getMorans_test,
#   mecr =mecrtest
# )


plotMetrics <- function(metrics_list = NULL, PlotAutocorr = T, PlotSparsity = T,
                        PlotEntropy = T,
                        PlotClusterSilhouette = T,
                        ncol = 14, nrow = 1,
                        rel_widths = c(0.75, 0.4, 0.5, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,0.3,0.75, 0.75, 0.8, 0.3, 0.8, 0.8)) {
  # PLOT
  p0 <- plotSampleLabel(metrics_list[['sample_meta']])
  #p1 <- plotPanelSize(metrics_list[["panel_size"]])
  p2 <- plotCellCount(metrics_list[["cell_count"]])
  p3 <- plotTxPerCell(metrics_list[["tx_per_cell"]])
  p4 <- plotTxPerArea(metrics_list[["tx_per_um2"]])
  p5 <- plotTxPerNuc(metrics_list[["tx_per_nuc"]])
  p6 <- plotTxPerCellNorm(metrics_list[["tx_per_cell_norm"]])
  p7 <- plotFractionTxInCell(metrics_list[["tx_fraction_in_cell"]])
  p8 <- plotSignalRatio(metrics_list[["signal_ratio"]])
  p9 <- plotMeanExpression(metrics_list[["mean_expression"]])
  p10 <- plotMECR(metrics_list[["mecr"]])
  if(PlotAutocorr == TRUE) {
    p11 <- plotMorans(metrics_list[["morans"]])
  }
  if(PlotSparsity == TRUE) {
    p13 <- plotSparsity(metrics_list[["sparsity"]])
  }
  if(PlotEntropy == TRUE) {
    p14 <- plotEntropy(metrics_list[["entropy"]])
  }
  if(PlotClusterSilhouette == TRUE) {
    p12 <- plotSilhouette(metrics_list[["silhouette"]])
  }
  p <- cowplot::plot_grid(p0, #p1,
                          p2,
                          p3, p4,
                          p5, p6,  p7, p8, p9,
                          p10,
                          p11, p12,p13, p14,
                          ncol=ncol, align='h',
                          rel_widths = rel_widths , nrow = nrow)

  return(p)
}

## QCreport Genereate

#'@description
#' Generate QC Report Table for Seurat Object
#' @details
#' Generate a QC report for a given Seurat object, allowing users to specify features and
#' output to a PDF file. It provides various statistical analyses including entropy, sparsity,
#' transcript counts per area, transcript counts per cell, and signal ratios. It also adjusts its behavior
#' based on the specified platform within the Seurat object (Xenium, CosMx, or Merscope).
#' @title genereateQCreport_table
#' @param seu_obj A Seurat object, which must have properties like `path`, `platform`, and specific assay data.
#'                Must inherit from "Seurat".
#' @param features Optional; A character vector of features to include in the analysis.
#'                 Defaults to all features in `seu_obj` if `NULL`.
#' @param pdfFile A character string specifying the path and name of the output PDF file.
#'                Defaults to "QCReportTable.pdf".
#' @param sample_id Identifier for the sample being processed.
#' @param path The file path from which to load the spatial transcriptomics data.
#' @param platform The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
#'        and 'Merscope'. Note: 'Merscope' is currently not supported.
#' @return A data frame with the QC metrics for the provided Seurat object, which is also printed to the console.
#'         Additionally, a PDF file is generated containing a plot and table of the metrics.
#'
#' @details Calculate various metrics such as entropy and sparsity using the `BioQC` and `coop` packages.
#'          It adjusts its processing logic based on the 'platform' attribute of the Seurat object to handle data from
#'          different technologies like Xenium, CosMx, or Merscope. It stops with an error if the object is not a Seurat object.
#'          The results are plotted using `ggplot2` and `gridExtra` for layout adjustments.
#' @import dplyr
#' @import ggplot2
#' @import gridExtra
#' @importFrom BioQC entropy
#' @importFrom coop sparsity
#' @export
generateQCreport_table <- function(seu_obj=NULL, features=NULL, pdfFile="QCReportTable.pdf", sample_id=NULL, path=NULL, platform=NULL) {
  if (!inherits(seu_obj, "Seurat") && !is.null(seu_obj)) {
    stop("Input must be a Seurat object or parameters must be provided to load one.")
  }

  # Check if a Seurat object needs to be loaded
  if (is.null(seu_obj)) {
    if (is.null(sample_id) || is.null(path) || is.null(platform)) {
      stop("Missing parameters to load spatial data. Please provide sample_id, path, and platform.")
    }
    seu_obj <- readSpatial(sample_id=sample_id, path=path, platform=platform, seurat=TRUE)
  }

  if (is.null(features)) {
    features <- rownames(seu_obj)
  }

  # Extract metrics
  entropy_value <- BioQC::entropy(as.matrix(seu_obj@assays$RNA$counts))
  sparsity_value <- coop::sparsity(as.matrix(seu_obj@assays$RNA$counts))
  ncell <- ncol(seu_obj)

  tx_means <- rowMeans(seu_obj[["RNA"]]$counts[features, ])
  neg_probe_means <- rowMeans(seu_obj[["ControlProbe"]]$counts)
  ratio <- log10(tx_means) - log10(mean(neg_probe_means))
  mean_signal_ratio <- mean(ratio)

  tx_count <- colSums(seu_obj[["RNA"]]$counts[features, ])
  mean_tx_norm <- mean(tx_count / seu_obj$cell_area)
  tx_perarea <- mean_tx_norm

  mean_tx <- mean(colSums(seu_obj[["RNA"]]$counts[features, ]))
  tx_percell <- mean_tx

  max_ratio <- log10(max(tx_means)) - log10(mean(neg_probe_means))

  # prepare result table
  res <- data.frame(
      sample_id = unique(seu_obj$sample_id),
      platform = unique(seu_obj$platform),
      ncell = ncell,
      entropy_value = round(entropy_value, digits = 3),
      sparsity_value = round(sparsity_value, digits = 3),
      tx_perarea =  tx_perarea,
      tx_percell = tx_percell,
      cell_tx_fraction = cell_tx_fraction ,
      mean_ratio = mean_signal_ratio,
      max_ratio = max_ratio
      #max_detection= max_detection,

    )

    Metric_table = data.frame(
      Metric = c("Entropy", "Sparsity", "Tx per Area","Cell tx Fraction","Mean Singal Ratio","Max Ratio"),
      Value = c(entropy_value, sparsity_value, tx_perarea,cell_tx_fraction,mean_signal_ratio,max_ratio )
    )
    Metric_table$Metric <- factor(Metric_table$Metric, levels = Metric_table$Metric)

    pdf(pdfFile, width = 12, height = 6)

    p <- ggplot( Metric_table, aes(x = Metric, y = Value)) +
      geom_point(size = 4) +
      theme_minimal() +
      labs(title = "Single Sample Metrics", x = "", y = "Value")
    res_table <- tableGrob(res, rows = NULL)
    Metric_table <- tableGrob(Metric_table)
    grid.arrange(p,  Metric_table, ncol = 2, heights = c(5/6, 1/6))
    grid.newpage()
    res_table <- tableGrob(res, rows = NULL)
    grid.draw(res_table)
    dev.off()

    print(res)
    return(res)
  }
