% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_update_final.R
\name{getPanelSize}
\alias{getPanelSize}
\title{getPanelSize
Get Panel Size from Expression Matrix or Seurat Object.}
\usage{
getPanelSize(seu_obj = NULL, expMat = "path_to_expMat", platform = NULL)
}
\arguments{
\item{seu_obj}{An optional Seurat object containing the expression data. If provided,
the function will calculate the panel size based on the number of features
present in the object.}

\item{expMat}{The path to the expression matrix file or directory containing
expression data files. This parameter is used if no Seurat object is provided.
The function supports both flat files and Hierarchical Data Format (HDF5) files.}

\item{platform}{The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
and 'Merscope'. Note: 'Merscope' is currently not supported.}
}
\value{
Returns an integer indicating the number of valid features (excluding controls
and system features) present in the expression data.
}
\description{
Calculate the size of the panel by counting the number of features
in an expression matrix, excluding negative or control features. It can handle
different types of expression data input, including direct paths to expression
matrices or Seurat object with embedded expression data.
}
