# Install 'devtools' if it's not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)
# Install BiocManager if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Function to install packages only if they are not already installed
install_if_missing <- function(package_list) {
  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

# Function to install Bioconductor packages using BiocManager
install_bioc_if_missing <- function(package_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  for (pkg in package_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
}


# # Function to install GitHub packages if they are not already installed

install_github_if_missing <- function(repo) {
  package_name <- unlist(strsplit(repo, "/"))[2]  # Extracts package name from repo string
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    devtools::install_github(repo)
  }
}



# List of required CRAN packages
required_cran_packages <- c("dplyr", "ggplot2", "shadowtext", "scales", "cowplot",
                            "data.table", "Matrix", "matrixStats", "fs", "gridExtra",
                            "magrittr", "tibble", "lsa")

# List of required Bioconductor packages
required_bioc_packages <- c("SingleCellExperiment", "SpatialExperiment",
                            "bluster", "BiocParallel", "BioQC",  "NMI", "coop","scater")

# Install CRAN packages if they are missing
install_if_missing(required_cran_packages)

# Install Bioconductor packages if they are missing
install_bioc_if_missing(required_bioc_packages)

# Install 'scater' with all dependencies
if (!requireNamespace("scater", quietly = TRUE)) {
  BiocManager::install("scater", dependencies = TRUE)
}

#options(repos = BiocManager::repositories(release = TRUE))
#BiocManager::install("SpatialFeatureExperiment")
# install.packages("devtools")

# Define the function to install GitHub packages if they are not already installed
install_github_if_missing <- function(repo, ref = NULL) {
  package_name <- unlist(strsplit(repo, "/"))[2]  # Extracts package name from repo string
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    # Check if a specific branch/reference is provided
    if (!is.null(ref)) {
      devtools::install_github(repo, ref = ref)
    } else {
      devtools::install_github(repo)
    }
  }
}

# Install 'SpatialFeatureExperiment' from GitHub if missing
install_github_if_missing("pachterlab/SpatialFeatureExperiment")

# Install 'Voyager' from GitHub's 'devel' branch if missing
install_github_if_missing("pachterlab/voyager", ref = "devel")
# Install `RCTD`
install_github_if_missing("dmcable/RCTD")

#library(SpatialExperiment)
#library(Voyager)




### Suggest for InSiteType
#
# if (!requireNamespace("InSituType", quietly = TRUE)) {
#   stop("The 'InSituType' package is required but not installed. ",
#        "Install it using devtools::install_github("https://github.com/Nanostring-Biostats/InSituType") ",
#        "if you wish to use the annotateData function.")
# }
#
## For MAC user


# gcc and gfortran are required for mac user to install `InSituType`

#/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# brew install gcc
# gcc --version
# mkdir -p ~/.R
# touch ~/.R/Makevars
# vi ~/.R/Makevars
  # CC=/opt/homebrew/bin/gcc-14
  # CXX=/opt/homebrew/bin/g++-14
  # CXX11=/opt/homebrew/bin/g++-14
  # CXX14=/opt/homebrew/bin/g++-14
  # CXX17=/opt/homebrew/bin/g++-14
  # FC=/opt/homebrew/bin/gfortran-14
  # F77=/opt/homebrew/bin/gfortran-14


#library(InSituType)
#library(spacexr)
remotes::install_github('https://github.com/Center-for-Spatial-OMICs/SpatialQC')


## still has warnings


# Warning message:
#   replacing previous import ‘dplyr::combine’ by ‘gridExtra::combine’ when loading ‘SpatialQC’
