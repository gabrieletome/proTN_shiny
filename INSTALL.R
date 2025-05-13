################################################################################
# ProTN: an integrative pipeline for complete analysis of proteomics           # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# Issue at: https://github.com/TebaldiLab/ProTN/issues                         #
# PI: Dr. Toma Tebaldi, PhD                                                    #
################################################################################
#Install packages for Shiny App Launcher
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("shiny", quietly = TRUE))
  install.packages("shiny", dependencies = T)
if (!require("shinydashboard", quietly = TRUE))
  install.packages("shinydashboard", dependencies = T)
if (!require("shinydashboardPlus", quietly = TRUE))
  install.packages("shinydashboardPlus", dependencies = T)
if (!require("shinymaterial", quietly = TRUE))
  install.packages("shinymaterial", dependencies = T)
if (!require("shinyjs", quietly = TRUE))
  install.packages("shinyjs", dependencies = T)
if (!require("shinyBS", quietly = TRUE))
  install.packages("shinyBS", dependencies = T)
if (!require("rmdformats", quietly = TRUE))
  install.packages("rmdformats", dependencies = T)
if (!require("extrafont", quietly = TRUE))
  install.packages("extrafont", dependencies = T)

#Install packages for ProTN
if (!require("devtools", quietly = TRUE))
  install.packages("devtools", dependencies = T)
if (!require("wesanderson", quietly = TRUE))
  install.packages("wesanderson", dependencies = T)
if (!require("svgPanZoom", quietly = TRUE))
  install.packages("svgPanZoom", dependencies = T)

if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")
if (!require("impute", quietly = TRUE))
  BiocManager::install("impute")
if (!require("pvca", quietly = TRUE))
  BiocManager::install("pvca")
if (!require("sva", quietly = TRUE))
  BiocManager::install("sva")

if (!require("SummarizedExperiment", quietly = TRUE))
  BiocManager::install("SummarizedExperiment")
if (!require("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb")
if (!require("DEqMS", quietly = TRUE))
  BiocManager::install("DEqMS")
if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")
if (!require("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")
if (!require("sva", quietly = TRUE))

if (!require("proBatch", quietly = TRUE))
  devtools::install_github("symbioticMe/proBatch", dependencies = T)
if (!require("PhosR", quietly = TRUE))
  devtools::install_github("PYangLab/PhosR", dependencies = T)
if (!require("proTN", quietly = TRUE))
  devtools::install_github("tomegabriele/proTN_package", dependencies = T)

#Install additional packages for PhosProTN
if (!require("shinyWidgets", quietly = TRUE))
  install.packages("shinyWidgets", dependencies = T)

# Install reticulate for kaleido of plotly
if (!require("reticulate", quietly = TRUE)){
  install.packages('reticulate')
  reticulate::install_miniconda()
  reticulate::conda_install('r-reticulate', 'python-kaleido')
  reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
  reticulate::use_miniconda('r-reticulate')
}

#Check installation package
list.of.packages <- c("BiocManager","shiny","markdown","knitr",
                      "shinydashboard","shinydashboardPlus","shinymaterial",
                      "shinyjs","shinyBS", "rmdformats","extrafont",
                      "devtools","wesanderson","svgPanZoom", "biomaRt", 
                      'SummarizedExperiment', 'STRINGdb', 'DEqMS', 'limma', 'preprocessCore',
                      "impute", "pvca","sva","proBatch","PhosR","proTN")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  message(paste0("ERROR: Some error occur durig the installation of the current package:\n\t",
                 paste0(new.packages, collapse = "\n\t"),
                 "\nPlease try to install them again."))}
