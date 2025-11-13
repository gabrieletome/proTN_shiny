# ProTN v0.3.1.1: an integrative pipeline for complete analysis of proteomics    # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# PI: Dr. Toma Tebaldi, PhD                                                    #
#                                                                              #
list.of.packages <- c("shiny","tidyverse","markdown","knitr","shinydashboard",
                      "shinydashboardPlus","shinymaterial","shinyjs","magrittr",
                      "dplyr","stringr","shinyBS","DT","bslib","readr",
                      "plotly","rhandsontable","shinyalert","ggplot2","webshot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) suppressMessages(suppressWarnings({install.packages(new.packages, dependencies = T)}))

#Load library
library(shiny)
library(tidyverse)
library(markdown)
library(knitr)
library(shinydashboard)
library(shinydashboardPlus)
library(shinymaterial)
library(shinyalert)
library(shinyjs)
library(magrittr)
library(ggplot2)
library(dplyr)
library(stringr)
library(stringi)
library(shinyBS)
library(bslib)
library(plotly)
library(readr)
library(DT)
library(limma)
library(proTN)
library(data.table)
library(rhandsontable)
library(webshot2)

tmpdir<-stri_replace_all_regex(tempdir(), pattern = "/Rtmp\\w+", replacement = "")
if (!dir.exists(file.path(tmpdir, "ProTN_shiny"))){
  dir.create(file.path(tmpdir, "ProTN_shiny"), showWarnings = FALSE)
}
Sys.setenv(TMPDIR=file.path(tmpdir, "ProTN_shiny"))
unlink(tempdir(), recursive = T)
tempdir(check = T)

jsCode_STRINGdb <- "
    shinyjs.loadStringData = function(params) {
      var defaultParams = {
        taxonomy : null,
        gene : ['TP53'],
        score_thr : 700
      };
      params = shinyjs.getParams(params, defaultParams);

      getSTRING('https://string-db.org', {
            'species': params.taxonomy,
            'identifiers': params.gene,
            'required_score': params.score_thr,
            'hide_disconnected_nodes': '1'})
    }"


# Define UI for application that draws a histogram
ui <- tagList(
  dashboardPage(
    skin = "black",
    title = "ProTN",
    # Application title
    header=dashboardHeader(
      title=tags$a(href="#shiny-tab-info",
                   tags$img(id="logo_protn"))
    ),
    #Sidebar menu
    sidebar=dashboardSidebar(
      sidebarMenu(
        id="tabs",
        menuItem("Info ProTN", tabName = "info", icon = icon("info-circle", lib="font-awesome")),
        menuItem("Run ProTN", tabName = "analysis_protn", icon = icon("rocket", "fa-regular")),
        menuItem("Info PhosProTN", tabName = "info_phos", icon = icon("info-circle", lib="font-awesome")),
        menuItem("Run PhosProTN", tabName = "analysis_protn_phos", icon = icon("rocket", "fa-regular")),
        menuItem("Info PhosProTN with prot.", tabName = "info_phos_protn", icon = icon("info-circle", lib="font-awesome")),
        menuItem("Run PhosProTN with prot.", tabName = "analysis_protn_phos_protn", icon = icon("rocket", "fa-regular")),
        menuItem("Info InteracTN", tabName = "info_interactn", icon = icon("info-circle", lib="font-awesome")),
        menuItem("Run InteracTN", tabName = "analysis_interactn", icon = icon("rocket", "fa-regular")),
        menuItem("Contacts", tabName = "contacts", icon = icon("comment", lib="font-awesome"))
      )
    ),
    #Main body web page
    body=dashboardBody(
      #Load associated file, CSS, JS, logo
      useShinyjs(),
      useShinyalert(),
      tags$meta(charset = "UTF-8"),
      tags$head(tags$link(rel="shortcut icon", href="images/logo_black.ico")),
      tags$head(tags$script(src="https://kit.fontawesome.com/5d5f342cf8.js")),
      tags$head(tags$script(src = "http://string-db.org/javascript/combined_embedded_network_v2.0.4.js")),
      tags$link(rel="stylesheet", href="https://fonts.googleapis.com/css?family=El+Messiri"),
      includeCSS("www/css/custom_theme.css"),
      includeCSS("www/css/materialize.css"),
      includeScript("www/js/materialize.js"),
      includeScript("www/js/full_screen_plot.js"),
      extendShinyjs(text = jsCode_STRINGdb, functions = c("loadStringData")),

      #Busy panel when app is running
      conditionalPanel(
        condition = "$(\'html\').hasClass(\'shiny-busy\')",
        tags$div(class = "loader"),
        tags$div(class = "prevent_click")
      ),    
      tabItems(
        #INFO tab ProTN
        tabItem(
          tabName = "info",
          includeHTML("www/README.html")
        ),
        #Execution tab of ProTN
        tabItem(
          tabName = "analysis_protn",
          tagList(
            fluidRow(
              column(
                width = 3,
                fluidRow(
                  downloadButton("download_CS_proteome", "Download case study!"),
                ),
                fluidRow(
                  textInput("title_exp", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp", "Brief description", rows = 4),
                ),
                fluidRow(
                  selectInput("sw_analyzer", "Software Analyzer:",
                              choice = list("ProteomeDiscoverer" = "PD", 
                                            "MaxQuant by evidence.txt" = "MQ_ev", 
                                            "MaxQuant by peptides.txt and proteinGroups.txt" = "MQ_prot",
                                            "Spectronaut" = "SP",
                                            "FragPipe" = "FP"),
                                 selected = "PD", multiple = FALSE
                  ),
                ),
                uiOutput("input_proteome"),
                checkboxInput("batch_correction", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui"),
                checkboxInput("advance_filter", "Advance Filter", FALSE),
                uiOutput("advance_filter_ui"),
                actionButton("report_proteome", "Load data!"),
                tags$h3("Select what execute:"),
                checkboxInput("abundance_plot", "% missing values", TRUE),
                checkboxInput("peptide_distribution", "N° peptides per protein", TRUE),
                checkboxInput("raw_violin", "Distribution raw abundance", FALSE),
                checkboxInput("complexity_plot", "Complexity plot", FALSE),
                checkboxInput("protein_violin", "Distribution abundance proteins", FALSE),
                checkboxInput("peptide_violin", "Distribution abundance peptides", FALSE),
                checkboxInput("mds_protein", "MDS based on protein", FALSE),
                checkboxInput("mds_peptide", "MDS based on peptide", FALSE),
                checkboxInput("pca_protein", "PCA based on protein", TRUE),
                checkboxInput("pca_peptide", "PCA based on peptide", FALSE),
                checkboxInput("boxplot_protein", "Boxplot selected proteins", FALSE),
                checkboxInput("heatmap_protein", "Heatmap selected proteins", FALSE),
                uiOutput("list_protein_ui"),
                tags$h3("Differential Analysis:"),
                checkboxInput("differential_analysis_checkbox", "Execute differential analysis", FALSE),
                uiOutput("differential_params_ui")
              ),
              column(
                id="panel_results",
                width = 9,
                tags$br(),
                # textOutput("messagge_read"),
                uiOutput("protn_results_ui"),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_abundance_plot")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_distribution")
                #   )
                # ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_abundance_plot")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_distribution")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_raw_violin")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_complexity_plot")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_protein_violin")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_violin")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_protein")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_protein")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide")
                  )
                ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_protein_violin")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_violin")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_protein")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_peptide")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_protein")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_peptide")
                #   )
                # ),
                uiOutput("render_protein_boxplot"),
                uiOutput("render_protein_heatmap"),
                uiOutput("render_differential_analysis"),
                uiOutput("render_protein_diff_table"),
                uiOutput("render_peptide_diff_table"),
                uiOutput("render_protein_diff_barplot"),
                uiOutput("render_peptide_diff_barplot"),
                uiOutput("render_protein_upset"),
                uiOutput("render_peptide_upset"),
                uiOutput("render_protein_ma_plot"),
                uiOutput("render_peptide_ma_plot"),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_protein_vulcano")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_vulcano")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_protein_diff")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_diff")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_protein_diff")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_diff")
                  )
                ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_protein_vulcano")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_vulcano")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_protein_diff")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_peptide_diff")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_protein_diff")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_peptide_diff")
                #   )
                # ),
                uiOutput("render_enrichement_analysis"),
                uiOutput("render_stringdb")
              )
            ),
            
            
            # Fullscreen overlay for any clicked plot
            tags$div(
              id = "fullscreen-container",
              tags$div(class = "close-btn", onclick = "hideFullscreenPlot()", "×"),
              plotOutput("fullscreen_plot", height = "60%", width = "80%")
            )
          )
        ),
        #INFO tab PhosProTN
        tabItem(
          tabName = "info_phos",
          includeHTML("www/README_phos.html")
        ),
        #Execution tab of PhosProTN
        tabItem(
          tabName = "analysis_protn_phos",
          tagList(
            fluidRow(
              column(
                width = 3,
                fluidRow(
                  downloadButton("download_CS_phos", "Download case study!"),
                ),
                fluidRow(
                  textInput("title_exp_phos", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp_phos", "Brief description", rows = 4),
                ),
                fluidRow(
                  selectInput("sw_analyzer_phos", "Software Analyzer:",
                              choice = list("ProteomeDiscoverer" = "PD", 
                                            "MaxQuant by evidence.txt" = "MQ"),
                              selected = "PD", multiple = FALSE
                  ),
                ),
                uiOutput("input_proteome_phos"),
                checkboxInput("batch_correction_phos", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui_phos"),
                sliderInput("phos_thr", "Phosphorylation threshold", 0, 100, step = 5, value = 75),
                checkboxInput("advance_filter_phos", "Advance Filter", FALSE),
                uiOutput("advance_filter_ui_phos"),
                actionButton("report_proteome_phos", "Load data!"),
                tags$h3("Select what execute:"),
                checkboxInput("phospho_percentage_plot_phos", "% phosphorylated site", TRUE),
                checkboxInput("abundance_plot_phos", "% missing values", TRUE),
                checkboxInput("peptide_distribution_phos", "N° peptides per protein", TRUE),
                checkboxInput("raw_violin_phos", "Distribution raw abundance", FALSE),
                checkboxInput("complexity_plot_phos", "Complexity plot", FALSE),
                checkboxInput("peptide_violin_phos", "Distribution abundance phospho-peptides", FALSE),
                checkboxInput("mds_peptide_phos", "MDS based on phospho-peptide", FALSE),
                checkboxInput("pca_peptide_phos", "PCA based on phospho-peptide", TRUE),
                checkboxInput("boxplot_protein_phos", "Boxplot selected proteins", FALSE),
                checkboxInput("heatmap_protein_phos", "Heatmap selected proteins", FALSE),
                uiOutput("list_protein_ui_phos"),
                tags$h3("Differential Analysis:"),
                checkboxInput("differential_analysis_checkbox_phos", "Execute differential analysis", FALSE),
                uiOutput("differential_params_ui_phos")
              ),
              column(
                id="panel_results_phos",
                width = 9,
                tags$br(),
                # textOutput("messagge_read_phos"),
                uiOutput("protn_results_ui_phos"),
                uiOutput("render_phospho_percentage_plot_phos"),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_abundance_plot_phos")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_distribution_phos")
                #   )
                # ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_abundance_plot_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_distribution_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_raw_violin_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_complexity_plot_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_violin_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_phos")
                  )
                ),
                uiOutput("render_protein_boxplot_phos"),
                uiOutput("render_protein_heatmap_phos"),
                uiOutput("render_differential_analysis_phos"),
                uiOutput("render_peptide_diff_table_phos"),
                uiOutput("render_peptide_diff_barplot_phos"),
                uiOutput("render_peptide_upset_phos"),
                uiOutput("render_peptide_ma_plot_phos"),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_vulcano_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_diff_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_diff_phos")
                  )
                ),
                uiOutput("render_enrichement_analysis_phos"),
                uiOutput("render_stringdb_phos"),
                uiOutput("render_kinase_tree_phos")
              )
            ),
            
            
            # Fullscreen overlay for any clicked plot
            tags$div(
              id = "fullscreen-container-phos",
              tags$div(class = "close-btn", onclick = "hideFullscreenPlot_phos()", "×"),
              plotOutput("fullscreen_plot_phos", height = "60%", width = "80%")
            )
          )
        ),
        #INFO tab PhosProTN_with_prot
        tabItem(
          tabName = "info_phos_protn",
          includeHTML("www/README_phos_protn.html")
        ),
        #Execution tab of PhosProTN_with_prot
        tabItem(
          tabName = "analysis_protn_phos_protn",
          tagList(
            fluidRow(
              column(
                width = 3,
                fluidRow(
                  downloadButton("download_CS_phos_protn", "Download case study!"),
                ),
                fluidRow(
                  textInput("title_exp_phos_protn", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp_phos_protn", "Brief description", rows = 4),
                ),
                fluidRow(
                  selectInput("sw_analyzer_phos_protn", "Software Analyzer:",
                              choice = list("ProteomeDiscoverer" = "PD", 
                                            "MaxQuant by evidence.txt" = "MQ"),
                              selected = "PD", multiple = FALSE
                  ),
                ),
                uiOutput("input_proteome_phos_protn"),
                checkboxInput("batch_correction_phos_protn", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui_phos_protn"),
                sliderInput("phos_thr_phos_protn", "Phosphorylation threshold", 0, 100, step = 5, value = 75),
                checkboxInput("advance_filter_phos_protn", "Advance Filter", FALSE),
                uiOutput("advance_filter_ui_phos_protn"),
                actionButton("report_proteome_phos_protn", "Load data!"),
                tags$h3("Select what execute:"),
                checkboxInput("phospho_percentage_plot_phos_protn", "% phosphorylated site", TRUE),
                checkboxInput("abundance_plot_phos_protn", "% missing values", TRUE),
                checkboxInput("peptide_distribution_phos_protn", "N° peptides per protein", TRUE),
                checkboxInput("raw_violin_prot_phos_protn", "Distribution raw abundance proteome", FALSE),
                checkboxInput("raw_violin_pep_phos_protn", "Distribution raw abundance phospho-peptides", FALSE),
                checkboxInput("complexity_plot_phos_protn", "Complexity plot", FALSE),
                checkboxInput("protein_violin_phos_protn", "Distribution abundance proteins", FALSE),
                checkboxInput("peptide_violin_phos_protn", "Distribution abundance phospho-peptides", FALSE),
                checkboxInput("mds_protein_phos_protn", "MDS based on protein", FALSE),
                checkboxInput("mds_peptide_phos_protn", "MDS based on phospho-peptide", FALSE),
                checkboxInput("pca_protein_phos_protn", "PCA based on protein", FALSE),
                checkboxInput("pca_peptide_phos_protn", "PCA based on phospho-peptide", TRUE),
                checkboxInput("boxplot_protein_phos_protn", "Boxplot selected proteins", FALSE),
                checkboxInput("heatmap_protein_phos_protn", "Heatmap selected proteins", FALSE),
                uiOutput("list_protein_ui_phos_protn"),
                tags$h3("Differential Analysis:"),
                checkboxInput("differential_analysis_checkbox_phos_protn", "Execute differential analysis", FALSE),
                uiOutput("differential_params_ui_phos_protn")
              ),
              column(
                id="panel_results_phos_protn",
                width = 9,
                tags$br(),
                # textOutput("messagge_read_phos_protn"),
                uiOutput("protn_results_ui_phos_protn"),
                uiOutput("render_phospho_percentage_plot_phos_protn"),
                uiOutput("render_abundance_plot_phos_protn"),
                uiOutput("render_peptide_distribution_phos_protn"),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_raw_violin_prot_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_raw_violin_pep_phos_protn")
                  )
                ),
                uiOutput("render_complexity_plot_phos_protn"),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_violin_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_phos_protn")
                  )
                ),
                uiOutput("render_protein_boxplot_phos_protn"),
                uiOutput("render_protein_heatmap_phos_protn"),
                uiOutput("render_differential_analysis_phos_protn"),
                uiOutput("render_peptide_diff_table_phos_protn"),
                uiOutput("render_peptide_diff_barplot_phos_protn"),
                uiOutput("render_peptide_upset_phos_protn"),
                uiOutput("render_peptide_ma_plot_phos_protn"),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_vulcano_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_diff_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_diff_phos_protn")
                  )
                ),
                uiOutput("render_enrichement_analysis_phos_protn"),
                uiOutput("render_stringdb_phos_protn"),
                uiOutput("render_kinase_tree_phos_protn")
              )
            ),
            
            
            # Fullscreen overlay for any clicked plot
            tags$div(
              id = "fullscreen-container-phos-protn",
              tags$div(class = "close-btn", onclick = "hideFullscreenPlot_phos_protn()", "×"),
              plotOutput("fullscreen_plot_phos_protn", height = "60%", width = "80%")
            )
          )
        ),
        #INFO tab InteracTN
        tabItem(
          tabName = "info_interactn",
          includeHTML("www/README_interactn.html")
        ),
        #Execution tab of InteracTN
        tabItem(
          tabName = "analysis_interactn",
          tagList(
            fluidRow(
              column(
                width = 3,
                fluidRow(
                  downloadButton("download_CS_interactn", "Download case study!"),
                ),
                fluidRow(
                  textInput("title_exp_interactn", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp_interactn", "Brief description", rows = 4),
                ),
                fluidRow(
                  selectInput("sw_analyzer_interactn", "Software Analyzer:",
                              choice = list("ProteomeDiscoverer" = "PD", 
                                            "MaxQuant by evidence.txt" = "MQ_ev", 
                                            "MaxQuant by peptides.txt and proteinGroups.txt" = "MQ_prot",
                                            "Spectronaut" = "SP",
                                            "FragPipe" = "FP"),
                              selected = "PD", multiple = FALSE
                  ),
                ),
                uiOutput("input_proteome_interactn"),
                checkboxInput("batch_correction_interactn", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui_interactn"),
                checkboxInput("advance_filter_interactn", "Advance Filter", FALSE),
                uiOutput("advance_filter_ui_interactn"),
                actionButton("report_proteome_interactn", "Load data!"),
                tags$h3("Select what execute:"),
                checkboxInput("abundance_plot_interactn", "% missing values", TRUE),
                checkboxInput("peptide_distribution_interactn", "N° peptides per protein", TRUE),
                checkboxInput("raw_violin_interactn", "Distribution raw abundance", FALSE),
                checkboxInput("complexity_plot_interactn", "Complexity plot", FALSE),
                checkboxInput("protein_violin_interactn", "Distribution abundance proteins", FALSE),
                checkboxInput("peptide_violin_interactn", "Distribution abundance peptides", FALSE),
                checkboxInput("mds_protein_interactn", "MDS based on protein", FALSE),
                checkboxInput("mds_peptide_interactn", "MDS based on peptide", FALSE),
                checkboxInput("pca_protein_interactn", "PCA based on protein", TRUE),
                checkboxInput("pca_peptide_interactn", "PCA based on peptide", FALSE),
                checkboxInput("boxplot_protein_interactn", "Boxplot selected proteins", FALSE),
                checkboxInput("heatmap_protein_interactn", "Heatmap selected proteins", FALSE),
                uiOutput("list_protein_ui_interactn"),
                tags$h3("Differential Analysis:"),
                checkboxInput("differential_analysis_checkbox_interactn", "Execute differential analysis", FALSE),
                uiOutput("differential_params_ui_interactn")
              ),
              column(
                id="panel_results_interactn",
                width = 9,
                tags$br(),
                # textOutput("messagge_read"),
                uiOutput("protn_results_ui_interactn"),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_abundance_plot_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_distribution_interactn")
                #   )
                # ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_abundance_plot_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_distribution_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_raw_violin_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_complexity_plot_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_protein_violin_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_violin_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_protein_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_protein_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_interactn")
                  )
                ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_protein_violin_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_violin_interactn")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_protein_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_peptide_interactn")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_protein_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_peptide_interactn")
                #   )
                # ),
                uiOutput("render_protein_boxplot_interactn"),
                uiOutput("render_protein_heatmap_interactn"),
                uiOutput("render_differential_analysis_interactn"),
                uiOutput("render_protein_diff_table_interactn"),
                uiOutput("render_peptide_diff_table_interactn"),
                uiOutput("render_protein_diff_barplot_interactn"),
                uiOutput("render_peptide_diff_barplot_interactn"),
                uiOutput("render_protein_upset_interactn"),
                uiOutput("render_peptide_upset_interactn"),
                uiOutput("render_protein_ma_plot_interactn"),
                uiOutput("render_peptide_ma_plot_interactn"),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_protein_vulcano_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_peptide_vulcano_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_protein_diff_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_mds_peptide_diff_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_protein_diff_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 11,
                    uiOutput("render_pca_peptide_diff_interactn")
                  )
                ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_protein_vulcano_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_peptide_vulcano_interactn")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_protein_diff_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_mds_peptide_diff_interactn")
                #   )
                # ),
                # fluidRow(
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_protein_diff_interactn")
                #   ),
                #   column(
                #     width = 6,
                #     uiOutput("render_pca_peptide_diff_interactn")
                #   )
                # ),
                uiOutput("render_enrichement_analysis_interactn"),
                uiOutput("render_stringdb_interactn")
              )
            ),
            
            
            # Fullscreen overlay for any clicked plot
            tags$div(
              id = "fullscreen-container-interactn",
              tags$div(class = "close-btn", onclick = "hideFullscreenPlot_interactn()", "×"),
              plotOutput("fullscreen_plot_interactn", height = "60%", width = "80%")
            )
          )
        ),
        #Contact tab
        tabItem(
          tabName = "contacts",
          tagList(
            includeHTML("www/contacts.html")
          )
        )
      )
    )
  ),
  #Footer with lab, name and university
  tags$footer(
    tags$div(
      class="footer-row",
      tags$div(class="footer-col", style="text-align: left;",
               tags$p("This website is free and open to all users and there is no login requirement."))),
    align = "center")
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 10 * 1024^3)
  # PROTN: db variable analysis ----
  db_execution <- reactiveValues(session = session$token,
                                 parameter = list(),
                                 data_loaded = FALSE,
                                 dirOutput = "",
                                 proteome_data = list(),
                                 imputed_data = list(),
                                 normalized_data = list(),
                                 formule_contrast = list(),
                                 dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                 differential_results = list(),
                                 enrichment_results = list(),
                                 stringdb_res = list(),
                                 kinase_tree_res = list(),
                                 phospho_percentage = NULL,
                                 generate_abundance = NULL,
                                 generate_complexity = NULL, raw_abundance_distribution = NULL,
                                 generate_peptide_distribution = NULL,
                                 protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                 protein_MDS = NULL, peptide_MDS = NULL,
                                 protein_PCA = NULL, peptide_PCA = NULL,
                                 protein_boxplot = NULL, protein_heatmap = NULL,
                                 protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                 protein_upset_plot = NULL, peptide_upset_plot = NULL,
                                 protein_ma_plot = NULL, peptide_ma_plot = NULL,
                                 protein_vulcano = NULL, peptide_vulcano = NULL,
                                 protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                 protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  # PHOSPROTN: db variable analysis ----
  db_execution_phos <- reactiveValues(session = session$token,
                                      parameter = list(),
                                      data_loaded = FALSE,
                                      dirOutput = "",
                                      proteome_data = list(),
                                      imputed_data = list(),
                                      normalized_data = list(),
                                      formule_contrast = list(),
                                      dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                      differential_results = list(),
                                      enrichment_results = list(),
                                      stringdb_res = list(),
                                      kinase_tree_res = list(),
                                      phospho_percentage = NULL,
                                      generate_abundance = NULL,
                                      generate_peptide_distribution = NULL, 
                                      raw_abundance_distribution = NULL,
                                      generate_complexity = NULL,
                                      protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                      protein_MDS = NULL, peptide_MDS = NULL,
                                      protein_PCA = NULL, peptide_PCA = NULL,
                                      protein_boxplot = NULL, protein_heatmap = NULL,
                                      protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                      protein_upset_plot = NULL, peptide_upset_plot = NULL,
                                      protein_ma_plot = NULL, peptide_ma_plot = NULL,
                                      protein_vulcano = NULL, peptide_vulcano = NULL,
                                      protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                      protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  # PhosProTN_with_prot: db variable analysis ----
  db_execution_phos_protn <- reactiveValues(session = session$token,
                                            parameter = list(),
                                            data_loaded = FALSE,
                                            dirOutput = "",
                                            proteome_data = list(),
                                            imputed_data = list(),
                                            normalized_data = list(),
                                            formule_contrast = list(),
                                            dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                            differential_results = list(),
                                            enrichment_results = list(),
                                            stringdb_res = list(),
                                            kinase_tree_res = list(),
                                            phospho_percentage = NULL,
                                            generate_abundance = NULL,
                                            generate_peptide_distribution = NULL,
                                            raw_proteome_abundance_distribution = NULL, raw_phospho_abundance_distribution = NULL,
                                            generate_complexity = NULL,
                                            protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                            protein_MDS = NULL, peptide_MDS = NULL,
                                            protein_PCA = NULL, peptide_PCA = NULL,
                                            protein_boxplot = NULL, protein_heatmap = NULL,
                                            protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                            protein_upset_plot = NULL, peptide_upset_plot = NULL,
                                            protein_ma_plot = NULL, peptide_ma_plot = NULL,
                                            protein_vulcano = NULL, peptide_vulcano = NULL,
                                            protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                            protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  
  # InteracTN: db variable analysis ----
  db_execution_interactn <- reactiveValues(session = session$token,
                                           parameter = list(),
                                 data_loaded = FALSE,
                                 dirOutput = "",
                                 proteome_data = list(),
                                 imputed_data = list(),
                                 normalized_data = list(),
                                 formule_contrast = list(),
                                 dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                 differential_results = list(),
                                 enrichment_results = list(),
                                 stringdb_res = list(),
                                 kinase_tree_res = list(),
                                 phospho_percentage = NULL,
                                 generate_abundance = NULL,
                                 generate_peptide_distribution = NULL, raw_abundance_distribution = NULL,
                                 generate_complexity = NULL,
                                 protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                 protein_MDS = NULL, peptide_MDS = NULL,
                                 protein_PCA = NULL, peptide_PCA = NULL,
                                 protein_boxplot = NULL, protein_heatmap = NULL,
                                 protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                 protein_upset_plot = NULL, peptide_upset_plot = NULL,
                                 protein_ma_plot = NULL, peptide_ma_plot = NULL,
                                 protein_vulcano = NULL, peptide_vulcano = NULL,
                                 protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                 protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  ##############################################################################
  ### PROTN ----
  # Optional visibility based on the selection ----
  
  ## PROTN: Visibility of the proteomics files for ProTN ----
  output$input_proteome <- renderUI({
    message(input$sw_analyzer)
    if (input$sw_analyzer == "PD"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome", "Select the PEP file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_proteome", "Select the PROT file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer == "MQ_ev"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome", "Select the EVIDENCE file of the PROTEOMICS..."),
        )
      )
    } else if (input$sw_analyzer == "MQ_prot"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome", "Select the Peptides.txt file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_proteome", "Select the ProteinGroups.txt file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer == "SP"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome", "Select the Spectronaut report file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer == "FP"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome", "Select the combined_modified_peptide.tsv file of the PROTEOMICS..."),
        )
      )
    } else{
      tagList(
        tags$p("BACK")
      )
    }
  })
  
  ## PROTN: textbox for batch correction----
  output$batch_correction_ui <- renderUI({ 
    if(input$batch_correction){
      textInput("batch_correction_col", "Column in Annotation file with the batch:")
    } 
  })
  ## PROTN: advance filters----
  output$advance_filter_ui <- renderUI({ 
    if(input$advance_filter){
      tagList(
        numericInput("NA_allow_condition", "N° missing value allow per condition", value = 0, min = 0, max = 5),
        numericInput("min_peptide_protein", "Minimum peptide per protein", value = 1, min = 1),
        selectizeInput("impute_algorithm", "Select impute algorithm:",
                       choices = list("PhosR and normalization" = "phosr_norm", 
                                      "Gaussian estimation and normalization" = "gaussian_norm",
                                      "missForest and normalization" = "missForest_norm",
                                      "Pre-normalization and PhosR" = "norm_phosr", 
                                      "Pre-normalization and Gaussian estimation" = "norm_gaussian",
                                      "Pre-normalization and missForest" = "norm_missForest",
                                      "Pre-normalization and pcaMethods" = "norm_pcaMethods"),
                       selected = "norm_phosr", multiple = FALSE
        ),
        textInput("sample_column", "Column name with the sample name:", value = "Sample")
      )
    } 
  })
  
  ## PROTN: textbox for list proteins ----
  output$list_protein_ui <- renderUI({ 
    if(input$boxplot_protein | input$heatmap_protein){
      textInput("list_proteins", "List proteins to show (separate by: \",\"):")
    } 
  })
  
  ## PROTN: show parameter for differential analysis ----
  output$differential_params_ui <- renderUI({ 
    if(input$differential_analysis_checkbox){
      tagList(
        tags$label("Write in each line a different comparison"),
        tags$label("(right click to add row)"),
        rHandsontableOutput('render_formule_contrast_table'),
        textInput("FC_thr", "Fold change threshold for significance:",value = 0.5),
        radioButtons("pval_fdr", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_val"),
        textInput("pval_thr", "P.value threshold for significance:", value = 0.05),
        actionButton("execute_differential_analysis_btn", "Run!"),
        checkboxInput("protein_diff_table", "Proteins differentiated table", FALSE),
        checkboxInput("peptide_diff_table", "Peptides differentiated table", FALSE),
        checkboxInput("protein_diff_barplot", "Proteins differentiated barplot", TRUE),
        checkboxInput("peptide_diff_barplot", "Peptides differentiated barplot", FALSE),
        checkboxInput("protein_upset", "Proteins upset plot", FALSE),
        checkboxInput("peptide_upset", "Peptides upset plot", FALSE),
        checkboxInput("protein_vulcano", "Proteins vulcano plot", FALSE),
        checkboxInput("peptide_vulcano", "Peptides vulcano plot", FALSE),
        checkboxInput("protein_ma_plot", "Proteins MA plot", FALSE),
        checkboxInput("peptide_ma_plot", "Peptides MA plot", FALSE),
        checkboxInput("mds_diff_protein", "MDS based on diffential protein", FALSE),
        checkboxInput("mds_diff_peptide", "MDS based on diffential peptide", FALSE),
        checkboxInput("pca_diff_protein", "PCA based on diffential protein", FALSE),
        checkboxInput("pca_diff_peptide", "PCA based on diffential peptide", FALSE),
        tags$h3("Enrichment Analysis:"),
        checkboxInput("enrichment_analysis", "Execute enrichment analysis", FALSE),
        uiOutput("enrichment_params_ui"),
        tags$h3("STRINGdb network:"),
        checkboxInput("stringdb_analysis", "Execute STRINGdb", FALSE),
        uiOutput("stringdb_params_ui")
      )
    } else{
      # Reset UI elements
      updateCheckboxInput(session, "enrichment_analysis", value = FALSE)
      updateCheckboxInput(session, "stringdb_analysis", value = FALSE)
      
      db_execution$formule_contrast <- list()
      db_execution$dt_formule_contrast <- data.table("Name"=c("","","",""),"Formule"=c("","","",""))
      db_execution$differential_results <- list()
      
      updateCheckboxInput(session, "protein_diff_barplot", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_barplot",  value = FALSE)
      updateCheckboxInput(session, "protein_diff_table", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_table", value = FALSE)
      updateCheckboxInput(session, "protein_upset", value = FALSE)
      updateCheckboxInput(session, "peptide_upset", value = FALSE)
      updateCheckboxInput(session, "protein_vulcano",  value =  FALSE)
      updateCheckboxInput(session, "peptide_vulcano", value = FALSE)
      updateCheckboxInput(session, "protein_ma_plot", value = FALSE)
      updateCheckboxInput(session, "peptide_ma_plot", value = FALSE)
      updateCheckboxInput(session, "mds_diff_protein", value = FALSE)
      updateCheckboxInput(session, "mds_diff_peptide", value = FALSE)
      updateCheckboxInput(session, "pca_diff_protein", value = FALSE)
      updateCheckboxInput(session, "pca_diff_peptide", value = FALSE)
      

      db_execution$protein_differential_barplot <- NULL
      db_execution$peptide_differential_barplot <- NULL
      db_execution$protein_upset_plot <- NULL
      db_execution$peptide_upset_plot <- NULL
      db_execution$protein_ma_plot <- NULL
      db_execution$peptide_ma_plot <- NULL
      db_execution$protein_vulcano <- NULL
      db_execution$peptide_vulcano <- NULL
      db_execution$protein_differential_MDS <- NULL
      db_execution$peptide_differential_MDS <- NULL
      db_execution$protein_differential_PCA <- NULL
      db_execution$peptide_differential_PCA <- NULL
      

      output$render_differential_analysis <- renderUI({NULL})
      output$render_protein_diff_table <- renderUI({NULL})
      output$render_peptide_diff_table <- renderUI({NULL})
      output$render_protein_diff_barplot <- renderUI({NULL})
      output$render_peptide_diff_barplot <- renderUI({NULL})
      output$render_protein_upset <- renderUI({NULL})
      output$render_peptide_upset <- renderUI({NULL})
      output$render_protein_ma_plot <- renderUI({NULL})
      output$render_peptide_ma_plot <- renderUI({NULL})
      output$render_protein_vulcano <- renderUI({NULL})
      output$render_peptide_vulcano <- renderUI({NULL})
      output$render_mds_protein_diff <- renderUI({NULL})
      output$render_mds_peptide_diff <- renderUI({NULL})
      output$render_pca_protein_diff <- renderUI({NULL})
      output$render_pca_peptide_diff <- renderUI({NULL})
    }
  })
  
  output$render_formule_contrast_table <- renderRHandsontable({
    rhandsontable(db_execution$dt_formule_contrast, rowHeaders = NULL, stretchH = "all")
  })
  
  ## PROTN: show enrichment parameter ----
  output$enrichment_params_ui <- renderUI({ 
    if(input$enrichment_analysis){
      tagList(
        selectizeInput("DB_enrichment", "DB to analyse:",
                       choices = lapply(split(read_tsv("data/dbs_enrichR.txt", col_names = FALSE)$X1,
                                              read_tsv("data/dbs_enrichR.txt", col_names = FALSE)[,2]), as.list),
                       selected = NULL, multiple = TRUE
        ),
        textInput("terms_enrich", "Terms to search (separated by \",\"):"),
        radioButtons("pval_fdr_enrich", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_adj"),
        textInput("pvalue_enrich", "P.value threshold for significance:", value = 0.05),
        sliderInput("os_enrich", "Overlap size thr for enrichment", 1, 30, step = 1, value = 5),
        checkboxInput("enrich_with_background", "Enrichment with background", FALSE),
        actionButton("execute_enrichment_analysis_btn", "Run!")
      )
    } else{
      db_execution$enrichment_results <- list()
      output$render_enrichement_analysis <- renderUI({NULL})
    }
  })
  
  ## PROTN: show stringdb parameter ----
  output$stringdb_params_ui <- renderUI({
    if(input$stringdb_analysis){
      tagList(
        selectizeInput("taxonomy", "NCBI Taxonomy ID", 
                       choice = data.table::fread("data/subset_tax.csv", select = "name"), 
                       selected = "Homo sapiens", multiple = F),
        sliderInput("score_thr_stringdb", "Score thr for STRINGdb", 500, 1000, step = 10, value = 700),
        actionButton("execute_stringdb_analysis_btn", "Run!"),
        tags$br()
      )
    } else{
      db_execution$stringdb_res <- list()
      output$render_stringdb <- renderUI({NULL})
    }
  })
  
  ## PROTN: function genereting plot ----
  generate_abundance <- reactive({
    req(input$abundance_plot)
    if(input$abundance_plot){
      generate_abundance_fig <- generate_abundance_plot(proteome_data = db_execution$proteome_data)$plot
      db_execution$generate_abundance = generate_abundance_fig
      generate_abundance_fig
    } else{
      db_execution$generate_abundance = NULL
    }
  })
  
  generate_peptide_distribution <- reactive({
    req(input$peptide_distribution)
    if(input$peptide_distribution){
      peptide_distribution_fig <- generate_peptide_distribution_plot(proteome_data = db_execution$proteome_data)$plot
      db_execution$generate_peptide_distribution = peptide_distribution_fig
      peptide_distribution_fig
    } else{
      db_execution$generate_peptide_distribution = NULL
    }
  })
  
  generate_raw_violin <- reactive({
    req(input$raw_violin)
    if(input$raw_violin){
      raw_abundance_distribution_fig <- plot_raw_abundance_distribution(proteome_data = db_execution$proteome_data,
                                                                        type = "protein")$plot
      db_execution$raw_abundance_distribution = raw_abundance_distribution_fig
      raw_abundance_distribution_fig
    } else{
      db_execution$raw_abundance_distribution = NULL
    }
  })
  
  generate_complexity <- reactive({
    req(input$complexity_plot)
    if(input$complexity_plot){
      generate_complexity_fig <- complexity_plot(proteome_data = db_execution$proteome_data)$plot
      db_execution$generate_complexity = generate_complexity_fig
      generate_complexity_fig
    } else{
      db_execution$generate_complexity = NULL
    }
  })
  
  generate_protein_violin <- reactive({
    req(input$protein_violin)
    if(input$protein_violin){
      protein_abundance_distribution_fig <- plot_abundance_distribution(proteome_data = db_execution$normalized_data,
                                                                        type = "protein")$plot
      db_execution$protein_abundance_distribution = protein_abundance_distribution_fig
      protein_abundance_distribution_fig
    } else{
      db_execution$protein_abundance_distribution = NULL
    }
  })
  
  generate_peptide_violin <- reactive({
    req(input$peptide_violin)
    if(input$peptide_violin){
      peptide_abundance_distirbution_fig <- plot_abundance_distribution(proteome_data = db_execution$normalized_data,
                                                                        type = "peptide")$plot
      db_execution$peptide_abundance_distirbution = peptide_abundance_distirbution_fig
      peptide_abundance_distirbution_fig
    } else{
      db_execution$peptide_abundance_distirbution = NULL
    }
  })
  
  generate_mds_protein <- reactive({
    req(input$mds_protein)
    if(input$mds_protein){
      mds_protein_fig <- mds_plot(proteome_data = db_execution$normalized_data,
                                  type = "protein")$plot
      db_execution$protein_MDS = mds_protein_fig
      mds_protein_fig
    } else{
      db_execution$protein_MDS = NULL
    }
  })
  
  generate_mds_peptide <- reactive({
    req(input$mds_peptide)
    if(input$mds_peptide){
      mds_peptide_fig <- mds_plot(proteome_data = db_execution$normalized_data,
                                  type = "peptide")$plot
      db_execution$peptide_MDS = mds_peptide_fig
      mds_peptide_fig
    } else{
      db_execution$peptide_MDS = NULL
    }
  })
  
  generate_pca_protein <- reactive({
    req(input$pca_protein)
    if(input$pca_protein){
      pca_protein_fig <- pca_plot(proteome_data = db_execution$normalized_data,
                                  type = "protein")$plot
      db_execution$protein_PCA = pca_protein_fig
      pca_protein_fig
    } else{
      db_execution$protein_PCA = NULL
    }
  })
  
  generate_pca_peptide <- reactive({
    req(input$pca_peptide)
    if(input$pca_peptide){
      pca_peptide_fig <- pca_plot(proteome_data = db_execution$normalized_data,
                                  type = "peptide")$plot
      db_execution$peptide_PCA = pca_peptide_fig
      pca_peptide_fig
    } else{
      db_execution$peptide_PCA = NULL
    }
  })
  
  generate_protein_boxplot <- reactive({
    req(input$boxplot_protein)
    if(input$boxplot_protein){
      req(input$list_proteins)
      db_execution$parameter <- c(db_execution$parameter, "List proteins boxplot abundance: "=input$list_proteins)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins), regex=",")
      boxplot_protein_fig <- plot_selected_proteins(proteome_data = db_execution$normalized_data,
                                             list_protein = unlist(list_proteins))$plot
      db_execution$protein_boxplot = boxplot_protein_fig
      boxplot_protein_fig
    } else{
      db_execution$protein_boxplot = NULL
    }
  })
  
  generate_protein_heatmap <- reactive({
    req(input$heatmap_protein)
    if(input$heatmap_protein){
      req(input$list_proteins)
      db_execution$parameter <- c(db_execution$parameter, "List proteins heatmap abundance: "=input$list_proteins)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins), regex=",")
      heatmap_protein_fig <- heatmap_selected_proteins(proteome_data = db_execution$normalized_data, list_protein = unlist(list_proteins))$plot
      db_execution$protein_heatmap = heatmap_protein_fig
      heatmap_protein_fig
    } else{
      db_execution$protein_heatmap = NULL
    }
  })
  
  generate_mds_protein_diff <- reactive({
    req(input$mds_diff_protein)
    if(input$mds_diff_protein){
      mds_protein_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                             proteome_data = db_execution$normalized_data,
                                                             type = "protein")$plot
      db_execution$protein_differential_MDS = mds_protein_diff_fig
      mds_protein_diff_fig
    } else{
      db_execution$protein_differential_MDS = NULL
    }
  })
  
  generate_mds_peptide_diff <- reactive({
    req(input$mds_diff_peptide)
    if(input$mds_diff_peptide){
      mds_peptide_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                             proteome_data = db_execution$normalized_data,
                                                             type = "peptide")$plot
      db_execution$peptide_differential_MDS = mds_peptide_diff_fig
      mds_peptide_diff_fig
    } else{
      db_execution$peptide_differential_MDS = NULL
    }
  })
  
  generate_pca_protein_diff <- reactive({
    req(input$pca_diff_protein)
    if(input$pca_diff_protein){
      pca_protein_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                             proteome_data = db_execution$normalized_data,
                                                             type = "protein")$plot
      db_execution$protein_differential_PCA = pca_protein_diff_fig
      pca_protein_diff_fig
    } else{
      db_execution$protein_differential_PCA = NULL
    }
  })
  
  generate_pca_peptide_diff <- reactive({
    req(input$pca_diff_peptide)
    if(input$pca_diff_peptide){
      pca_peptide_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                             proteome_data = db_execution$normalized_data,
                                                             type = "peptide")$plot
      db_execution$peptide_differential_PCA = pca_peptide_diff_fig
      pca_peptide_diff_fig
    } else{
      db_execution$peptide_differential_PCA = NULL
    }
  })
  
  generate_protein_diff_barplot <- reactive(function(size_text){
    req(input$protein_diff_barplot)
    if(input$protein_diff_barplot){
      ploft_diff_number <- generate_differential_barplots(db_execution$differential_results,
                                                          data_type="protein", size_text=size_text)$plot
      db_execution$protein_differential_barplot = ploft_diff_number
      ploft_diff_number
    } else{
      db_execution$protein_differential_barplot = NULL
    }
  })
  
  generate_peptide_diff_barplot <- reactive(function(size_text){
    req(input$peptide_diff_barplot)
    if(input$peptide_diff_barplot){
      ploft_diff_number_pep <- generate_differential_barplots(db_execution$differential_results,
                                                              data_type="peptide", size_text=size_text)$plot
      db_execution$peptide_differential_barplot = ploft_diff_number_pep
      ploft_diff_number_pep
    } else{
      db_execution$peptide_differential_barplot = NULL
    }
  })
  
  generate_protein_upset <- reactive({
    # req(input$protein_upset)
    if(input$protein_upset){
      message("test")
      ploft_diff_number <- generate_upset_plot(db_execution$differential_results,
                                               type="protein", 
                                               DE_class = "all")$plot
      db_execution$protein_upset_plot = ploft_diff_number
      ploft_diff_number
    } else{
      message("strange")
      db_execution$protein_upset_plot = NULL
    }
  })
  
  generate_peptide_upset <- reactive({
    req(input$peptide_upset)
    if(input$peptide_upset){
      ploft_diff_number_pep <- generate_upset_plot(db_execution$differential_results,
                                                   type="peptide", 
                                                   DE_class = "all")$plot
      db_execution$peptide_upset_plot = ploft_diff_number_pep
      ploft_diff_number_pep
    } else{
      db_execution$peptide_upset_plot = NULL
    }
  })
  
  # PROTN: Execution pipeline ----
  observeEvent(input$report_proteome, {
    
    output$protn_results_ui <- renderUI({
      isolate({
        
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
              # Reset other analysis
              db_execution$parameter <- list()
              updateCheckboxInput(session, "differential_analysis_checkbox", value = FALSE)

              
              message(session$token)
              message(tempdir())
              #Creation directory for the results
              dirOutput_2 <- tempdir()
              currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
              dirOutput_1 <- paste("/", currentTime, "/", sep = "")
              dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
              dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
              message(dirOutput_Server)
              db_execution$dirOutput <- dirOutput_Server
              #Save folder for the download
              readr::write_csv(data.frame("session"=session$token,
                                          "outdir"=dirOutput_Server),
                               file = paste0(tempdir(),"/outdir_log_ProTN.log"), append = T)
              
              
              #Read parameter and execution
              software <- input$sw_analyzer
              file_input_proteome = input$input_file_proteome$name
              file_prot_proteome = if(software%in%c("PD","MQ_prot")){input$prot_file_proteome$name}else{NA}
              file_pep_proteome = input$pep_file_proteome$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_protn"), showWarnings = FALSE)
              dir_input <- paste(dirOutput_Server, "input_protn", sep = "")
              file.copy(from = input$input_file_proteome$datapath, to = paste0(dir_input,'/ANNOTATION_',file_input_proteome)) 
              if(software%in%c("PD","MQ_prot")){file.copy(from = input$prot_file_proteome$datapath, to =paste0(dir_input,'/PROT_',file_prot_proteome))} 
              file.copy(from = input$pep_file_proteome$datapath, to = paste0(dir_input,'/PEP_',file_pep_proteome)) 
              
              # If advance filter
              if(input$advance_filter){
                NA_allow_condition <- input$NA_allow_condition
                min_peptide_protein <- input$min_peptide_protein
                impute_algorithm <- unlist(tstrsplit(input$impute_algorithm, "_"))
                if(input$sample_column == "Sample"){
                  sample_column <- input$sample_column
                } else{
                  if(software=="PD"){
                    sample_column <- "File Name"
                  } else{
                    sample_column <- "Sample"
                  }
                }
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
                impute_algorithm <- c("norm","phosr")
                if(software=="PD"){
                  sample_column <- "File Name"
                } else{
                  sample_column <- "Sample"
                }
              }
              message("SAMPLE COLUMN:")
              message(as.character(sample_column))
              # If to batch corrected read column
              if(input$batch_correction){
                batch_corr <- TRUE
                batch_correction_col <- input$batch_correction_col
              } else{
                batch_corr <- FALSE
                batch_correction_col <- "batch"
              }
              
              message(software)
              progress=0
              msg_read_function <- NULL
              withCallingHandlers(
                {
                  shinyjs::html("text", "")
                  if(software == "PD"){
                    db_execution$proteome_data <- read_proteomics(software = "PD",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_",
                                                                  proteinGroup_filename = "PROT_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ_ev"){
                    db_execution$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ_prot"){
                    db_execution$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  proteinGroup_filename = "PROT_", 
                                                                  sample_col = sample_column,
                                                                  use_proteinGroups_MQ = TRUE,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "SP"){
                    db_execution$proteome_data <- read_proteomics(software = "SP",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "FP"){
                    db_execution$proteome_data <- read_proteomics(software = "FP",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  }
                },
                message = function(m) {
                  msg_read_function <<- append(msg_read_function, conditionMessage(m))
                  progress=progress+0.05
                  setProgress(value = progress)
                }
              )
              
              write_lines(msg_read_function, file = paste0(db_execution$dirOutput,"log_filter_read_function.txt"))
              db_execution$data_loaded <- TRUE
              
              if(impute_algorithm[1] != "norm"){
                message("Doing before imputation")
                message(impute_algorithm[1])
                db_execution$imputed_data <- impute_intensity(proteome_data = db_execution$proteome_data, type = impute_algorithm[1])
                db_execution$normalized_data <- normalization_ProTN(proteome_data = db_execution$imputed_data)
              } else{
                message("Doing before normalization")
                message(impute_algorithm[2])
                db_execution$normalized_data <- normalization_ProTN(proteome_data = db_execution$proteome_data)
                db_execution$normalized_data <- impute_intensity(proteome_data = db_execution$normalized_data, type = impute_algorithm[2])
              }
              
              if(batch_corr){
                message("Executing batch correction...")
                db_execution$normalized_data <- batch_correction(proteome_data = db_execution$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              db_execution$parameter<-list("Imputation and normalization algorithm: " = ifelse(impute_algorithm[1] != "norm", impute_algorithm[1], impute_algorithm[2]), 
                                           "Sample column in annotation file: " = sample_column, 
                                           "Batch correction: " = ifelse(batch_corr, batch_correction_col, "FALSE"), 
                                           "N° missing value allow per condition: " = NA_allow_condition, 
                                           "Minimum peptide per protein: " = min_peptide_protein)
              
              output$c_anno <- DT::renderDT(db_execution$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome", "Download results (ZIP file)", width = "240px")
                ),
                tags$h3("Statistics:"),
                tags$h4(paste0("Number of proteins: ", uniqueN(db_execution$normalized_data$dat_gene$GeneName))),
                if("dat_pep" %in% names(db_execution$normalized_data)){
                  tags$h4(paste0("Number of peptides: ", uniqueN(db_execution$normalized_data$dat_pep$ID_peptide)))
                },
                tags$h3("Annotation table"),
                DT::DTOutput("c_anno")
              )
            })
          },
          error = function(e) {
            #Create error report and reactivate the click in the page
            showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
            html_text<-str_replace(read_file("R/error.html"), 
                                   pattern = "The page you’re looking for doesn’t exist.</p>", 
                                   replacement = paste0("Description:", e, "</p>"))
            write_file(html_text, file = paste0(tempdir(), "/error.html"))
            tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
          }
        )
      })
      
    })
    
    output$render_abundance_plot <- renderUI({
      if (input$abundance_plot) {
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('abundance_plot')",
            plotOutput("small_abundance_plot")
          )
        )
      } else{
        db_execution$generate_abundance = NULL
      }
    })
    output$small_abundance_plot <- renderPlot({
      generate_abundance()
    })
    
    output$render_peptide_distribution <- renderUI({
      if (input$peptide_distribution) {
        tagList(
          tags$h3("N° peptides per proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('peptide_distribution_plot')",
            plotOutput("small_peptide_distribution")
          )
        )
      } else{
        db_execution$generate_peptide_distribution = NULL
      }
    })
    output$small_peptide_distribution <- renderPlot({
      generate_peptide_distribution()
    })
    
    
    output$render_raw_violin <- renderUI({
      if (input$raw_violin) {
        tagList(
          tags$h3("Distribution raw abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('raw_violin_plot')",
            plotOutput("small_raw_violin")
          )
        )
      } else{
        db_execution$raw_abundance_distribution = NULL
      }
    })
    output$small_raw_violin <- renderPlot({
      generate_raw_violin()
    })
    
    
    output$render_complexity_plot <- renderUI({
      if (input$complexity_plot) {
        tagList(
          tags$h3("Complexity plot of raw abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('complexity_plot')",
            plotOutput("small_complexity_plot")
          )
        )
      } else{
        db_execution$generate_complexity = NULL
      }
    })
    output$small_complexity_plot <- renderPlot({
      generate_complexity()
    })
    
    
    output$render_protein_violin <- renderUI({
      if (input$protein_violin) {
        tagList(
          tags$h3("Distribution protein abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('protein_violin_plot')",
            plotOutput("small_protein_violin")
          )
        )
      } else{
        db_execution$protein_abundance_distribution = NULL
      }
    })
    output$small_protein_violin <- renderPlot({
      generate_protein_violin()
    })
    
    output$render_peptide_violin <- renderUI({
      if (input$peptide_violin) {
        tagList(
          tags$h3("Distribution peptide abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('peptide_violin_plot')",
            plotOutput("small_peptide_violin")
          )
        )
      } else{
        db_execution$peptide_abundance_distirbution = NULL
      }
    })
    output$small_peptide_violin <- renderPlot({
      generate_peptide_violin()
    })
    
    output$render_mds_protein <- renderUI({
      if (input$mds_protein) {
        tagList(
          tags$h3("MDS based on proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('mds_protein')",
            plotOutput("small_mds_protein")
          )
        )
      } else{
        db_execution$protein_MDS = NULL
      }
    })
    output$small_mds_protein <- renderPlot({
      generate_mds_protein()
    })
    
    output$render_mds_peptide <- renderUI({
      if (input$mds_peptide) {
        tagList(
          tags$h3("MDS based on peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('mds_peptide')",
            plotOutput("small_mds_peptide")
          )
        )
      } else{
        db_execution$peptide_MDS = NULL
      }
    })
    output$small_mds_peptide <- renderPlot({
      generate_mds_peptide()
    })
    
    output$render_pca_protein <- renderUI({
      if (input$pca_protein) {
        tagList(
          tags$h3("PCA based on proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('pca_protein')",
            plotOutput("small_pca_protein")
          )
        )
      } else{
        db_execution$protein_PCA = NULL
      }
    })
    output$small_pca_protein <- renderPlot({
      generate_pca_protein()
    })
    
    output$render_pca_peptide <- renderUI({
      if (input$pca_peptide) {
        tagList(
          tags$h3("PCA based on peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('pca_peptide')",
            plotOutput("small_pca_peptide")
          )
        )
      } else{
        db_execution$peptide_PCA = NULL
      }
    })
    output$small_pca_peptide <- renderPlot({
      generate_pca_peptide()
    })
    
    output$render_protein_boxplot <- renderUI({
      if (input$boxplot_protein) {
        req(input$list_proteins)
        tagList(
          tags$h3("Boxplot selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('protein_boxplot')",
            plotOutput("small_protein_boxplot")
          )
        )
      } else{
        db_execution$protein_boxplot = NULL
      }
    })
    output$small_protein_boxplot <- renderPlot({
      generate_protein_boxplot()
    })
    
    output$render_protein_heatmap <- renderUI({
      if (input$heatmap_protein) {
        req(input$list_proteins)
        tagList(
          tags$h3("Heatmap selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('protein_heatmap')",
            plotOutput("small_protein_heatmap")
          )
        )
      } else{
        db_execution$protein_heatmap = NULL
      }
    })
    output$small_protein_heatmap <- renderPlot({
      generate_protein_heatmap()
    })
    
  })
  
  ## PROTN: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn, {
    output$render_differential_analysis <- renderUI({
      isolate({
        updateCheckboxInput(session, "enrichment_analysis", value = FALSE)
        updateCheckboxInput(session, "stringdb_analysis", value = FALSE)
        
        db_execution$dt_formule_contrast <- as.data.table(hot_to_r(input$render_formule_contrast_table))
        db_execution$dt_formule_contrast <- db_execution$dt_formule_contrast[Formule!=""]
        print(db_execution$dt_formule_contrast)
        formule_diff <- as.list(db_execution$dt_formule_contrast$Formule)
        names(formule_diff) <- stri_replace_all(db_execution$dt_formule_contrast$Name, replacement = "_", regex = "-")
        
        names(formule_diff) <- lapply(1:length(formule_diff), function(x){
          if(names(formule_diff)[x] == ""){
            stri_replace_all(formule_diff[[x]], replacement = "_VS_", regex = "-")
          } else{
            names(formule_diff)[x]
          }
        })
        db_execution$formule_contrast <- formule_diff
        message(db_execution$formule_contrast)

        withProgress(message = "Differential analysis in process, please wait!", {
          message(session$token)
          message(tempdir())

          db_execution$differential_results <- differential_analysis(proteome_data = db_execution$normalized_data,
                                                                     formule_contrast = db_execution$formule_contrast,
                                                                     fc_thr=as.double(input$FC_thr),
                                                                     pval_fdr = input$pval_fdr,
                                                                     pval_thr=as.double(input$pval_thr),
                                                                     signal_thr=0)
          db_execution$formule_contrast <- db_execution$formule_contrast[unique(union(db_execution$differential_results$protein_results_long$comp, 
                                                                                      db_execution$differential_results$peptide_results_long$comp))]
        })
        db_execution$parameter<-c(db_execution$parameter,
                                      "Fold change threshold for significance: "=input$FC_thr,
                                      "P.value type used: "=input$pval_fdr,
                                      "P.value threshold for significance: "=input$pval_thr)
        
        
        tags$h2("Differential Analysis")
      })
    })
    
    output$render_protein_diff_table <- renderUI({
      if(input$protein_diff_table){
        output$protein_results_long <- DT::renderDT(db_execution$differential_results$protein_results_long)
        DT::DTOutput("protein_results_long")
      }
    })
    
    output$render_peptide_diff_table <- renderUI({
      if(input$peptide_diff_table){
        output$peptide_results_long <- DT::renderDT(db_execution$differential_results$peptide_results_long)
        DT::DTOutput("peptide_results_long")
      }
    })
    
    output$render_protein_diff_barplot <- renderUI({
      if (input$protein_diff_barplot) {
        tagList(
          tags$h3("N° differential proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('protein_diff_barplot')",
            plotOutput("small_protein_diff_barplot")
          )
        )
      } else{
        db_execution$protein_differential_barplot = NULL
      }
    })
    output$small_protein_diff_barplot <- renderPlot({
      generate_protein_diff_barplot()(3)
    })
    
    output$render_peptide_diff_barplot <- renderUI({
      if (input$peptide_diff_barplot) {
        tagList(
          tags$h3("N° differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('peptide_diff_barplot')",
            plotOutput("small_peptide_diff_barplot")
          )
        )
      } else{
        db_execution$peptide_differential_barplot = NULL
      }
    })
    output$small_peptide_diff_barplot <- renderPlot({
      generate_peptide_diff_barplot()(3)
    })
    
    output$render_protein_upset <- renderUI({
      if (input$protein_upset) {
        message("plotting protein upset")
        tagList(
          tags$h3("Differential proteins upset plot"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('protein_upset_l')",
            plotOutput("small_protein_upset")
          )
        )
      } else{
        db_execution$protein_upset_plot = NULL
      }
    })
    output$small_protein_upset <- renderPlot({
      generate_protein_upset()
    })
    
    output$render_peptide_upset <- renderUI({
      if (input$peptide_upset) {
        tagList(
          tags$h3("Differential peptides upset plot"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('peptide_upset_l')",
            plotOutput("small_peptide_upset")
          )
        )
      } else{
        db_execution$peptide_upset_plot = NULL
      }
    })
    output$small_peptide_upset <- renderPlot({
      generate_peptide_upset()
    })
    
    output$render_protein_ma_plot <- renderUI({
      if (input$protein_ma_plot) {
        c_anno <- db_execution$proteome_data$c_anno
        generate_ma_plots_protein <- list()
        for(comp in names(db_execution$formule_contrast)){
          message(comp)
          design <- model.matrix(~0 + c_anno$condition)
          colnames(design) <- levels(as.factor(c_anno$condition))
          rownames(design) <- c_anno$sample
          
          conds <- as.data.table(makeContrasts(contrasts = db_execution$formule_contrast[[comp]], levels = design), keep.rownames = T)
          conds <- conds[as.vector(conds[,2]!=0), rn]
          message(conds)
          
          generate_ma_plots_protein[[comp]] <- ma_plot(differential_results = db_execution$differential_results, 
                                                    proteome_data = db_execution$normalized_data,
                                                    type="protein", comparison = comp, condition = conds)$plot
        }
        db_execution$protein_ma_plot = lapply(generate_ma_plots_protein, function(x){ggplotly(x, tooltip = c("text"))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_ma_plots_protein)) {
          plot_id <- paste0(names(generate_ma_plots_protein)[i], "_ma_prot")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            # output[[my_plot_id]] <- renderPlot(generate_ma_plots_protein[[names(generate_ma_plots_protein)[my_i]]])
            output[[my_plot_id]] <- renderPlotly(ggplotly(generate_ma_plots_protein[[names(generate_ma_plots_protein)[my_i]]], tooltip = c("text")))
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_ma_plots_protein)[i]),
            # plotOutput(plot_id)
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("MA Plot differential proteins"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_ma_protein"), tabs))
        )
      } else{
        db_execution$protein_ma_plot = NULL
      }
    })
    
    output$render_peptide_ma_plot <- renderUI({
      if (input$peptide_ma_plot) {
        c_anno <- db_execution$proteome_data$c_anno
        generate_ma_plots_peptide <- list()
        for(comp in names(db_execution$formule_contrast)){
          message(comp)
          design <- model.matrix(~0 + c_anno$condition)
          colnames(design) <- levels(as.factor(c_anno$condition))
          rownames(design) <- c_anno$sample
          
          conds <- as.data.table(makeContrasts(contrasts = db_execution$formule_contrast[[comp]], levels = design), keep.rownames = T)
          conds <- conds[as.vector(conds[,2]!=0), rn]
          message(conds)
          
          generate_ma_plots_peptide[[comp]] <- ma_plot(differential_results = db_execution$differential_results, 
                                                       proteome_data = db_execution$normalized_data,
                                                       type="peptide", comparison = comp, condition = conds)$plot
        }
        db_execution$peptide_ma_plot = lapply(generate_ma_plots_peptide, function(x){ggplotly(x, tooltip = c("text"))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_ma_plots_peptide)) {
          plot_id <- paste0(names(generate_ma_plots_peptide)[i], "_ma_pep")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(ggplotly(generate_ma_plots_peptide[[names(generate_ma_plots_peptide)[my_i]]], tooltip = "text"))
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_ma_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("MA Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_ma_peptide"), tabs))
        )
      } else{
        db_execution$peptide_ma_plot = NULL
      }
    })
  
    output$render_protein_vulcano <- renderUI({
      if(input$protein_vulcano){
        generate_volcano_plots_protein <- list()
        for(comp in names(db_execution$formule_contrast)){
          generate_volcano_plots_protein<-c(generate_volcano_plots_protein,
                                            generate_volcano_plots(db_execution$differential_results,
                                                                 data_type="protein",
                                                                 comparison=comp,
                                                                 fc_thr=as.double(input$FC_thr),
                                                                 pval_fdr = input$pval_fdr,
                                                                 pval_thr=as.double(input$pval_thr)))
        }
        db_execution$protein_vulcano = generate_volcano_plots_protein
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_volcano_plots_protein)) {
          plot_id <- paste0(names(generate_volcano_plots_protein)[i], "_prot")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(generate_volcano_plots_protein[[names(generate_volcano_plots_protein)[my_i]]])
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_volcano_plots_protein)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential proteins"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_protein"), tabs))
        )
        
      }else{
        db_execution$protein_vulcano = NULL
      }
    })
    
    output$render_peptide_vulcano <- renderUI({
      if(input$peptide_vulcano){
        generate_volcano_plots_peptide <- list()
        for(comp in names(db_execution$formule_contrast)){
          generate_volcano_plots_peptide<-c(generate_volcano_plots_peptide,
                                            generate_volcano_plots(db_execution$differential_results,
                                                                   data_type="peptide",
                                                                   comparison=comp,
                                                                   fc_thr=as.double(input$FC_thr),
                                                                   pval_fdr = input$pval_fdr,
                                                                   pval_thr=as.double(input$pval_thr)))
        }
        db_execution$peptide_vulcano = generate_volcano_plots_peptide
        # Generate tabPanels in a for loop
        tabs_pep_vulcano <- list()
        for (i in seq_along(generate_volcano_plots_peptide)) {
          plot_id <- paste0(names(generate_volcano_plots_peptide)[i], "_pep")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(generate_volcano_plots_peptide[[names(generate_volcano_plots_peptide)[my_i]]])
          })
          
          tabs_pep_vulcano[[i]] <- tabPanel(
            title = paste(names(generate_volcano_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_peptide"), tabs_pep_vulcano))
        )
      }else{
        db_execution$peptide_vulcano = NULL
      }
    })
    
    output$render_mds_protein_diff <- renderUI({
      if (input$mds_diff_protein) {
        tagList(
          tags$h3("MDS based on differential proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('mds_protein_diff')",
            plotOutput("small_mds_protein_diff")
          )
        )
      } else{
        db_execution$protein_differential_MDS = NULL
      }
    })
    output$small_mds_protein_diff <- renderPlot({
      generate_mds_protein_diff()
    })
    
    output$render_mds_peptide_diff <- renderUI({
      if (input$mds_diff_peptide) {
        tagList(
          tags$h3("MDS based on differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('mds_peptide_diff')",
            plotOutput("small_mds_peptide_diff")
          )
        )
      } else{
        db_execution$peptide_differential_MDS = NULL
      }
    })
    output$small_mds_peptide_diff <- renderPlot({
      generate_mds_peptide_diff()
    })
    
    output$render_pca_protein_diff <- renderUI({
      if (input$pca_diff_protein) {
        tagList(
          tags$h3("PCA based on differential proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('pca_protein_diff')",
            plotOutput("small_pca_protein_diff")
          )
        )
      } else{
        db_execution$protein_differential_PCA = NULL
      }
    })
    output$small_pca_protein_diff <- renderPlot({
      generate_pca_protein_diff()
    })
    
    output$render_pca_peptide_diff <- renderUI({
      if (input$pca_diff_peptide) {
        tagList(
          tags$h3("PCA based on differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot('pca_peptide_diff')",
            plotOutput("small_pca_peptide_diff")
          )
        )
      } else{
        db_execution$peptide_differential_PCA = NULL
      }
    })
    output$small_pca_peptide_diff <- renderPlot({
      generate_pca_peptide_diff()
    })
    
  })
  
  ## PROTN: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn, {
    output$render_enrichement_analysis <- renderUI({
      isolate({
        db_execution$enrichment_results <- perform_enrichment_analysis(differential_results = db_execution$differential_results,
                                                          enrichR_custom_DB = T,
                                                          enrich_filter_DBs=input$DB_enrichment,    
                                                          overlap_size_enrich_thr=as.double(input$os_enrich),
                                                          pval_fdr_enrich = input$pval_fdr,
                                                          pval_enrich_thr=as.double(input$pval_thr),
                                                          dirOutput=db_execution$dirOutput, 
                                                          with_background = input$enrich_with_background)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich), regex=","))
        
        db_execution$parameter <- c(db_execution$parameter,
                                       "Enrichment databases selected: "=paste(input$DB_enrichment, collapse = ", "),
                                       "P.value type used for enrichment: "=input$pval_fdr,
                                       "P.value threshold for enrichment significance: "=input$pval_thr,
                                       "Overlap size threshold for enrichment significance: "=input$os_enrich,
                                       "Enrichment filter terms: "=if(length(terms_enrich)>0){paste(terms_enrich, collapse = ", ")}else{"None"},
                                       "Enrichment with background: "=input$enrich_with_background)
        
        plots_down <- enrichment_figure(enr_df = db_execution$enrichment_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution$enrichment_results$anno_class), input$DB_enrichment))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(plots_down)) {
          plot_id <- names(plots_down)[i]
          height_id <- max(min(20, length(unique(plots_down[[names(plots_down)[i]]]$data$y_col))*0.4),3)*96
          message(paste0("Height for ",names(plots_down)[i], ": ", height_id))
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlot({
              plots_down[[names(plots_down)[my_i]]]
            }, height = height_id)
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(plots_down)[i]),
            plotOutput(plot_id, height = height_id)
          )
        }
        
        tagList(
          tags$h2("Enrichment Analysis"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_enrichment"), tabs))
        )
        
      })
    })
  })
  ## PROTN: stringdb analysis ----
  observeEvent(input$execute_stringdb_analysis_btn, {
    output$render_stringdb <- renderUI({
      isolate({
        withProgress(message = "STRINGdb analysis in process, please wait!", {
          
          db_execution$stringdb_res <- STRINGdb_network(differential_results = db_execution$differential_results,
                                                        species=input$taxonomy, 
                                                        dirOutput=db_execution$dirOutput, 
                                                        score_thr=input$score_thr_stringdb,
                                                        shiny = T)
          db_execution$parameter <- c(db_execution$parameter,
                                         "STRINGdb taxonomy: "=input$taxonomy,
                                         "STRINGdb score threshold: "=input$score_thr_stringdb)
            
          tagList(
            tags$h2("STRINGdb analysis"),
            fluidRow(
              selectInput("stringdb_show", label = "Select StringDB to show: (click on STRING logo to open the results on stringDB website)", 
                          choices = names(db_execution$stringdb_res), width = "15%"),
              actionButton("stringdb_selected", "Select!", width = "10%")  
            ),
            tags$div(id = "stringEmbedded")
          )
        })
      })
    })
  })
  
  observeEvent(input$stringdb_selected, {
    js$loadStringData(input$taxonomy, db_execution$stringdb_res[[input$stringdb_show]], input$score_thr_stringdb)
  })
  # PROTN: download results ----
  output$download_proteome <- downloadHandler(
    filename = "results_proTN.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            message(db_execution$dirOutput)
            setProgress(value = 0.01)
            
            # Generate report
            params <- list(
              doc_title = input$title_exp,
              description = input$description_exp,
              readPD_files = if (input$sw_analyzer == "PD") {TRUE} else {FALSE},
              readMQ_files = if (input$sw_analyzer == "MQ_ev") {TRUE} else {FALSE},
              readMQ_prot_files = if (input$sw_analyzer == "MQ_prot") {TRUE} else {FALSE},
              impute_algorithm = if(input$advance_filter){input$impute_algorithm} else {"norm_phosr"},
              db_execution = reactiveValuesToList(db_execution),
              file_input = paste(db_execution$dirOutput, "input_protn", sep = ""),
              batch_corr_exe = if(input$batch_correction){input$batch_correction_col}else{NULL},
              prot_boxplot = if(input$boxplot_protein | input$heatmap_protein){input$list_proteins}else{NULL},
              fc_thr = if(is.null(input$FC_thr)){"0.75"}else{input$FC_thr},
              pval_fdr = input$pval_fdr,
              pval_thr = if(is.null(input$pval_thr)){"0.05"}else{input$pval_thr},
              pval_fdr_enrich = input$pval_fdr_enrich,
              pval_enrich_thr = if(is.null(input$pvalue_enrich)){"0.05"}else{input$pvalue_enrich},
              overlap_size_enrich_thr = if(is.null(input$os_enrich)){as.integer(5)}else{input$os_enrich},
              enrich_filter_term = input$terms_enrich,
              enrich_filter_DBs = input$DB_enrichment,
              taxonomy=input$taxonomy, 
              score_thr=input$score_thr_stringdb,
              dirOutput = db_execution$dirOutput
            )
            
            # Render in background the report
            p = callr::r_bg(
              func = function(db_execution, params, dirOutput, env) {
                rmarkdown::render("R/protn_report.Rmd",
                                  output_file = "protn_report.html",
                                  output_dir = dirOutput,
                                  params = params,
                                  envir = env
                )
              },
              args = list(db_execution, params, db_execution$dirOutput, new.env(parent = globalenv())),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            # Render in background the report
            # p_pdf = callr::r_bg(
            #   func = function(db_execution, params, dirOutput, env) {
            #     rmarkdown::render("R/protn_report_pdf.Rmd",
            #                       output_file = "protn_report.pdf",
            #                       output_dir = dirOutput,
            #                       params = params,
            #                       envir = env
            #     )
            #   },
            #   args = list(db_execution, params, db_execution$dirOutput, new.env(parent = globalenv())),
            #   stdout = "|",
            #   stderr = "|",
            #   error = getOption("callr.error", "error")
            # )
            
            # Saving RData in background
            db_results = reactiveValuesToList(db_execution)
            db_results <- db_results[!(unlist(lapply(db_results, is.null)))]
            p_rdata = callr::r_bg(
              func = function(db_results, dirOutput) {
                save(db_results, file = paste0(dirOutput,"db_results.RData"))
              },
              args = list(db_results, db_results$dirOutput),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            
            # Prepare file for the download
            if(length(db_execution$normalized_data)>0){
              save_abundance_tables(proteome_data = db_execution$normalized_data, 
                                    dirOutput = db_execution$dirOutput)
            }
            setProgress(value = 0.1)
            
            if(length(db_execution$differential_results)>0){
              save_differential_analysis_table(proteome_data = db_execution$normalized_data,
                                               differential_results = db_execution$differential_results,
                                               dirOutput=db_execution$dirOutput)
            }
            setProgress(value = 0.2)
            
            if(input$abundance_plot & !is.null(db_execution$generate_abundance)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/missing_available_abundance.pdf"), 
                     plot = db_execution$generate_abundance, 
                     create.dir = T, width = 7, height = 5)
            } else if("missing_available_abundance.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/missing_available_abundance.pdf"))
            }
            setProgress(value = 0.25)
            
            if(input$peptide_distribution & !is.null(db_execution$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_per_protein.pdf"), 
                     plot = db_execution$generate_peptide_distribution, 
                     create.dir = T, width = 7, height = 5)
            }  else if("peptide_per_protein.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_per_protein.pdf"))
            }
            setProgress(value = 0.30)
            
            if(input$raw_violin & !is.null(db_execution$raw_abundance_distribution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/raw_abundance_distribution.pdf"), 
                     plot = db_execution$raw_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("raw_abundance_distribution.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/raw_abundance_distribution.pdf"))
            }
            
            if(input$complexity_plot & !is.null(db_execution$generate_complexity)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/complexity_plot.pdf"), 
                     plot = db_execution$generate_complexity, 
                     create.dir = T, width = 10, height = 8)
            } else if("complexity_plot.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/complexity_plot.pdf"))
            }
            setProgress(value = 0.33)
            
            if(input$protein_violin & !is.null(db_execution$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_abundance_distribution.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_abundance_distribution.pdf"))
            }
            setProgress(value = 0.35)
            
            if(input$peptide_violin & !is.null(db_execution$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_abundance_distribution.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_abundance_distribution.pdf"))
            }
            setProgress(value = 0.40)
            
            if(input$mds_protein & !is.null(db_execution$protein_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_MDS.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_MDS.pdf"))
            }
            setProgress(value = 0.43)
            
            if(input$mds_peptide & !is.null(db_execution$peptide_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_MDS.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_MDS.pdf"))
            }
            setProgress(value = 0.45)
            
            if(input$pca_protein & !is.null(db_execution$protein_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_PCA.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_PCA.pdf"))
            }
            setProgress(value = 0.47)
            
            if(input$pca_peptide & !is.null(db_execution$peptide_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_PCA.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_PCA.pdf"))
            }
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein & !is.null(db_execution$protein_boxplot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_boxplot.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_boxplot.pdf"))
            }
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein & !is.null(db_execution$protein_heatmap)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_heatmap.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_heatmap.pdf"))
            }
            setProgress(value = 0.55)
            
            if(!is.null(db_execution$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution$protein_differential_barplot, 
                     create.dir = T,width = 17, height = 6)
            } else if("protein_differential_barplot.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_differential_barplot.pdf"))
            }
            setProgress(value = 0.58)
            
            if(!is.null(db_execution$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution$peptide_differential_barplot, 
                     create.dir = T, width = 17, height = 6)
            } else if("peptide_differential_barplot.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_differential_barplot.pdf"))
            }
            setProgress(value = 0.60)
            
            if(!is.null(db_execution$protein_upset_plot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_upset_plot.pdf"), 
                     plot = db_execution$protein_upset_plot, 
                     create.dir = T,width = 12, height = 6)
            } else if("protein_upset_plot.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_upset_plot.pdf"))
            }
            setProgress(value = 0.62)
            
            if(!is.null(db_execution$peptide_upset_plot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_upset_plot.pdf"), 
                     plot = db_execution$peptide_upset_plot, 
                     create.dir = T, width = 12, height = 6)
            } else if("peptide_upset_plot.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_upset_plot.pdf"))
            }
            setProgress(value = 0.63)
            
            if(!is.null(db_execution$protein_ma_plot)){
              dir.create(file.path(paste0(db_execution$dirOutput,"pics/"), "protein_ma_plot"), showWarnings = FALSE)
              for(comp in names(db_execution$protein_ma_plot)){
                htmlwidgets::saveWidget(db_execution$protein_ma_plot[[comp]], 
                                        file = paste0(db_execution$dirOutput,"pics/protein_ma_plot/",comp,"_protein_ma_plot.html"))
                webshot2::webshot(url = paste0(db_execution$dirOutput,"pics/protein_ma_plot/",comp,"_protein_ma_plot.html"), 
                                  file = paste0(db_execution$dirOutput,"pics/protein_ma_plot/",comp,"_protein_ma_plot.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution$dirOutput,"pics/protein_ma_plot"))
            }
            setProgress(value = 0.64)
            
            
            if(!is.null(db_execution$peptide_ma_plot)){
              dir.create(file.path(paste0(db_execution$dirOutput,"pics/"), "peptide_ma_plot"), showWarnings = FALSE)
              for(comp in names(db_execution$peptide_ma_plot)){
                htmlwidgets::saveWidget(db_execution$peptide_ma_plot[[comp]], 
                                        file = paste0(db_execution$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"))
                webshot2::webshot(url = paste0(db_execution$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"), 
                                  file = paste0(db_execution$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution$dirOutput,"pics/peptide_ma_plot"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution$protein_vulcano)){
              dir.create(file.path(paste0(db_execution$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution$protein_vulcano)){
                # plotly::save_image(db_execution$protein_vulcano[[comp]], 
                #                         file = paste0(str_replace_all(db_execution$dirOutput, pattern="\\\\", replacement="/"),"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution$protein_vulcano[[comp]], 
                                   file = paste0(db_execution$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"), 
                                  file = paste0(db_execution$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution$dirOutput,"pics/protein_vulcano"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution$peptide_vulcano)){
                # plotly::save_image(db_execution$peptide_vulcano[[comp]], 
                #              file = paste0(str_replace_all(db_execution$dirOutput, pattern="\\\\", replacement="/"),"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution$peptide_vulcano[[comp]], 
                                   file = paste0(db_execution$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"), 
                                  file = paste0(db_execution$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution$dirOutput,"pics/peptide_vulcano"))
            }
            setProgress(value = 0.68)
            
            if(!is.null(db_execution$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_differential_MDS.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_differential_MDS.pdf"))
            }
            setProgress(value = 0.69)
            
            if(!is.null(db_execution$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_MDS.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_differential_MDS.pdf"))
            }
            setProgress(value = 0.70)
            
            if(!is.null(db_execution$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_differential_PCA.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/protein_differential_PCA.pdf"))
            }
            setProgress(value = 0.72)
            
            if(!is.null(db_execution$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_PCA.pdf" %in% list.files(paste0(db_execution$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution$dirOutput,"pics/peptide_differential_PCA.pdf"))
            }
            setProgress(value = 0.75)
            
            if(length(db_execution$enrichment_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution$enrichment_results,
                                              category = c("down","up"), 
                                              enrich_filter_term = terms_enrich,
                                              save=T, 
                                              dirOutput = db_execution$dirOutput)
            } 
            setProgress(value = 0.82)
            
            if(length(db_execution$stringdb_res)>0){
              tmp_res <- STRINGdb_network(differential_results = db_execution$differential_results,
                                                            species=input$taxonomy, 
                                                            dirOutput=db_execution$dirOutput,
                                                            score_thr=input$score_thr_stringdb,
                                                            shiny = F)
              
            } 
            setProgress(value = 0.95)
            
            # Write tsv file with parameter
            params <- data.table("Parameter" = names(db_execution$parameter),
                                 "Value" = unlist(db_execution$parameter))
            fwrite(params, paste0(db_execution$dirOutput,"parameters_used.txt"), sep = "\t", col.names = F)
            
            #Get results Report
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p$read_output()
            p$wait(30000)
            for (i in 1:15) {
              p$read_output()
              p$wait(1000*60)  
            }
            if(p$is_alive() | is.null(p$get_result())){
              p$kill()
              print("\n ERROR: An error occur during the report rendering. \n ")
            } else{
              report<-p$get_result()
              p$kill()
              message("Render report DONE.")
            }
            
            #Get results Report
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            # hide_res<-p_pdf$read_output()
            # p_pdf$wait(30000)
            # for (i in 1:15) {
            #   p_pdf$read_output()
            #   p_pdf$wait(1000*60)  
            # }
            # if(p_pdf$is_alive() | is.null(p_pdf$get_result())){
            #   p_pdf$kill()
            #   print("\n ERROR: An error occur during the report rendering. \n ")
            # } else{
            #   report<-p_pdf$get_result()
            #   p_pdf$kill()
            #   message("Render report DONE.")
            # }
            
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p_rdata$read_output()
            p_rdata$wait(30000)
            for (i in 1:15) {
              p_rdata$read_output()
              p_rdata$wait(1000*60)
            }
            if(p_rdata$is_alive() | is.null(p_rdata$get_result())){
              p_rdata$kill()
              print("\n ERROR: An error occur during the RData saving. \n ")
            } else{
              report<-p_rdata$get_result()
              p_rdata$kill()
              message("RData saved.")
            }
            
            # Save RData db_execution
            # setProgress(value = 0.95, message = "Saving RData...")
            # db_results_proTN = reactiveValuesToList(db_execution)
            # db_results_proTN <- db_results_proTN[!(unlist(lapply(db_results_proTN, is.null)))]
            # save(db_results_proTN, file = paste0(db_results_proTN$dirOutput,"db_results_proTN.RData"))
            # 
            #Save folder for the download
            oldwd <- getwd()
            message(db_execution$dirOutput)
            setwd(db_execution$dirOutput)
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  ## PROTN: full screen trigger ----
  
  # ReactiveVal for currently selected plot to fullscreen
  selected_plot <- reactiveVal(NULL)
  
  # Update selected_plot when JS sends fullscreen_trigger id
  observeEvent(input$fullscreen_trigger, {
    selected_plot(input$fullscreen_trigger)
  })
  
  # Render fullscreen plot dynamically based on selected_plot()
  output$fullscreen_plot <- renderPlot({
    req(selected_plot())
    message(selected_plot())
    switch(selected_plot(),
           "abundance_plot" = generate_abundance() + ggtitle("Percentage missing values respect detected abundance")+theme(text=element_text(size=25)),
           "peptide_distribution_plot" = generate_peptide_distribution() + ggtitle("N° peptides per proteins")+theme(text=element_text(size=25)),
           "raw_violin_plot" = generate_raw_violin() + ggtitle("Raw abundance distribution")+theme(text=element_text(size=25)),
           "complexity_plot" = generate_complexity() + ggtitle("Complexity plot of raw abundance")+theme(text=element_text(size=25)),
           "protein_violin_plot" = generate_protein_violin() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "peptide_violin_plot" = generate_peptide_violin() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "mds_protein" = generate_mds_protein() + ggtitle("MDS based on protein")+theme(text=element_text(size=25)),
           "mds_peptide" = generate_mds_peptide() + ggtitle("MDS based on peptides")+theme(text=element_text(size=25)),
           "pca_protein" = generate_pca_protein() + ggtitle("PCA based on protein")+theme(text=element_text(size=25)),
           "pca_peptide" = generate_pca_peptide() + ggtitle("PCA based on peptides")+theme(text=element_text(size=25)),
           "protein_boxplot" = generate_protein_boxplot() + ggtitle("Boxplot selected proteins")+theme(text=element_text(size=25)),
           "protein_heatmap" = generate_protein_heatmap() + ggtitle("Heatmap selected proteins")+theme(text=element_text(size=25)),
           "protein_diff_barplot" = generate_protein_diff_barplot()(8) + ggtitle("N° differential proteins")+theme(text=element_text(size=25)),
           "peptide_diff_barplot" = generate_peptide_diff_barplot()(8) + ggtitle("N° differential peptides")+theme(text=element_text(size=25)),
           "protein_upset_l" = generate_protein_upset(),
           "peptide_upset_l" = generate_peptide_upset(),
           "protein_ma_plot" = generate_protein_ma_plot() + ggtitle("Differential proteins MA plot")+theme(text=element_text(size=25)),
           "peptide_ma_plot" = generate_peptide_ma_plot() + ggtitle("Differential peptides MA plot")+theme(text=element_text(size=25)),
           "mds_protein_diff" = generate_mds_protein_diff() + ggtitle("MDS based on differential protein")+theme(text=element_text(size=25)),
           "mds_peptide_diff" = generate_mds_peptide_diff() + ggtitle("MDS based on differential peptides")+theme(text=element_text(size=25)),
           "pca_protein_diff" = generate_pca_protein_diff() + ggtitle("PCA based on differential protein")+theme(text=element_text(size=25)),
           "pca_peptide_diff" = generate_pca_peptide_diff() + ggtitle("PCA based on differential peptides")+theme(text=element_text(size=25)),
           # default fallback:
           NULL
    )
  })
  
  ##############################################################################
  ### PHOSPROTN ----
  # Optional visibility based on the selection ----
  
  ## PHOSPROTN: Visibility of the proteomics files for PHOSPROTN ----
  output$input_proteome_phos <- renderUI({
    if (input$sw_analyzer_phos == "PD"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_phos", "Select the SAMPLE_ANNOTATION file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_phos", "Select the PEP file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_proteome_phos", "Select the PROT file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("psm_file_proteome_phos", "Select the PSM file of the PHOSPHO-PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer_phos == "MQ"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_phos", "Select the SAMPLE_ANNOTATION file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_phos", "Select the EVIDENCE file of the PHOSPHO-PROTEOMICS..."),
        )
      )
    } else{
      tagList(
        tags$p("BACK")
      )
    }
  })
  
  ## PHOSPROTN: textbox for batch correction----
  output$batch_correction_ui_phos <- renderUI({ 
    if(input$batch_correction_phos){
      textInput("batch_correction_col_phos", "Column in Annotation file with the batch:")
    } 
  })
  
  ## PHOSPROTN: advance filters----
  output$advance_filter_ui_phos <- renderUI({ 
    if(input$advance_filter_phos){
      tagList(
        numericInput("NA_allow_condition_phos", "N° missing value allow per condition", value = 0, min = 0, max = 5),
        numericInput("min_peptide_protein_phos", "Minimum peptide per protein", value = 1, min = 1),
        selectizeInput("impute_algorithm_phos", "Select impute algorithm:",
                       choices = list("PhosR and normalization" = "phosr_norm", 
                                      "Gaussian estimation and normalization" = "gaussian_norm",
                                      "missForest and normalization" = "missForest_norm",
                                      "Pre-normalization and PhosR" = "norm_phosr", 
                                      "Pre-normalization and Gaussian estimation" = "norm_gaussian",
                                      "Pre-normalization and missForest" = "norm_missForest",
                                      "Pre-normalization and pcaMethods" = "norm_pcaMethods"),
                       selected = "norm_phosr", multiple = FALSE
        ),
        textInput("sample_column_phos", "Column name with the sample name:", value = "Sample")
      )
    } 
  })
  ## PHOSPROTN: textbox for list proteins ----
  output$list_protein_ui_phos <- renderUI({ 
    if(input$boxplot_protein_phos | input$heatmap_protein_phos){
      textInput("list_proteins_phos", "List proteins to show (separate by: \",\"):")
    } 
  })
  
  ## PHOSPROTN: show parameter for differential analysis ----
  output$differential_params_ui_phos <- renderUI({ 
    if(input$differential_analysis_checkbox_phos){
      tagList(
        tags$label("Write in each line a different comparison"),
        tags$label("(right click to add row)"),
        rHandsontableOutput('render_formule_contrast_table_phos'),
        textInput("FC_thr_phos", "Fold change threshold for significance:",value = 0.5),
        radioButtons("pval_fdr_phos", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_val"),
        textInput("pval_thr_phos", "P.value threshold for significance:", value = 0.05),
        actionButton("execute_differential_analysis_btn_phos", "Run!"),
        checkboxInput("peptide_diff_table_phos", "Peptides differentiated table", FALSE),
        checkboxInput("peptide_diff_barplot_phos", "Peptides differentiated barplot", TRUE),
        checkboxInput("peptide_upset_phos", "Peptides upset plot", FALSE),
        checkboxInput("peptide_ma_plot_phos", "Peptides MA plot", FALSE),
        checkboxInput("peptide_vulcano_phos", "Peptides vulcano plot", FALSE),
        checkboxInput("mds_diff_peptide_phos", "MDS based on diffential peptide", FALSE),
        checkboxInput("pca_diff_peptide_phos", "PCA based on diffential peptide", FALSE),
        tags$h3("Enrichment Analysis:"),
        checkboxInput("enrichment_analysis_phos", "Execute enrichment analysis", FALSE),
        uiOutput("enrichment_params_ui_phos"),
        tags$h3("STRINGdb network:"),
        checkboxInput("stringdb_analysis_phos", "Execute STRINGdb", FALSE),
        uiOutput("stringdb_params_ui_phos"),
        tags$h3("Kinase tree:"),
        checkboxInput("kinase_tree_analysis_phos", "Execute PhosR kinase tree", FALSE),
        uiOutput("kinase_tree_params_ui_phos")
      )
    } else{
      # Reset UI elements
      updateCheckboxInput(session, "enrichment_analysis_phos", value = FALSE)
      updateCheckboxInput(session, "stringdb_analysis_phos", value = FALSE)
      updateCheckboxInput(session, "kinase_tree_analysis_phos", value = FALSE)
      
      db_execution_phos$formule_contrast <- list()
      db_execution_phos$dt_formule_contrast <- data.table("Name"=c("","","",""),"Formule"=c("","","",""))
      db_execution_phos$differential_results <- list()
      
      updateCheckboxInput(session, "protein_diff_barplot_phos", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_barplot_phos",  value = FALSE)
      updateCheckboxInput(session, "protein_diff_table_phos", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_table_phos", value = FALSE)
      updateCheckboxInput(session, "protein_upset_phos", value = FALSE)
      updateCheckboxInput(session, "peptide_upset_phos", value = FALSE)
      updateCheckboxInput(session, "protein_vulcano_phos",  value =  FALSE)
      updateCheckboxInput(session, "peptide_vulcano_phos", value = FALSE)
      updateCheckboxInput(session, "protein_ma_plot_phos", value = FALSE)
      updateCheckboxInput(session, "peptide_ma_plot_phos", value = FALSE)
      updateCheckboxInput(session, "mds_diff_protein_phos", value = FALSE)
      updateCheckboxInput(session, "mds_diff_peptide_phos", value = FALSE)
      updateCheckboxInput(session, "pca_diff_protein_phos", value = FALSE)
      updateCheckboxInput(session, "pca_diff_peptide_phos", value = FALSE)
      
      db_execution_phos$protein_differential_barplot <- NULL
      db_execution_phos$peptide_differential_barplot <- NULL
      db_execution_phos$protein_upset_plot <- NULL
      db_execution_phos$peptide_upset_plot <- NULL
      db_execution_phos$protein_ma_plot <- NULL
      db_execution_phos$peptide_ma_plot <- NULL
      db_execution_phos$protein_vulcano <- NULL
      db_execution_phos$peptide_vulcano <- NULL
      db_execution_phos$protein_differential_MDS <- NULL
      db_execution_phos$peptide_differential_MDS <- NULL
      db_execution_phos$protein_differential_PCA <- NULL
      db_execution_phos$peptide_differential_PCA <- NULL
      
      output$render_differential_analysis_phos <- renderUI({NULL})
      output$render_protein_diff_table_phos <- renderUI({NULL})
      output$render_peptide_diff_table_phos <- renderUI({NULL})
      output$render_protein_diff_barplot_phos <- renderUI({NULL})
      output$render_peptide_diff_barplot_phos <- renderUI({NULL})
      output$render_protein_upset_phos <- renderUI({NULL})
      output$render_peptide_upset_phos <- renderUI({NULL})
      output$render_protein_ma_plot_phos <- renderUI({NULL})
      output$render_peptide_ma_plot_phos <- renderUI({NULL})
      output$render_protein_vulcano_phos <- renderUI({NULL})
      output$render_peptide_vulcano_phos <- renderUI({NULL})
      output$render_mds_protein_diff_phos <- renderUI({NULL})
      output$render_mds_peptide_diff_phos <- renderUI({NULL})
      output$render_pca_protein_diff_phos <- renderUI({NULL})
      output$render_pca_peptide_diff_phos <- renderUI({NULL})
    }
  })
  
  output$render_formule_contrast_table_phos <- renderRHandsontable({
    rhandsontable(db_execution_phos$dt_formule_contrast, rowHeaders = NULL, stretchH = "all")
  })
  
  ## PHOSPROTN: show enrichment parameter ----
  output$enrichment_params_ui_phos <- renderUI({ 
    if(input$enrichment_analysis_phos){
      tagList(
        # radioButtons("enrichR_universe", "Execute enrichment of the whole Universe", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
        selectizeInput("DB_enrichment_phos", "DB to analyse:",
                       choices = lapply(split(read_tsv("data/dbs_enrichR.txt", col_names = FALSE)$X1,
                                              read_tsv("data/dbs_enrichR.txt", col_names = FALSE)[,2]), as.list),
                       selected = NULL, multiple = TRUE
        ),
        textInput("terms_enrich_phos", "Terms to search (separated by \",\"):"),
        radioButtons("pval_fdr_enrich_phos", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_adj"),
        textInput("pvalue_enrich_phos", "P.value threshold for significance:", value = 0.05),
        sliderInput("os_enrich_phos", "Overlap size thr for enrichment", 1, 30, step = 1, value = 5),
        checkboxInput("enrich_with_background_phos", "Enrichment with background", FALSE),
        actionButton("execute_enrichment_analysis_btn_phos", "Run!")
      )
    } else{
      db_execution_phos$enrichment_results <- list()
      output$render_enrichement_analysis_phos <- renderUI({NULL})
    }
  })
  
  ## PHOSPROTN: show stringdb parameter ----
  output$stringdb_params_ui_phos <- renderUI({
    if(input$stringdb_analysis_phos){
      tagList(
        selectizeInput("taxonomy_phos", "NCBI Taxonomy ID", 
                       choice = data.table::fread("data/subset_tax.csv", select = "name"), 
                       selected = "Homo sapiens", multiple = F),
        sliderInput("score_thr_stringdb_phos", "Score thr for STRINGdb", 500, 1000, step = 10, value = 700),
        actionButton("execute_stringdb_analysis_btn_phos", "Run!"),
        tags$br()
      )
    } else{
      db_execution_phos$stringdb_res <- list()
      output$render_stringdb_phos <- renderUI({NULL})
    }
  })
  
  ## PHOSPROTN: show kinase_tree parameter ----
  output$kinase_tree_params_ui_phos <- renderUI({
    if(input$kinase_tree_analysis_phos){
      tax_sel <- if(is.null(input$taxonomy_phos)){
        character(0)
      } else if(input$taxonomy_phos == "Homo sapiens"){
        "Homo sapiens"
      } else if(input$taxonomy_phos == "Mus musculus"){
        "Mus musculus"
      } else{
        NULL
      }
      
      if(is.null(tax_sel)){
        shinyalert::shinyalert("Kinase analysis", 
                               "Kinase analysis can be performed only for Homo Sapiens (Human) or Mus Musculus (Mouse)", 
                               type = "info")
        updateCheckboxInput(session, "kinase_tree_analysis_phos", value = FALSE)
      } else{
        tagList(
          radioButtons("taxonomy_kinase_phos", "Select species (CORAL tree will be print only for Homo sapiens)", 
                       choiceNames = c("Homo Sapiens", "Mus Musculus"),
                       choiceValues = c("Homo sapiens","Mus musculus"), inline = TRUE, 
                       selected = tax_sel),
          sliderInput("score_thr_phosr_phos", "Score thr for PhosR", 0, 1, step = 0.05, value = 0.7),
          actionButton("execute_kinase_tree_analysis_btn_phos", "Run!"),
          tags$br()
        )
      }
    }else{
      db_execution_phos$kinase_tree_res <- list()
      output$render_kinase_tree_phos <- renderUI({NULL})
    }
  })
  
  ## PHOSPROTN: function genereting plot ----
  generate_phospho_percentage_plot_phos <- reactive(function(size_text,zoom) {
    req(input$phospho_percentage_plot_phos)
    if(input$phospho_percentage_plot_phos){
      generate_abundance_fig <- create_phosphosite_plot(proteome_data = db_execution_phos$proteome_data, 
                                                        software = input$sw_analyzer_phos, 
                                                        size_text = size_text)$plot
      if(!zoom){
        db_execution_phos$phospho_percentage = generate_abundance_fig
      }
      generate_abundance_fig
    } else{
      db_execution_phos$phospho_percentage = NULL
    }
  })
  
  generate_abundance_phos <- reactive({
    req(input$abundance_plot_phos)
    if(input$abundance_plot_phos){
      generate_abundance_fig <- generate_abundance_plot(proteome_data = db_execution_phos$proteome_data)$plot
      db_execution_phos$generate_abundance = generate_abundance_fig
      generate_abundance_fig
    } else{
      db_execution_phos$generate_abundance = NULL
    }
  })
  
  generate_raw_violin_phos <- reactive({
    req(input$raw_violin_phos)
    if(input$raw_violin_phos){
      raw_abundance_distribution_fig <- plot_raw_abundance_distribution(proteome_data = db_execution_phos$proteome_data,
                                                                        type = "protein")$plot
      db_execution_phos$raw_abundance_distribution = raw_abundance_distribution_fig
      raw_abundance_distribution_fig
    } else{
      db_execution_phos$raw_abundance_distribution = NULL
    }
  })
  
  generate_complexity_phos <- reactive({
    req(input$complexity_plot_phos)
    if(input$complexity_plot_phos){
      generate_complexity_fig <- complexity_plot(proteome_data = db_execution_phos$proteome_data)$plot
      db_execution_phos$generate_complexity = generate_complexity_fig
      generate_complexity_fig
    } else{
      db_execution_phos$generate_complexity = NULL
    }
  })
  
  generate_peptide_distribution_phos <- reactive({
    req(input$peptide_distribution_phos)
    if(input$peptide_distribution_phos){
      peptide_distribution_fig <- generate_peptide_distribution_plot(proteome_data = db_execution_phos$proteome_data)$plot
      db_execution_phos$generate_peptide_distribution = peptide_distribution_fig
      peptide_distribution_fig
    } else{
      db_execution_phos$generate_peptide_distribution = NULL
    }
  })
  
  generate_peptide_violin_phos <- reactive({
    req(input$peptide_violin_phos)
    if(input$peptide_violin_phos){
      peptide_abundance_distirbution_fig <- plot_abundance_distribution(proteome_data = db_execution_phos$normalized_data,
                                                                        type = "peptide")$plot
      db_execution_phos$peptide_abundance_distirbution = peptide_abundance_distirbution_fig
      peptide_abundance_distirbution_fig
    } else{
      db_execution_phos$peptide_abundance_distirbution = NULL
    }
  })
  
  generate_mds_peptide_phos <- reactive({
    req(input$mds_peptide_phos)
    if(input$mds_peptide_phos){
      mds_peptide_fig <- mds_plot(proteome_data = db_execution_phos$normalized_data,
                                  type = "peptide")$plot
      db_execution_phos$peptide_MDS = mds_peptide_fig
      mds_peptide_fig
    } else{
      db_execution_phos$peptide_MDS = NULL
    }
  })
  
  generate_pca_peptide_phos <- reactive({
    req(input$pca_peptide_phos)
    if(input$pca_peptide_phos){
      pca_peptide_fig <- pca_plot(proteome_data = db_execution_phos$normalized_data,
                                  type = "peptide")$plot
      db_execution_phos$peptide_PCA = pca_peptide_fig
      pca_peptide_fig
    } else{
      db_execution_phos$peptide_PCA = NULL
    }
  })
  
  generate_protein_boxplot_phos <- reactive({
    req(input$boxplot_protein_phos)
    if(input$boxplot_protein_phos){
      req(input$list_proteins_phos)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos), regex=",")
      db_execution_phos$parameter <- c(db_execution_phos$parameter, "List proteins boxplot abundance: "=input$list_proteins_phos)
      boxplot_protein_fig <- plot_selected_proteins(proteome_data = db_execution_phos$normalized_data,
                                                    list_protein = unlist(list_proteins))$plot
      db_execution_phos$protein_boxplot = boxplot_protein_fig
      boxplot_protein_fig
    } else{
      db_execution_phos$protein_boxplot = NULL
    }
  })
  
  generate_protein_heatmap_phos <- reactive({
    req(input$heatmap_protein_phos)
    if(input$heatmap_protein_phos){
      req(input$list_proteins_phos)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos), regex=",")
      db_execution_phos$parameter <- c(db_execution_phos$parameter, "List proteins heatmap abundance: "=input$list_proteins_phos)
      heatmap_protein_fig <- heatmap_selected_proteins(proteome_data = db_execution_phos$normalized_data, list_protein = unlist(list_proteins))$plot
      db_execution_phos$protein_heatmap = heatmap_protein_fig
      heatmap_protein_fig
    } else{
      db_execution_phos$protein_heatmap = NULL
    }
  })
  
  # generate_mds_protein_diff_phos <- reactive({
  #   req(input$mds_diff_protein_phos)
  #   if(input$mds_diff_protein_phos){
  #     mds_protein_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
  #                                                            proteome_data = db_execution_phos$normalized_data,
  #                                                            type = "protein")$plot
  #     db_execution_phos$protein_differential_MDS = mds_protein_diff_fig
  #     mds_protein_diff_fig
  #   } else{
  #     db_execution_phos$protein_differential_MDS = NULL
  #   }
  # })
  # 
  generate_mds_peptide_diff_phos <- reactive({
    req(input$mds_diff_peptide_phos)
    if(input$mds_diff_peptide_phos){
      mds_peptide_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
                                                             proteome_data = db_execution_phos$normalized_data,
                                                             type = "peptide")$plot
      db_execution_phos$peptide_differential_MDS = mds_peptide_diff_fig
      mds_peptide_diff_fig
    } else{
      db_execution_phos$peptide_differential_MDS = NULL
    }
  })
  
  # generate_pca_protein_diff_phos <- reactive({
  #   req(input$pca_diff_protein_phos)
  #   if(input$pca_diff_protein_phos){
  #     pca_protein_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
  #                                                            proteome_data = db_execution_phos$normalized_data,
  #                                                            type = "protein")$plot
  #     db_execution_phos$protein_differential_PCA = pca_protein_diff_fig
  #     pca_protein_diff_fig
  #   } else{
  #     db_execution_phos$protein_differential_PCA = NULL
  #   }
  # })
  # 
  generate_pca_peptide_diff_phos <- reactive({
    req(input$pca_diff_peptide_phos)
    if(input$pca_diff_peptide_phos){
      pca_peptide_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
                                                             proteome_data = db_execution_phos$normalized_data,
                                                             type = "peptide")$plot
      db_execution_phos$peptide_differential_PCA = pca_peptide_diff_fig
      pca_peptide_diff_fig
    } else{
      db_execution_phos$peptide_differential_PCA = NULL
    }
  })
  
  
  # generate_protein_diff_barplot_phos <- reactive(function(size_text){
  #   req(input$protein_diff_barplot_phos)
  #   if(input$protein_diff_barplot_phos){
  #     ploft_diff_number <- generate_differential_barplots(db_execution_phos$differential_results,
  #                                                         data_type="protein", size_text=size_text)$plot
  #     db_execution_phos$protein_differential_barplot = ploft_diff_number
  #     ploft_diff_number
  #   } else{
  #     db_execution_phos$protein_differential_barplot = NULL
  #   }
  # })
  # 
  generate_peptide_diff_barplot_phos <- reactive(function(size_text, zoom){
    req(input$peptide_diff_barplot_phos)
    if(input$peptide_diff_barplot_phos){
      ploft_diff_number_pep <- generate_differential_barplots(db_execution_phos$differential_results,
                                                              data_type="peptide", size_text=size_text)$plot
      if(!zoom){
        db_execution_phos$peptide_differential_barplot = ploft_diff_number_pep
      }
      ploft_diff_number_pep
    } else{
      db_execution_phos$peptide_differential_barplot = NULL
    }
  })
  
  generate_peptide_upset_phos <- reactive({
    req(input$peptide_upset_phos)
    if(input$peptide_upset_phos){
      ploft_diff_number_pep <- generate_upset_plot(db_execution_phos$differential_results,
                                                   type="peptide", 
                                                   DE_class = "all")$plot
      db_execution_phos$peptide_upset_plot = ploft_diff_number_pep
      ploft_diff_number_pep
    } else{
      db_execution_phos$peptide_upset_plot = NULL
    }
  })
  
  # PHOSPROTN: Execution pipeline ----
  observeEvent(input$report_proteome_phos, {
    
    output$protn_results_ui_phos <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
              # Reset other analysis
              db_execution_phos$parameter <- list()
              updateCheckboxInput(session, "differential_analysis_checkbox_phos", value = FALSE)
              
              message(session$token)
              message(tempdir())
              #Creation directory for the results
              dirOutput_2 <- tempdir()
              currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
              dirOutput_1 <- paste("/", currentTime, "/", sep = "")
              dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
              dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
              message(dirOutput_Server)
              db_execution_phos$dirOutput <- dirOutput_Server
              #Save folder for the download
              readr::write_csv(data.frame("session"=session$token,
                                          "outdir"=dirOutput_Server),
                               file = paste0(tempdir(),"/outdir_log_PhosProTN.log"), append = T)
              
              
              #Read parameter and execution
              software <- input$sw_analyzer_phos
              file_input_proteome = input$input_file_proteome_phos$name
              file_prot_proteome = if(software=="PD"){input$prot_file_proteome_phos$name}else{NA}
              file_psm_proteome = if(software=="PD"){input$psm_file_proteome_phos$name}else{NA}
              file_pep_proteome = input$pep_file_proteome_phos$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_phosprotn"), showWarnings = FALSE)
              dir_input <- paste(dirOutput_Server, "input_phosprotn", sep = "")
              file.copy(from = input$input_file_proteome_phos$datapath, to = paste0(dir_input,'/ANNOTATION_',file_input_proteome)) 
              if(software=="PD"){file.copy(from = input$prot_file_proteome_phos$datapath, to =paste0(dir_input,'/PROT_',file_prot_proteome))} 
              if(software=="PD"){file.copy(from = input$psm_file_proteome_phos$datapath, to =paste0(dir_input,'/PSM_',file_psm_proteome))} 
              file.copy(from = input$pep_file_proteome_phos$datapath, to = paste0(dir_input,'/PEP_',file_pep_proteome)) 
              
              # If to batch corrected read column
              if(input$batch_correction_phos){
                batch_corr <- TRUE
                batch_correction_col <- input$batch_correction_col_phos
              } else{
                batch_corr <- FALSE
                batch_correction_col <- "batch"
              }
              
              # If advance filter
              if(input$advance_filter_phos){
                NA_allow_condition <- input$NA_allow_condition_phos
                min_peptide_protein <- input$min_peptide_protein_phos
                impute_algorithm <- unlist(tstrsplit(input$impute_algorithm_phos, "_"))
                if(input$sample_column_phos == "Sample"){
                  sample_column <- input$sample_column_phos
                } else{
                  if(software=="PD"){
                    sample_column <- "File Name"
                  } else{
                    sample_column <- "Sample"
                  }
                }
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
                impute_algorithm <- c("norm","phosr")
                if(software=="PD"){
                  sample_column <- "File Name"
                } else{
                  sample_column <- "Sample"
                }
              }
              
              message(software)
              progress=0
              msg_read_function <-NULL
              withCallingHandlers(
                {
                  shinyjs::html("text", "")
                  if(software == "PD"){
                    db_execution_phos$proteome_data <-read_phosphoproteomics(software = "PD",
                                                                        folder = dir_input,
                                                                        peptide_filename = "PEP_",
                                                                        annotation_filename = "ANNOTATION_",
                                                                        proteinGroup_filename = "PROT_", 
                                                                        psm_filename = "PSM_",
                                                                        sample_col = sample_column,
                                                                        batch_corr_exe = batch_corr, 
                                                                        batch_col = batch_correction_col,
                                                                        phospho_thr = input$phos_thr/100, 
                                                                        filt_absent_value = NA_allow_condition, 
                                                                        min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ"){
                    db_execution_phos$proteome_data <- read_phosphoproteomics(software = "MQ",
                                                                         folder = dir_input,
                                                                         peptide_filename = "PEP_",
                                                                         annotation_filename = "ANNOTATION_", 
                                                                         sample_col = sample_column,
                                                                         batch_corr_exe = batch_corr, 
                                                                         batch_col = batch_correction_col,
                                                                         phospho_thr = input$phos_thr/100, 
                                                                         filt_absent_value = NA_allow_condition, 
                                                                         min_peptide_protein = min_peptide_protein)
                  }
                },
                message = function(m) {
                  msg_read_function <<- append(msg_read_function, conditionMessage(m))
                  # shinyjs::html(id = "messagge_read_phos_protn", html = paste0("<p>",m$message,"</p>"), add = TRUE)
                  progress=progress+0.05
                  setProgress(value = progress)
                }
              )
              
              write_lines(msg_read_function, file = paste0(db_execution_phos$dirOutput,"log_filter_read_function.txt"))
              
              db_execution_phos$data_loaded <- TRUE
              
              if(impute_algorithm[1] != "norm"){
                message("Doing before imputation")
                message(impute_algorithm[1])
                db_execution_phos$imputed_data <- impute_intensity(proteome_data = db_execution_phos$proteome_data, type = impute_algorithm[1])
                db_execution_phos$normalized_data <- normalization_ProTN(proteome_data = db_execution_phos$imputed_data)
              } else{
                message("Doing before normalization")
                message(impute_algorithm[2])
                db_execution_phos$normalized_data <- normalization_ProTN(proteome_data = db_execution_phos$proteome_data)
                db_execution_phos$normalized_data <- impute_intensity(proteome_data = db_execution_phos$normalized_data, type = impute_algorithm[2])
              }
              
              if(batch_corr){
                message("Executing batch correction...")
                db_execution_phos$normalized_data <- batch_correction(proteome_data = db_execution_phos$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              db_execution_phos$parameter<-list("Imputation and normalization algorithm: " = ifelse(impute_algorithm[1] != "norm", impute_algorithm[1], impute_algorithm[2]), 
                                           "Sample column in annotation file: " = sample_column, 
                                           "Batch correction: " = ifelse(batch_corr, batch_correction_col, "FALSE"), 
                                           "N° missing value allow per condition: " = NA_allow_condition, 
                                           "Minimum peptide per protein: " = min_peptide_protein,
                                           "Phospho site threshold: " = input$phos_thr)
              
              output$c_anno_phos <- DT::renderDT(db_execution_phos$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome_phos", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
                tags$h3("Statistics:"),
                tags$h4(paste0("Number of phospho-proteins: ", uniqueN(db_execution_phos$normalized_data$dat_gene$GeneName))),
                tags$h4(paste0("Number of phospho-peptides: ", uniqueN(db_execution_phos$normalized_data$dat_pep$ID_peptide))),
                tags$h3("Annotation table"),
                DT::DTOutput("c_anno_phos")
              )
            })
          },
          error = function(e) {
            #Create error report and reactivate the click in the page
            showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
            html_text<-str_replace(read_file("R/error.html"), 
                                   pattern = "The page you’re looking for doesn’t exist.</p>", 
                                   replacement = paste0("Description:", e, "</p>"))
            write_file(html_text, file = paste0(tempdir(), "/error.html"))
            tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
          }
        )
      })
      
    })
    
    output$render_phospho_percentage_plot_phos <- renderUI({
      if (input$phospho_percentage_plot_phos) {
        tagList(
          tags$h3("Percentage of phosphosite residue"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('phospho_percentage_plot_phos')",
            plotOutput("small_phospho_percentage_plot_phos")
          )
        )
      } else{
        db_execution_phos$phospho_percentage = NULL
      }
    })
    output$small_phospho_percentage_plot_phos <- renderPlot({
      generate_phospho_percentage_plot_phos()(4,zoom=F)
    })
    
    output$render_abundance_plot_phos <- renderUI({
      if (input$abundance_plot_phos) {
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('abundance_plot_phos')",
            plotOutput("small_abundance_plot_phos")
          )
        )
      } else{
        db_execution_phos$generate_abundance = NULL
      }
    })
    output$small_abundance_plot_phos <- renderPlot({
      generate_abundance_phos()
    })
    
    output$render_raw_violin_phos <- renderUI({
      if (input$raw_violin_phos) {
        tagList(
          tags$h3("Distribution raw abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('raw_violin_plot_phos')",
            plotOutput("small_raw_violin_phos")
          )
        )
      } else{
        db_execution_phos$raw_abundance_distribution = NULL
      }
    })
    output$small_raw_violin_phos <- renderPlot({
      generate_raw_violin_phos()
    })
    
    output$render_complexity_plot_phos <- renderUI({
      if (input$complexity_plot_phos) {
        tagList(
          tags$h3("Complexity plot of raw abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('complexity_plot_phos')",
            plotOutput("small_complexity_plot_phos")
          )
        )
      } else{
        db_execution_phos$generate_complexity = NULL
      }
    })
    output$small_complexity_plot_phos <- renderPlot({
      generate_complexity_phos()
    })
    
    output$render_peptide_distribution_phos <- renderUI({
      if (input$peptide_distribution_phos) {
        tagList(
          tags$h3("N° peptides per proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('peptide_distribution_plot_phos')",
            plotOutput("small_peptide_distribution_phos")
          )
        )
      } else{
        db_execution_phos$generate_peptide_distribution = NULL
      }
    })
    output$small_peptide_distribution_phos <- renderPlot({
      generate_peptide_distribution_phos()
    })
    
    output$render_peptide_violin_phos <- renderUI({
      if (input$peptide_violin_phos) {
        tagList(
          tags$h3("Distribution phospho-peptide abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('peptide_violin_plot_phos')",
            plotOutput("small_peptide_violin_phos")
          )
        )
      } else{
        db_execution_phos$peptide_abundance_distirbution = NULL
      }
    })
    output$small_peptide_violin_phos <- renderPlot({
      generate_peptide_violin_phos()
    })
    
    output$render_mds_peptide_phos <- renderUI({
      if (input$mds_peptide_phos) {
        tagList(
          tags$h3("MDS based on phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('mds_peptide_phos')",
            plotOutput("small_mds_peptide_phos")
          )
        )
      } else{
        db_execution_phos$peptide_MDS = NULL
      }
    })
    output$small_mds_peptide_phos <- renderPlot({
      generate_mds_peptide_phos()
    })

    output$render_pca_peptide_phos <- renderUI({
      if (input$pca_peptide_phos) {
        tagList(
          tags$h3("PCA based on phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('pca_peptide_phos')",
            plotOutput("small_pca_peptide_phos")
          )
        )
      } else{
        db_execution_phos$peptide_PCA = NULL
      }
    })
    output$small_pca_peptide_phos <- renderPlot({
      generate_pca_peptide_phos()
    })
    
    output$render_protein_boxplot_phos <- renderUI({
      if (input$boxplot_protein_phos) {
        req(input$list_proteins_phos)
        tagList(
          tags$h3("Boxplot selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('protein_boxplot_phos')",
            plotOutput("small_protein_boxplot_phos")
          )
        )
      } else{
        db_execution_phos$protein_boxplot = NULL
      }
    })
    output$small_protein_boxplot_phos <- renderPlot({
      generate_protein_boxplot_phos()
    })
    
    output$render_protein_heatmap_phos <- renderUI({
      if (input$heatmap_protein_phos) {
        req(input$list_proteins_phos)
        tagList(
          tags$h3("Heatmap selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('protein_heatmap_phos')",
            plotOutput("small_protein_heatmap_phos")
          )
        )
      } else{
        db_execution_phos$protein_heatmap = NULL
      }
    })
    output$small_protein_heatmap_phos <- renderPlot({
      generate_protein_heatmap_phos()
    })
  })
  
  ## PHOSPROTN: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn_phos, {
    output$render_differential_analysis_phos <- renderUI({
      isolate({
        updateCheckboxInput(session, "enrichment_analysis_phos", value = FALSE)
        updateCheckboxInput(session, "stringdb_analysis_phos", value = FALSE)
        updateCheckboxInput(session, "kinase_tree_analysis_phos", value = FALSE)
        
        db_execution_phos$dt_formule_contrast <- as.data.table(hot_to_r(input$render_formule_contrast_table_phos))
        db_execution_phos$dt_formule_contrast <- db_execution_phos$dt_formule_contrast[Formule!=""]
        print(db_execution_phos$dt_formule_contrast)
        formule_diff <- as.list(db_execution_phos$dt_formule_contrast$Formule)
        names(formule_diff) <- stri_replace_all(db_execution_phos$dt_formule_contrast$Name, replacement = "_", regex = "-")
        
        names(formule_diff) <- lapply(1:length(formule_diff), function(x){
          if(names(formule_diff)[x] == ""){
            stri_replace_all(formule_diff[[x]], replacement = "_VS_", regex = "-")
          } else{
            names(formule_diff)[x]
          }
        })
        db_execution_phos$formule_contrast <- formule_diff
        message(db_execution_phos$formule_contrast)
        
        withProgress(message = "Differential analysis in process, please wait!", {
          message(session$token)
          message(tempdir())
          
          db_execution_phos$differential_results <- differential_analysis(proteome_data = db_execution_phos$normalized_data,
                                                                     formule_contrast = db_execution_phos$formule_contrast,
                                                                     fc_thr=as.double(input$FC_thr_phos),
                                                                     pval_fdr = input$pval_fdr_phos,
                                                                     pval_thr=as.double(input$pval_thr_phos),
                                                                     signal_thr=0)
          db_execution_phos$formule_contrast <- db_execution_phos$formule_contrast[unique(union(db_execution_phos$differential_results$protein_results_long$comp, 
                                                                                                db_execution_phos$differential_results$peptide_results_long$comp))]
          db_execution_phos$parameter <- c(db_execution_phos$parameter,
                                           "Fold change threshold for significance: "=input$FC_thr_phos,
                                           "P.value type used: "=input$pval_fdr_phos,
                                           "P.value threshold for significance: "=input$pval_thr_phos)
        })
        
        tags$h2("Differential Analysis")
      })
    })
    
    output$render_peptide_diff_table_phos <- renderUI({
      if(input$peptide_diff_table_phos){
        output$peptide_results_long_phos <- DT::renderDT(db_execution_phos$differential_results$peptide_results_long)
        DT::DTOutput("peptide_results_long_phos")
      }
    })
    
    output$render_peptide_diff_barplot_phos <- renderUI({
      if (input$peptide_diff_barplot_phos) {
        tagList(
          tags$h3("N° differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('peptide_diff_barplot_phos')",
            plotOutput("small_peptide_diff_barplot_phos")
          )
        )
      } else{
        db_execution_phos$peptide_differential_barplot = NULL
      }
    })
    output$small_peptide_diff_barplot_phos <- renderPlot({
      generate_peptide_diff_barplot_phos()(3,zoom=F)
    })
    
    output$render_peptide_upset_phos <- renderUI({
      if (input$peptide_upset_phos) {
        tagList(
          tags$h3("Differential peptides upset plot"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('peptide_upset_phos')",
            plotOutput("small_peptide_upset_phos")
          )
        )
      } else{
        db_execution_phos$peptide_upset_plot = NULL
      }
    })
    output$small_peptide_upset_phos <- renderPlot({
      generate_peptide_upset_phos()
    })
    
    output$render_peptide_ma_plot_phos <- renderUI({
      if (input$peptide_ma_plot_phos) {
        c_anno <- db_execution_phos$proteome_data$c_anno
        generate_ma_plots_peptide <- list()
        for(comp in names(db_execution_phos$formule_contrast)){
          message(comp)
          design <- model.matrix(~0 + c_anno$condition)
          colnames(design) <- levels(as.factor(c_anno$condition))
          rownames(design) <- c_anno$sample
          
          conds <- as.data.table(makeContrasts(contrasts = db_execution_phos$formule_contrast[[comp]], levels = design), keep.rownames = T)
          conds <- conds[as.vector(conds[,2]!=0), rn]
          message(conds)
          
          generate_ma_plots_peptide[[comp]] <- ma_plot(differential_results = db_execution_phos$differential_results, 
                                                       proteome_data = db_execution_phos$normalized_data,
                                                       type="peptide", comparison = comp, condition = conds)$plot
        }
        db_execution_phos$peptide_ma_plot = lapply(generate_ma_plots_peptide, function(x){ggplotly(x, tooltip = c("text"))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_ma_plots_peptide)) {
          plot_id <- paste0(names(generate_ma_plots_peptide)[i], "_ma_pep_phos")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(ggplotly(generate_ma_plots_peptide[[names(generate_ma_plots_peptide)[my_i]]], tooltip = "text"))
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_ma_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("MA Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_ma_peptide_phos"), tabs))
        )
      } else{
        db_execution_phos$peptide_ma_plot = NULL
      }
    })
    
    output$render_peptide_vulcano_phos <- renderUI({
      if(input$peptide_vulcano_phos){
        generate_volcano_plots_peptide <- list()
        for(comp in names(db_execution_phos$formule_contrast)){
          generate_volcano_plots_peptide<-c(generate_volcano_plots_peptide,
                                            generate_volcano_plots(db_execution_phos$differential_results,
                                                                   data_type="peptide",
                                                                   comparison=comp,
                                                                   fc_thr=as.double(input$FC_thr_phos),
                                                                   pval_fdr = input$pval_fdr_phos,
                                                                   pval_thr=as.double(input$pval_thr_phos)))
        }
        db_execution_phos$peptide_vulcano = generate_volcano_plots_peptide
        # Generate tabPanels in a for loop
        tabs_pep_vulcano_phos <- list()
        for (i_phos in seq_along(generate_volcano_plots_peptide)) {
          plot_id_phos <- paste0(names(generate_volcano_plots_peptide)[i_phos], "_pep_phos")
          # Create an output slot for each plot
          local({
            my_i_phos <- i_phos
            my_plot_id_phos <- plot_id_phos
            output[[my_plot_id_phos]] <- renderPlotly(generate_volcano_plots_peptide[[names(generate_volcano_plots_peptide)[my_i_phos]]])
          })
          
          tabs_pep_vulcano_phos[[i_phos]] <- tabPanel(
            title = paste(names(generate_volcano_plots_peptide)[i_phos]),
            plotlyOutput(plot_id_phos, width = "99%")
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_peptide_phos"), tabs_pep_vulcano_phos))
        )
      } else{
        db_execution_phos$peptide_vulcano = NULL
      }
    })
    
    output$render_mds_peptide_diff_phos <- renderUI({
      if (input$mds_diff_peptide_phos) {
        tagList(
          tags$h3("MDS based on differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('mds_peptide_diff_phos')",
            plotOutput("small_mds_peptide_diff_phos")
          )
        )
      } else{
        db_execution_phos$peptide_differential_MDS = NULL
      }
    })
    output$small_mds_peptide_diff_phos <- renderPlot({
      generate_mds_peptide_diff_phos()
    })
    
    output$render_pca_peptide_diff_phos <- renderUI({
      if (input$pca_diff_peptide_phos) {
        tagList(
          tags$h3("PCA based on differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos('pca_peptide_diff_phos')",
            plotOutput("small_pca_peptide_diff_phos")
          )
        )
      } else{
        db_execution_phos$peptide_differential_PCA = NULL
      }
    })
    output$small_pca_peptide_diff_phos <- renderPlot({
      generate_pca_peptide_diff_phos()
    })
    
  })
  
  ## PHOSPROTN: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn_phos, {
    output$render_enrichement_analysis_phos <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution_phos$enrichment_results <- perform_enrichment_analysis(differential_results = db_execution_phos$differential_results,
                                                                        enrichR_custom_DB = T,
                                                                        enrich_filter_DBs=input$DB_enrichment_phos,    
                                                                        overlap_size_enrich_thr=as.double(input$os_enrich_phos),
                                                                        pval_fdr_enrich = input$pval_fdr_phos,
                                                                        pval_enrich_thr=as.double(input$pval_thr_phos),
                                                                        dirOutput=db_execution_phos$dirOutput, 
                                                                        with_background = input$enrich_with_background_phos)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich_phos), regex=","))
        
        db_execution_phos$parameter <- c(db_execution_phos$parameter,
                                    "Enrichment databases selected: "=paste(input$DB_enrichment_phos, collapse = ", "),
                                    "P.value type used for enrichment: "=input$pval_fdr_phos,
                                    "P.value threshold for enrichment significance: "=input$pval_thr_phos,
                                    "Overlap size threshold for enrichment significance: "=input$os_enrich_phos,
                                    "Enrichment filter terms: "=if(length(terms_enrich)>0){paste(terms_enrich, collapse = ", ")}else{"None"},
                                    "Enrichment with background: "=input$enrich_with_background_phos)
        
        
        plots_down <- enrichment_figure(enr_df = db_execution_phos$enrichment_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution_phos$enrichment_results$anno_class), input$DB_enrichment_phos))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(plots_down)) {
          plot_id <- names(plots_down)[i]
          height_id <- max(min(20, length(unique(plots_down[[names(plots_down)[i]]]$data$y_col))*0.4),3)*96
          message(paste0("Height for ",names(plots_down)[i], ": ", height_id))
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlot({
              plots_down[[names(plots_down)[my_i]]]
            }, height = height_id)
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(plots_down)[i]),
            plotOutput(plot_id, height = height_id)
          )
        }
        
        tagList(
          tags$h2("Enrichment Analysis"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_enrichment_phos"), tabs))
        )
        
      })
    })
  })
  ## PHOSPROTN: stringdb analysis ----
  observeEvent(input$execute_stringdb_analysis_btn_phos, {
    output$render_stringdb_phos <- renderUI({
      isolate({
        withProgress(message = "STRINGdb analysis in process, please wait!", {
          
          db_execution_phos$stringdb_res <- STRINGdb_network(differential_results = db_execution_phos$differential_results,
                                                        species=input$taxonomy_phos, 
                                                        dirOutput=db_execution_phos$dirOutput, 
                                                        score_thr=input$score_thr_stringdb_phos,
                                                        shiny = T)
          db_execution_phos$parameter <- c(db_execution_phos$parameter,
                                      "STRINGdb taxonomy: "=input$taxonomy_phos,
                                      "STRINGdb score threshold: "=input$score_thr_stringdb_phos)
          
          tagList(
            tags$h2("STRINGdb analysis"),
            fluidRow(
              selectInput("stringdb_show_phos", label = "Select StringDB to show: (click on STRING logo to open the results on stringDB website)", 
                          choices = names(db_execution_phos$stringdb_res), width = "15%"),
              actionButton("stringdb_selected_phos", "Select!", width = "10%")  
            ),
            tags$div(id = "stringEmbedded")
          )
        })
      })
    })
  })
  
  observeEvent(input$stringdb_selected_phos, {
    js$loadStringData(input$taxonomy_phos, db_execution_phos$stringdb_res[[input$stringdb_show_phos]], input$score_thr_stringdb_phos)
  })
  
  ## PHOSPROTN: kinase tree analysis ----
  observeEvent(input$execute_kinase_tree_analysis_btn_phos, {
    output$render_kinase_tree_phos <- renderUI({
      isolate({
        withProgress(message = "Kinase Tree analysis in process, please wait! (Can take several minutes)", {
          
          db_execution_phos$kinase_tree_res <- kinase_tree(proteome_data = db_execution_phos$normalized_data, 
                                                      differential_results = db_execution_phos$differential_results, 
                                                      formule_CORAL = db_execution_phos$formule_contrast, 
                                                      dirOutput=db_execution_phos$dirOutput, 
                                                      phosR_thr = input$score_thr_phosr_phos, 
                                                      species = input$taxonomy_kinase_phos)
          db_execution_phos$parameter <- c(db_execution_phos$parameter,
                                      "Kinase tree taxonomy: "=input$taxonomy_kinase_phos,
                                      "Score thr for PhosR: "=input$score_thr_phosr_phos)
          
          if(input$taxonomy_kinase_phos == "Homo sapiens"){
            tagList(
              tags$h2("Kinase Tree analysis"),
              fluidRow(
                selectInput("kinase_tree_show", label = "Select Kinase Tree to show:", 
                            choices = names(db_execution_phos$kinase_tree_res), width = "15%"),
                actionButton("kinase_tree_selected", "Select!", width = "10%")  
              ),
              imageOutput("render_kin_tree", height = "auto")
            )
          } else{
            tagList(
              tags$h2("Kinase Tree analysis"),
              tags$h4("For Mouse the graphical representation of the kinome tree is not done. The results can be downloaded."),
              tags$hr()
            )
          }
        })
      })
    })
  })
  
  observeEvent(input$kinase_tree_selected, {
    output$render_kin_tree <- renderImage({
      isolate({
        list(src = paste0(db_execution_phos$dirOutput, "pics/kinaseTree/",input$kinase_tree_show,"_kinase_Tree_CORAL.svg"),
             alt = "Kinase Tree"
        )
      })
    }, deleteFile = FALSE)
  })
  
  # PHOSPROTN: download results ----
  output$download_proteome_phos <- downloadHandler(
    filename = "results_PhosProTN.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            message(db_execution_phos$dirOutput)
            setProgress(value = 0.01)
            
            # Generate report
            params <- list(
              doc_title = input$title_exp_phos,
              description = input$description_exp_phos,
              readPD_files = if (input$sw_analyzer_phos == "PD") {TRUE} else {FALSE},
              readMQ_files = if (input$sw_analyzer_phos == "MQ") {TRUE} else {FALSE},
              impute_algorithm = if(input$advance_filter_phos){input$impute_algorithm_phos} else {"norm_phosr"},
              db_execution = reactiveValuesToList(db_execution_phos),
              file_input = paste(db_execution_phos$dirOutput, "input_phosprotn", sep = ""),
              batch_corr_exe = if(input$batch_correction_phos){input$batch_correction_col_phos}else{NULL},
              prot_boxplot = if(input$boxplot_protein_phos | input$heatmap_protein_phos){input$list_proteins_phos}else{NULL},
              fc_thr = if(is.null(input$FC_thr_phos)){"0.75"}else{input$FC_thr_phos},
              pval_fdr = input$pval_fdr_phos,
              pval_thr = if(is.null(input$pval_thr_phos)){"0.05"}else{input$pval_thr_phos},
              pval_fdr_enrich = input$pval_fdr_enrich_phos,
              pval_enrich_thr = if(is.null(input$pvalue_enrich_phos)){"0.05"}else{input$pvalue_enrich_phos},
              overlap_size_enrich_thr = if(is.null(input$os_enrich_phos)){as.integer(5)}else{input$os_enrich_phos},
              enrich_filter_term = input$terms_enrich_phos,
              enrich_filter_DBs = input$DB_enrichment_phos,
              taxonomy=input$taxonomy_phos, 
              score_thr=input$score_thr_stringdb_phos,
              dirOutput = db_execution_phos$dirOutput
            )
            
            # Render in background the report
            p = callr::r_bg(
              func = function(db_execution_phos, params, dirOutput, env) {
                rmarkdown::render("R/phosprotn_report.Rmd",
                                  output_file = "phosprotn_report.html",
                                  output_dir = dirOutput,
                                  params = params,
                                  envir = env
                )
              },
              args = list(db_execution_phos, params, db_execution_phos$dirOutput, new.env(parent = globalenv())),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            # Saving RData in background
            db_results_phos = reactiveValuesToList(db_execution_phos)
            db_results_phos <- db_results_phos[!(unlist(lapply(db_results_phos, is.null)))]
            p_rdata = callr::r_bg(
              func = function(db_results_phos, dirOutput) {
                save(db_results_phos, file = paste0(dirOutput,"db_results_PhosProTN.RData"))
              },
              args = list(db_results_phos, db_results_phos$dirOutput),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            
            # Prepare file for the download
            if(length(db_execution_phos$normalized_data)>0){
              save_abundance_tables(proteome_data = db_execution_phos$normalized_data, 
                                    dirOutput = db_execution_phos$dirOutput)
            }
            setProgress(value = 0.1)
            
            if(length(db_execution_phos$differential_results)>0){
              save_differential_analysis_table(proteome_data = db_execution_phos$normalized_data,
                                               differential_results = db_execution_phos$differential_results,
                                               dirOutput=db_execution_phos$dirOutput, phospho_analysis = TRUE)
            }
            setProgress(value = 0.2)
            
            if(input$phospho_percentage_plot_phos & !is.null(db_execution_phos$phospho_percentage)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/phospho_percentage.pdf"), 
                     plot = db_execution_phos$phospho_percentage, 
                     create.dir = T, width = 7, height = 3)
            } else if("phospho_percentage.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/phospho_percentage.pdf"))
            }
            
            if(input$abundance_plot_phos & !is.null(db_execution_phos$generate_abundance)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/missing_available_abundance.pdf"), 
                     plot = db_execution_phos$generate_abundance, 
                     create.dir = T, width = 7, height = 5)
            } else if("missing_available_abundance.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/missing_available_abundance.pdf"))
            }
            setProgress(value = 0.25)
            
            if(input$raw_violin_phos & !is.null(db_execution_phos$raw_abundance_distribution)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/raw_abundance_distribution.pdf"), 
                     plot = db_execution_phos$raw_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("raw_abundance_distribution.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/raw_abundance_distribution.pdf"))
            }
            
            
            if(input$peptide_distribution_phos & !is.null(db_execution_phos$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_per_protein.pdf"), 
                     plot = db_execution_phos$generate_peptide_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_per_protein.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_per_protein.pdf"))
            }
            setProgress(value = 0.30)
            
            if(input$complexity_plot_phos & !is.null(db_execution_phos$generate_complexity)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/complexity_plot.pdf"), 
                     plot = db_execution_phos$generate_complexity, 
                     create.dir = T, width = 10, height = 8)
            } else if("complexity_plot.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/complexity_plot.pdf"))
            }
            setProgress(value = 0.33)
            
            setProgress(value = 0.35)
            
            if(input$peptide_violin_phos & !is.null(db_execution_phos$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution_phos$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_abundance_distribution.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_abundance_distribution.pdf"))
            }
            setProgress(value = 0.40)
            
            setProgress(value = 0.43)
            
            if(input$mds_peptide_phos & !is.null(db_execution_phos$peptide_MDS)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution_phos$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_MDS.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_MDS.pdf"))
            }
            setProgress(value = 0.45)
            
            setProgress(value = 0.47)
            
            if(input$pca_peptide_phos & !is.null(db_execution_phos$peptide_PCA)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution_phos$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_PCA.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_PCA.pdf"))
            }
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein_phos & !is.null(db_execution_phos$protein_boxplot)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution_phos$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_boxplot.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/protein_boxplot.pdf"))
            }
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein_phos & !is.null(db_execution_phos$protein_heatmap)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution_phos$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_heatmap.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/protein_heatmap.pdf"))
            }
            setProgress(value = 0.55)
            
            setProgress(value = 0.58)
            
            if(!is.null(db_execution_phos$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution_phos$peptide_differential_barplot, 
                     create.dir = T, width = 17, height = 6)
            } else if("peptide_differential_barplot.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_differential_barplot.pdf"))
            }
            setProgress(value = 0.60)
            
            if(!is.null(db_execution_phos$peptide_upset_plot)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_upset_plot.pdf"), 
                     plot = db_execution_phos$peptide_upset_plot, 
                     create.dir = T, width = 12, height = 6)
            } else if("peptide_upset_plot.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_upset_plot.pdf"))
            }
            setProgress(value = 0.63)
            
            
            if(!is.null(db_execution_phos$peptide_ma_plot)){
              dir.create(file.path(paste0(db_execution_phos$dirOutput,"pics/"), "peptide_ma_plot"), showWarnings = FALSE)
              for(comp in names(db_execution_phos$peptide_ma_plot)){
                
                htmlwidgets::saveWidget(db_execution_phos$peptide_ma_plot[[comp]], 
                                        file = paste0(db_execution_phos$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"))
                webshot2::webshot(url = paste0(db_execution_phos$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"), 
                                  file = paste0(db_execution_phos$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_phos$dirOutput,"pics/peptide_ma_plot"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_phos$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution_phos$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos$peptide_vulcano)){
                # plotly::save_image(db_execution_phos$peptide_vulcano[[comp]], 
                #                    file = paste0(str_replace_all(db_execution_phos$dirOutput, pattern="\\\\", replacement="/"),"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos$peptide_vulcano[[comp]], 
                                        file = paste0(db_execution_phos$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution_phos$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"), 
                                  file = paste0(db_execution_phos$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_phos$dirOutput,"pics/peptide_vulcano"))
            }
            setProgress(value = 0.68)
            
            if(!is.null(db_execution_phos$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution_phos$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_MDS.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_differential_MDS.pdf"))
            }
            setProgress(value = 0.70)
            
            if(!is.null(db_execution_phos$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution_phos$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_PCA.pdf" %in% list.files(paste0(db_execution_phos$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos$dirOutput,"pics/peptide_differential_PCA.pdf"))
            }
            setProgress(value = 0.75)
            
            if(length(db_execution_phos$enrichment_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich_phos), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution_phos$enrichment_results,
                                              category = c("down","up"), 
                                              enrich_filter_term = terms_enrich,
                                              save=T, 
                                              dirOutput = db_execution_phos$dirOutput)
            } 
            setProgress(value = 0.82)
            
            if(length(db_execution_phos$stringdb_res)>0){
              tmp_res <- STRINGdb_network(differential_results = db_execution_phos$differential_results,
                                          species=input$taxonomy_phos, 
                                          dirOutput=db_execution_phos$dirOutput,
                                          score_thr=input$score_thr_stringdb_phos,
                                          shiny = F)
              
            } 
            setProgress(value = 0.95)
            
            # Write tsv file with parameter
            params <- data.table("Parameter" = names(db_execution_phos$parameter),
                                 "Value" = unlist(db_execution_phos$parameter))
            fwrite(params, paste0(db_execution_phos$dirOutput,"parameters_used.txt"), sep = "\t", col.names = F)
            
            #Get results Report
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p$read_output()
            p$wait(30000)
            for (i in 1:15) {
              p$read_output()
              p$wait(1000*60)  
            }
            if(p$is_alive() | is.null(p$get_result())){
              p$kill()
              print("\n ERROR: An error occur during the report rendering. \n ")
            } else{
              report<-p$get_result()
              p$kill()
              message("Render report DONE.")
            }
            
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p_rdata$read_output()
            p_rdata$wait(30000)
            for (i in 1:15) {
              p_rdata$read_output()
              p_rdata$wait(1000*60)
            }
            if(p_rdata$is_alive() | is.null(p_rdata$get_result())){
              p_rdata$kill()
              print("\n ERROR: An error occur during the RData saving. \n ")
            } else{
              report<-p_rdata$get_result()
              p_rdata$kill()
              message("RData saved.")
            }
            
            # Save RData db_execution_phos
            # setProgress(value = 0.95, message = "Saving RData...")
            # db_results_PhosProTN = reactiveValuesToList(db_execution_phos)
            # db_results_PhosProTN <- db_results_PhosProTN[!(unlist(lapply(db_results_PhosProTN, is.null)))]
            # save(db_results_PhosProTN, file = paste0(db_results_PhosProTN$dirOutput,"db_results_PhosProTN.RData"))
            
            #Save folder for the download
            oldwd <- getwd()
            message(db_execution_phos$dirOutput)
            setwd(db_execution_phos$dirOutput)
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  ## PHOSPROTN: full screen trigger ----
  
  # ReactiveVal for currently selected plot to fullscreen
  selected_plot_phos <- reactiveVal(NULL)
  
  # Update selected_plot when JS sends fullscreen_trigger id
  observeEvent(input$fullscreen_trigger_phos, {
    selected_plot_phos(input$fullscreen_trigger_phos)
  })
  
  # Render fullscreen plot dynamically based on selected_plot()
  output$fullscreen_plot_phos <- renderPlot({
    req(selected_plot_phos())
    switch(selected_plot_phos(),
           "phospho_percentage_plot_phos" = generate_phospho_percentage_plot_phos()(size_text=8,zoom=T) + ggtitle("Percentage of phosphosite residue")+theme(text=element_text(size=25), axis.text.y = element_text(size = 25)),
           "abundance_plot_phos" = generate_abundance_phos() + ggtitle("Percentage missing values respect detected abundance")+theme(text=element_text(size=25)),
           "raw_violin_plot_phos" = generate_raw_violin_phos() + ggtitle("Raw abundance distribution")+theme(text=element_text(size=25)),
           "peptide_distribution_plot_phos" = generate_peptide_distribution_phos() + ggtitle("N° peptides per proteins")+theme(text=element_text(size=25)),
           "complexity_plot_phos" = generate_complexity_phos() + ggtitle("Complexity plot of raw abundance")+theme(text=element_text(size=25)),
           "peptide_violin_plot_phos" = generate_peptide_violin_phos() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "mds_peptide_phos" = generate_mds_peptide_phos() + ggtitle("MDS based on peptides")+theme(text=element_text(size=25)),
           "pca_peptide_phos" = generate_pca_peptide_phos() + ggtitle("PCA based on peptides")+theme(text=element_text(size=25)),
           "protein_boxplot_phos" = generate_protein_boxplot_phos() + ggtitle("Boxplot selected proteins")+theme(text=element_text(size=25)),
           "protein_heatmap_phos" = generate_protein_heatmap_phos() + ggtitle("Heatmap selected proteins")+theme(text=element_text(size=25)),
           "peptide_diff_barplot_phos" = generate_peptide_diff_barplot_phos()(8,zoom=T) + ggtitle("N° differential peptides")+theme(text=element_text(size=25)),
           "peptide_upset_phos" = generate_peptide_upset_phos(),
           "mds_peptide_diff_phos" = generate_mds_peptide_diff_phos() + ggtitle("MDS based on differential peptides")+theme(text=element_text(size=25)),
           "pca_peptide_diff_phos" = generate_pca_peptide_diff_phos() + ggtitle("PCA based on differential peptides")+theme(text=element_text(size=25)),
           # default fallback:
           NULL
    )
  })
  
  ##############################################################################
  ### PhosProTN_with_prot ----
  # Optional visibility based on the selection ----
  
  ## PhosProTN_with_prot: Visibility of the proteomics files for PhosProTN_with_prot ----
  output$input_proteome_phos_protn <- renderUI({
    if (input$sw_analyzer_phos_protn == "PD"){
      tagList(
        fluidRow(
          fileInput("input_file_phospho_phos_protn", "Select the SAMPLE_ANNOTATION file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_phospho_phos_protn", "Select the PEP file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_phospho_phos_protn", "Select the PROT file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("psm_file_phospho_phos_protn", "Select the PSM file of the PHOSPHO-PROTEOMICS..."),
        ),
        tags$br(),
        fluidRow(
          fileInput("input_file_proteome_phos_protn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_phos_protn", "Select the PEP file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_proteome_phos_protn", "Select the PROT file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer_phos_protn == "MQ"){
      tagList(
        fluidRow(
          fileInput("input_file_phospho_phos_protn", "Select the SAMPLE_ANNOTATION file of the PHOSPHO-PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_phospho_phos_protn", "Select the EVIDENCE file of the PHOSPHO-PROTEOMICS..."),
        ),
        tags$br(),
        fluidRow(
          fileInput("input_file_proteome_phos_protn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_phos_protn", "Select the EVIDENCE file of the PROTEOMICS..."),
        )
      )
    } else{
      tagList(
        tags$p("BACK")
      )
    }
  })
  
  ## PhosProTN_with_prot: textbox for batch correction----
  output$batch_correction_ui_phos_protn <- renderUI({ 
    if(input$batch_correction_phos_protn){
      textInput("batch_correction_col_phos_protn", "Column in Annotation file with the batch:")
    } 
  })
  
  ## PhosProTN_with_prot: advance filters----
  output$advance_filter_ui_phos_protn <- renderUI({ 
    if(input$advance_filter_phos_protn){
      tagList(
        numericInput("NA_allow_condition_phos_protn", "N° missing value allow per condition", value = 0, min = 0, max = 5),
        numericInput("min_peptide_protein_phos_protn", "Minimum peptide per protein", value = 1, min = 1),
        selectizeInput("impute_algorithm_phos_protn", "Select impute algorithm:",
                       choices = list("PhosR and normalization" = "phosr_norm", 
                                      "Gaussian estimation and normalization" = "gaussian_norm",
                                      "missForest and normalization" = "missForest_norm",
                                      "Pre-normalization and PhosR" = "norm_phosr", 
                                      "Pre-normalization and Gaussian estimation" = "norm_gaussian",
                                      "Pre-normalization and missForest" = "norm_missForest",
                                      "Pre-normalization and pcaMethods" = "norm_pcaMethods"),
                       selected = "norm_phosr", multiple = FALSE
        ),
        textInput("sample_column_phos_protn", "Column name with the sample name:", value = "Sample")
      )
    } 
  })
  ## PhosProTN_with_prot: textbox for list proteins ----
  output$list_protein_ui_phos_protn <- renderUI({ 
    if(input$boxplot_protein_phos_protn | input$heatmap_protein_phos_protn){
      textInput("list_proteins_phos_protn", "List proteins to show (separate by: \",\"):")
    } 
  })
  
  ## PhosProTN_with_prot: show parameter for differential analysis ----
  output$differential_params_ui_phos_protn <- renderUI({ 
    if(input$differential_analysis_checkbox_phos_protn){
      tagList(
        tags$label("Write in each line a different comparison"),
        tags$label("(right click to add row)"),
        rHandsontableOutput('render_formule_contrast_table_phos_protn'),
        # textAreaInput("formule_contrast", "Write in each line a different comparison", rows = 4),
        textInput("FC_thr_phos_protn", "Fold change threshold for significance:",value = 0.5),
        radioButtons("pval_fdr_phos_protn", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_val"),
        textInput("pval_thr_phos_protn", "P.value threshold for significance:", value = 0.05),
        actionButton("execute_differential_analysis_btn_phos_protn", "Run!"),
        checkboxInput("peptide_diff_table_phos_protn", "Phospho-peptides differentiated table", FALSE),
        checkboxInput("peptide_diff_barplot_phos_protn", "Phospho-peptides differentiated barplot", TRUE),
        checkboxInput("peptide_upset_phos_protn", "Peptides upset plot", FALSE),
        checkboxInput("peptide_ma_plot_phos_protn", "Peptides MA plot", FALSE),
        checkboxInput("peptide_vulcano_phos_protn", "Phospho-peptides vulcano plot", FALSE),
        checkboxInput("mds_diff_peptide_phos_protn", "MDS based on diffential phospho-peptide", FALSE),
        checkboxInput("pca_diff_peptide_phos_protn", "PCA based on diffential phospho-peptide", FALSE),
        tags$h3("Enrichment Analysis:"),
        checkboxInput("enrichment_analysis_phos_protn", "Execute enrichment analysis", FALSE),
        uiOutput("enrichment_params_ui_phos_protn"),
        tags$h3("STRINGdb network:"),
        checkboxInput("stringdb_analysis_phos_protn", "Execute STRINGdb", FALSE),
        uiOutput("stringdb_params_ui_phos_protn"),
        tags$h3("Kinase tree:"),
        checkboxInput("kinase_tree_analysis_phos_protn", "Execute PhosR kinase tree", FALSE),
        uiOutput("kinase_tree_params_ui_phos_protn")
      )
    } else{
      # Reset UI elements
      updateCheckboxInput(session, "enrichment_analysis_phos_protn", value = FALSE)
      updateCheckboxInput(session, "stringdb_analysis_phos_protn", value = FALSE)
      updateCheckboxInput(session, "kinase_tree_analysis_phos_protn", value = FALSE)
      
      db_execution_phos_protn$formule_contrast <- list()
      db_execution_phos_protn$dt_formule_contrast <- data.table("Name"=c("","","",""),"Formule"=c("","","",""))
      db_execution_phos_protn$differential_results <- list()
      
      updateCheckboxInput(session, "protein_diff_barplot_phos_protn", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_barplot_phos_protn",  value = FALSE)
      updateCheckboxInput(session, "protein_diff_table_phos_protn", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_table_phos_protn", value = FALSE)
      updateCheckboxInput(session, "protein_upset_phos_protn", value = FALSE)
      updateCheckboxInput(session, "peptide_upset_phos_protn", value = FALSE)
      updateCheckboxInput(session, "protein_vulcano_phos_protn",  value =  FALSE)
      updateCheckboxInput(session, "peptide_vulcano_phos_protn", value = FALSE)
      updateCheckboxInput(session, "protein_ma_plot_phos_protn", value = FALSE)
      updateCheckboxInput(session, "peptide_ma_plot_phos_protn", value = FALSE)
      updateCheckboxInput(session, "mds_diff_protein_phos_protn", value = FALSE)
      updateCheckboxInput(session, "mds_diff_peptide_phos_protn", value = FALSE)
      updateCheckboxInput(session, "pca_diff_protein_phos_protn", value = FALSE)
      updateCheckboxInput(session, "pca_diff_peptide_phos_protn", value = FALSE)
      
      db_execution_phos_protn$protein_differential_barplot <- NULL
      db_execution_phos_protn$peptide_differential_barplot <- NULL
      db_execution_phos_protn$protein_upset_plot <- NULL
      db_execution_phos_protn$peptide_upset_plot <- NULL
      db_execution_phos_protn$protein_ma_plot <- NULL
      db_execution_phos_protn$peptide_ma_plot <- NULL
      db_execution_phos_protn$protein_vulcano <- NULL
      db_execution_phos_protn$peptide_vulcano <- NULL
      db_execution_phos_protn$protein_differential_MDS <- NULL
      db_execution_phos_protn$peptide_differential_MDS <- NULL
      db_execution_phos_protn$protein_differential_PCA <- NULL
      db_execution_phos_protn$peptide_differential_PCA <- NULL
      
      output$render_differential_analysis_phos_protn <- renderUI({NULL})
      output$render_protein_diff_table_phos_protn <- renderUI({NULL})
      output$render_peptide_diff_table_phos_protn <- renderUI({NULL})
      output$render_protein_diff_barplot_phos_protn <- renderUI({NULL})
      output$render_peptide_diff_barplot_phos_protn <- renderUI({NULL})
      output$render_protein_upset_phos_protn <- renderUI({NULL})
      output$render_peptide_upset_phos_protn <- renderUI({NULL})
      output$render_protein_ma_plot_phos_protn <- renderUI({NULL})
      output$render_peptide_ma_plot_phos_protn <- renderUI({NULL})
      output$render_protein_vulcano_phos_protn <- renderUI({NULL})
      output$render_peptide_vulcano_phos_protn <- renderUI({NULL})
      output$render_mds_protein_diff_phos_protn <- renderUI({NULL})
      output$render_mds_peptide_diff_phos_protn <- renderUI({NULL})
      output$render_pca_protein_diff_phos_protn <- renderUI({NULL})
      output$render_pca_peptide_diff_phos_protn <- renderUI({NULL})
    }
  })
  
  output$render_formule_contrast_table_phos_protn <- renderRHandsontable({
    rhandsontable(db_execution_phos_protn$dt_formule_contrast, rowHeaders = NULL, stretchH = "all")
  })
  
  ## PhosProTN_with_prot: show enrichment parameter ----
  output$enrichment_params_ui_phos_protn <- renderUI({ 
    if(input$enrichment_analysis_phos_protn){
      tagList(
        # radioButtons("enrichR_universe", "Execute enrichment of the whole Universe", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
        selectizeInput("DB_enrichment_phos_protn", "DB to analyse:",
                       choices = lapply(split(read_tsv("data/dbs_enrichR.txt", col_names = FALSE)$X1,
                                              read_tsv("data/dbs_enrichR.txt", col_names = FALSE)[,2]), as.list),
                       selected = NULL, multiple = TRUE
        ),
        textInput("terms_enrich_phos_protn", "Terms to search (separated by \",\"):"),
        radioButtons("pval_fdr_enrich_phos_protn", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_adj"),
        textInput("pvalue_enrich_phos_protn", "P.value threshold for significance:", value = 0.05),
        sliderInput("os_enrich_phos_protn", "Overlap size thr for enrichment", 1, 30, step = 1, value = 5),
        checkboxInput("enrich_with_background_phos_protn", "Enrichment with background", FALSE),
        actionButton("execute_enrichment_analysis_btn_phos_protn", "Run!")
      )
    } else{
      db_execution_phos_protn$enrichment_results <- list()
      output$render_enrichement_analysis_phos_protn <- renderUI({NULL})
    }
  })
  
  ## PhosProTN_with_prot: show stringdb parameter ----
  output$stringdb_params_ui_phos_protn <- renderUI({
    if(input$stringdb_analysis_phos_protn){
      tagList(
        selectizeInput("taxonomy_phos_protn", "NCBI Taxonomy ID", 
                       choice = data.table::fread("data/subset_tax.csv", select = "name"), 
                       selected = "Homo sapiens", multiple = F),
        sliderInput("score_thr_stringdb_phos_protn", "Score thr for STRINGdb", 500, 1000, step = 10, value = 700),
        actionButton("execute_stringdb_analysis_btn_phos_protn", "Run!"),
        tags$br()
      )
    } else{
      db_execution_phos_protn$stringdb_res <- list()
      output$render_stringdb_phos_protn <- renderUI({NULL})
    }
  })
  ## PhosProTN_with_prot: show kinase_tree parameter ----
  output$kinase_tree_params_ui_phos_protn <- renderUI({
    if(input$kinase_tree_analysis_phos_protn){
      tax_sel <- if(is.null(input$taxonomy_phos_protn)){
        character(0)
      } else if(input$taxonomy_phos_protn == "Homo sapiens"){
        "Homo sapiens"
      } else if(input$taxonomy_phos_protn == "Mus musculus"){
        "Mus musculus"
      } else{
        NULL
      }
      
      if(is.null(tax_sel)){
        shinyalert::shinyalert("Kinase analysis", 
                               "Kinase analysis can be performed only for Homo Sapiens (Human) or Mus Musculus (Mouse)", 
                               type = "info")
        updateCheckboxInput(session, "kinase_tree_analysis_phos_protn", value = FALSE)
      } else{
        tagList(
          radioButtons("taxonomy_kinase_phos_protn", "Select species (CORAL tree will be print only for Homo sapiens)", 
                       choiceNames = c("Homo Sapiens", "Mus Musculus"),
                       choiceValues = c("Homo sapiens","Mus musculus"), inline = TRUE, 
                       selected = tax_sel),
          sliderInput("score_thr_phosr_phos_protn", "Score thr for PhosR", 0, 1, step = 0.05, value = 0.7),
          actionButton("execute_kinase_tree_analysis_btn_phos_protn", "Run!"),
          tags$br()
        )
      } else{
        db_execution_phos_protn$kinase_tree_res <- list()
        output$render_kinase_tree_phos_protn <- renderUI({NULL})
      }
    }
  })
  
  ## PhosProTN_with_prot: function genereting plot ----
  generate_phospho_percentage_plot_phos_protn <- reactive(function(size_text, zoom) {
    req(input$phospho_percentage_plot_phos_protn)
    if(input$phospho_percentage_plot_phos_protn){
      generate_abundance_fig <- create_phosphosite_plot(proteome_data = db_execution_phos_protn$proteome_data, 
                                                        software = input$sw_analyzer_phos_protn, 
                                                        size_text = size_text)$plot
      if(!zoom){
        db_execution_phos_protn$phospho_percentage = generate_abundance_fig
      }
      generate_abundance_fig
    } else{
      db_execution_phos_protn$phospho_percentage = NULL
    }
  })
  
  generate_abundance_phos_protn_prot <- reactive({
    req(input$abundance_plot_phos_protn)
    if(input$abundance_plot_phos_protn){
      generate_abundance_fig <- generate_abundance_plot(proteome_data = db_execution_phos_protn$proteome_data)
      db_execution_phos_protn$generate_abundance = generate_abundance_fig
      generate_abundance_fig$proteome_plot
    } else{
      db_execution_phos_protn$generate_abundance = NULL
    }
  })
  generate_abundance_phos_protn_phos <- reactive({
    req(input$abundance_plot_phos_protn)
    if(input$abundance_plot_phos_protn){
      generate_abundance_fig <- generate_abundance_plot(proteome_data = db_execution_phos_protn$proteome_data)
      db_execution_phos_protn$generate_abundance = generate_abundance_fig
      generate_abundance_fig$phospho_plot
    } else{
      db_execution_phos_protn$generate_abundance = NULL
    }
  })
  
  generate_peptide_distribution_phos_protn_prot <- reactive({
    req(input$peptide_distribution_phos_protn)
    if(input$peptide_distribution_phos_protn){
      peptide_distribution_fig <- generate_peptide_distribution_plot(proteome_data = db_execution_phos_protn$proteome_data)
      db_execution_phos_protn$generate_peptide_distribution = peptide_distribution_fig
      peptide_distribution_fig$proteome_plot
    } else{
      db_execution_phos_protn$generate_peptide_distribution = NULL
    }
  })
  generate_peptide_distribution_phos_protn_phos <- reactive({
    req(input$peptide_distribution_phos_protn)
    if(input$peptide_distribution_phos_protn){
      peptide_distribution_fig <- generate_peptide_distribution_plot(proteome_data = db_execution_phos_protn$proteome_data)
      db_execution_phos_protn$generate_peptide_distribution = peptide_distribution_fig
      peptide_distribution_fig$phospho_plot
    } else{
      db_execution_phos_protn$generate_peptide_distribution = NULL
    }
  })
  
  generate_raw_violin_prot_phos_protn <- reactive({
    req(input$raw_violin_prot_phos_protn)
    if(input$raw_violin_prot_phos_protn){
      raw_abundance_distribution_fig <- plot_raw_abundance_distribution(proteome_data = db_execution_phos_protn$proteome_data,
                                                                        type = "protein")$plot
      db_execution_phos_protn$raw_proteome_abundance_distribution = raw_abundance_distribution_fig
      raw_abundance_distribution_fig
    } else{
      db_execution_phos_protn$raw_proteome_abundance_distribution = NULL
    }
  })
  generate_raw_violin_pep_phos_protn <- reactive({
    req(input$raw_violin_pep_phos_protn)
    if(input$raw_violin_pep_phos_protn){
      raw_abundance_distribution_fig <- plot_raw_abundance_distribution(proteome_data = db_execution_phos_protn$proteome_data,
                                                                        type = "peptide")$plot
      db_execution_phos_protn$raw_phospho_abundance_distribution = raw_abundance_distribution_fig
      raw_abundance_distribution_fig
    } else{
      db_execution_phos_protn$raw_phospho_abundance_distribution = NULL
    }
  })
  
  generate_complexity_phos_protn_prot <- reactive({
    req(input$complexity_plot_phos_protn)
    if(input$complexity_plot_phos_protn){
      generate_complexity_fig <- complexity_plot(proteome_data = db_execution_phos_protn$proteome_data, phospho_with_proteome = T)
      db_execution_phos_protn$generate_complexity = generate_complexity_fig
      generate_complexity_fig$proteome_plot
    } else{
      db_execution_phos_protn$generate_complexity = NULL
    }
  })
  generate_complexity_phos_protn_phos <- reactive({
    req(input$complexity_plot_phos_protn)
    if(input$complexity_plot_phos_protn){
      generate_complexity_fig <- complexity_plot(proteome_data = db_execution_phos_protn$proteome_data, phospho_with_proteome = T)
      db_execution_phos_protn$generate_complexity = generate_complexity_fig
      generate_complexity_fig$phospho_plot
    } else{
      db_execution_phos_protn$generate_complexity = NULL
    }
  })
  
  generate_protein_violin_phos_protn <- reactive({
    req(input$protein_violin_phos_protn)
    if(input$protein_violin_phos_protn){
      protein_abundance_distribution_fig <- plot_abundance_distribution(proteome_data = db_execution_phos_protn$normalized_data,
                                                                        type = "protein")$plot
      db_execution_phos_protn$protein_abundance_distribution = protein_abundance_distribution_fig
      protein_abundance_distribution_fig
    } else{
      db_execution_phos_protn$protein_abundance_distribution = NULL
    }
  })
  
  generate_peptide_violin_phos_protn <- reactive({
    req(input$peptide_violin_phos_protn)
    if(input$peptide_violin_phos_protn){
      peptide_abundance_distirbution_fig <- plot_abundance_distribution(proteome_data = db_execution_phos_protn$normalized_data,
                                                                        type = "peptide")$plot
      db_execution_phos_protn$peptide_abundance_distirbution = peptide_abundance_distirbution_fig
      peptide_abundance_distirbution_fig
    } else{
      db_execution_phos_protn$peptide_abundance_distirbution = NULL
    }
  })
  
  generate_mds_protein_phos_protn <- reactive({
    req(input$mds_protein_phos_protn)
    if(input$mds_protein_phos_protn){
      mds_protein_fig <- mds_plot(proteome_data = db_execution_phos_protn$normalized_data,
                                  type = "protein")$plot
      db_execution_phos_protn$protein_MDS = mds_protein_fig
      mds_protein_fig
    } else{
      db_execution_phos_protn$protein_MDS = NULL
    }
  })
  
  generate_mds_peptide_phos_protn <- reactive({
    req(input$mds_peptide_phos_protn)
    if(input$mds_peptide_phos_protn){
      mds_peptide_fig <- mds_plot(proteome_data = db_execution_phos_protn$normalized_data,
                                  type = "peptide")$plot
      db_execution_phos_protn$peptide_MDS = mds_peptide_fig
      mds_peptide_fig
    } else{
      db_execution_phos_protn$peptide_MDS = NULL
    }
  })
  
  generate_pca_protein_phos_protn <- reactive({
    req(input$pca_protein_phos_protn)
    if(input$pca_protein_phos_protn){
      pca_protein_fig <- pca_plot(proteome_data = db_execution_phos_protn$normalized_data,
                                  type = "protein")$plot
      db_execution_phos_protn$protein_PCA = pca_protein_fig
      pca_protein_fig
    } else{
      db_execution_phos_protn$protein_PCA = NULL
    }
  })
  
  generate_pca_peptide_phos_protn <- reactive({
    req(input$pca_peptide_phos_protn)
    if(input$pca_peptide_phos_protn){
      pca_peptide_fig <- pca_plot(proteome_data = db_execution_phos_protn$normalized_data,
                                  type = "peptide")$plot
      db_execution_phos_protn$peptide_PCA = pca_peptide_fig
      pca_peptide_fig
    } else{
      db_execution_phos_protn$peptide_PCA = NULL
    }
  })
  
  generate_protein_boxplot_phos_protn <- reactive({
    req(input$boxplot_protein_phos_protn)
    if(input$boxplot_protein_phos_protn){
      req(input$list_proteins_phos_protn)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos_protn), regex=",")
      db_execution_phos_protn$parameter <- c(db_execution_phos_protn$parameter, "List proteins boxplot abundance: "=input$list_proteins_phos_protn)
      boxplot_protein_fig <- plot_selected_proteins(proteome_data = db_execution_phos_protn$normalized_data,
                                                    list_protein = unlist(list_proteins))$plot
      db_execution_phos_protn$protein_boxplot = boxplot_protein_fig
      boxplot_protein_fig
    } else{
      db_execution_phos_protn$protein_boxplot = NULL
    }
  })
  
  generate_protein_heatmap_phos_protn <- reactive({
    req(input$heatmap_protein_phos_protn)
    if(input$heatmap_protein_phos_protn){
      req(input$list_proteins_phos_protn)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos_protn), regex=",")
      db_execution_phos_protn$parameter <- c(db_execution_phos_protn$parameter, "List proteins heatmap abundance: "=input$list_proteins_phos_protn)
      heatmap_protein_fig <- heatmap_selected_proteins(proteome_data = db_execution_phos_protn$normalized_data, list_protein = unlist(list_proteins))$plot
      db_execution_phos_protn$protein_heatmap = heatmap_protein_fig
      heatmap_protein_fig
    } else{
      db_execution_phos_protn$protein_heatmap = NULL
    }
  })
  
  generate_mds_peptide_diff_phos_protn <- reactive({
    req(input$mds_diff_peptide_phos_protn)
    if(input$mds_diff_peptide_phos_protn){
      mds_peptide_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution_phos_protn$differential_results,
                                                             proteome_data = db_execution_phos_protn$normalized_data,
                                                             type = "peptide")$plot
      db_execution_phos_protn$peptide_differential_MDS = mds_peptide_diff_fig
      mds_peptide_diff_fig
    } else{
      db_execution_phos_protn$peptide_differential_MDS = NULL
    }
  })
  
  generate_pca_peptide_diff_phos_protn <- reactive({
    req(input$pca_diff_peptide_phos_protn)
    if(input$pca_diff_peptide_phos_protn){
      pca_peptide_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution_phos_protn$differential_results,
                                                             proteome_data = db_execution_phos_protn$normalized_data,
                                                             type = "peptide")$plot
      db_execution_phos_protn$peptide_differential_PCA = pca_peptide_diff_fig
      pca_peptide_diff_fig
    } else{
      db_execution_phos_protn$peptide_differential_PCA = NULL
    }
  })
  
  generate_peptide_diff_barplot_phos_protn <- reactive(function(size_text, zoom){
    req(input$peptide_diff_barplot_phos_protn)
    if(input$peptide_diff_barplot_phos_protn){
      ploft_diff_number_pep <- generate_differential_barplots(db_execution_phos_protn$differential_results,
                                                              data_type="peptide", size_text = size_text)$plot
      if(!zoom){
        db_execution_phos_protn$peptide_differential_barplot = ploft_diff_number_pep
      }
      ploft_diff_number_pep
    } else{
      db_execution_phos_protn$peptide_differential_barplot = NULL
    }
  })
  
  generate_peptide_upset_phos_protn <- reactive({
    req(input$peptide_upset_phos_protn)
    if(input$peptide_upset_phos_protn){
      ploft_diff_number_pep <- generate_upset_plot(db_execution_phos_protn$differential_results,
                                                   type="peptide", 
                                                   DE_class = "all")$plot
      db_execution_phos_protn$peptide_upset_plot = ploft_diff_number_pep
      ploft_diff_number_pep
    } else{
      db_execution_phos_protn$peptide_upset_plot = NULL
    }
  })
  
  # PhosProTN_with_prot: Execution pipeline ----
  observeEvent(input$report_proteome_phos_protn, {
    
    output$protn_results_ui_phos_protn <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
              # Reset other analysis
              db_execution_phos$parameter <- list()
              updateCheckboxInput(session, "differential_analysis_checkbox_phos_protn", value = FALSE)
              
              message(session$token)
              message(tempdir())
              #Creation directory for the results
              dirOutput_2 <- tempdir()
              currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
              dirOutput_1 <- paste("/", currentTime, "/", sep = "")
              dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
              dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
              message(dirOutput_Server)
              db_execution_phos_protn$dirOutput <- dirOutput_Server
              #Save folder for the download
              readr::write_csv(data.frame("session"=session$token,
                                          "outdir"=dirOutput_Server),
                               file = paste0(tempdir(),"/outdir_log_PhosProTN_proteome.log"), append = T)
              
              ### PHOSPHO ----
              #Read parameter and execution
              software <- input$sw_analyzer_phos_protn
              file_input_phospho = input$input_file_phospho_phos_protn$name
              file_prot_phospho = if(software=="PD"){input$prot_file_phospho_phos_protn$name}else{NA}
              file_psm_phospho = if(software=="PD"){input$psm_file_phospho_phos_protn$name}else{NA}
              file_pep_phospho = input$pep_file_phospho_phos_protn$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_phospho"), showWarnings = FALSE)
              dir_input_phospho <- paste(dirOutput_Server, "input_phospho", sep = "")
              file.copy(from = input$input_file_phospho_phos_protn$datapath, to = paste0(dir_input_phospho,'/ANNOTATION_',file_input_phospho)) 
              if(software=="PD"){file.copy(from = input$prot_file_phospho_phos_protn$datapath, to =paste0(dir_input_phospho,'/PROT_',file_prot_phospho))} 
              if(software=="PD"){file.copy(from = input$psm_file_phospho_phos_protn$datapath, to =paste0(dir_input_phospho,'/PSM_',file_psm_phospho))} 
              file.copy(from = input$pep_file_phospho_phos_protn$datapath, to = paste0(dir_input_phospho,'/PEP_',file_pep_phospho)) 
              
              ### PROTEOME ----
              #Read parameter and execution
              file_input_proteome = input$input_file_proteome_phos_protn$name
              file_prot_proteome = if(software=="PD"){input$prot_file_proteome_phos_protn$name}else{NA}
              file_pep_proteome = input$pep_file_proteome_phos_protn$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_proteome"), showWarnings = FALSE)
              dir_input_proteome <- paste(dirOutput_Server, "input_proteome", sep = "")
              file.copy(from = input$input_file_proteome_phos_protn$datapath, to = paste0(dir_input_proteome,'/ANNOTATION_',file_input_proteome)) 
              if(software=="PD"){file.copy(from = input$prot_file_proteome_phos_protn$datapath, to =paste0(dir_input_proteome,'/PROT_',file_prot_proteome))} 
              file.copy(from = input$pep_file_proteome_phos_protn$datapath, to = paste0(dir_input_proteome,'/PEP_',file_pep_proteome)) 
              
              ### ----
              # If advance filter
              if(input$advance_filter_phos_protn){
                NA_allow_condition <- input$NA_allow_condition_phos_protn
                min_peptide_protein <- input$min_peptide_protein_phos_protn
                impute_algorithm <- unlist(tstrsplit(input$impute_algorithm_phos_protn, "_"))

                if(input$sample_column_phos_protn == "Sample"){
                  sample_column <- input$sample_column_phos_protn
                } else{
                  if(software=="PD"){
                    sample_column <- "File Name"
                  } else{
                    sample_column <- "Sample"
                  }
                }
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
                impute_algorithm <- c("norm","phosr")
                if(software=="PD"){
                  sample_column <- "File Name"
                } else{
                  sample_column <- "Sample"
                }
              }
              message(input$impute_algorithm_phos_protn)
              message(impute_algorithm)
              
              # If to batch corrected read column
              if(input$batch_correction_phos_protn){
                batch_corr <- TRUE
                batch_correction_col <- input$batch_correction_col_phos_protn
              } else{
                batch_corr <- FALSE
                batch_correction_col <- "batch"
              }
              
              message(software)
              progress=0
              msg_read_function <- NULL
              withCallingHandlers(
                {
                  shinyjs::html("text", "")
                  if(software == "PD"){
                    db_execution_phos_protn$proteome_data <- read_phospho_proteome_proteomics(software = "PD", 
                                                                                   folder_proteome = dir_input_proteome,
                                                                                   folder_phospho = dir_input_phospho,
                                                                                   peptide_proteome_filename = "PEP_", 
                                                                                   peptide_phospho_filename = "PEP_", 
                                                                                   annotation_proteome_filename = "ANNOTATION_",
                                                                                   proteinGroup_proteome_filename = "PROT_", 
                                                                                   annotation_phospho_filename = "ANNOTATION_",
                                                                                   proteinGroup_phospho_filename = "PROT_", 
                                                                                   psm_phospho_filename = "PSM_", 
                                                                                   sample_proteome_col = sample_column, 
                                                                                   sample_phospho_col = sample_column,
                                                                                   batch_corr_exe = batch_corr, 
                                                                                   batch_col = batch_correction_col,
                                                                                   phospho_thr = input$phos_thr_phos_protn/100, 
                                                                                   filt_absent_value = NA_allow_condition, 
                                                                                   min_peptide_protein = min_peptide_protein)

                  } else if(software == "MQ"){
                    db_execution_phos_protn$proteome_data <- read_phospho_proteome_proteomics(software = "MQ", 
                                                                                   folder_proteome = dir_input_proteome,
                                                                                   folder_phospho = dir_input_phospho,
                                                                                   peptide_proteome_filename = "PEP_", 
                                                                                   peptide_phospho_filename = "PEP_", 
                                                                                   annotation_proteome_filename = "ANNOTATION_",
                                                                                   annotation_phospho_filename = "ANNOTATION_",
                                                                                   sample_proteome_col = sample_column, 
                                                                                   sample_phospho_col = sample_column,
                                                                                   batch_corr_exe = batch_corr, 
                                                                                   batch_col = batch_correction_col,
                                                                                   phospho_thr = input$phos_thr_phos_protn/100, 
                                                                                   filt_absent_value = NA_allow_condition, 
                                                                                   min_peptide_protein = min_peptide_protein)
                  }
                },
                message = function(m) {
                  msg_read_function <<- append(msg_read_function, conditionMessage(m))
                  # shinyjs::html(id = "messagge_read_phos_protn", html = paste0("<p>",m$message,"</p>"), add = TRUE)
                  progress=progress+0.05
                  setProgress(value = progress)
                }
              )
              
              write_lines(msg_read_function, file = paste0(db_execution_phos_protn$dirOutput,"log_filter_read_function.txt"))
              
              db_execution_phos_protn$data_loaded <- TRUE
              
              
              if(impute_algorithm[1] != "norm"){
                message("Doing before imputation")
                message(impute_algorithm[1])
                db_execution_phos_protn$imputed_data <- impute_intensity(proteome_data = db_execution_phos_protn$proteome_data, type = impute_algorithm[1])
                db_execution_phos_protn$normalized_data <- normalization_ProTN(proteome_data = db_execution_phos_protn$imputed_data)
              } else{
                message("Doing before normalization")
                message(impute_algorithm[2])
                db_execution_phos_protn$normalized_data <- normalization_ProTN(proteome_data = db_execution_phos_protn$proteome_data)
                db_execution_phos_protn$normalized_data <- impute_intensity(proteome_data = db_execution_phos_protn$normalized_data, type = impute_algorithm[2])
              }
              
              if(batch_corr){
                message("Executing batch correction...")
                db_execution_phos_protn$normalized_data <- batch_correction(proteome_data = db_execution_phos_protn$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              db_execution_phos_protn$parameter<-list("Imputation and normalization algorithm: " = ifelse(impute_algorithm[1] != "norm", impute_algorithm[1], impute_algorithm[2]), 
                                                "Sample column in annotation file: " = sample_column, 
                                                "Batch correction: " = ifelse(batch_corr, batch_correction_col, "FALSE"), 
                                                "N° missing value allow per condition: " = NA_allow_condition, 
                                                "Minimum peptide per protein: " = min_peptide_protein,
                                                "Phospho site threshold: " = input$phos_thr_phos_protn)
              
              
              output$c_anno_phos_protn <- DT::renderDT(db_execution_phos_protn$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome_phos_protn", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
                tags$h3("Statistics:"),
                tags$h4(paste0("Number of proteins from proteome dataset: ", uniqueN(db_execution_phos_protn$normalized_data$dat_gene$GeneName))),
                tags$h4(paste0("Number of phospho-peptides from phospho-proteome dataset: ", uniqueN(db_execution_phos_protn$normalized_data$dat_pep$ID_peptide))),
                tags$h3("Annotation table"),
                DT::DTOutput("c_anno_phos_protn")
              )
            })
          },
          error = function(e) {
            #Create error report and reactivate the click in the page
            showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
            html_text<-str_replace(read_file("R/error.html"), 
                                   pattern = "The page you’re looking for doesn’t exist.</p>", 
                                   replacement = paste0("Description:", e, "</p>"))
            write_file(html_text, file = paste0(tempdir(), "/error.html"))
            tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
          }
        )
      })
      
    })
    
    output$render_phospho_percentage_plot_phos_protn <- renderUI({
      if (input$phospho_percentage_plot_phos_protn) {
        tagList(
          tags$h3("Percentage of phosphosite residue"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('phospho_percentage_plot_phos_protn')",
            plotOutput("small_phospho_percentage_plot_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$phospho_percentage = NULL
      }
    })
    output$small_phospho_percentage_plot_phos_protn <- renderPlot({
      generate_phospho_percentage_plot_phos_protn()(4, zoom=F)
    })
    
    output$render_abundance_plot_phos_protn <- renderUI({
      if (input$abundance_plot_phos_protn) {
        tagList(
          fluidRow(
            column(
              width = 6,
              tags$h3("Percentage missing values respect detected abundance - Proteomics"),
              tags$div(
                style = "cursor:pointer;",
                onclick = "showFullscreenPlot_phos_protn('abundance_plot_phos_protn_prot')",
                plotOutput("small_abundance_plot_phos_protn_prot")
              )
            ),
            column(
              width = 6,
              tags$h3("Percentage missing values respect detected abundance - Phospho-proteomics"),
              tags$div(
                style = "cursor:pointer;",
                onclick = "showFullscreenPlot_phos_protn('abundance_plot_phos_protn_phos')",
                plotOutput("small_abundance_plot_phos_protn_phos")
              )
            )
          )
        )
      } else{
        db_execution_phos_protn$generate_abundance = NULL
      }
    })
    output$small_abundance_plot_phos_protn_prot <- renderPlot({
      generate_abundance_phos_protn_prot()
    })
    output$small_abundance_plot_phos_protn_phos <- renderPlot({
      generate_abundance_phos_protn_phos()
    })
    
    output$render_peptide_distribution_phos_protn <- renderUI({
      if (input$peptide_distribution_phos_protn) {
        tagList(
          
          fluidRow(
            column(
              width = 6,
              tags$h3("N° peptides per proteins - Proteomics"),
              tags$div(
                style = "cursor:pointer;",
                onclick = "showFullscreenPlot_phos_protn('peptide_distribution_plot_phos_protn_prot')",
                plotOutput("small_peptide_distribution_phos_protn_prot")
              )
            ),
            column(
              width = 6,
              tags$h3("N° peptides per proteins - Phospho-proteomics"),
              tags$div(
                style = "cursor:pointer;",
                onclick = "showFullscreenPlot_phos_protn('peptide_distribution_plot_phos_protn_phos')",
                plotOutput("small_peptide_distribution_phos_protn_phos")
              )
            )
          )
        )
      } else{
        db_execution_phos_protn$generate_peptide_distribution = NULL
      }
    })
    output$small_peptide_distribution_phos_protn_prot <- renderPlot({
      generate_peptide_distribution_phos_protn_prot()
    })
    output$small_peptide_distribution_phos_protn_phos <- renderPlot({
      generate_peptide_distribution_phos_protn_phos()
    })
    
    output$render_raw_violin_pep_phos_protn <- renderUI({
      if (input$raw_violin_pep_phos_protn) {
        tagList(
          tags$h3("Distribution raw abundance - Phospho-proteome"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('raw_violin_plot_pep_phos_protn')",
            plotOutput("small_raw_violin_pep_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$raw_phospho_abundance_distribution = NULL
      }
    })
    output$small_raw_violin_pep_phos_protn <- renderPlot({
      generate_raw_violin_pep_phos_protn()
    })
    output$render_raw_violin_prot_phos_protn <- renderUI({
      if (input$raw_violin_prot_phos_protn) {
        tagList(
          tags$h3("Distribution raw abundance - Proteome"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('raw_violin_plot_prot_phos_protn')",
            plotOutput("small_raw_violin_prot_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$raw_proteome_abundance_distribution = NULL
      }
    })
    output$small_raw_violin_prot_phos_protn <- renderPlot({
      generate_raw_violin_prot_phos_protn()
    })
    
    output$render_complexity_plot_phos_protn <- renderUI({
      if (input$complexity_plot_phos_protn) {
        tagList(
          fluidRow(
            column(
              width = 6,
              tags$h3("Complexity plot of raw abundance - Proteomics"),
              tags$div(
                style = "cursor:pointer;",
                onclick = "showFullscreenPlot_phos_protn('complexity_plot_phos_protn_prot')",
                plotOutput("small_complexity_plot_phos_protn_prot")
              )
            ),
            column(
              width = 6,
              tags$h3("Complexity plot of raw abundance - Phospho-proteomics"),
              tags$div(
                style = "cursor:pointer;",
                onclick = "showFullscreenPlot_phos_protn('complexity_plot_phos_protn_phos')",
                plotOutput("small_complexity_plot_phos_protn_phos")
              )
            )
          )
        )
      } else{
        db_execution_phos_protn$generate_complexity = NULL
      }
    })
    output$small_complexity_plot_phos_protn_prot <- renderPlot({
      generate_complexity_phos_protn_prot()
    })
    output$small_complexity_plot_phos_protn_phos <- renderPlot({
      generate_complexity_phos_protn_phos()
    })
    
    output$render_protein_violin_phos_protn <- renderUI({
      if (input$protein_violin_phos_protn) {
        tagList(
          tags$h3("Distribution protein abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('protein_violin_plot_phos_protn')",
            plotOutput("small_protein_violin_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$protein_abundance_distribution = NULL
      }
    })
    output$small_protein_violin_phos_protn <- renderPlot({
      generate_protein_violin_phos_protn()
    })
    
    output$render_peptide_violin_phos_protn <- renderUI({
      if (input$peptide_violin_phos_protn) {
        tagList(
          tags$h3("Distribution phospho-peptide abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('peptide_violin_plot_phos_protn')",
            plotOutput("small_peptide_violin_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_abundance_distirbution = NULL
      }
    })
    output$small_peptide_violin_phos_protn <- renderPlot({
      generate_peptide_violin_phos_protn()
    })
    
    output$render_mds_protein_phos_protn <- renderUI({
      if (input$mds_protein_phos_protn) {
        tagList(
          tags$h3("MDS based on proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('mds_protein_phos_protn')",
            plotOutput("small_mds_protein_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$protein_MDS = NULL
      }
    })
    output$small_mds_protein_phos_protn <- renderPlot({
      generate_mds_protein_phos_protn()
    })
    
    output$render_mds_peptide_phos_protn <- renderUI({
      if (input$mds_peptide_phos_protn) {
        tagList(
          tags$h3("MDS based on phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('mds_peptide_phos_protn')",
            plotOutput("small_mds_peptide_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_MDS = NULL
      }
    })
    output$small_mds_peptide_phos_protn <- renderPlot({
      generate_mds_peptide_phos_protn()
    })
    
    output$render_pca_protein_phos_protn <- renderUI({
      if (input$pca_protein_phos_protn) {
        tagList(
          tags$h3("PCA based on proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('pca_protein_phos_protn')",
            plotOutput("small_pca_protein_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$protein_PCA = NULL
      }
    })
    output$small_pca_protein_phos_protn <- renderPlot({
      generate_pca_protein_phos_protn()
    })
    
    output$render_pca_peptide_phos_protn <- renderUI({
      if (input$pca_peptide_phos_protn) {
        tagList(
          tags$h3("PCA based on phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('pca_peptide_phos_protn')",
            plotOutput("small_pca_peptide_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_PCA = NULL
      }
    })
    output$small_pca_peptide_phos_protn <- renderPlot({
      generate_pca_peptide_phos_protn()
    })
    
    output$render_protein_boxplot_phos_protn <- renderUI({
      if (input$boxplot_protein_phos_protn) {
        req(input$list_proteins_phos_protn)
        tagList(
          tags$h3("Boxplot selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('protein_boxplot_phos_protn')",
            plotOutput("small_protein_boxplot_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$protein_boxplot = NULL
      }
    })
    output$small_protein_boxplot_phos_protn <- renderPlot({
      generate_protein_boxplot_phos_protn()
    })
    
    output$render_protein_heatmap_phos_protn <- renderUI({
      if (input$heatmap_protein_phos_protn) {
        req(input$list_proteins_phos_protn)
        tagList(
          tags$h3("Heatmap selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('protein_heatmap_phos_protn')",
            plotOutput("small_protein_heatmap_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$protein_heatmap = NULL
      }
    })
    output$small_protein_heatmap_phos_protn <- renderPlot({
      generate_protein_heatmap_phos_protn()
    })
  })
  
  ## PhosProTN_with_prot: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn_phos_protn, {
    output$render_differential_analysis_phos_protn <- renderUI({
      isolate({
        db_execution_phos_protn$dt_formule_contrast <- as.data.table(hot_to_r(input$render_formule_contrast_table_phos_protn))
        db_execution_phos_protn$dt_formule_contrast <- db_execution_phos_protn$dt_formule_contrast[Formule!=""]
        print(db_execution_phos_protn$dt_formule_contrast)
        formule_diff <- as.list(db_execution_phos_protn$dt_formule_contrast$Formule)
        names(formule_diff) <- stri_replace_all(db_execution_phos_protn$dt_formule_contrast$Name, replacement = "_", regex = "-")
        
        names(formule_diff) <- lapply(1:length(formule_diff), function(x){
          if(names(formule_diff)[x] == ""){
            stri_replace_all(formule_diff[[x]], replacement = "_VS_", regex = "-")
          } else{
            names(formule_diff)[x]
          }
        })
        db_execution_phos_protn$formule_contrast <- formule_diff
        message(db_execution_phos_protn$formule_contrast)
        
        withProgress(message = "Differential analysis in process, please wait!", {
          message(session$token)
          message(tempdir())
          
          db_execution_phos_protn$differential_results <- differential_analysis(proteome_data = db_execution_phos_protn$normalized_data,
                                                                     formule_contrast = db_execution_phos_protn$formule_contrast,
                                                                     fc_thr=as.double(input$FC_thr_phos_protn),
                                                                     pval_fdr = input$pval_fdr_phos_protn,
                                                                     pval_thr=as.double(input$pval_thr_phos_protn),
                                                                     signal_thr=0)
          ll<-reactiveValuesToList(db_execution_phos_protn)
          db_execution_phos_protn$formule_contrast <- db_execution_phos_protn$formule_contrast[names(db_execution_phos_protn$formule_contrast) %in% unique((db_execution_phos_protn$differential_results$peptide_results_long$comp))]
          
          db_execution_phos_protn$parameter <- c(db_execution_phos_protn$parameter,
                                           "Fold change threshold for significance: "=input$FC_thr_phos_protn,
                                           "P.value type used: "=input$pval_fdr_phos_protn,
                                           "P.value threshold for significance: "=input$pval_thr_phos_protn)
        })
        
        tags$h2("Differential Analysis")
      })
    })
    
    output$render_peptide_diff_table_phos_protn <- renderUI({
      if(input$peptide_diff_table_phos_protn){
        output$peptide_results_long_phos_protn <- DT::renderDT(db_execution_phos_protn$differential_results$peptide_results_long)
        DT::DTOutput("peptide_results_long_phos_protn")
      }
    })
    
    output$render_peptide_diff_barplot_phos_protn <- renderUI({
      if (input$peptide_diff_barplot_phos_protn) {
        tagList(
          tags$h3("N° differential phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('peptide_diff_barplot_phos_protn')",
            plotOutput("small_peptide_diff_barplot_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_differential_barplot = NULL
      }
    })
    output$small_peptide_diff_barplot_phos_protn <- renderPlot({
      generate_peptide_diff_barplot_phos_protn()(3, zoom=F)
    })
    
    output$render_peptide_upset_phos_protn <- renderUI({
      if (input$peptide_upset_phos_protn) {
        tagList(
          tags$h3("Differential peptides upset plot"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('peptide_upset_phos_protn')",
            plotOutput("small_peptide_upset_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_upset_plot = NULL
      }
    })
    output$small_peptide_upset_phos_protn <- renderPlot({
      generate_peptide_upset_phos_protn()
    })
    
    output$render_peptide_ma_plot_phos_protn <- renderUI({
      if (input$peptide_ma_plot_phos_protn) {
        c_anno <- db_execution_phos_protn$proteome_data$c_anno_phospho
        generate_ma_plots_peptide <- list()
        for(comp in names(db_execution_phos_protn$formule_contrast)){
          message(comp)
          design <- model.matrix(~0 + c_anno$condition)
          colnames(design) <- levels(as.factor(c_anno$condition))
          rownames(design) <- c_anno$sample
          
          conds <- as.data.table(makeContrasts(contrasts = db_execution_phos_protn$formule_contrast[[comp]], levels = design), keep.rownames = T)
          conds <- conds[as.vector(conds[,2]!=0), rn]
          message(conds)
          
          generate_ma_plots_peptide[[comp]] <- ma_plot(differential_results = db_execution_phos_protn$differential_results, 
                                                       proteome_data = db_execution_phos_protn$normalized_data,
                                                       type="peptide", comparison = comp, condition = conds)$plot
        }
        db_execution_phos_protn$peptide_ma_plot = lapply(generate_ma_plots_peptide, function(x){ggplotly(x, tooltip = c("text"))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_ma_plots_peptide)) {
          plot_id <- paste0(names(generate_ma_plots_peptide)[i], "_ma_pep_phos_protn")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(ggplotly(generate_ma_plots_peptide[[names(generate_ma_plots_peptide)[my_i]]], tooltip = "text"))
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_ma_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("MA Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_ma_peptide_phos_protn"), tabs))
        )
      } else{
        db_execution_phos_protn$peptide_ma_plot = NULL
      }
    })
    
    output$render_peptide_vulcano_phos_protn <- renderUI({
      if(input$peptide_vulcano_phos_protn){
        
        generate_volcano_plots_peptide <- list()
        for(comp in names(db_execution_phos_protn$formule_contrast)){
          message(comp)
          generate_volcano_plots_peptide<-c(generate_volcano_plots_peptide,
                                            generate_volcano_plots(db_execution_phos_protn$differential_results, 
                                                                   data_type="peptide",
                                                                   comparison=comp,
                                                                   fc_thr=as.double(input$FC_thr_phos_protn),
                                                                   pval_fdr = input$pval_fdr_phos_protn,
                                                                   pval_thr=as.double(input$pval_thr_phos_protn)))
        }
        db_execution_phos_protn$peptide_vulcano = generate_volcano_plots_peptide
        # Generate tabPanels in a for loop
        tabs_pep_vulcano_phos_protn <- list()
        for (i in seq_along(generate_volcano_plots_peptide)) {
          plot_id <- paste0(names(generate_volcano_plots_peptide)[i], "_pep_phos_protn")
          message(plot_id)
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(generate_volcano_plots_peptide[[names(generate_volcano_plots_peptide)[my_i]]])
          })
          
          tabs_pep_vulcano_phos_protn[[i]] <- tabPanel(
            title = paste(names(generate_volcano_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential phospho-peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_peptide_phos_protn"), tabs_pep_vulcano_phos_protn))
        )
      } else{
        db_execution_phos_protn$peptide_vulcano = NULL
      }
    })
    
    output$render_mds_peptide_diff_phos_protn <- renderUI({
      if (input$mds_diff_peptide_phos_protn) {
        tagList(
          tags$h3("MDS based on differential phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('mds_peptide_diff_phos_protn')",
            plotOutput("small_mds_peptide_diff_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_differential_MDS = NULL
      }
    })
    output$small_mds_peptide_diff_phos_protn <- renderPlot({
      generate_mds_peptide_diff_phos_protn()
    })
    
    output$render_pca_peptide_diff_phos_protn <- renderUI({
      if (input$pca_diff_peptide_phos_protn) {
        tagList(
          tags$h3("PCA based on differential phospho-peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_phos_protn('pca_peptide_diff_phos_protn')",
            plotOutput("small_pca_peptide_diff_phos_protn")
          )
        )
      } else{
        db_execution_phos_protn$peptide_differential_PCA = NULL
      }
    })
    output$small_pca_peptide_diff_phos_protn <- renderPlot({
      generate_pca_peptide_diff_phos_protn()
    })
    
  })
  
  ## PhosProTN_with_prot: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn_phos_protn, {
    output$render_enrichement_analysis_phos_protn <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution_phos_protn$enrichment_results <- perform_enrichment_analysis(differential_results = db_execution_phos_protn$differential_results,
                                                                        enrichR_custom_DB = T,
                                                                        enrich_filter_DBs=input$DB_enrichment_phos_protn,    
                                                                        overlap_size_enrich_thr=as.double(input$os_enrich_phos_protn),
                                                                        pval_fdr_enrich = input$pval_fdr_phos_protn,
                                                                        pval_enrich_thr=as.double(input$pval_thr_phos_protn),
                                                                        dirOutput=db_execution_phos_protn$dirOutput, 
                                                                        with_background = input$enrich_with_background_phos_protn)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich_phos_protn), regex=","))
        
        db_execution_phos_protn$parameter <- c(db_execution_phos_protn$parameter,
                                         "Enrichment databases selected: "=paste(input$DB_enrichment_phos_protn, collapse = ", "),
                                         "P.value type used for enrichment: "=input$pval_fdr_phos_protn,
                                         "P.value threshold for enrichment significance: "=input$pval_thr_phos_protn,
                                         "Overlap size threshold for enrichment significance: "=input$os_enrich_phos_protn,
                                         "Enrichment filter terms: "=if(length(terms_enrich)>0){paste(terms_enrich, collapse = ", ")}else{"None"},
                                         "Enrichment with background: "=input$enrich_with_background_phos_protn)
        
        
        plots_down <- enrichment_figure(enr_df = db_execution_phos_protn$enrichment_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution_phos_protn$enrichment_results$anno_class), input$DB_enrichment_phos_protn))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(plots_down)) {
          plot_id <- names(plots_down)[i]
          height_id <- max(min(20, length(unique(plots_down[[names(plots_down)[i]]]$data$y_col))*0.4),3)*96
          message(paste0("Height for ",names(plots_down)[i], ": ", height_id))
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlot({
              plots_down[[names(plots_down)[my_i]]]
            }, height = height_id)
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(plots_down)[i]),
            plotOutput(plot_id, height = height_id)
          )
        }
        
        tagList(
          tags$h2("Enrichment Analysis"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_enrichment_phos_protn"), tabs))
        )
        
      })
    })
  })
  ## PhosProTN_with_prot: stringdb analysis ----
  observeEvent(input$execute_stringdb_analysis_btn_phos_protn, {
    output$render_stringdb_phos_protn <- renderUI({
      isolate({
        withProgress(message = "STRINGdb analysis in process, please wait!", {
          
          db_execution_phos_protn$stringdb_res <- STRINGdb_network(differential_results = db_execution_phos_protn$differential_results,
                                                        species=input$taxonomy_phos_protn, 
                                                        dirOutput=db_execution_phos_protn$dirOutput, 
                                                        score_thr=input$score_thr_stringdb_phos_protn,
                                                        shiny = T)
          db_execution_phos_protn$parameter <- c(db_execution_phos_protn$parameter,
                                           "STRINGdb taxonomy: "=input$taxonomy_phos_protn,
                                           "STRINGdb score threshold: "=input$score_thr_stringdb_phos_protn)
          tagList(
            tags$h2("STRINGdb analysis"),
            fluidRow(
              selectInput("stringdb_show_phos_protn", label = "Select StringDB to show: (click on STRING logo to open the results on stringDB website)", 
                          choices = names(db_execution_phos_protn$stringdb_res), width = "15%"),
              actionButton("stringdb_selected_phos_protn", "Select!", width = "10%")  
            ),
            tags$div(id = "stringEmbedded")
          )
        })
      })
    })
  })
  
  observeEvent(input$stringdb_selected_phos_protn, {
    js$loadStringData(input$taxonomy_phos_protn, db_execution_phos_protn$stringdb_res[[input$stringdb_show_phos_protn]], input$score_thr_stringdb_phos_protn)
  })
  
  ## PhosProTN_with_prot: kinase tree analysis ----
  observeEvent(input$execute_kinase_tree_analysis_btn_phos_protn, {
    output$render_kinase_tree_phos_protn <- renderUI({
      isolate({
        withProgress(message = "Kinase Tree analysis in process, please wait!", {
          
          db_execution_phos_protn$kinase_tree_res <- kinase_tree(proteome_data = db_execution_phos_protn$normalized_data, 
                                                      differential_results = db_execution_phos_protn$differential_results, 
                                                      formule_CORAL = db_execution_phos_protn$formule_contrast, 
                                                      dirOutput=db_execution_phos_protn$dirOutput, 
                                                      phosR_thr = input$score_thr_phosr_phos_protn, 
                                                      species = input$taxonomy_kinase_phos_protn)
          db_execution_phos_protn$parameter <- c(db_execution_phos_protn$parameter,
                                           "Kinase tree taxonomy: "=input$taxonomy_kinase_phos_protn,
                                           "Score thr for PhosR: "=input$score_thr_phosr_phos_protn)
          if(input$taxonomy_kinase_phos_protn == "Homo sapiens"){
            tagList(
              tags$h2("Kinase Tree analysis"),
              fluidRow(
                selectInput("kinase_tree_show_phos_protn", label = "Select Kinase Tree to show:", 
                            choices = names(db_execution_phos_protn$kinase_tree_res), width = "15%"),
                actionButton("kinase_tree_selected_phos_protn", "Select!", width = "10%")  
              ),
              imageOutput("render_kin_tree_phos_protn", height = "auto")
            )
          } else{
            tagList(
              tags$h2("Kinase Tree analysis"),
              tags$h4("For Mouse the graphical representation of the kinome tree is not done. The results can be downloaded."),
              tags$hr()
            )
          }
        })
      })
    })
  })
  
  observeEvent(input$kinase_tree_selected_phos_protn, {
    output$render_kin_tree_phos_protn <- renderImage({
      isolate({
        list(src = paste0(db_execution_phos_protn$dirOutput, "pics/kinaseTree/",input$kinase_tree_show_phos_protn,"_kinase_Tree_CORAL.svg"),
             alt = "Kinase Tree"
        )
      })
    }, deleteFile = FALSE)
  })
  
  # PhosProTN_with_prot: download results ----
  output$download_proteome_phos_protn <- downloadHandler(
    filename = "results_PhosProTN_with_proteome.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            message(db_execution_phos_protn$dirOutput)
            setProgress(value = 0.01)
            
            # Generate report
            params <- list(
              doc_title = input$title_exp_phos_protn,
              description = input$description_exp_phos_protn,
              readPD_files = if (input$sw_analyzer_phos_protn == "PD") {TRUE} else {FALSE},
              readMQ_files = if (input$sw_analyzer_phos_protn == "MQ") {TRUE} else {FALSE},
              impute_algorithm = if(input$advance_filter_phos_protn){input$impute_algorithm_phos_protn} else {"norm_phosr"},
              db_execution = reactiveValuesToList(db_execution_phos_protn),
              file_input_phospho = paste(db_execution_phos_protn$dirOutput, "input_phospho", sep = ""),
              file_input_proteome = paste(db_execution_phos_protn$dirOutput, "input_proteome", sep = ""),
              batch_corr_exe = if(input$batch_correction_phos_protn){input$batch_correction_col_phos_protn}else{NULL},
              prot_boxplot = if(input$boxplot_protein_phos_protn | input$heatmap_protein_phos_protn){input$list_proteins_phos_protn}else{NULL},
              fc_thr = if(is.null(input$FC_thr_phos_protn)){"0.75"}else{input$FC_thr_phos_protn},
              pval_fdr = input$pval_fdr_phos_protn,
              pval_thr = if(is.null(input$pval_thr_phos_protn)){"0.05"}else{input$pval_thr_phos_protn},
              pval_fdr_enrich = input$pval_fdr_enrich_phos_protn,
              pval_enrich_thr = if(is.null(input$pvalue_enrich_phos_protn)){"0.05"}else{input$pvalue_enrich_phos_protn},
              overlap_size_enrich_thr = if(is.null(input$os_enrich_phos_protn)){as.integer(5)}else{input$os_enrich_phos_protn},
              enrich_filter_term = input$terms_enrich_phos_protn,
              enrich_filter_DBs = input$DB_enrichment_phos_protn,
              taxonomy=input$taxonomy_phos_protn, 
              score_thr=input$score_thr_stringdb_phos_protn,
              dirOutput = db_execution_phos_protn$dirOutput
            )
            
            # Render in background the report
            p = callr::r_bg(
              func = function(db_execution_phos_protn, params, dirOutput, env) {
                rmarkdown::render("R/phosprotn_protn_report.Rmd",
                                  output_file = "phosprotn_with_proteome_report.html",
                                  output_dir = dirOutput,
                                  params = params,
                                  envir = env
                )
              },
              args = list(db_execution_phos_protn, params, db_execution_phos_protn$dirOutput, new.env(parent = globalenv())),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            # Saving RData in background
            db_results_PhosProTN_with_proteome = reactiveValuesToList(db_execution_phos_protn)
            db_results_PhosProTN_with_proteome <- db_results_PhosProTN_with_proteome[!(unlist(lapply(db_results_PhosProTN_with_proteome, is.null)))]
            p_rdata = callr::r_bg(
              func = function(db_results_PhosProTN_with_proteome, dirOutput) {
                save(db_results_PhosProTN_with_proteome, file = paste0(dirOutput,"db_results_PhosProTN_with_proteome.RData"))
              },
              args = list(db_results_PhosProTN_with_proteome, db_results_PhosProTN_with_proteome$dirOutput),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            
            # Prepare file for the download
            if(length(db_execution_phos_protn$normalized_data)>0){
              save_abundance_tables(proteome_data = db_execution_phos_protn$normalized_data, 
                                    dirOutput = db_execution_phos_protn$dirOutput)
            }
            setProgress(value = 0.1)
            
            if(length(db_execution_phos_protn$differential_results)>0){
              save_differential_analysis_table(proteome_data = db_execution_phos_protn$normalized_data,
                                               differential_results = db_execution_phos_protn$differential_results,
                                               dirOutput=db_execution_phos_protn$dirOutput)
            }
            setProgress(value = 0.2)
            
            
            if(input$phospho_percentage_plot_phos & !is.null(db_execution_phos_protn$phospho_percentage)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/phospho_percentage.pdf"), 
                     plot = db_execution_phos_protn$phospho_percentage, 
                     create.dir = T, width = 7, height = 3)
            } else if("phospho_percentage.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/phospho_percentage.pdf"))
            }
            
            if(input$abundance_plot_phos_protn & !is.null(db_execution_phos_protn$generate_abundance)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/missing_available_abundance_proteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_abundance$proteome_plot, 
                     create.dir = T, width = 7, height = 5)
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/missing_available_abundance_phosphoproteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_abundance$phospho_plot, 
                     create.dir = T, width = 7, height = 5)
            } else if("missing_available_abundance_proteomics.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/missing_available_abundance_*"))
            }
            setProgress(value = 0.25)
            
            if(input$peptide_distribution_phos_protn & !is.null(db_execution_phos_protn$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_per_protein_proteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_peptide_distribution$proteome_plot, 
                     create.dir = T, width = 7, height = 5)
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_per_protein_phosphoproteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_peptide_distribution$phospho_plot, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_per_protein_proteomics.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_per_protein_*"))
            }
            setProgress(value = 0.30)
            
            if(input$raw_violin_prot_phos_protn & !is.null(db_execution_phos_protn$raw_proteome_abundance_distribution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/raw_abundance_distribution_proteome.pdf"), 
                     plot = db_execution_phos_protn$raw_proteome_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("raw_abundance_distribution.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/raw_abundance_distribution_proteome.pdf"))
            }
            
            if(input$raw_violin_pep_phos_protn & !is.null(db_execution_phos_protn$raw_phospho_abundance_distribution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/raw_abundance_distribution_phosphoproteomics.pdf"), 
                     plot = db_execution_phos_protn$raw_phospho_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("raw_abundance_distribution.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/raw_abundance_distribution_phosphoproteomics.pdf"))
            }
            
            if(input$complexity_plot_phos_protn & !is.null(db_execution_phos_protn$generate_complexity)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/complexity_proteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_complexity$proteome_plot, 
                     create.dir = T, width = 7, height = 5)
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/complexity_phosphoproteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_complexity$phospho_plot, 
                     create.dir = T, width = 7, height = 5)
            } else if("complexity_proteomics.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/complexity_*"))
            }
            setProgress(value = 0.33)
            
            if(input$protein_violin_phos_protn & !is.null(db_execution_phos_protn$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution_phos_protn$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_abundance_distribution.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_abundance_distribution.pdf"))
            }
            setProgress(value = 0.35)
            
            if(input$peptide_violin_phos_protn & !is.null(db_execution_phos_protn$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution_phos_protn$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_abundance_distribution.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_abundance_distribution.pdf"))
            }
            setProgress(value = 0.40)
            
            if(input$mds_protein_phos_protn & !is.null(db_execution_phos_protn$protein_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution_phos_protn$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_MDS.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_MDS.pdf"))
            }
            setProgress(value = 0.43)
            
            if(input$mds_peptide_phos_protn & !is.null(db_execution_phos_protn$peptide_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution_phos_protn$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_MDS.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_MDS.pdf"))
            }
            setProgress(value = 0.45)
            
            if(input$pca_protein_phos_protn & !is.null(db_execution_phos_protn$protein_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution_phos_protn$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_PCA.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_PCA.pdf"))
            }
            setProgress(value = 0.47)
            
            if(input$pca_peptide_phos_protn & !is.null(db_execution_phos_protn$peptide_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution_phos_protn$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_PCA.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_PCA.pdf"))
            }
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein_phos_protn & !is.null(db_execution_phos_protn$protein_boxplot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution_phos_protn$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_boxplot.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_boxplot.pdf"))
            }
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein_phos_protn & !is.null(db_execution_phos_protn$protein_heatmap)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution_phos_protn$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_heatmap.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_heatmap.pdf"))
            }
            setProgress(value = 0.55)
            
            if(!is.null(db_execution_phos_protn$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution_phos_protn$protein_differential_barplot, 
                     create.dir = T,width = 17, height = 6)
            } else if("protein_differential_barplot.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_differential_barplot.pdf"))
            }
            setProgress(value = 0.58)
            
            if(!is.null(db_execution_phos_protn$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution_phos_protn$peptide_differential_barplot, 
                     create.dir = T, width = 17, height = 6)
            } else if("peptide_differential_barplot.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_differential_barplot.pdf"))
            }
            setProgress(value = 0.60)
            
            if(!is.null(db_execution_phos_protn$peptide_upset_plot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_upset_plot.pdf"), 
                     plot = db_execution_phos_protn$peptide_upset_plot, 
                     create.dir = T, width = 12, height = 6)
            } else if("peptide_upset_plot.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_upset_plot.pdf"))
            }
            setProgress(value = 0.63)
            
            if(!is.null(db_execution_phos_protn$peptide_ma_plot)){
              dir.create(file.path(paste0(db_execution_phos_protn$dirOutput,"pics/"), "peptide_ma_plot"), showWarnings = FALSE)
              for(comp in names(db_execution_phos_protn$peptide_ma_plot)){
                htmlwidgets::saveWidget(db_execution_phos_protn$peptide_ma_plot[[comp]], 
                                        file = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"))
                webshot2::webshot(url = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"), 
                                  file = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_phos_protn$dirOutput,"pics/peptide_ma_plot"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_phos_protn$protein_vulcano)){
              dir.create(file.path(paste0(db_execution_phos_protn$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos_protn$protein_vulcano)){
                # plotly::save_image(db_execution_phos_protn$protein_vulcano[[comp]], 
                #                    file = paste0(str_replace_all(db_execution_phos_protn$dirOutput, pattern="\\\\", replacement="/"),"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos_protn$protein_vulcano[[comp]], 
                                        file = paste0(db_execution_phos_protn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution_phos_protn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"), 
                                  file = paste0(db_execution_phos_protn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_phos_protn$dirOutput,"pics/protein_vulcano"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_phos_protn$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution_phos_protn$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos_protn$peptide_vulcano)){
                # plotly::save_image(db_execution_phos_protn$peptide_vulcano[[comp]], 
                #                    file = paste0(str_replace_all(db_execution_phos_protn$dirOutput, pattern="\\\\", replacement="/"),"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos_protn$peptide_vulcano[[comp]], 
                                        file = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"), 
                                  file = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_phos_protn$dirOutput,"pics/peptide_vulcano"))
            }
            setProgress(value = 0.68)
            
            if(!is.null(db_execution_phos_protn$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution_phos_protn$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_differential_MDS.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_differential_MDS.pdf"))
            }
            setProgress(value = 0.69)
            
            if(!is.null(db_execution_phos_protn$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution_phos_protn$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_MDS.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_differential_MDS.pdf"))
            }
            setProgress(value = 0.70)
            
            if(!is.null(db_execution_phos_protn$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution_phos_protn$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_differential_PCA.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/protein_differential_PCA.pdf"))
            }
            setProgress(value = 0.72)
            
            if(!is.null(db_execution_phos_protn$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution_phos_protn$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_PCA.pdf" %in% list.files(paste0(db_execution_phos_protn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_phos_protn$dirOutput,"pics/peptide_differential_PCA.pdf"))
            }
            setProgress(value = 0.75)
            
            if(length(db_execution_phos_protn$enrichment_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich_phos_protn), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution_phos_protn$enrichment_results,
                                              category = c("down","up"), 
                                              enrich_filter_term = terms_enrich,
                                              save=T, 
                                              dirOutput = db_execution_phos_protn$dirOutput)
            } 
            setProgress(value = 0.82)
            
            if(length(db_execution_phos_protn$stringdb_res)>0){
              tmp_res <- STRINGdb_network(differential_results = db_execution_phos_protn$differential_results,
                                          species=input$taxonomy_phos_protn, 
                                          dirOutput=db_execution_phos_protn$dirOutput,
                                          score_thr=input$score_thr_stringdb_phos_protn,
                                          shiny = F)
              
            } 
            setProgress(value = 0.95)
            
            # Write tsv file with parameter
            params <- data.table("Parameter" = names(db_execution_phos_protn$parameter),
                                 "Value" = unlist(db_execution_phos_protn$parameter))
            fwrite(params, paste0(db_execution_phos_protn$dirOutput,"parameters_used.txt"), sep = "\t", col.names = F)
            
            
            #Get results Report
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p$read_output()
            p$wait(30000)
            for (i in 1:15) {
              p$read_output()
              p$wait(1000*60)  
            }
            if(p$is_alive() | is.null(p$get_result())){
              p$kill()
              print("\n ERROR: An error occur during the report rendering. \n ")
            } else{
              report<-p$get_result()
              p$kill()
              message("Render report DONE.")
            }
            
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p_rdata$read_output()
            p_rdata$wait(30000)
            for (i in 1:15) {
              p_rdata$read_output()
              p_rdata$wait(1000*60)
            }
            if(p_rdata$is_alive() | is.null(p_rdata$get_result())){
              p_rdata$kill()
              print("\n ERROR: An error occur during the RData saving. \n ")
            } else{
              report<-p_rdata$get_result()
              p_rdata$kill()
              message("RData saved.")
            }
            
            # Save RData db_execution_phos_protn
            # setProgress(value = 0.95, message = "Saving RData...")
            # db_results_PhosProTN_with_proteome = reactiveValuesToList(db_execution_phos_protn)
            # db_results_PhosProTN_with_proteome <- db_results_PhosProTN_with_proteome[!(unlist(lapply(db_results_PhosProTN_with_proteome, is.null)))]
            # save(db_results_PhosProTN_with_proteome, file = paste0(db_results_PhosProTN_with_proteome$dirOutput,"db_results_PhosProTN_with_proteome.RData"))
            
            #Save folder for the download
            oldwd <- getwd()
            message(db_execution_phos_protn$dirOutput)
            setwd(db_execution_phos_protn$dirOutput)
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  ## PhosProTN_with_prot: full screen trigger ----
  
  # ReactiveVal for currently selected plot to fullscreen
  selected_plot_phos_protn <- reactiveVal(NULL)
  
  # Update selected_plot when JS sends fullscreen_trigger id
  observeEvent(input$fullscreen_trigger_phos_protn, {
    selected_plot_phos_protn(input$fullscreen_trigger_phos_protn)
  })
  
  # Render fullscreen plot dynamically based on selected_plot()
  output$fullscreen_plot_phos_protn <- renderPlot({
    req(selected_plot_phos_protn())
    switch(selected_plot_phos_protn(),
           "phospho_percentage_plot_phos_protn" = generate_phospho_percentage_plot_phos_protn()(size_text=8, zoom=T) + ggtitle("Percentage of phosphosite residue")+theme(text=element_text(size=25), axis.text.y = element_text(size = 25)),
           "abundance_plot_phos_protn_prot" = generate_abundance_phos_protn_prot() + ggtitle("Percentage missing values respect detected abundance - Proteomics")+theme(text=element_text(size=25)),
           "abundance_plot_phos_protn_phos" = generate_abundance_phos_protn_phos() + ggtitle("Percentage missing values respect detected abundance - Phospho-proteomics")+theme(text=element_text(size=25)),
           "peptide_distribution_plot_phos_protn_prot" = generate_peptide_distribution_phos_protn_prot() + ggtitle("N° peptides per proteins - Proteomics")+theme(text=element_text(size=25)),
           "peptide_distribution_plot_phos_protn_phos" = generate_peptide_distribution_phos_protn_phos() + ggtitle("N° peptides per proteins - Phospho-proteomics")+theme(text=element_text(size=25)),
           "raw_violin_plot_prot_phos_protn" = generate_raw_violin_prot_phos_protn() + ggtitle("Raw abundance distribution - Proteome")+theme(text=element_text(size=25)),
           "raw_violin_plot_pep_phos_protn" = generate_raw_violin_pep_phos_protn() + ggtitle("Raw abundance distribution - Phospho-proteome")+theme(text=element_text(size=25)),
           "complexity_plot_phos_protn_prot" = generate_complexity_phos_protn_prot() + ggtitle("Complexity plot of raw abundance - Proteomics")+theme(text=element_text(size=25)),
           "complexity_plot_phos_protn_phos" = generate_complexity_phos_protn_prot() + ggtitle("Complexity plot of raw abundance - Phospho-proteomics")+theme(text=element_text(size=25)),
           "protein_violin_plot_phos_protn" = generate_protein_violin_phos_protn() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "peptide_violin_plot_phos_protn" = generate_peptide_violin_phos_protn() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "mds_protein_phos_protn" = generate_mds_protein_phos_protn() + ggtitle("MDS based on protein")+theme(text=element_text(size=25)),
           "mds_peptide_phos_protn" = generate_mds_peptide_phos_protn() + ggtitle("MDS based on peptides")+theme(text=element_text(size=25)),
           "pca_protein_phos_protn" = generate_pca_protein_phos_protn() + ggtitle("PCA based on protein")+theme(text=element_text(size=25)),
           "pca_peptide_phos_protn" = generate_pca_peptide_phos_protn() + ggtitle("PCA based on peptides")+theme(text=element_text(size=25)),
           "protein_boxplot_phos_protn" = generate_protein_boxplot_phos_protn() + ggtitle("Boxplot selected proteins")+theme(text=element_text(size=25)),
           "protein_heatmap_phos_protn" = generate_protein_heatmap_phos_protn() + ggtitle("Heatmap selected proteins")+theme(text=element_text(size=25)),
           "peptide_diff_barplot_phos_protn" = generate_peptide_diff_barplot_phos_protn()(8, zoom=T) + ggtitle("N° differential phospho-peptides")+theme(text=element_text(size=25)),
           "peptide_upset_phos_protn" = generate_peptide_upset_phos_protn(),
           "mds_peptide_diff_phos_protn" = generate_mds_peptide_diff_phos_protn() + ggtitle("MDS based on differential phospho-peptides")+theme(text=element_text(size=25)),
           "pca_peptide_diff_phos_protn" = generate_pca_peptide_diff_phos_protn() + ggtitle("PCA based on differential phospho-peptides")+theme(text=element_text(size=25)),
           # default fallback:
           NULL
    )
  })
  
  ##############################################################################
  
  ### InteracTN ----
  # Optional visibility based on the selection ----
  
  ## InteracTN: Visibility of the proteomics files for ProTN ----
  output$input_proteome_interactn <- renderUI({
    if (input$sw_analyzer_interactn == "PD"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_interactn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_interactn", "Select the PEP file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_proteome_interactn", "Select the PROT file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer_interactn == "MQ_ev"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_interactn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_interactn", "Select the EVIDENCE file of the PROTEOMICS..."),
        )
      )
    } else if (input$sw_analyzer_interactn == "MQ_prot"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_interactn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_interactn", "Select the Peptides.txt file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("prot_file_proteome_interactn", "Select the ProteinGroups.txt file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer_interactn == "SP"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_interactn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_interactn", "Select the Spectronaut report file of the PROTEOMICS..."),
        )
      )
    } else if(input$sw_analyzer_interactn == "FP"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_interactn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_interactn", "Select the combined_modified_peptide.tsv file of the PROTEOMICS..."),
        )
      )
    } else{
      tagList(
        tags$p("BACK")
      )
    }
  })
  
  ## InteracTN: textbox for batch correction----
  output$batch_correction_ui_interactn <- renderUI({ 
    if(input$batch_correction_interactn){
      textInput("batch_correction_col_interactn", "Column in Annotation file with the batch:")
    } 
  })
  ## InteracTN: advance filters----
  output$advance_filter_ui_interactn <- renderUI({ 
    if(input$advance_filter_interactn){
      tagList(
        numericInput("NA_allow_condition_interactn", "N° missing value allow per condition", value = 0, min = 0, max = 5),
        numericInput("min_peptide_protein_interactn", "Minimum peptide per protein", value = 1, min = 1),
        selectizeInput("impute_algorithm_interactn", "Select impute algorithm:",
                       choices = list("PhosR and normalization" = "phosr_norm", 
                                      "Gaussian estimation and normalization" = "gaussian_norm",
                                      "missForest and normalization" = "missForest_norm",
                                      "Pre-normalization and PhosR" = "norm_phosr", 
                                      "Pre-normalization and Gaussian estimation" = "norm_gaussian",
                                      "Pre-normalization and missForest" = "norm_missForest",
                                      "Pre-normalization and pcaMethods" = "norm_pcaMethods"),
                       selected = "norm_phosr", multiple = FALSE
        ),
        textInput("sample_column_interactn", "Column name with the sample name:", value = "Sample")
      )
    } 
  })
  
  ## InteracTN: textbox for list proteins ----
  output$list_protein_ui_interactn <- renderUI({ 
    if(input$boxplot_protein_interactn | input$heatmap_protein_interactn){
      textInput("list_proteins_interactn", "List proteins to show (separate by: \",\"):")
    } 
  })
  
  ## InteracTN: show parameter for differential analysis ----
  output$differential_params_ui_interactn <- renderUI({ 
    if(input$differential_analysis_checkbox_interactn){
      tagList(
        tags$label("Write in each line a different comparison"),
        tags$label("(right click to add row)"),
        rHandsontableOutput('render_formule_contrast_table_interactn'),
        # textAreaInput("formule_contrast", "Write in each line a different comparison", rows = 4),
        textInput("FC_thr_interactn", "Fold change threshold for significance:",value = 0.5),
        radioButtons("pval_fdr_interactn", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_val"),
        textInput("pval_thr_interactn", "P.value threshold for significance:", value = 0.05),
        actionButton("execute_differential_analysis_btn_interactn", "Run!"),
        checkboxInput("protein_diff_table_interactn", "Proteins differentiated table", FALSE),
        checkboxInput("peptide_diff_table_interactn", "Peptides differentiated table", FALSE),
        checkboxInput("protein_diff_barplot_interactn", "Proteins differentiated barplot", TRUE),
        checkboxInput("peptide_diff_barplot_interactn", "Peptides differentiated barplot", FALSE),
        checkboxInput("protein_upset_interactn", "Proteins upset plot", FALSE),
        checkboxInput("peptide_upset_interactn", "Peptides upset plot", FALSE),
        checkboxInput("protein_ma_plot_interactn", "Proteins MA plot", FALSE),
        checkboxInput("peptide_ma_plot_interactn", "Peptides MA plot", FALSE),
        checkboxInput("protein_vulcano_interactn", "Proteins vulcano plot", FALSE),
        checkboxInput("peptide_vulcano_interactn", "Peptides vulcano plot", FALSE),
        checkboxInput("mds_diff_protein_interactn", "MDS based on diffential protein", FALSE),
        checkboxInput("mds_diff_peptide_interactn", "MDS based on diffential peptide", FALSE),
        checkboxInput("pca_diff_protein_interactn", "PCA based on diffential protein", FALSE),
        checkboxInput("pca_diff_peptide_interactn", "PCA based on diffential peptide", FALSE),
        tags$h3("Enrichment Analysis:"),
        checkboxInput("enrichment_analysis_interactn", "Execute enrichment analysis", FALSE),
        uiOutput("enrichment_params_ui_interactn"),
        tags$h3("STRINGdb network:"),
        checkboxInput("stringdb_analysis_interactn", "Execute STRINGdb", FALSE),
        uiOutput("stringdb_params_ui_interactn")
      )
    } else{
      # Reset UI elements
      updateCheckboxInput(session, "enrichment_analysis_interactn", value = FALSE)
      updateCheckboxInput(session, "stringdb_analysis_interactn", value = FALSE)
      
      db_execution_interactn$formule_contrast <- list()
      db_execution_interactn$dt_formule_contrast <- data.table("Name"=c("","","",""),"Formule"=c("","","",""))
      db_execution_interactn$differential_results <- list()
      
      updateCheckboxInput(session, "protein_diff_barplot_interactn", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_barplot_interactn",  value = FALSE)
      updateCheckboxInput(session, "protein_diff_table_interactn", value = FALSE)
      updateCheckboxInput(session, "peptide_diff_table_interactn", value = FALSE)
      updateCheckboxInput(session, "protein_upset_interactn", value = FALSE)
      updateCheckboxInput(session, "peptide_upset_interactn", value = FALSE)
      updateCheckboxInput(session, "protein_vulcano_interactn",  value =  FALSE)
      updateCheckboxInput(session, "peptide_vulcano_interactn", value = FALSE)
      updateCheckboxInput(session, "protein_ma_plot_interactn", value = FALSE)
      updateCheckboxInput(session, "peptide_ma_plot_interactn", value = FALSE)
      updateCheckboxInput(session, "mds_diff_protein_interactn", value = FALSE)
      updateCheckboxInput(session, "mds_diff_peptide_interactn", value = FALSE)
      updateCheckboxInput(session, "pca_diff_protein_interactn", value = FALSE)
      updateCheckboxInput(session, "pca_diff_peptide_interactn", value = FALSE)
      
      db_execution_interactn$protein_differential_barplot <- NULL
      db_execution_interactn$peptide_differential_barplot <- NULL
      db_execution_interactn$protein_upset_plot <- NULL
      db_execution_interactn$peptide_upset_plot <- NULL
      db_execution_interactn$protein_ma_plot <- NULL
      db_execution_interactn$peptide_ma_plot <- NULL
      db_execution_interactn$protein_vulcano <- NULL
      db_execution_interactn$peptide_vulcano <- NULL
      db_execution_interactn$protein_differential_MDS <- NULL
      db_execution_interactn$peptide_differential_MDS <- NULL
      db_execution_interactn$protein_differential_PCA <- NULL
      db_execution_interactn$peptide_differential_PCA <- NULL
      
      output$render_differential_analysis_interactn <- renderUI({NULL})
      output$render_protein_diff_table_interactn <- renderUI({NULL})
      output$render_peptide_diff_table_interactn <- renderUI({NULL})
      output$render_protein_diff_barplot_interactn <- renderUI({NULL})
      output$render_peptide_diff_barplot_interactn <- renderUI({NULL})
      output$render_protein_upset_interactn <- renderUI({NULL})
      output$render_peptide_upset_interactn <- renderUI({NULL})
      output$render_protein_ma_plot_interactn <- renderUI({NULL})
      output$render_peptide_ma_plot_interactn <- renderUI({NULL})
      output$render_protein_vulcano_interactn <- renderUI({NULL})
      output$render_peptide_vulcano_interactn <- renderUI({NULL})
      output$render_mds_protein_diff_interactn <- renderUI({NULL})
      output$render_mds_peptide_diff_interactn <- renderUI({NULL})
      output$render_pca_protein_diff_interactn <- renderUI({NULL})
      output$render_pca_peptide_diff_interactn <- renderUI({NULL})
    }
  })
  
  output$render_formule_contrast_table_interactn <- renderRHandsontable({
    rhandsontable(db_execution_interactn$dt_formule_contrast, rowHeaders = NULL, stretchH = "all")
  })
  
  ## InteracTN: show enrichment parameter ----
  output$enrichment_params_ui_interactn <- renderUI({ 
    if(input$enrichment_analysis_interactn){
      tagList(
        # radioButtons("enrichR_universe", "Execute enrichment of the whole Universe", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
        selectizeInput("DB_enrichment_interactn", "DB to analyse:",
                       choices = lapply(split(read_tsv("data/dbs_enrichR.txt", col_names = FALSE)$X1,
                                              read_tsv("data/dbs_enrichR.txt", col_names = FALSE)[,2]), as.list),
                       selected = NULL, multiple = TRUE
        ),
        textInput("terms_enrich_interactn", "Terms to search (separated by \",\"):"),
        radioButtons("pval_fdr_enrich_interactn", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_adj"),
        textInput("pvalue_enrich_interactn", "P.value threshold for significance:", value = 0.05),
        sliderInput("os_enrich_interactn", "Overlap size thr for enrichment", 1, 30, step = 1, value = 5),
        checkboxInput("enrich_with_background_interactn", "Enrichment with background", FALSE),
        actionButton("execute_enrichment_analysis_btn_interactn", "Run!")
      )
    } else{
      db_execution_interactn$enrichment_results <- list()
      output$render_enrichement_analysis_interactn <- renderUI({NULL})
    }
  })
  
  ## InteracTN: show stringdb parameter ----
  output$stringdb_params_ui_interactn <- renderUI({
    if(input$stringdb_analysis_interactn){
      tagList(
        selectizeInput("taxonomy_interactn", "NCBI Taxonomy ID", 
                       choice = data.table::fread("data/subset_tax.csv", select = "name"), 
                       selected = "Homo sapiens", multiple = F),
        sliderInput("score_thr_stringdb_interactn", "Score thr for STRINGdb", 500, 1000, step = 10, value = 700),
        actionButton("execute_stringdb_analysis_btn_interactn", "Run!"),
        tags$br()
      )
    } else{
      db_execution_interactn$stringdb_res <- list()
      output$render_stringdb_interactn <- renderUI({NULL})
    }
  })
  
  ## InteracTN: function genereting plot ----
  generate_abundance_interactn <- reactive({
    req(input$abundance_plot_interactn)
    if(input$abundance_plot_interactn){
      generate_abundance_fig <- generate_abundance_plot(proteome_data = db_execution_interactn$proteome_data)$plot
      db_execution_interactn$generate_abundance = generate_abundance_fig
      generate_abundance_fig
    } else{
      db_execution_interactn$generate_abundance = NULL
    }
  })
  
  generate_peptide_distribution_interactn <- reactive({
    req(input$peptide_distribution_interactn)
    if(input$peptide_distribution_interactn){
      peptide_distribution_fig <- generate_peptide_distribution_plot(proteome_data = db_execution_interactn$proteome_data)$plot
      db_execution_interactn$generate_peptide_distribution = peptide_distribution_fig
      peptide_distribution_fig
    } else{
      db_execution_interactn$generate_peptide_distribution = NULL
    }
  })
  
  generate_raw_violin_interactn <- reactive({
    req(input$raw_violin_interactn)
    if(input$raw_violin_interactn){
      raw_abundance_distribution_fig <- plot_raw_abundance_distribution(proteome_data = db_execution_interactn$proteome_data,
                                                                        type = "protein")$plot
      db_execution_interactn$raw_abundance_distribution = raw_abundance_distribution_fig
      raw_abundance_distribution_fig
    } else{
      db_execution_interactn$raw_abundance_distribution = NULL
    }
  })
  
  generate_complexity_interactn <- reactive({
    req(input$complexity_plot_interactn)
    if(input$complexity_plot_interactn){
      generate_complexity_fig <- complexity_plot(proteome_data = db_execution_interactn$proteome_data)$plot
      db_execution_interactn$generate_complexity = generate_complexity_fig
      generate_complexity_fig
    } else{
      db_execution_interactn$generate_complexity = NULL
    }
  })
  
  generate_protein_violin_interactn <- reactive({
    req(input$protein_violin_interactn)
    if(input$protein_violin_interactn){
      protein_abundance_distribution_fig <- plot_abundance_distribution(proteome_data = db_execution_interactn$normalized_data,
                                                                        type = "protein")$plot
      db_execution_interactn$protein_abundance_distribution = protein_abundance_distribution_fig
      protein_abundance_distribution_fig
    } else{
      db_execution_interactn$protein_abundance_distribution = NULL
    }
  })
  
  generate_peptide_violin_interactn <- reactive({
    req(input$peptide_violin_interactn)
    if(input$peptide_violin_interactn){
      peptide_abundance_distirbution_fig <- plot_abundance_distribution(proteome_data = db_execution_interactn$normalized_data,
                                                                        type = "peptide")$plot
      db_execution_interactn$peptide_abundance_distirbution = peptide_abundance_distirbution_fig
      peptide_abundance_distirbution_fig
    } else{
      db_execution_interactn$peptide_abundance_distirbution = NULL
    }
  })
  
  generate_mds_protein_interactn <- reactive({
    req(input$mds_protein_interactn)
    if(input$mds_protein_interactn){
      mds_protein_fig <- mds_plot(proteome_data = db_execution_interactn$normalized_data,
                                  type = "protein")$plot
      db_execution_interactn$protein_MDS = mds_protein_fig
      mds_protein_fig
    } else{
      db_execution_interactn$protein_MDS = NULL
    }
  })
  
  generate_mds_peptide_interactn <- reactive({
    req(input$mds_peptide_interactn)
    if(input$mds_peptide_interactn){
      mds_peptide_fig <- mds_plot(proteome_data = db_execution_interactn$normalized_data,
                                  type = "peptide")$plot
      db_execution_interactn$peptide_MDS = mds_peptide_fig
      mds_peptide_fig
    } else{
      db_execution_interactn$peptide_MDS = NULL
    }
  })
  
  generate_pca_protein_interactn <- reactive({
    req(input$pca_protein_interactn)
    if(input$pca_protein_interactn){
      pca_protein_fig <- pca_plot(proteome_data = db_execution_interactn$normalized_data,
                                  type = "protein")$plot
      db_execution_interactn$protein_PCA = pca_protein_fig
      pca_protein_fig
    } else{
      db_execution_interactn$protein_PCA = NULL
    }
  })
  
  generate_pca_peptide_interactn <- reactive({
    req(input$pca_peptide_interactn)
    if(input$pca_peptide_interactn){
      pca_peptide_fig <- pca_plot(proteome_data = db_execution_interactn$normalized_data,
                                  type = "peptide")$plot
      db_execution_interactn$peptide_PCA = pca_peptide_fig
      pca_peptide_fig
    } else{
      db_execution_interactn$peptide_PCA = NULL
    }
  })
  
  generate_protein_boxplot_interactn <- reactive({
    req(input$boxplot_protein_interactn)
    if(input$boxplot_protein_interactn){
      req(input$list_proteins_interactn)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_interactn), regex=",")
      db_execution_interactn$parameter <- c(db_execution_interactn$parameter, "List proteins boxplot abundance: "=input$list_proteins_interactn)
      boxplot_protein_fig <- plot_selected_proteins(proteome_data = db_execution_interactn$normalized_data,
                                                    list_protein = unlist(list_proteins))$plot
      db_execution_interactn$protein_boxplot = boxplot_protein_fig
      boxplot_protein_fig
    } else{
      db_execution_interactn$protein_boxplot = NULL
    }
  })
  
  generate_protein_heatmap_interactn <- reactive({
    req(input$heatmap_protein_interactn)
    if(input$heatmap_protein_interactn){
      req(input$list_proteins_interactn)
      list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_interactn), regex=",")
      db_execution_interactn$parameter <- c(db_execution_interactn$parameter, "List proteins heatmap abundance: "=input$list_proteins_interactn)
      heatmap_protein_fig <- heatmap_selected_proteins(proteome_data = db_execution_interactn$normalized_data, list_protein = unlist(list_proteins))$plot
      db_execution_interactn$protein_heatmap = heatmap_protein_fig
      heatmap_protein_fig
    } else{
      db_execution_interactn$protein_heatmap = NULL
    }
  })
  
  generate_mds_protein_diff_interactn <- reactive({
    req(input$mds_diff_protein_interactn)
    if(input$mds_diff_protein_interactn){
      mds_protein_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                             proteome_data = db_execution_interactn$normalized_data,
                                                             type = "protein")$plot
      db_execution_interactn$protein_differential_MDS = mds_protein_diff_fig
      mds_protein_diff_fig
    } else{
      db_execution_interactn$protein_differential_MDS = NULL
    }
  })
  
  generate_mds_peptide_diff_interactn <- reactive({
    req(input$mds_diff_peptide_interactn)
    if(input$mds_diff_peptide_interactn){
      mds_peptide_diff_fig <- mds_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                             proteome_data = db_execution_interactn$normalized_data,
                                                             type = "peptide")$plot
      db_execution_interactn$peptide_differential_MDS = mds_peptide_diff_fig
      mds_peptide_diff_fig
    } else{
      db_execution_interactn$peptide_differential_MDS = NULL
    }
  })
  
  generate_pca_protein_diff_interactn <- reactive({
    req(input$pca_diff_protein_interactn)
    if(input$pca_diff_protein_interactn){
      pca_protein_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                             proteome_data = db_execution_interactn$normalized_data,
                                                             type = "protein")$plot
      db_execution_interactn$protein_differential_PCA = pca_protein_diff_fig
      pca_protein_diff_fig
    } else{
      db_execution_interactn$protein_differential_PCA = NULL
    }
  })
  
  generate_pca_peptide_diff_interactn <- reactive({
    req(input$pca_diff_peptide_interactn)
    if(input$pca_diff_peptide_interactn){
      pca_peptide_diff_fig <- pca_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                             proteome_data = db_execution_interactn$normalized_data,
                                                             type = "peptide")$plot
      db_execution_interactn$peptide_differential_PCA = pca_peptide_diff_fig
      pca_peptide_diff_fig
    } else{
      db_execution_interactn$peptide_differential_PCA = NULL
    }
  })
  
  generate_protein_diff_barplot_interactn <- reactive(function(size_text){
    req(input$protein_diff_barplot_interactn)
    if(input$protein_diff_barplot_interactn){
      ploft_diff_number <- generate_differential_barplots(db_execution_interactn$differential_results,
                                                          data_type="protein", size_text=size_text)$plot
      db_execution_interactn$protein_differential_barplot = ploft_diff_number
      ploft_diff_number
    } else{
      db_execution_interactn$protein_differential_barplot = NULL
    }
  })
  
  generate_peptide_diff_barplot_interactn <- reactive(function(size_text){
    req(input$peptide_diff_barplot_interactn)
    if(input$peptide_diff_barplot_interactn){
      ploft_diff_number_pep <- generate_differential_barplots(db_execution_interactn$differential_results,
                                                              data_type="peptide", size_text=size_text)$plot
      db_execution_interactn$peptide_differential_barplot = ploft_diff_number_pep
      ploft_diff_number_pep
    } else{
      db_execution_interactn$peptide_differential_barplot = NULL
    }
  })
  
  generate_protein_upset_interactn <- reactive({
    req(input$protein_upset_interactn)
    if(input$protein_upset_interactn){
      ploft_diff_number <- generate_upset_plot(db_execution_interactn$differential_results,
                                               type="protein", 
                                               DE_class = "all")$plot
      db_execution_interactn$protein_upset_plot = ploft_diff_number
      ploft_diff_number
    } else{
      db_execution_interactn$protein_upset_plot = NULL
    }
  })
  
  generate_peptide_upset_interactn <- reactive({
    req(input$peptide_upset_interactn)
    if(input$peptide_upset_interactn){
      ploft_diff_number_pep <- generate_upset_plot(db_execution_interactn$differential_results,
                                                   type="peptide", 
                                                   DE_class = "all")$plot
      db_execution_interactn$peptide_upset_plot = ploft_diff_number_pep
      ploft_diff_number_pep
    } else{
      db_execution_interactn$peptide_upset_plot = NULL
    }
  })
  
  # InteracTN: Execution pipeline ----
  observeEvent(input$report_proteome_interactn, {
    
    output$protn_results_ui_interactn <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
              # Reset other analysis
              db_execution_interactn$parameter <- list()
              updateCheckboxInput(session, "differential_analysis_checkbox_interactn", value = FALSE)
              
              message(session$token)
              message(tempdir())
              #Creation directory for the results
              dirOutput_2 <- tempdir()
              currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
              dirOutput_1 <- paste("/", currentTime, "/", sep = "")
              dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
              dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
              message(dirOutput_Server)
              db_execution_interactn$dirOutput <- dirOutput_Server
              #Save folder for the download
              readr::write_csv(data.frame("session"=session$token,
                                          "outdir"=dirOutput_Server),
                               file = paste0(tempdir(),"/outdir_log_InteracTN.log"), append = T)
              
              
              #Read parameter and execution
              software <- input$sw_analyzer_interactn
              file_input_proteome = input$input_file_proteome_interactn$name
              file_prot_proteome = if(software%in%c("PD","MQ_prot")){input$prot_file_proteome_interactn$name}else{NA}
              file_pep_proteome = input$pep_file_proteome_interactn$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_protn"), showWarnings = FALSE)
              dir_input <- paste(dirOutput_Server, "input_protn", sep = "")
              file.copy(from = input$input_file_proteome_interactn$datapath, to = paste0(dir_input,'/ANNOTATION_',file_input_proteome)) 
              if(software%in%c("PD","MQ_prot")){file.copy(from = input$prot_file_proteome_interactn$datapath, to =paste0(dir_input,'/PROT_',file_prot_proteome))} 
              file.copy(from = input$pep_file_proteome_interactn$datapath, to = paste0(dir_input,'/PEP_',file_pep_proteome)) 
              
              # If advance filter
              if(input$advance_filter_interactn){
                NA_allow_condition <- input$NA_allow_condition_interactn
                min_peptide_protein <- input$min_peptide_protein_interactn
                impute_algorithm <- unlist(tstrsplit(input$impute_algorithm_interactn, "_"))
                if(input$sample_column_interactn == "Sample"){
                  sample_column <- input$sample_column_interactn
                } else{
                  if(software=="PD"){
                    sample_column <- "File Name"
                  } else{
                    sample_column <- "Sample"
                  }
                }
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
                impute_algorithm <- c("norm","phosr")
                if(software=="PD"){
                  sample_column <- "File Name"
                } else{
                  sample_column <- "Sample"
                }
              }
              
              # If to batch corrected read column
              if(input$batch_correction_interactn){
                batch_corr <- TRUE
                batch_correction_col <- input$batch_correction_col_interactn
              } else{
                batch_corr <- FALSE
                batch_correction_col <- "batch"
              }
              
              message(software)
              progress=0
              msg_read_function <- NULL
              withCallingHandlers(
                {
                  shinyjs::html("text", "")
                  if(software == "PD"){
                    db_execution_interactn$proteome_data <- read_proteomics(software = "PD",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_",
                                                                  proteinGroup_filename = "PROT_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ_ev"){
                    db_execution_interactn$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ_prot"){
                    db_execution_interactn$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  proteinGroup_filename = "PROT_", 
                                                                  sample_col = sample_column,
                                                                  use_proteinGroups_MQ = TRUE,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "SP"){
                    db_execution_interactn$proteome_data <- read_proteomics(software = "SP",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  sample_col = sample_column,
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "FP"){
                    db_execution_interactn$proteome_data <- read_proteomics(software = "FP",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  sample_col = sample_column,
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  }
                },
                message = function(m) {
                  msg_read_function <<- append(msg_read_function, conditionMessage(m))
                  # shinyjs::html(id = "messagge_read_phos_protn", html = paste0("<p>",m$message,"</p>"), add = TRUE)
                  progress=progress+0.05
                  setProgress(value = progress)
                }
              )
              
              write_lines(msg_read_function, file = paste0(db_execution_interactn$dirOutput,"log_filter_read_function.txt"))
              
              db_execution_interactn$data_loaded <- TRUE
              
              if(impute_algorithm[1] != "norm"){
                message("Doing before imputation")
                message(impute_algorithm[1])
                db_execution_interactn$imputed_data <- impute_intensity(proteome_data = db_execution_interactn$proteome_data, type = impute_algorithm[1])
                db_execution_interactn$normalized_data <- normalization_ProTN(proteome_data = db_execution_interactn$imputed_data)
              } else{
                message("Doing before normalization")
                message(impute_algorithm[2])
                db_execution_interactn$normalized_data <- normalization_ProTN(proteome_data = db_execution_interactn$proteome_data)
                db_execution_interactn$normalized_data <- impute_intensity(proteome_data = db_execution_interactn$normalized_data, type = impute_algorithm[2])
              }
              
              if(batch_corr){
                message("Executing batch correction...")
                db_execution_interactn$normalized_data <- batch_correction(proteome_data = db_execution_interactn$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              db_execution_interactn$parameter<-list("Imputation and normalization algorithm: " = ifelse(impute_algorithm[1] != "norm", impute_algorithm[1], impute_algorithm[2]), 
                                           "Sample column in annotation file: " = sample_column, 
                                           "Batch correction: " = ifelse(batch_corr, batch_correction_col, "FALSE"), 
                                           "N° missing value allow per condition: " = NA_allow_condition, 
                                           "Minimum peptide per protein: " = min_peptide_protein)
              
              output$c_anno_interactn <- DT::renderDT(db_execution_interactn$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome_interactn", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
                tags$h3("Statistics:"),
                tags$h4(paste0("Number of proteins: ", uniqueN(db_execution_interactn$normalized_data$dat_gene$GeneName))),
                tags$h4(paste0("Number of peptides: ", uniqueN(db_execution_interactn$normalized_data$dat_pep$ID_peptide))),
                tags$h3("Annotation table"),
                DT::DTOutput("c_anno_interactn")
              )
            })
          },
          error = function(e) {
            #Create error report and reactivate the click in the page
            showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
            html_text<-str_replace(read_file("R/error.html"), 
                                   pattern = "The page you’re looking for doesn’t exist.</p>", 
                                   replacement = paste0("Description:", e, "</p>"))
            write_file(html_text, file = paste0(tempdir(), "/error.html"))
            tags$iframe(src = "basedir/error.html", height = "100%", width = "100%", scrolling = "yes")
          }
        )
      })
      
    })
    
    output$render_abundance_plot_interactn <- renderUI({
      if (input$abundance_plot_interactn) {
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('abundance_plot_interactn')",
            plotOutput("small_abundance_plot_interactn")
          )
        )
      } else{
        db_execution_interactn$generate_abundance = NULL
      }
    })
    output$small_abundance_plot_interactn <- renderPlot({
      generate_abundance_interactn()
    })
    
    output$render_peptide_distribution_interactn <- renderUI({
      if (input$peptide_distribution_interactn) {
        tagList(
          tags$h3("N° peptides per proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('peptide_distribution_plot_interactn')",
            plotOutput("small_peptide_distribution_interactn")
          )
        )
      } else{
        db_execution_interactn$generate_peptide_distribution = NULL
      }
    })
    output$small_peptide_distribution_interactn <- renderPlot({
      generate_peptide_distribution_interactn()
    })
    
    
    output$render_raw_violin_interactn <- renderUI({
      if (input$raw_violin_interactn) {
        tagList(
          tags$h3("Distribution raw abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('raw_violin_plot_interactn')",
            plotOutput("small_raw_violin_interactn")
          )
        )
      } else{
        db_execution_interactn$raw_abundance_distribution = NULL
      }
    })
    output$small_raw_violin_interactn <- renderPlot({
      generate_raw_violin_interactn()
    })
    
    
    
    output$render_complexity_plot_interactn <- renderUI({
      if (input$complexity_plot_interactn) {
        tagList(
          tags$h3("Complexity plot of raw abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('complexity_plot_interactn')",
            plotOutput("small_complexity_plot_interactn")
          )
        )
      } else{
        db_execution_interactn$generate_complexity = NULL
      }
    })
    output$small_complexity_plot_interactn <- renderPlot({
      generate_complexity_interactn()
    })
    
    
    output$render_protein_violin_interactn <- renderUI({
      if (input$protein_violin_interactn) {
        tagList(
          tags$h3("Distribution protein abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('protein_violin_plot_interactn')",
            plotOutput("small_protein_violin_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_abundance_distribution = NULL
      }
    })
    output$small_protein_violin_interactn <- renderPlot({
      generate_protein_violin_interactn()
    })
    
    output$render_peptide_violin_interactn <- renderUI({
      if (input$peptide_violin_interactn) {
        tagList(
          tags$h3("Distribution peptide abundance"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('peptide_violin_plot_interactn')",
            plotOutput("small_peptide_violin_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_abundance_distirbution = NULL
      }
    })
    output$small_peptide_violin_interactn <- renderPlot({
      generate_peptide_violin_interactn()
    })
    
    output$render_mds_protein_interactn <- renderUI({
      if (input$mds_protein_interactn) {
        tagList(
          tags$h3("MDS based on proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('mds_protein_interactn')",
            plotOutput("small_mds_protein_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_MDS = NULL
      }
    })
    output$small_mds_protein_interactn <- renderPlot({
      generate_mds_protein_interactn()
    })
    
    output$render_mds_peptide_interactn <- renderUI({
      if (input$mds_peptide_interactn) {
        tagList(
          tags$h3("MDS based on peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('mds_peptide_interactn')",
            plotOutput("small_mds_peptide_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_MDS = NULL
      }
    })
    output$small_mds_peptide_interactn <- renderPlot({
      generate_mds_peptide_interactn()
    })
    
    output$render_pca_protein_interactn <- renderUI({
      if (input$pca_protein_interactn) {
        tagList(
          tags$h3("PCA based on proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('pca_protein_interactn')",
            plotOutput("small_pca_protein_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_PCA = NULL
      }
    })
    output$small_pca_protein_interactn <- renderPlot({
      generate_pca_protein_interactn()
    })
    
    output$render_pca_peptide_interactn <- renderUI({
      if (input$pca_peptide_interactn) {
        tagList(
          tags$h3("PCA based on peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('pca_peptide_interactn')",
            plotOutput("small_pca_peptide_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_PCA = NULL
      }
    })
    output$small_pca_peptide_interactn <- renderPlot({
      generate_pca_peptide_interactn()
    })
    
    output$render_protein_boxplot_interactn <- renderUI({
      if (input$boxplot_protein_interactn) {
        req(input$list_proteins_interactn)
        tagList(
          tags$h3("Boxplot selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('protein_boxplot_interactn')",
            plotOutput("small_protein_boxplot_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_boxplot = NULL
      }
    })
    output$small_protein_boxplot_interactn <- renderPlot({
      generate_protein_boxplot_interactn()
    })
    
    output$render_protein_heatmap_interactn <- renderUI({
      if (input$heatmap_protein_interactn) {
        req(input$list_proteins_interactn)
        tagList(
          tags$h3("Heatmap selected proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('protein_heatmap_interactn')",
            plotOutput("small_protein_heatmap_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_heatmap = NULL
      }
    })
    output$small_protein_heatmap_interactn <- renderPlot({
      generate_protein_heatmap_interactn()
    })
  })
  
  ## InteracTN: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn_interactn, {
    output$render_differential_analysis_interactn <- renderUI({
      isolate({
        updateCheckboxInput(session, "enrichment_analysis_interactn", value = FALSE)
        updateCheckboxInput(session, "stringdb_analysis_interactn", value = FALSE)
        
        db_execution_interactn$dt_formule_contrast <- as.data.table(hot_to_r(input$render_formule_contrast_table_interactn))
        db_execution_interactn$dt_formule_contrast <- db_execution_interactn$dt_formule_contrast[Formule!=""]
        print(db_execution_interactn$dt_formule_contrast)
        formule_diff <- as.list(db_execution_interactn$dt_formule_contrast$Formule)
        names(formule_diff) <- stri_replace_all(db_execution_interactn$dt_formule_contrast$Name, replacement = "_", regex = "-")
        
        names(formule_diff) <- lapply(1:length(formule_diff), function(x){
          if(names(formule_diff)[x] == ""){
            stri_replace_all(formule_diff[[x]], replacement = "_VS_", regex = "-")
          } else{
            names(formule_diff)[x]
          }
        })
        db_execution_interactn$formule_contrast <- formule_diff
        message(db_execution_interactn$formule_contrast)
        
        withProgress(message = "Differential analysis in process, please wait!", {
          message(session$token)
          message(tempdir())
          
          db_execution_interactn$differential_results <- differential_analysis(proteome_data = db_execution_interactn$normalized_data,
                                                                     formule_contrast = db_execution_interactn$formule_contrast,
                                                                     fc_thr=as.double(input$FC_thr_interactn),
                                                                     pval_fdr = input$pval_fdr_interactn,
                                                                     pval_thr=as.double(input$pval_thr_interactn),
                                                                     signal_thr=0, 
                                                                     interactomics = TRUE)
          db_execution_interactn$formule_contrast <- db_execution_interactn$formule_contrast[unique(union(db_execution_interactn$differential_results$protein_results_long$comp, 
                                                                                                          db_execution_interactn$differential_results$peptide_results_long$comp))]
          
          db_execution_interactn$parameter<-c(db_execution_interactn$parameter,
                                    "Fold change threshold for significance: "=input$FC_thr_interactn,
                                    "P.value type used: "=input$pval_fdr_interactn,
                                    "P.value threshold for significance: "=input$pval_thr_interactn)
          
        })
        
        tags$h2("Differential Analysis")
      })
    })
    
    output$render_protein_diff_table_interactn <- renderUI({
      if(input$protein_diff_table_interactn){
        output$protein_results_long_interactn <- DT::renderDT(db_execution_interactn$differential_results$protein_results_long)
        DT::DTOutput("protein_results_long_interactn")
      }
    })
    
    output$render_peptide_diff_table_interactn <- renderUI({
      if(input$peptide_diff_table_interactn){
        output$peptide_results_long_interactn <- DT::renderDT(db_execution_interactn$differential_results$peptide_results_long)
        DT::DTOutput("peptide_results_long_interactn")
      }
    })
    
    output$render_protein_diff_barplot_interactn <- renderUI({
      if (input$protein_diff_barplot_interactn) {
        tagList(
          tags$h3("N° differential proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('protein_diff_barplot_interactn')",
            plotOutput("small_protein_diff_barplot_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_differential_barplot = NULL
      }
    })
    output$small_protein_diff_barplot_interactn <- renderPlot({
      generate_protein_diff_barplot_interactn()(4)
    })
    
    output$render_peptide_diff_barplot_interactn <- renderUI({
      if (input$peptide_diff_barplot_interactn) {
        tagList(
          tags$h3("N° differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('peptide_diff_barplot_interactn')",
            plotOutput("small_peptide_diff_barplot_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_differential_barplot = NULL
      }
    })
    output$small_peptide_diff_barplot_interactn <- renderPlot({
      generate_peptide_diff_barplot_interactn()(4)
    })
    
    output$render_protein_upset_interactn <- renderUI({
      if (input$protein_upset_interactn) {
        tagList(
          tags$h3("Differential proteins upset plot"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('protein_upset_interactn')",
            plotOutput("small_protein_upset_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_upset_plot = NULL
      }
    })
    output$small_protein_upset_interactn <- renderPlot({
      generate_protein_upset_interactn()
    })
    
    output$render_peptide_upset_interactn <- renderUI({
      if (input$peptide_upset_interactn) {
        tagList(
          tags$h3("Differential peptides upset plot"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('peptide_upset_interactn')",
            plotOutput("small_peptide_upset_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_upset_plot = NULL
      }
    })
    output$small_peptide_upset_interactn <- renderPlot({
      generate_peptide_upset_interactn()
    })
    
    output$render_protein_ma_plot_interactn <- renderUI({
      if (input$protein_ma_plot_interactn) {
        c_anno <- db_execution_interactn$proteome_data$c_anno
        generate_ma_plots_protein <- list()
        for(comp in names(db_execution_interactn$formule_contrast)){
          message(comp)
          design <- model.matrix(~0 + c_anno$condition)
          colnames(design) <- levels(as.factor(c_anno$condition))
          rownames(design) <- c_anno$sample
          
          conds <- as.data.table(makeContrasts(contrasts = db_execution_interactn$formule_contrast[[comp]], levels = design), keep.rownames = T)
          conds <- conds[as.vector(conds[,2]!=0), rn]
          message(conds)
          
          generate_ma_plots_protein[[comp]] <- ma_plot(differential_results = db_execution_interactn$differential_results, 
                                                       proteome_data = db_execution_interactn$normalized_data,
                                                       type="protein", comparison = comp, condition = conds)$plot
        }
        db_execution_interactn$protein_ma_plot = lapply(generate_ma_plots_protein, function(x){ggplotly(x, tooltip = c("text"))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_ma_plots_protein)) {
          plot_id <- paste0(names(generate_ma_plots_protein)[i], "_ma_prot_interactn")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            # output[[my_plot_id]] <- renderPlot(generate_ma_plots_protein[[names(generate_ma_plots_protein)[my_i]]])
            output[[my_plot_id]] <- renderPlotly(ggplotly(generate_ma_plots_protein[[names(generate_ma_plots_protein)[my_i]]], tooltip = c("text")))
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_ma_plots_protein)[i]),
            # plotOutput(plot_id)
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("MA Plot differential proteins"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_ma_protein_interactn"), tabs))
        )
      } else{
        db_execution_interactn$protein_ma_plot = NULL
      }
    })
    
    output$render_peptide_ma_plot_interactn <- renderUI({
      if (input$peptide_ma_plot_interactn) {
        c_anno <- db_execution_interactn$proteome_data$c_anno
        generate_ma_plots_peptide <- list()
        for(comp in names(db_execution_interactn$formule_contrast)){
          message(comp)
          design <- model.matrix(~0 + c_anno$condition)
          colnames(design) <- levels(as.factor(c_anno$condition))
          rownames(design) <- c_anno$sample
          
          conds <- as.data.table(makeContrasts(contrasts = db_execution_interactn$formule_contrast[[comp]], levels = design), keep.rownames = T)
          conds <- conds[as.vector(conds[,2]!=0), rn]
          message(conds)
          
          generate_ma_plots_peptide[[comp]] <- ma_plot(differential_results = db_execution_interactn$differential_results, 
                                                       proteome_data = db_execution_interactn$normalized_data,
                                                       type="peptide", comparison = comp, condition = conds)$plot
        }
        db_execution_interactn$peptide_ma_plot = lapply(generate_ma_plots_peptide, function(x){ggplotly(x, tooltip = c("text"))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_ma_plots_peptide)) {
          plot_id <- paste0(names(generate_ma_plots_peptide)[i], "_ma_pep_interactn")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(ggplotly(generate_ma_plots_peptide[[names(generate_ma_plots_peptide)[my_i]]], tooltip = "text"))
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_ma_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("MA Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_ma_peptide_interactn"), tabs))
        )
      } else{
        db_execution_interactn$peptide_ma_plot = NULL
      }
    })
    
    output$render_protein_vulcano_interactn <- renderUI({
      if(input$protein_vulcano_interactn){
        generate_volcano_plots_protein <- list()
        for(comp in names(db_execution_interactn$formule_contrast)){
          generate_volcano_plots_protein<-c(generate_volcano_plots_protein,
                                            generate_volcano_plots(db_execution_interactn$differential_results,
                                                                   data_type="protein",
                                                                   comparison=comp,
                                                                   fc_thr=as.double(input$FC_thr_interactn),
                                                                   pval_fdr = input$pval_fdr_interactn,
                                                                   pval_thr=as.double(input$pval_thr_interactn), 
                                                                   interactomics = TRUE))
        }
        db_execution_interactn$protein_vulcano = generate_volcano_plots_protein
        # Generate tabPanels in a for loop

        tabs <- list()
        for (i in seq_along(generate_volcano_plots_protein)) {
          plot_id <- paste0(names(generate_volcano_plots_protein)[i], "_prot_interactn")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(generate_volcano_plots_protein[[names(generate_volcano_plots_protein)[my_i]]])
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_volcano_plots_protein)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential proteins"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_protein_interactn"), tabs))
          # renderPlotly(generate_volcano_plots_protein[[names(db_execution_interactn$formule_contrast)[[1]]]])
        )
        
      }else{
        db_execution_interactn$protein_vulcano = NULL
      }
    })
    
    output$render_peptide_vulcano_interactn <- renderUI({
      if(input$peptide_vulcano_interactn){
        generate_volcano_plots_peptide <- list()
        for(comp in names(db_execution_interactn$formule_contrast)){
          generate_volcano_plots_peptide<-c(generate_volcano_plots_peptide,
                                            generate_volcano_plots(db_execution_interactn$differential_results,
                                                                   data_type="peptide",
                                                                   comparison=comp,
                                                                   fc_thr=as.double(input$FC_thr_interactn),
                                                                   pval_fdr = input$pval_fdr_interactn,
                                                                   pval_thr=as.double(input$pval_thr_interactn), 
                                                                   interactomics = TRUE))
        }
        db_execution_interactn$peptide_vulcano = generate_volcano_plots_peptide
        # Generate tabPanels in a for loop
        tabs_pep_vulcano <- list()
        for (i in seq_along(generate_volcano_plots_peptide)) {
          plot_id <- paste0(names(generate_volcano_plots_peptide)[i], "_pep_interactn")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(generate_volcano_plots_peptide[[names(generate_volcano_plots_peptide)[my_i]]])
          })
          
          tabs_pep_vulcano[[i]] <- tabPanel(
            title = paste(names(generate_volcano_plots_peptide)[i]),
            plotlyOutput(plot_id)
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_peptide_interactn"), tabs_pep_vulcano))
        )
      }else{
        db_execution_interactn$peptide_vulcano = NULL
      }
    })
    
    output$render_mds_protein_diff_interactn <- renderUI({
      if (input$mds_diff_protein_interactn) {
        tagList(
          tags$h3("MDS based on differential proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('mds_protein_diff_interactn')",
            plotOutput("small_mds_protein_diff_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_differential_MDS = NULL
      }
    })
    output$small_mds_protein_diff_interactn <- renderPlot({
      generate_mds_protein_diff_interactn()
    })
    
    output$render_mds_peptide_diff_interactn <- renderUI({
      if (input$mds_diff_peptide_interactn) {
        tagList(
          tags$h3("MDS based on differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('mds_peptide_diff_interactn')",
            plotOutput("small_mds_peptide_diff_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_differential_MDS = NULL
      }
    })
    output$small_mds_peptide_diff_interactn <- renderPlot({
      generate_mds_peptide_diff_interactn()
    })
    
    output$render_pca_protein_diff_interactn <- renderUI({
      if (input$pca_diff_protein_interactn) {
        tagList(
          tags$h3("PCA based on differential proteins"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('pca_protein_diff_interactn')",
            plotOutput("small_pca_protein_diff_interactn")
          )
        )
      } else{
        db_execution_interactn$protein_differential_PCA = NULL
      }
    })
    output$small_pca_protein_diff_interactn <- renderPlot({
      generate_pca_protein_diff_interactn()
    })
    
    output$render_pca_peptide_diff_interactn <- renderUI({
      if (input$pca_diff_peptide_interactn) {
        tagList(
          tags$h3("PCA based on differential peptides"),
          tags$div(
            style = "cursor:pointer;",
            onclick = "showFullscreenPlot_interactn('pca_peptide_diff_interactn')",
            plotOutput("small_pca_peptide_diff_interactn")
          )
        )
      } else{
        db_execution_interactn$peptide_differential_PCA = NULL
      }
    })
    output$small_pca_peptide_diff_interactn <- renderPlot({
      generate_pca_peptide_diff_interactn()
    })
    
  })
  
  ## InteracTN: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn_interactn, {
    output$render_enrichement_analysis_interactn <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution_interactn$enrichment_results <- perform_enrichment_analysis(differential_results = db_execution_interactn$differential_results,
                                                                        enrichR_custom_DB = T,
                                                                        enrich_filter_DBs=input$DB_enrichment_interactn,    
                                                                        overlap_size_enrich_thr=as.double(input$os_enrich_interactn),
                                                                        pval_fdr_enrich = input$pval_fdr_interactn,
                                                                        pval_enrich_thr=as.double(input$pval_thr_interactn),
                                                                        dirOutput=db_execution_interactn$dirOutput, 
                                                                        with_background = input$enrich_with_background_interactn)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich_interactn), regex=","))
        
        db_execution_interactn$parameter <- c(db_execution_interactn$parameter,
                                    "Enrichment databases selected: "=paste(input$DB_enrichment_interactn, collapse = ", "),
                                    "P.value type used for enrichment: "=input$pval_fdr_interactn,
                                    "P.value threshold for enrichment significance: "=input$pval_thr_interactn,
                                    "Overlap size threshold for enrichment significance: "=input$os_enrich_interactn,
                                    "Enrichment filter terms: "=if(length(terms_enrich)>0){paste(terms_enrich, collapse = ", ")}else{"None"},
                                    "Enrichment with background: "=input$enrich_with_background_interactn)
        
        plots_down <- enrichment_figure(enr_df = db_execution_interactn$enrichment_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution_interactn$enrichment_results$anno_class), input$DB_enrichment_interactn))})
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(plots_down)) {
          plot_id <- names(plots_down)[i]
          height_id <- max(min(20, length(unique(plots_down[[names(plots_down)[i]]]$data$y_col))*0.4),3)*96
          message(paste0("Height for ",names(plots_down)[i], ": ", height_id))
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlot({
              plots_down[[names(plots_down)[my_i]]]
            }, height = height_id)
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(plots_down)[i]),
            plotOutput(plot_id, height = height_id)
          )
        }
        
        tagList(
          tags$h2("Enrichment Analysis"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_enrichment_interactn"), tabs))
        )
        
      })
    })
  })
  ## InteracTN: stringdb analysis ----
  observeEvent(input$execute_stringdb_analysis_btn_interactn, {
    output$render_stringdb_interactn <- renderUI({
      isolate({
        withProgress(message = "STRINGdb analysis in process, please wait!", {
          
          db_execution_interactn$stringdb_res <- STRINGdb_network(differential_results = db_execution_interactn$differential_results,
                                                        species=input$taxonomy_interactn, 
                                                        dirOutput=db_execution_interactn$dirOutput, 
                                                        score_thr=input$score_thr_stringdb_interactn,
                                                        shiny = T)
          db_execution_interactn$parameter <- c(db_execution_interactn$parameter,
                                      "STRINGdb taxonomy: "=input$taxonomy_interactn,
                                      "STRINGdb score threshold: "=input$score_thr_stringdb_interactn)
          
          tagList(
            tags$h2("STRINGdb analysis"),
            fluidRow(
              selectInput("stringdb_show_interactn", label = "Select StringDB to show: (click on STRING logo to open the results on stringDB website)", 
                          choices = names(db_execution_interactn$stringdb_res), width = "15%"),
              actionButton("stringdb_selected_interactn", "Select!", width = "10%")  
            ),
            tags$div(id = "stringEmbedded")
          )
        })
      })
    })
  })
  
  observeEvent(input$stringdb_selected_interactn, {
    js$loadStringData(input$taxonomy_interactn, db_execution_interactn$stringdb_res[[input$stringdb_show_interactn]], input$score_thr_stringdb_interactn)
  })
  # InteracTN: download results ----
  output$download_proteome_interactn <- downloadHandler(
    filename = "results_InteracTN.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            message(db_execution_interactn$dirOutput)
            setProgress(value = 0.01)
            
            # Generate report
            params <- list(
              doc_title = input$title_exp_interactn,
              description = input$description_exp_interactn,
              readPD_files = if (input$sw_analyzer_interactn == "PD") {TRUE} else {FALSE},
              readMQ_files = if (input$sw_analyzer_interactn == "MQ") {TRUE} else {FALSE},
              impute_algorithm = if(input$advance_filter_interactn){input$impute_algorithm_interactn} else {"norm_phosr"},
              db_execution = reactiveValuesToList(db_execution_interactn),
              file_input = paste(db_execution_interactn$dirOutput, "input_protn", sep = ""),
              batch_corr_exe = if(input$batch_correction_interactn){input$batch_correction_col_interactn}else{NULL},
              prot_boxplot = if(input$boxplot_protein_interactn | input$heatmap_protein_interactn){input$list_proteins_interactn}else{NULL},
              fc_thr = if(is.null(input$FC_thr_interactn)){"0.75"}else{input$FC_thr_interactn},
              pval_fdr = input$pval_fdr_interactn,
              pval_thr = if(is.null(input$pval_thr_interactn)){"0.05"}else{input$pval_thr_interactn},
              pval_fdr_enrich = input$pval_fdr_enrich_interactn,
              pval_enrich_thr = if(is.null(input$pvalue_enrich_interactn)){"0.05"}else{input$pvalue_enrich_interactn},
              overlap_size_enrich_thr = if(is.null(input$os_enrich_interactn)){as.integer(5)}else{input$os_enrich_interactn},
              enrich_filter_term = input$terms_enrich_interactn,
              enrich_filter_DBs = input$DB_enrichment_interactn,
              taxonomy=input$taxonomy_interactn, 
              score_thr=input$score_thr_stringdb_interactn,
              dirOutput = db_execution_interactn$dirOutput
            )
            
            # Render in background the report
            p = callr::r_bg(
              func = function(db_execution_interactn, params, dirOutput, env) {
                rmarkdown::render("R/interactn_report.Rmd",
                                  output_file = "interactn_report.html",
                                  output_dir = dirOutput,
                                  params = params,
                                  envir = env
                )
              },
              args = list(db_execution_interactn, params, db_execution_interactn$dirOutput, new.env(parent = globalenv())),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            
            
            # Saving RData in background
            setProgress(value = 0.95, message = "Saving RData...")
            db_results_interacTN = reactiveValuesToList(db_execution_interactn)
            db_results_interacTN <- db_results_interacTN[!(unlist(lapply(db_results_interacTN, is.null)))]
            p_rdata = callr::r_bg(
              func = function(db_results_interacTN, dirOutput) {
                save(db_results_interacTN, file = paste0(dirOutput,"db_results_InteracTN.RData"))
              },
              args = list(db_results_interacTN, db_results_interacTN$dirOutput),
              stdout = "|",
              stderr = "|",
              error = getOption("callr.error", "error")
            )
            
            # Prepare file for the download
            if(length(db_execution_interactn$normalized_data)>0){
              save_abundance_tables(proteome_data = db_execution_interactn$normalized_data, 
                                    dirOutput = db_execution_interactn$dirOutput)
            }
            setProgress(value = 0.1)
            
            if(length(db_execution_interactn$differential_results)>0){
              save_differential_analysis_table(proteome_data = db_execution_interactn$normalized_data,
                                               differential_results = db_execution_interactn$differential_results,
                                               dirOutput=db_execution_interactn$dirOutput)
            }
            setProgress(value = 0.2)
            
            if(input$abundance_plot_interactn & !is.null(db_execution_interactn$generate_abundance)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/missing_available_abundance.pdf"), 
                     plot = db_execution_interactn$generate_abundance, 
                     create.dir = T, width = 7, height = 5)
            } else if("missing_available_abundance.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/missing_available_abundance.pdf"))
            }
            setProgress(value = 0.25)
            
            if(input$raw_violin_interactn & !is.null(db_execution_interactn$raw_abundance_distribution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/raw_abundance_distribution.pdf"), 
                     plot = db_execution_interactn$raw_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("raw_abundance_distribution.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/raw_abundance_distribution.pdf"))
            }
            
            
            if(input$peptide_distribution_interactn & !is.null(db_execution_interactn$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_per_protein.pdf"), 
                     plot = db_execution_interactn$generate_peptide_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_per_protein.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_per_protein.pdf"))
            }
            setProgress(value = 0.30)
            
            if(input$complexity_plot_interactn & !is.null(db_execution_interactn$generate_complexity)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/complexity_plot.pdf"), 
                     plot = db_execution_interactn$generate_complexity, 
                     create.dir = T, width = 10, height = 8)
            } else if("complexity_plot.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/complexity_plot.pdf"))
            }
            setProgress(value = 0.33)
            
            if(input$protein_violin_interactn & !is.null(db_execution_interactn$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution_interactn$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_abundance_distribution.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_abundance_distribution.pdf"))
            }
            setProgress(value = 0.35)
            
            if(input$peptide_violin_interactn & !is.null(db_execution_interactn$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution_interactn$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_abundance_distribution.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_abundance_distribution.pdf"))
            }
            setProgress(value = 0.40)
            
            if(input$mds_protein_interactn & !is.null(db_execution_interactn$protein_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution_interactn$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_MDS.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_MDS.pdf"))
            }
            setProgress(value = 0.43)
            
            if(input$mds_peptide_interactn & !is.null(db_execution_interactn$peptide_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution_interactn$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_MDS.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_MDS.pdf"))
            }
            setProgress(value = 0.45)
            
            if(input$pca_protein_interactn & !is.null(db_execution_interactn$protein_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution_interactn$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_PCA.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_PCA.pdf"))
            }
            setProgress(value = 0.47)
            
            if(input$pca_peptide_interactn & !is.null(db_execution_interactn$peptide_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution_interactn$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_PCA.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_PCA.pdf"))
            }
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein_interactn & !is.null(db_execution_interactn$protein_boxplot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution_interactn$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_boxplot.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_boxplot.pdf"))
            }
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein_interactn & !is.null(db_execution_interactn$protein_heatmap)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution_interactn$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } else if("protein_heatmap.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_heatmap.pdf"))
            }
            setProgress(value = 0.55)
            
            if(!is.null(db_execution_interactn$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution_interactn$protein_differential_barplot, 
                     create.dir = T, width = 17, height = 6)
            } else if("protein_differential_barplot.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_differential_barplot.pdf"))
            }
            setProgress(value = 0.58)
            
            if(!is.null(db_execution_interactn$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution_interactn$peptide_differential_barplot, 
                     create.dir = T, width = 17, height = 6)
            } else if("peptide_differential_barplot.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_differential_barplot.pdf"))
            }
            setProgress(value = 0.60)
            
            if(!is.null(db_execution_interactn$protein_upset_plot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_upset_plot.pdf"), 
                     plot = db_execution_interactn$protein_upset_plot, 
                     create.dir = T,width = 12, height = 6)
            } else if("protein_upset_plot.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_upset_plot.pdf"))
            }
            setProgress(value = 0.62)
            
            if(!is.null(db_execution_interactn$peptide_upset_plot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_upset_plot.pdf"), 
                     plot = db_execution_interactn$peptide_upset_plot, 
                     create.dir = T, width = 12, height = 6)
            } else if("peptide_upset_plot.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_upset_plot.pdf"))
            }
            setProgress(value = 0.63)
            
            if(!is.null(db_execution_interactn$protein_ma_plot)){
              dir.create(file.path(paste0(db_execution_interactn$dirOutput,"pics/"), "protein_ma_plot"), showWarnings = FALSE)
              for(comp in names(db_execution_interactn$protein_ma_plot)){
                htmlwidgets::saveWidget(db_execution_interactn$protein_ma_plot[[comp]], 
                                        file = paste0(db_execution_interactn$dirOutput,"pics/protein_ma_plot/",comp,"_protein_ma_plot.html"))
                webshot2::webshot(url = paste0(db_execution_interactn$dirOutput,"pics/protein_ma_plot/",comp,"_protein_ma_plot.html"), 
                                  file = paste0(db_execution_interactn$dirOutput,"pics/protein_ma_plot/",comp,"_protein_ma_plot.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_interactn$dirOutput,"pics/protein_ma_plot"))
            }
            setProgress(value = 0.64)
            
            
            if(!is.null(db_execution_interactn$peptide_ma_plot)){
              dir.create(file.path(paste0(db_execution_interactn$dirOutput,"pics/"), "peptide_ma_plot"), showWarnings = FALSE)
              for(comp in names(db_execution_interactn$peptide_ma_plot)){
                htmlwidgets::saveWidget(db_execution_interactn$peptide_ma_plot[[comp]], 
                                        file = paste0(db_execution_interactn$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"))
                webshot2::webshot(url = paste0(db_execution_interactn$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.html"), 
                                  file = paste0(db_execution_interactn$dirOutput,"pics/peptide_ma_plot/",comp,"_peptide_ma_plot.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_interactn$dirOutput,"pics/peptide_ma_plot"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_interactn$protein_vulcano)){
              dir.create(file.path(paste0(db_execution_interactn$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_interactn$protein_vulcano)){
                # plotly::save_image(db_execution_interactn$protein_vulcano[[comp]], 
                #                    file = paste0(str_replace_all(db_execution_interactn$dirOutput, pattern="\\\\", replacement="/"),"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_interactn$protein_vulcano[[comp]], 
                                        file = paste0(db_execution_interactn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution_interactn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"), 
                                  file = paste0(db_execution_interactn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_interactn$dirOutput,"pics/protein_vulcano"))
            }
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_interactn$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution_interactn$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_interactn$peptide_vulcano)){
                # plotly::save_image(db_execution_interactn$peptide_vulcano[[comp]], 
                #                    file = paste0(str_replace_all(db_execution_interactn$dirOutput, pattern="\\\\", replacement="/"),"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_interactn$peptide_vulcano[[comp]], 
                                        file = paste0(db_execution_interactn$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"))
                webshot2::webshot(url = paste0(db_execution_interactn$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.html"), 
                                  file = paste0(db_execution_interactn$dirOutput,"pics/peptide_vulcano/",comp,"_peptide_vulcano.png"), delay = 1, zoom = 4)
              }
            } else{
              message("Removing old rendered plot")
              system(paste0("rm -r ",db_execution_interactn$dirOutput,"pics/peptide_vulcano"))
            }
            setProgress(value = 0.68)
            
            if(!is.null(db_execution_interactn$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution_interactn$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_differential_MDS.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_differential_MDS.pdf"))
            }
            setProgress(value = 0.69)
            
            if(!is.null(db_execution_interactn$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution_interactn$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_MDS.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_differential_MDS.pdf"))
            }
            setProgress(value = 0.70)
            
            if(!is.null(db_execution_interactn$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution_interactn$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("protein_differential_PCA.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/protein_differential_PCA.pdf"))
            }
            setProgress(value = 0.72)
            
            if(!is.null(db_execution_interactn$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution_interactn$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } else if("peptide_differential_PCA.pdf" %in% list.files(paste0(db_execution_interactn$dirOutput,"pics"))){
              message("Removing old rendered plot")
              system(paste0("rm ",db_execution_interactn$dirOutput,"pics/peptide_differential_PCA.pdf"))
            }
            setProgress(value = 0.75)
            
            if(length(db_execution_interactn$enrichment_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich_interactn), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution_interactn$enrichment_results,
                                              category = c("down","up"), 
                                              enrich_filter_term = terms_enrich,
                                              save=T, 
                                              dirOutput = db_execution_interactn$dirOutput)
            } 
            setProgress(value = 0.82)
            
            if(length(db_execution_interactn$stringdb_res)>0){
              tmp_res <- STRINGdb_network(differential_results = db_execution_interactn$differential_results,
                                          species=input$taxonomy_interactn, 
                                          dirOutput=db_execution_interactn$dirOutput,
                                          score_thr=input$score_thr_stringdb_interactn,
                                          shiny = F)
              
            } 
            setProgress(value = 0.95)
            
            # Write tsv file with parameter
            params <- data.table("Parameter" = names(db_execution_interactn$parameter),
                                 "Value" = unlist(db_execution_interactn$parameter))
            fwrite(params, paste0(db_execution_interactn$dirOutput,"parameters_used.txt"), sep = "\t", col.names = F)
            
            
            #Get results Report
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p$read_output()
            p$wait(30000)
            for (i in 1:15) {
              p$read_output()
              p$wait(1000*60)  
            }
            if(p$is_alive() | is.null(p$get_result())){
              p$kill()
              print("\n ERROR: An error occur during the report rendering. \n ")
            } else{
              report<-p$get_result()
              p$kill()
              message("Render report DONE.")
            }
            
            #Wait 10 minutes. If do not end in 10 minutes, kill the process
            hide_res<-p_rdata$read_output()
            p_rdata$wait(30000)
            for (i in 1:15) {
              p_rdata$read_output()
              p_rdata$wait(1000*60)
            }
            if(p_rdata$is_alive() | is.null(p_rdata$get_result())){
              p_rdata$kill()
              print("\n ERROR: An error occur during the RData saving. \n ")
            } else{
              report<-p_rdata$get_result()
              p_rdata$kill()
              message("RData saved.")
            }

            # Save RData db_execution_interactn
            
            #Save folder for the download
            oldwd <- getwd()
            message(db_execution_interactn$dirOutput)
            setwd(db_execution_interactn$dirOutput)
            files2zip <- list.files("./", recursive = TRUE)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  ## InteracTN: full screen trigger ----
  
  # ReactiveVal for currently selected plot to fullscreen
  selected_plot_interactn <- reactiveVal(NULL)
  
  # Update selected_plot when JS sends fullscreen_trigger id
  observeEvent(input$fullscreen_trigger_interactn, {
    selected_plot_interactn(input$fullscreen_trigger_interactn)
  })
  
  # Render fullscreen plot dynamically based on selected_plot()
  output$fullscreen_plot_interactn <- renderPlot({
    req(selected_plot_interactn())
    switch(selected_plot_interactn(),
           "abundance_plot_interactn" = generate_abundance_interactn() + ggtitle("Percentage missing values respect detected abundance")+theme(text=element_text(size=25)),
           "peptide_distribution_plot_interactn" = generate_peptide_distribution_interactn() + ggtitle("N° peptides per proteins")+theme(text=element_text(size=25)),
           "raw_violin_plot_interactn" = generate_raw_violin_interactn() + ggtitle("Raw abundance distribution")+theme(text=element_text(size=25)),
           "complexity_plot_interactn" = generate_complexity_interactn() + ggtitle("Complexity plot of raw abundance")+theme(text=element_text(size=25)),
           "protein_violin_plot_interactn" = generate_protein_violin_interactn() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "peptide_violin_plot_interactn" = generate_peptide_violin_interactn() + ggtitle("Distribution peptide abundance")+theme(text=element_text(size=25)),
           "mds_protein_interactn" = generate_mds_protein_interactn() + ggtitle("MDS based on protein")+theme(text=element_text(size=25)),
           "mds_peptide_interactn" = generate_mds_peptide_interactn() + ggtitle("MDS based on peptides")+theme(text=element_text(size=25)),
           "pca_protein_interactn" = generate_pca_protein_interactn() + ggtitle("PCA based on protein")+theme(text=element_text(size=25)),
           "pca_peptide_interactn" = generate_pca_peptide_interactn() + ggtitle("PCA based on peptides")+theme(text=element_text(size=25)),
           "protein_boxplot_interactn" = generate_protein_boxplot_interactn() + ggtitle("Boxplot selected proteins")+theme(text=element_text(size=25)),
           "protein_heatmap_interactn" = generate_protein_heatmap_interactn() + ggtitle("Heatmap selected proteins")+theme(text=element_text(size=25)),
           "protein_diff_barplot_interactn" = generate_protein_diff_barplot_interactn()(8) + ggtitle("N° differential proteins")+theme(text=element_text(size=25)),
           "peptide_diff_barplot_interactn" = generate_peptide_diff_barplot_interactn()(8) + ggtitle("N° differential peptides")+theme(text=element_text(size=25)),
           "protein_upset_interactn" = generate_protein_upset_interactn(),
           "peptide_upset_interactn" = generate_peptide_upset_interactn(),
           "mds_protein_diff_interactn" = generate_mds_protein_diff_interactn() + ggtitle("MDS based on differential protein")+theme(text=element_text(size=25)),
           "mds_peptide_diff_interactn" = generate_mds_peptide_diff_interactn() + ggtitle("MDS based on differential peptides")+theme(text=element_text(size=25)),
           "pca_protein_diff_interactn" = generate_pca_protein_diff_interactn() + ggtitle("PCA based on differential protein")+theme(text=element_text(size=25)),
           "pca_peptide_diff_interactn" = generate_pca_peptide_diff_interactn() + ggtitle("PCA based on differential peptides")+theme(text=element_text(size=25)),
           # default fallback:
           NULL
    )
  })
  
  ##############################################################################
  # Download case study Proteome ----
  output$download_CS_proteome <- downloadHandler(
    filename = "case_study_proteomics.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            setProgress(value = 0.01)
            
            dirOutput_2 <- tempdir()
            currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
            dirOutput_1 <- paste("/", currentTime, "/", sep = "")
            dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
            dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
            message(dirOutput_Server)
            setProgress(value = 0.1)
            
            extract_example(path_proteome = dirOutput_Server)
            setProgress(value = 0.5)
            #Save folder for the download
            oldwd <- getwd()
            message(dirOutput_Server)
            setwd(dirOutput_Server)
            files2zip <- list.files("./", recursive = TRUE)
            setProgress(value = 0.9)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  # Download case study Phospho ----
  output$download_CS_phos <- downloadHandler(
    filename = "case_study_phospho.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            setProgress(value = 0.01)
            
            dirOutput_2 <- tempdir()
            currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
            dirOutput_1 <- paste("/", currentTime, "/", sep = "")
            dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
            dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
            message(dirOutput_Server)
            setProgress(value = 0.1)
            
            extract_example(path_phospho = dirOutput_Server)
            setProgress(value = 0.5)
            #Save folder for the download
            oldwd <- getwd()
            message(dirOutput_Server)
            setwd(dirOutput_Server)
            files2zip <- list.files("./", recursive = TRUE)
            setProgress(value = 0.9)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  # Download case study Phospho proteome ----
  output$download_CS_phos_protn <- downloadHandler(
    filename = "case_study.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            setProgress(value = 0.01)
            
            dirOutput_2 <- tempdir()
            currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
            dirOutput_1 <- paste("/", currentTime, "/", sep = "")
            dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
            dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
            dir.create(file.path(dirOutput_Server, "/proteome/"), showWarnings = FALSE)
            dirOutput_Server_protn <- paste(dirOutput_Server, "/proteome/", sep = "")
            dir.create(file.path(dirOutput_Server, "/phospho/"), showWarnings = FALSE)
            dirOutput_Server_phos <- paste(dirOutput_Server, "/phospho/", sep = "")
            message(dirOutput_Server)
            setProgress(value = 0.1)
            
            extract_example(path_proteome = dirOutput_Server_protn, path_phospho = dirOutput_Server_phos)
            setProgress(value = 0.5)
            #Save folder for the download
            oldwd <- getwd()
            message(dirOutput_Server)
            setwd(dirOutput_Server)
            files2zip <- list.files("./", recursive = TRUE)
            setProgress(value = 0.9)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  # Download case study Proteome ----
  output$download_CS_interactn <- downloadHandler(
    filename = "case_study_interactn.zip",
    content = function(file) {
      tryCatch(
        {
          withProgress(message = "Preparing files to download, please wait!", {
            #Zip the dir resutls
            message(session$token)
            setProgress(value = 0.01)
            
            dirOutput_2 <- tempdir()
            currentTime <- gsub(".*?([0-9]+).*?", "\\1", Sys.time())
            dirOutput_1 <- paste("/", currentTime, "/", sep = "")
            dir.create(file.path(dirOutput_2, dirOutput_1), showWarnings = FALSE)
            dirOutput_Server <- paste(dirOutput_2, dirOutput_1, sep = "")
            message(dirOutput_Server)
            setProgress(value = 0.1)
            
            extract_example(path_interactn = dirOutput_Server)
            setProgress(value = 0.5)
            #Save folder for the download
            oldwd <- getwd()
            message(dirOutput_Server)
            setwd(dirOutput_Server)
            files2zip <- list.files("./", recursive = TRUE)
            setProgress(value = 0.9)
            zip(zipfile = file, files = files2zip, extra = "-r")
            setwd(oldwd)
            
          })
        },
        error = function(e) {
          #Create error report and reactivate the click in the page
          showNotification(paste0("ERROR: ", e), type = "error", duration = 30)
          html_text<-str_replace(read_file("R/error.html"), 
                                 pattern = "The page you’re looking for doesn’t exist.</p>", 
                                 replacement = paste0("Description:", e, "</p>"))
          write_file(html_text, file = paste0(tempdir(), "/error.html"))
          zip(zipfile = file, files = paste0(tempdir(), "/error.html"), extra = "-j")
        }
      )
    }
  )
  
  
  
  ##############################################################################
  # -- DELETE TEMP FILES WHEN SESSION ENDS ----
  session$onSessionEnded(function() {
    message(session$token)
    if (file.exists(paste0(tempdir(),"/outdir_log_ProTN.log"))){
      log_clean <- fread(paste0(tempdir(),"/outdir_log_ProTN.log"), header = FALSE)
      dir_to_remove <- log_clean[V1==session$token, V2]
      unlink(dir_to_remove, recursive = T)
    }
    
    if (file.exists(paste0(tempdir(),"/outdir_log_PhosProTN.log"))){
      log_clean <- fread(paste0(tempdir(),"/outdir_log_PhosProTN.log"), header = FALSE)
      dir_to_remove <- log_clean[V1==session$token, V2]
      unlink(dir_to_remove, recursive = T)
    }
    
    if (file.exists(paste0(tempdir(),"/outdir_log_PhosProTN_proteome.log"))){
      log_clean <- fread(paste0(tempdir(),"/outdir_log_PhosProTN_proteome.log"), header = FALSE)
      dir_to_remove <- log_clean[V1==session$token, V2]
      unlink(dir_to_remove, recursive = T)
    }
    
    if (file.exists(paste0(tempdir(),"/outdir_log_InteracTN.log"))){
      log_clean <- fread(paste0(tempdir(),"/outdir_log_InteracTN.log"), header = FALSE)
      dir_to_remove <- log_clean[V1==session$token, V2]
      unlink(dir_to_remove, recursive = T)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 8100))
