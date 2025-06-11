# ProTN v0.1.2: an integrative pipeline for complete analysis of proteomics    # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# PI: Dr. Toma Tebaldi, PhD                                                    #
#                                                                              #
list.of.packages <- c("shiny","tidyverse","markdown","knitr","shinydashboard",
                      "shinydashboardPlus","shinymaterial","shinyjs","magrittr",
                      "dplyr","stringr","shinyBS","DT","bslib","readr",
                      "plotly","rhandsontable","shinyalert","ggplot2")
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
library(proTN)
library(data.table)
library(rhandsontable)

tmpdir<-stri_replace_all_regex(tempdir(), pattern = "/Rtmp\\w+", replacement = "")
if (!dir.exists(file.path(tmpdir, "ProTN_shiny"))){
  dir.create(file.path(tmpdir, "ProTN_shiny"), showWarnings = FALSE)
}
Sys.setenv(TMPDIR=file.path(tmpdir, "ProTN_shiny"))
unlink(tempdir(), recursive = T)
tempdir(check = T)
#Javascript code to disable click when is running the app
# jsCode <- "shinyjs.pageDisable = function(params){
#               $('body').css('pointer-events', params);
#             };"

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
      extendShinyjs(text = jsCode_STRINGdb, functions = c("loadStringData")),
      # extendShinyjs(text = jsCode, functions = c("pageDisable")),
      
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
                  textInput("title_exp", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp", "Brief description", rows = 4),
                ),
                fluidRow(
                  selectInput("sw_analyzer", "Software Analyzer:",
                              choice = list("ProteomeDiscoverer" = "PD", 
                                            "MaxQuant by evidence.txt" = "MQ_ev", 
                                            "MaxQuant by peptides.txt and proteinGroups.txt" = "MQ_prot"),
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
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_abundance_plot")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_distribution")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_violin")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_violin")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide")
                  )
                ),
                uiOutput("render_protein_boxplot"),
                uiOutput("render_protein_heatmap"),
                uiOutput("render_differential_analysis"),
                uiOutput("render_protein_diff_table"),
                uiOutput("render_peptide_diff_table"),
                uiOutput("render_protein_diff_barplot"),
                uiOutput("render_peptide_diff_barplot"),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_vulcano")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_vulcano")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_diff")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_diff")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_diff")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_diff")
                  )
                ),
                uiOutput("render_enrichement_analysis"),
                uiOutput("render_stringdb")
              )
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
                  textInput("title_exp_phos", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp_phos", "Brief description", rows = 4),
                ),
                fluidRow(
                  radioButtons("sw_analyzer_phos", "Software Analyzer", 
                               choiceNames = c("ProteomeDiscoverer", "MaxQuant", "TMT_PD"),
                               choiceValues = c("PD","MQ","TMT_PD"), inline = TRUE),
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
                checkboxInput("protein_violin_phos", "Distribution abundance proteins", FALSE),
                checkboxInput("peptide_violin_phos", "Distribution abundance peptides", FALSE),
                checkboxInput("mds_protein_phos", "MDS based on protein", FALSE),
                checkboxInput("mds_peptide_phos", "MDS based on peptide", FALSE),
                checkboxInput("pca_protein_phos", "PCA based on protein", FALSE),
                checkboxInput("pca_peptide_phos", "PCA based on peptide", TRUE),
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
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_abundance_plot_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_distribution_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_violin_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_violin_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_phos")
                  )
                ),
                uiOutput("render_protein_boxplot_phos"),
                uiOutput("render_protein_heatmap_phos"),
                uiOutput("render_differential_analysis_phos"),
                uiOutput("render_protein_diff_table_phos"),
                uiOutput("render_peptide_diff_table_phos"),
                uiOutput("render_protein_diff_barplot_phos"),
                uiOutput("render_peptide_diff_barplot_phos"),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_vulcano_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_vulcano_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_diff_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_diff_phos")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_diff_phos")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_diff_phos")
                  )
                ),
                uiOutput("render_enrichement_analysis_phos"),
                uiOutput("render_stringdb_phos"),
                uiOutput("render_kinase_tree_phos")
              )
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
                  textInput("title_exp_phos_protn", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp_phos_protn", "Brief description", rows = 4),
                ),
                fluidRow(
                  radioButtons("sw_analyzer_phos_protn", "Software Analyzer", 
                               choiceNames = c("ProteomeDiscoverer", "MaxQuant", "TMT_PD"),
                               choiceValues = c("PD","MQ","TMT_PD"), inline = TRUE),
                ),
                uiOutput("input_proteome_phos_protn"),
                checkboxInput("batch_correction_phos_protn", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui_phos_protn"),
                sliderInput("phos_thr", "Phosphorylation threshold", 0, 100, step = 5, value = 75),
                checkboxInput("advance_filter_phos_protn", "Advance Filter", FALSE),
                uiOutput("advance_filter_ui_phos_protn"),
                actionButton("report_proteome_phos_protn", "Load data!"),
                tags$h3("Select what execute:"),
                checkboxInput("phospho_percentage_plot_phos_protn", "% phosphorylated site", TRUE),
                checkboxInput("abundance_plot_phos_protn", "% missing values", TRUE),
                checkboxInput("peptide_distribution_phos_protn", "N° peptides per protein", TRUE),
                checkboxInput("protein_violin_phos_protn", "Distribution abundance proteins", FALSE),
                checkboxInput("peptide_violin_phos_protn", "Distribution abundance peptides", FALSE),
                checkboxInput("mds_protein_phos_protn", "MDS based on protein", FALSE),
                checkboxInput("mds_peptide_phos_protn", "MDS based on peptide", FALSE),
                checkboxInput("pca_protein_phos_protn", "PCA based on protein", FALSE),
                checkboxInput("pca_peptide_phos_protn", "PCA based on peptide", TRUE),
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
                    width = 6,
                    uiOutput("render_protein_violin_phos_protn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_violin_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_phos_protn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_phos_protn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_phos_protn")
                  )
                ),
                uiOutput("render_protein_boxplot_phos_protn"),
                uiOutput("render_protein_heatmap_phos_protn"),
                uiOutput("render_differential_analysis_phos_protn"),
                uiOutput("render_protein_diff_table_phos_protn"),
                uiOutput("render_peptide_diff_table_phos_protn"),
                uiOutput("render_protein_diff_barplot_phos_protn"),
                uiOutput("render_peptide_diff_barplot_phos_protn"),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_vulcano_phos_protn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_vulcano_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_diff_phos_protn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_diff_phos_protn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_diff_phos_protn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_diff_phos_protn")
                  )
                ),
                uiOutput("render_enrichement_analysis_phos_protn"),
                uiOutput("render_stringdb_phos_protn"),
                uiOutput("render_kinase_tree_phos_protn")
              )
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
                  textInput("title_exp_interactn", "Title of the analysis"),
                ),
                fluidRow(
                  textAreaInput("description_exp_interactn", "Brief description", rows = 4),
                ),
                fluidRow(
                  radioButtons("sw_analyzer_interactn", "Software Analyzer", 
                               choiceNames = c("ProteomeDiscoverer", "MaxQuant", "TMT_PD"),
                               choiceValues = c("PD","MQ","TMT_PD"), inline = TRUE),
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
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_abundance_plot_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_distribution_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_violin_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_violin_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_interactn")
                  )
                ),
                uiOutput("render_protein_boxplot_interactn"),
                uiOutput("render_protein_heatmap_interactn"),
                uiOutput("render_differential_analysis_interactn"),
                uiOutput("render_protein_diff_table_interactn"),
                uiOutput("render_peptide_diff_table_interactn"),
                uiOutput("render_protein_diff_barplot_interactn"),
                uiOutput("render_peptide_diff_barplot_interactn"),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_protein_vulcano_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_peptide_vulcano_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_mds_protein_diff_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_mds_peptide_diff_interactn")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    uiOutput("render_pca_protein_diff_interactn")
                  ),
                  column(
                    width = 6,
                    uiOutput("render_pca_peptide_diff_interactn")
                  )
                ),
                uiOutput("render_enrichement_analysis_interactn"),
                uiOutput("render_stringdb_interactn")
              )
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
                                 data_loaded = FALSE,
                                 dirOutput = "",
                                 proteome_data = list(),
                                 imputed_data = list(),
                                 normalized_data = list(),
                                 formule_contrast = list(),
                                 dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                 differential_results = list(),
                                 enrichmnent_results = list(),
                                 stringdb_res = list(),
                                 kinase_tree_res = list(),
                                 phospho_percentage = NULL,
                                 generate_abundance = NULL,
                                 generate_peptide_distribution = NULL,
                                 protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                 protein_MDS = NULL, peptide_MDS = NULL,
                                 protein_PCA = NULL, peptide_PCA = NULL,
                                 protein_boxplot = NULL, protein_heatmap = NULL,
                                 protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                 protein_vulcano = NULL, peptide_vulcano = NULL,
                                 protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                 protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  # PHOSPROTN: db variable analysis ----
  db_execution_phos <- reactiveValues(session = session$token,
                                      data_loaded = FALSE,
                                      dirOutput = "",
                                      proteome_data = list(),
                                      imputed_data = list(),
                                      normalized_data = list(),
                                      formule_contrast = list(),
                                      dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                      differential_results = list(),
                                      enrichmnent_results = list(),
                                      stringdb_res = list(),
                                      kinase_tree_res = list(),
                                      phospho_percentage = NULL,
                                      generate_abundance = NULL,
                                      generate_peptide_distribution = NULL,
                                      protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                      protein_MDS = NULL, peptide_MDS = NULL,
                                      protein_PCA = NULL, peptide_PCA = NULL,
                                      protein_boxplot = NULL, protein_heatmap = NULL,
                                      protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                      protein_vulcano = NULL, peptide_vulcano = NULL,
                                      protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                      protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  # PhosProTN_with_prot: db variable analysis ----
  db_execution_phos_protn <- reactiveValues(session = session$token,
                                            data_loaded = FALSE,
                                            dirOutput = "",
                                            proteome_data = list(),
                                            imputed_data = list(),
                                            normalized_data = list(),
                                            formule_contrast = list(),
                                            dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                            differential_results = list(),
                                            enrichmnent_results = list(),
                                            stringdb_res = list(),
                                            kinase_tree_res = list(),
                                            phospho_percentage = NULL,
                                            generate_abundance = NULL,
                                            generate_peptide_distribution = NULL,
                                            protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                            protein_MDS = NULL, peptide_MDS = NULL,
                                            protein_PCA = NULL, peptide_PCA = NULL,
                                            protein_boxplot = NULL, protein_heatmap = NULL,
                                            protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                            protein_vulcano = NULL, peptide_vulcano = NULL,
                                            protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                            protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
  
  # InteracTN: db variable analysis ----
  db_execution_interactn <- reactiveValues(session = session$token,
                                 data_loaded = FALSE,
                                 dirOutput = "",
                                 proteome_data = list(),
                                 imputed_data = list(),
                                 normalized_data = list(),
                                 formule_contrast = list(),
                                 dt_formule_contrast = data.table("Name"=c("","","",""),"Formule"=c("","","","")),
                                 differential_results = list(),
                                 enrichmnent_results = list(),
                                 stringdb_res = list(),
                                 kinase_tree_res = list(),
                                 phospho_percentage = NULL,
                                 generate_abundance = NULL,
                                 generate_peptide_distribution = NULL,
                                 protein_abundance_distribution = NULL, peptide_abundance_distirbution = NULL,
                                 protein_MDS = NULL, peptide_MDS = NULL,
                                 protein_PCA = NULL, peptide_PCA = NULL,
                                 protein_boxplot = NULL, protein_heatmap = NULL,
                                 protein_differential_barplot = NULL, peptide_differential_barplot = NULL,
                                 protein_vulcano = NULL, peptide_vulcano = NULL,
                                 protein_differential_MDS = NULL, peptide_differential_MDS = NULL,
                                 protein_differential_PCA = NULL, peptide_differential_PCA = NULL)
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
        numericInput("min_peptide_protein", "Minimum peptide per protein", value = 1, min = 1)
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
        # textAreaInput("formule_contrast", "Write in each line a different comparison", rows = 4),
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
        checkboxInput("protein_vulcano", "Proteins vulcano plot", FALSE),
        checkboxInput("peptide_vulcano", "Peptides vulcano plot", FALSE),
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
    } 
  })
  
  output$render_formule_contrast_table <- renderRHandsontable({
    rhandsontable(db_execution$dt_formule_contrast, rowHeaders = NULL, stretchH = "all")
  })
  
  ## PROTN: show enrichment parameter ----
  output$enrichment_params_ui <- renderUI({ 
    if(input$enrichment_analysis){
      tagList(
        # radioButtons("enrichR_universe", "Execute enrichment of the whole Universe", c("TRUE", "FALSE"), inline = TRUE, selected = FALSE),
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
    }
  })
  
  # PROTN: Execution pipeline ----
  observeEvent(input$report_proteome, {
    
    output$protn_results_ui <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
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
              file_prot_proteome = if(software!="MQ_ev"){input$prot_file_proteome$name}else{NA}
              file_pep_proteome = input$pep_file_proteome$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_protn"), showWarnings = FALSE)
              dir_input <- paste(dirOutput_Server, "input_protn", sep = "")
              file.copy(from = input$input_file_proteome$datapath, to = paste0(dir_input,'/ANNOTATION_',file_input_proteome)) 
              if(software!="MQ_ev"){file.copy(from = input$prot_file_proteome$datapath, to =paste0(dir_input,'/PROT_',file_prot_proteome))} 
              file.copy(from = input$pep_file_proteome$datapath, to = paste0(dir_input,'/PEP_',file_pep_proteome)) 
              
              # If advance filter
              if(input$advance_filter){
                NA_allow_condition <- input$NA_allow_condition
                min_peptide_protein <- input$min_peptide_protein
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
              }
              
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
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ_ev"){
                    db_execution$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
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
                                                                  use_proteinGroups_MQ = TRUE,
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
              
              write_lines(msg_read_function, file = paste0(db_execution$dirOutput,"log_filter_read_function.txt"))
              
              db_execution$data_loaded <- TRUE
              db_execution$imputed_data <- impute_intensity(proteome_data = db_execution$proteome_data)
              db_execution$normalized_data <- normalization_ProTN(proteome_data = db_execution$imputed_data)
              if(batch_corr){
                message("Executing batch correction...")
                db_execution$normalized_data <- batch_correction(proteome_data = db_execution$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              output$c_anno <- DT::renderDT(db_execution$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
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
      if(input$abundance_plot){
        generate_abundance <- generate_abundance_plot(proteome_data = db_execution$proteome_data)
        db_execution$generate_abundance = generate_abundance$plot
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          renderPlot(generate_abundance$plot)
        )
      } else{
        db_execution$generate_abundance = NULL
      }
    })
    
    output$render_peptide_distribution <- renderUI({ 
      if(input$peptide_distribution){
        generate_peptide_distribution <- generate_peptide_distribution_plot(proteome_data = db_execution$proteome_data)
        db_execution$generate_peptide_distribution = generate_peptide_distribution$plot
          tagList(
          tags$h3("N° peptides per proteins"),
          renderPlot(generate_peptide_distribution$plot)
        )
      } else{
        db_execution$generate_peptide_distribution = NULL
      }
    })
    
    output$render_protein_violin <- renderUI({ 
      if(input$protein_violin){
        generate_protein_violin <- plot_abundance_distribution(proteome_data = db_execution$normalized_data,
                                                               type = "protein")
        db_execution$protein_abundance_distribution = generate_protein_violin$plot
        tagList(
          tags$h3("Distribution protein abundance"),
          renderPlot(generate_protein_violin$plot)
        )
      } else{
        db_execution$protein_abundance_distribution = NULL
      }
    })
    
    output$render_peptide_violin <- renderUI({ 
      if(input$peptide_violin){
        generate_peptide_violin <- plot_abundance_distribution(proteome_data = db_execution$normalized_data,
                                                               type = "peptide")
        db_execution$peptide_abundance_distirbution = generate_peptide_violin$plot
        tagList(
          tags$h3("Distribution peptide abundace"),
          renderPlot(generate_peptide_violin$plot)
        )
      } else{
        db_execution$peptide_abundance_distirbution = NULL
      }
    })
    
    output$render_mds_protein <- renderUI({ 
      if(input$mds_protein){
        res_plot <- mds_plot(proteome_data = db_execution$normalized_data,
                             type = "protein")
        db_execution$protein_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution$protein_MDS = NULL
      }
    })
    
    output$render_mds_peptide <- renderUI({ 
      if(input$mds_peptide){
        res_plot <- mds_plot(proteome_data = db_execution$normalized_data,
                             type = "peptide")
        db_execution$peptide_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on peptides"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution$peptide_MDS = NULL
      }
    })
    
    output$render_pca_protein <- renderUI({ 
      if(input$pca_protein){
        res_plot <- pca_plot(proteome_data = db_execution$normalized_data,
                             type = "protein")
        db_execution$protein_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution$protein_PCA = NULL
      }
    })
    
    output$render_pca_peptide <- renderUI({ 
      if(input$pca_peptide){
        res_plot <- pca_plot(proteome_data = db_execution$normalized_data,
                             type = "peptide")
        db_execution$peptide_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on peptides"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution$peptide_PCA = NULL
      }
    })
    
    output$render_protein_boxplot <- renderUI({
      if(input$boxplot_protein){
        req(input$list_proteins)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins), regex=",")
        prot_boxplot <- plot_selected_proteins(proteome_data = db_execution$normalized_data,
                                               list_protein = unlist(list_proteins))
        db_execution$protein_boxplot = prot_boxplot$plot
        
        tagList(
          tags$h3("Boxplot selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      }else{
        db_execution$protein_boxplot = NULL
      }
    })
    
    output$render_protein_heatmap <- renderUI({
      if(input$heatmap_protein){
        req(input$list_proteins)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins), regex=",")
        prot_boxplot <- heatmap_selected_proteins(proteome_data = db_execution$normalized_data,
                                                  list_protein = unlist(list_proteins))
        db_execution$protein_heatmap = prot_boxplot$plot
        
        tagList(
          tags$h3("Heatmap selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      }else{
        db_execution$protein_heatmap = NULL
      }
    })
  })
  
  ## PROTN: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn, {
    output$render_differential_analysis <- renderUI({
      isolate({
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
        })

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
      if(input$protein_diff_barplot){
        ploft_diff_number <- generate_differential_barplots(db_execution$differential_results,
                                                            data_type="protein")
        db_execution$protein_differential_barplot = ploft_diff_number$plot
        tagList(
          tags$h3("N° differential proteins"),
          renderPlot(ploft_diff_number$plot)
        )
      }else{
        db_execution$protein_differential_barplot = NULL
      }
    })
    
    output$render_peptide_diff_barplot <- renderUI({
      if(input$peptide_diff_barplot){
        ploft_diff_number_pep <- generate_differential_barplots(db_execution$differential_results,
                                                                data_type="peptide")
        db_execution$peptide_differential_barplot = ploft_diff_number_pep$plot
        tagList(
          tags$h3("N° differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution$peptide_differential_barplot = NULL
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
          # renderPlotly(generate_volcano_plots_protein[[names(db_execution$formule_contrast)[[1]]]])
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
      if(input$mds_diff_protein){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "protein")
        db_execution$protein_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution$protein_differential_MDS = NULL
      }
    })
    
    output$render_mds_peptide_diff <- renderUI({
      if(input$mds_diff_peptide){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "peptide")
        db_execution$peptide_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution$peptide_differential_MDS = NULL
      }
    })
    
    output$render_pca_protein_diff <- renderUI({
      if(input$pca_diff_protein){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "protein")
        db_execution$protein_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution$protein_differential_PCA = NULL
      }
    })
    
    output$render_pca_peptide_diff <- renderUI({
      if(input$pca_diff_peptide){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "peptide")
        db_execution$peptide_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution$peptide_differential_PCA = NULL
      }
    })
    
  })
  
  ## PROTN: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn, {
    output$render_enrichement_analysis <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution$enrichmnent_results <- perform_enrichment_analysis(differential_results = db_execution$differential_results,
                                                          enrichR_custom_DB = T,
                                                          enrich_filter_DBs=input$DB_enrichment,    
                                                          overlap_size_enrich_thr=as.double(input$FC_thr),
                                                          pval_fdr_enrich = input$pval_fdr,
                                                          pval_enrich_thr=as.double(input$pval_thr),
                                                          dirOutput=db_execution$dirOutput, 
                                                          with_background = input$enrich_with_background)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich), regex=","))
        plots_down <- enrichment_figure(enr_df = db_execution$enrichmnent_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution$enrichmnent_results$anno_class), input$DB_enrichment))})
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
          withProgress(message = "Prepraring files to download, please wait!", {
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
            } 
            setProgress(value = 0.25)
            
            if(input$peptide_distribution & !is.null(db_execution$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_per_protein.pdf"), 
                     plot = db_execution$generate_peptide_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.30)
            
            if(input$protein_violin & !is.null(db_execution$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.35)
            
            if(input$peptide_violin & !is.null(db_execution$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.40)
            
            if(input$mds_protein & !is.null(db_execution$protein_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.43)
            
            if(input$mds_peptide & !is.null(db_execution$peptide_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.45)
            
            if(input$pca_protein & !is.null(db_execution$protein_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.47)
            
            if(input$pca_peptide & !is.null(db_execution$peptide_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein & !is.null(db_execution$protein_boxplot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein & !is.null(db_execution$protein_heatmap)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.55)
            
            if(!is.null(db_execution$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution$protein_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.58)
            
            if(!is.null(db_execution$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution$peptide_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.60)
            
            if(!is.null(db_execution$protein_vulcano)){
              dir.create(file.path(paste0(db_execution$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution$protein_vulcano)){
                plotly::save_image(db_execution$protein_vulcano[[comp]], 
                                        file = paste0(db_execution$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution$protein_vulcano[[comp]], 
                                   file = paste0(db_execution$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.64)
            
            if(!is.null(db_execution$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution$peptide_vulcano)){
                plotly::save_image(db_execution$peptide_vulcano[[comp]], 
                             file = paste0(db_execution$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution$peptide_vulcano[[comp]], 
                                   file = paste0(db_execution$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.68)
            
            if(!is.null(db_execution$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.69)
            
            if(!is.null(db_execution$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.70)
            
            if(!is.null(db_execution$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.72)
            
            if(!is.null(db_execution$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.75)
            
            if(length(db_execution$enrichmnent_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution$enrichmnent_results,
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
            
            # Save RData db_execution
            db_results_proTN = reactiveValuesToList(db_execution)
            db_results_proTN <- db_results_proTN[!(unlist(lapply(db_results_proTN, is.null)))]
            save(db_results_proTN, file = paste0(db_results_proTN$dirOutput,"db_results_proTN.RData"))
            
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
        numericInput("min_peptide_protein_phos", "Minimum peptide per protein", value = 1, min = 1)
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
        # textAreaInput("formule_contrast", "Write in each line a different comparison", rows = 4),
        textInput("FC_thr_phos", "Fold change threshold for significance:",value = 0.5),
        radioButtons("pval_fdr_phos", "Select which p.value use:", 
                     choiceNames = c("Adj.P.Val", "P.Val"),
                     choiceValues = c("p_adj","p_val"), inline = TRUE, selected = "p_val"),
        textInput("pval_thr_phos", "P.value threshold for significance:", value = 0.05),
        actionButton("execute_differential_analysis_btn_phos", "Run!"),
        checkboxInput("protein_diff_table_phos", "Proteins differentiated table", FALSE),
        checkboxInput("peptide_diff_table_phos", "Peptides differentiated table", FALSE),
        checkboxInput("protein_diff_barplot_phos", "Proteins differentiated barplot", TRUE),
        checkboxInput("peptide_diff_barplot_phos", "Peptides differentiated barplot", FALSE),
        checkboxInput("protein_vulcano_phos", "Proteins vulcano plot", FALSE),
        checkboxInput("peptide_vulcano_phos", "Peptides vulcano plot", FALSE),
        checkboxInput("mds_diff_protein_phos", "MDS based on diffential protein", FALSE),
        checkboxInput("mds_diff_peptide_phos", "MDS based on diffential peptide", FALSE),
        checkboxInput("pca_diff_protein_phos", "PCA based on diffential protein", FALSE),
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
    }
  })
  
  # PHOSPROTN: Execution pipeline ----
  observeEvent(input$report_proteome_phos, {
    
    output$protn_results_ui_phos <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
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
                               file = paste0(tempdir(),"/outdir_log_ProTN.log"), append = T)
              
              
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
              if(input$advance_filter){
                NA_allow_condition <- input$NA_allow_condition_phos
                min_peptide_protein <- input$min_peptide_protein_phos
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
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
              db_execution_phos$imputed_data <- impute_intensity(proteome_data = db_execution_phos$proteome_data)
              db_execution_phos$normalized_data <- normalization_ProTN(proteome_data = db_execution_phos$imputed_data)
              if(batch_corr){
                message("Executing batch correction...")
                db_execution_phos$normalized_data <- batch_correction(proteome_data = db_execution_phos$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              output$c_anno_phos <- DT::renderDT(db_execution_phos$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome_phos", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
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
      if(input$phospho_percentage_plot_phos){
        phospho_percentage <- create_phosphosite_plot(proteome_data = db_execution_phos$proteome_data, 
                                                      software = input$sw_analyzer_phos)
        db_execution_phos$phospho_percentage = phospho_percentage$plot
        tagList(
          tags$h3("Percentage of phosphosite residue"),
          renderPlot(phospho_percentage$plot, height = 250)
        )
      } else{
        db_execution_phos$phospho_percentage = NULL
      } 
    })
    
    output$render_abundance_plot_phos <- renderUI({ 
      if(input$abundance_plot_phos){
        generate_abundance <- generate_abundance_plot(proteome_data = db_execution_phos$proteome_data)
        db_execution_phos$generate_abundance = generate_abundance$plot
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          renderPlot(generate_abundance$plot)
        )
      } else{
        db_execution_phos$generate_abundance = NULL
      } 
    })
    
    output$render_peptide_distribution_phos <- renderUI({ 
      if(input$peptide_distribution_phos){
        generate_peptide_distribution <- generate_peptide_distribution_plot(proteome_data = db_execution_phos$proteome_data)
        db_execution_phos$generate_peptide_distribution = generate_peptide_distribution$plot
        tagList(
          tags$h3("N° peptides per proteins"),
          renderPlot(generate_peptide_distribution$plot)
        )
      } else{
        db_execution_phos$generate_peptide_distribution = NULL
      } 
    })
    
    output$render_protein_violin_phos <- renderUI({ 
      if(input$protein_violin_phos){
        generate_protein_violin <- plot_abundance_distribution(proteome_data = db_execution_phos$normalized_data,
                                                               type = "protein")
        db_execution_phos$protein_abundance_distribution = generate_protein_violin$plot
        tagList(
          tags$h3("Distribution protein abundance"),
          renderPlot(generate_protein_violin$plot)
        )
      } else{
        db_execution_phos$protein_abundance_distribution = NULL
      } 
    })
    
    output$render_peptide_violin_phos <- renderUI({ 
      if(input$peptide_violin_phos){
        generate_peptide_violin <- plot_abundance_distribution(proteome_data = db_execution_phos$normalized_data,
                                                               type = "peptide")
        db_execution_phos$peptide_abundance_distirbution = generate_peptide_violin$plot
        tagList(
          tags$h3("Distribution peptide abundace"),
          renderPlot(generate_peptide_violin$plot)
        )
      } else{
        db_execution_phos$peptide_abundance_distirbution = NULL
      } 
    })
    
    output$render_mds_protein_phos <- renderUI({ 
      if(input$mds_protein_phos){
        res_plot <- mds_plot(proteome_data = db_execution_phos$normalized_data,
                             type = "protein")
        db_execution_phos$protein_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_phos$protein_MDS = NULL
      } 
    })
    
    output$render_mds_peptide_phos <- renderUI({ 
      if(input$mds_peptide_phos){
        res_plot <- mds_plot(proteome_data = db_execution_phos$normalized_data,
                             type = "peptide")
        db_execution_phos$peptide_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on peptides"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_phos$peptide_MDS = NULL
      } 
    })
    
    output$render_pca_protein_phos <- renderUI({ 
      if(input$pca_protein_phos){
        res_plot <- pca_plot(proteome_data = db_execution_phos$normalized_data,
                             type = "protein")
        db_execution_phos$protein_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_phos$protein_PCA = NULL
      } 
    })
    
    output$render_pca_peptide_phos <- renderUI({ 
      if(input$pca_peptide_phos){
        res_plot <- pca_plot(proteome_data = db_execution_phos$normalized_data,
                             type = "peptide")
        db_execution_phos$peptide_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on peptides"),
          renderPlot(res_plot$plot)
        )
      }  else{
        db_execution_phos$peptide_PCA = NULL
      }
    })
    
    output$render_protein_boxplot_phos <- renderUI({
      if(input$boxplot_protein_phos){
        req(input$list_proteins_phos)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos), regex=",")
        prot_boxplot <- plot_selected_proteins(proteome_data = db_execution_phos$normalized_data,
                                               list_protein = unlist(list_proteins))
        db_execution_phos$protein_boxplot = prot_boxplot$plot
        
        tagList(
          tags$h3("Boxplot selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      } else{
        db_execution_phos$protein_boxplot = NULL
      }
    })
    
    output$render_protein_heatmap_phos <- renderUI({
      if(input$heatmap_protein_phos){
        req(input$list_proteins_phos)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos), regex=",")
        prot_boxplot <- heatmap_selected_proteins(proteome_data = db_execution_phos$normalized_data,
                                                  list_protein = unlist(list_proteins))
        db_execution_phos$protein_heatmap = prot_boxplot$plot
        
        tagList(
          tags$h3("Heatmap selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      } else{
        db_execution_phos$protein_heatmap = NULL
      }
    })
  })
  
  ## PHOSPROTN: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn_phos, {
    output$render_differential_analysis_phos <- renderUI({
      isolate({
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
        })
        
        tags$h2("Differential Analysis")
      })
    })
    
    output$render_protein_diff_table_phos <- renderUI({
      if(input$protein_diff_table_phos){
        output$protein_results_long_phos <- DT::renderDT(db_execution_phos$differential_results$protein_results_long)
        DT::DTOutput("protein_results_long_phos")
      }
    })
    
    output$render_peptide_diff_table_phos <- renderUI({
      if(input$peptide_diff_table_phos){
        output$peptide_results_long_phos <- DT::renderDT(db_execution_phos$differential_results$peptide_results_long)
        DT::DTOutput("peptide_results_long_phos")
      }
    })
    
    output$render_protein_diff_barplot_phos <- renderUI({
      if(input$protein_diff_barplot_phos){
        ploft_diff_number <- generate_differential_barplots(db_execution_phos$differential_results,
                                                            data_type="protein")
        db_execution_phos$protein_differential_barplot = ploft_diff_number$plot
        tagList(
          tags$h3("N° differential proteins"),
          renderPlot(ploft_diff_number$plot)
        )
      } else{
        db_execution_phos$protein_differential_barplot = NULL
      }
    })
    
    output$render_peptide_diff_barplot_phos <- renderUI({
      if(input$peptide_diff_barplot_phos){
        ploft_diff_number_pep <- generate_differential_barplots(db_execution_phos$differential_results,
                                                                data_type="peptide")
        db_execution_phos$peptide_differential_barplot = ploft_diff_number_pep$plot
        tagList(
          tags$h3("N° differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos$peptide_differential_barplot = NULL
      }
    })
    
    output$render_protein_vulcano_phos <- renderUI({
      if(input$protein_vulcano_phos){
        generate_volcano_plots_protein <- list()
        for(comp in names(db_execution_phos$formule_contrast)){
          generate_volcano_plots_protein<-c(generate_volcano_plots_protein,
                                            generate_volcano_plots(db_execution_phos$differential_results,
                                                                   data_type="protein",
                                                                   comparison=comp,
                                                                   fc_thr=as.double(input$FC_thr_phos),
                                                                   pval_fdr = input$pval_fdr_phos,
                                                                   pval_thr=as.double(input$pval_thr_phos)))
        }
        db_execution_phos$protein_vulcano = generate_volcano_plots_protein
        # Generate tabPanels in a for loop
        tabs <- list()
        for (i in seq_along(generate_volcano_plots_protein)) {
          plot_id <- paste0(names(generate_volcano_plots_protein)[i], "_prot_phos")
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlotly(generate_volcano_plots_protein[[names(generate_volcano_plots_protein)[my_i]]])
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(generate_volcano_plots_protein)[i]),
            plotlyOutput(plot_id, width = "99%")
          )
        }
        
        # Use do.call to unpack the tab list into tabsetPanel
        tagList(
          tags$h3("Vulcano Plot differential proteins"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_protein_phos"), tabs))
          # renderPlotly(generate_volcano_plots_protein[[names(db_execution_phos$formule_contrast)[[1]]]])
        )
        
      } else{
        db_execution_phos$protein_vulcano = NULL
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
    
    output$render_mds_protein_diff_phos <- renderUI({
      if(input$mds_diff_protein_phos){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
                                                                proteome_data = db_execution_phos$normalized_data,
                                                                type = "protein")
        db_execution_phos$protein_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos$protein_differential_MDS = NULL
      }
    })
    
    output$render_mds_peptide_diff_phos <- renderUI({
      if(input$mds_diff_peptide_phos){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
                                                                proteome_data = db_execution_phos$normalized_data,
                                                                type = "peptide")
        db_execution_phos$peptide_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos$peptide_differential_MDS = NULL
      }
    })
    
    output$render_pca_protein_diff_phos <- renderUI({
      if(input$pca_diff_protein_phos){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
                                                                proteome_data = db_execution_phos$normalized_data,
                                                                type = "protein")
        db_execution_phos$protein_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos$protein_differential_PCA = NULL
      }
    })
    
    output$render_pca_peptide_diff_phos <- renderUI({
      if(input$pca_diff_peptide_phos){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution_phos$differential_results,
                                                                proteome_data = db_execution_phos$normalized_data,
                                                                type = "peptide")
        db_execution_phos$peptide_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos$peptide_differential_PCA = NULL
      }
    })
    
  })
  
  ## PHOSPROTN: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn_phos, {
    output$render_enrichement_analysis_phos <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution_phos$enrichmnent_results <- perform_enrichment_analysis(differential_results = db_execution_phos$differential_results,
                                                                        enrichR_custom_DB = T,
                                                                        enrich_filter_DBs=input$DB_enrichment_phos,    
                                                                        overlap_size_enrich_thr=as.double(input$FC_thr_phos),
                                                                        pval_fdr_enrich = input$pval_fdr_phos,
                                                                        pval_enrich_thr=as.double(input$pval_thr_phos),
                                                                        dirOutput=db_execution_phos$dirOutput, 
                                                                        with_background = input$enrich_with_background_phos)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich_phos), regex=","))
        plots_down <- enrichment_figure(enr_df = db_execution_phos$enrichmnent_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution_phos$enrichmnent_results$anno_class), input$DB_enrichment_phos))})
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
          withProgress(message = "Prepraring files to download, please wait!", {
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
            
            
            
            # Prepare file for the download
            if(length(db_execution_phos$normalized_data)>0){
              save_abundance_tables(proteome_data = db_execution_phos$normalized_data, 
                                    dirOutput = db_execution_phos$dirOutput)
            }
            setProgress(value = 0.1)
            
            if(length(db_execution_phos$differential_results)>0){
              save_differential_analysis_table(proteome_data = db_execution_phos$normalized_data,
                                               differential_results = db_execution_phos$differential_results,
                                               dirOutput=db_execution_phos$dirOutput)
            }
            setProgress(value = 0.2)
            
            if(input$phospho_percentage_plot_phos & !is.null(db_execution_phos$phospho_percentage)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/phospho_percentage.pdf"), 
                     plot = db_execution_phos$phospho_percentage, 
                     create.dir = T, width = 7, height = 3)
            } 
            
            if(input$abundance_plot_phos & !is.null(db_execution_phos$generate_abundance)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/missing_available_abundance.pdf"), 
                     plot = db_execution_phos$generate_abundance, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.25)
            
            if(input$peptide_distribution_phos & !is.null(db_execution_phos$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_per_protein.pdf"), 
                     plot = db_execution_phos$generate_peptide_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.30)
            
            if(input$protein_violin_phos & !is.null(db_execution_phos$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution_phos$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.35)
            
            if(input$peptide_violin_phos & !is.null(db_execution_phos$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution_phos$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.40)
            
            if(input$mds_protein_phos & !is.null(db_execution_phos$protein_MDS)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution_phos$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.43)
            
            if(input$mds_peptide_phos & !is.null(db_execution_phos$peptide_MDS)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution_phos$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.45)
            
            if(input$pca_protein_phos & !is.null(db_execution_phos$protein_PCA)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution_phos$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.47)
            
            if(input$pca_peptide_phos & !is.null(db_execution_phos$peptide_PCA)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution_phos$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein_phos & !is.null(db_execution_phos$protein_boxplot)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution_phos$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein_phos & !is.null(db_execution_phos$protein_heatmap)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution_phos$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.55)
            
            if(!is.null(db_execution_phos$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution_phos$protein_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.58)
            
            if(!is.null(db_execution_phos$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution_phos$peptide_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.60)
            
            if(!is.null(db_execution_phos$protein_vulcano)){
              dir.create(file.path(paste0(db_execution_phos$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos$protein_vulcano)){
                plotly::save_image(db_execution_phos$protein_vulcano[[comp]], 
                                   file = paste0(db_execution_phos$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos$protein_vulcano[[comp]], 
                                        file = paste0(db_execution_phos$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_phos$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution_phos$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos$peptide_vulcano)){
                plotly::save_image(db_execution_phos$peptide_vulcano[[comp]], 
                                   file = paste0(db_execution_phos$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos$peptide_vulcano[[comp]], 
                                        file = paste0(db_execution_phos$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.68)
            
            if(!is.null(db_execution_phos$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution_phos$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.69)
            
            if(!is.null(db_execution_phos$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution_phos$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.70)
            
            if(!is.null(db_execution_phos$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution_phos$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.72)
            
            if(!is.null(db_execution_phos$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution_phos$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.75)
            
            if(length(db_execution_phos$enrichmnent_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich_phos), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution_phos$enrichmnent_results,
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
            
            # Save RData db_execution_phos
            db_results_PhosProTN = reactiveValuesToList(db_execution_phos)
            db_results_PhosProTN <- db_results_PhosProTN[!(unlist(lapply(db_results_PhosProTN, is.null)))]
            save(db_results_PhosProTN, file = paste0(db_results_PhosProTN$dirOutput,"db_results_PhosProTN.RData"))
            
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
  
  ## PHOSPROTN: modal for zoom image ----
  
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
          fileInput("psm_file_phospho_phos_protn", "Select the PROT file of the PHOSPHO-PROTEOMICS..."),
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
        numericInput("min_peptide_protein_phos_protn", "Minimum peptide per protein", value = 1, min = 1)
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
      }
    }
  })
  
  # PhosProTN_with_prot: Execution pipeline ----
  observeEvent(input$report_proteome_phos_protn, {
    
    output$protn_results_ui_phos_protn <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
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
              if(input$advance_filter){
                NA_allow_condition <- input$NA_allow_condition_phos_protn
                min_peptide_protein <- input$min_peptide_protein_phos_protn
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
              }
              
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
                                                                                   psm_proteome_filename = "PSM_",
                                                                                   annotation_phospho_filename = "ANNOTATION_",
                                                                                   proteinGroup_phospho_filename = "PROT_", 
                                                                                   psm_phospho_filename = "PSM_",
                                                                                   batch_corr_exe = batch_corr, 
                                                                                   batch_col = batch_correction_col,
                                                                                   phospho_thr = input$phos_thr/100, 
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
              
              write_lines(msg_read_function, file = paste0(db_execution_phos_protn$dirOutput,"log_filter_read_function.txt"))
              
              db_execution_phos_protn$data_loaded <- TRUE
              db_execution_phos_protn$imputed_data <- impute_intensity(proteome_data = db_execution_phos_protn$proteome_data)
              db_execution_phos_protn$normalized_data <- normalization_ProTN(proteome_data = db_execution_phos_protn$imputed_data)
              if(batch_corr){
                message("Executing batch correction...")
                db_execution_phos_protn$normalized_data <- batch_correction(proteome_data = db_execution_phos_protn$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              output$c_anno_phos_protn <- DT::renderDT(db_execution_phos_protn$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome_phos_protn", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
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
      if(input$phospho_percentage_plot_phos_protn){
        phospho_percentage <- create_phosphosite_plot(proteome_data = db_execution_phos_protn$proteome_data, 
                                                      software = input$sw_analyzer_phos_protn)
        db_execution_phos_protn$phospho_percentage = phospho_percentage$plot
        tagList(
          tags$h3("Percentage of phosphosite residue"),
          renderPlot(phospho_percentage$plot, height = 250)
        )
      } else{
        db_execution_phos_protn$phospho_percentage = NULL
      } 
    })
    
    output$render_abundance_plot_phos_protn <- renderUI({ 
      if(input$abundance_plot_phos_protn){
        generate_abundance <- generate_abundance_plot(proteome_data = db_execution_phos_protn$proteome_data, 
                                                      phospho_with_proteome = TRUE)
        db_execution_phos_protn$generate_abundance = generate_abundance
        tagList(
          fluidRow(
            column(
              width = 6,
              tags$h3("Percentage missing values respect detected abundance - Proteomics"),
              renderPlot(generate_abundance$proteome_plot)
            ),
            column(
              width = 6,
              tags$h3("Percentage missing values respect detected abundance - Phospho-proteomics"),
              renderPlot(generate_abundance$phospho_plot)
            )
          )
        )
      } else{
        db_execution_phos_protn$generate_abundance = NULL
      } 
    })
    
    output$render_peptide_distribution_phos_protn <- renderUI({ 
      if(input$peptide_distribution_phos_protn){
        generate_peptide_distribution <- generate_peptide_distribution_plot(proteome_data = db_execution_phos_protn$proteome_data, 
                                                                            phospho_with_proteome = TRUE)
        db_execution_phos_protn$generate_peptide_distribution = generate_peptide_distribution
        tagList(
          fluidRow(
            column(
              width = 6,
              tags$h3("N° peptides per proteins - Proteomics"),
              renderPlot(generate_peptide_distribution$proteome_plot)
            ),
            column(
              width = 6,
              tags$h3("N° peptides per proteins - Phospho-proteomics"),
              renderPlot(generate_peptide_distribution$phospho_plot)
            )
          )
        )
      } else{
        db_execution_phos_protn$generate_peptide_distribution = NULL
      } 
    })
    
    output$render_protein_violin_phos_protn <- renderUI({ 
      if(input$protein_violin_phos_protn){
        generate_protein_violin <- plot_abundance_distribution(proteome_data = db_execution_phos_protn$normalized_data,
                                                               type = "protein")
        db_execution_phos_protn$protein_abundance_distribution = generate_protein_violin$plot
        tagList(
          tags$h3("Distribution protein abundance"),
          renderPlot(generate_protein_violin$plot)
        )
      } else{
        db_execution_phos_protn$protein_abundance_distribution = NULL
      } 
    })
    
    output$render_peptide_violin_phos_protn <- renderUI({ 
      if(input$peptide_violin_phos_protn){
        generate_peptide_violin <- plot_abundance_distribution(proteome_data = db_execution_phos_protn$normalized_data,
                                                               type = "peptide")
        db_execution_phos_protn$peptide_abundance_distirbution = generate_peptide_violin$plot
        tagList(
          tags$h3("Distribution peptide abundace"),
          renderPlot(generate_peptide_violin$plot)
        )
      } else{
        db_execution_phos_protn$peptide_abundance_distirbution = NULL
      } 
    })
    
    output$render_mds_protein_phos_protn <- renderUI({ 
      if(input$mds_protein_phos_protn){
        res_plot <- mds_plot(proteome_data = db_execution_phos_protn$normalized_data,
                             type = "protein")
        db_execution_phos_protn$protein_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_phos_protn$protein_MDS = NULL
      } 
    })
    
    output$render_mds_peptide_phos_protn <- renderUI({ 
      if(input$mds_peptide_phos_protn){
        res_plot <- mds_plot(proteome_data = db_execution_phos_protn$normalized_data,
                             type = "peptide")
        db_execution_phos_protn$peptide_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on peptides"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_phos_protn$peptide_MDS = NULL
      } 
    })
    
    output$render_pca_protein_phos_protn <- renderUI({ 
      if(input$pca_protein_phos_protn){
        res_plot <- pca_plot(proteome_data = db_execution_phos_protn$normalized_data,
                             type = "protein")
        db_execution_phos_protn$protein_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_phos_protn$protein_PCA = NULL
      } 
    })
    
    output$render_pca_peptide_phos_protn <- renderUI({ 
      if(input$pca_peptide_phos_protn){
        res_plot <- pca_plot(proteome_data = db_execution_phos_protn$normalized_data,
                             type = "peptide")
        db_execution_phos_protn$peptide_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on peptides"),
          renderPlot(res_plot$plot)
        )
      }  else{
        db_execution_phos_protn$peptide_PCA = NULL
      }
    })
    
    output$render_protein_boxplot_phos_protn <- renderUI({
      if(input$boxplot_protein_phos_protn){
        req(input$list_proteins_phos_protn)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos_protn), regex=",")
        prot_boxplot <- plot_selected_proteins(proteome_data = db_execution_phos_protn$normalized_data,
                                               list_protein = unlist(list_proteins))
        db_execution_phos_protn$protein_boxplot = prot_boxplot$plot
        
        tagList(
          tags$h3("Boxplot selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      } else{
        db_execution_phos_protn$protein_boxplot = NULL
      }
    })
    
    output$render_protein_heatmap_phos_protn <- renderUI({
      if(input$heatmap_protein_phos_protn){
        req(input$list_proteins_phos_protn)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_phos_protn), regex=",")
        prot_boxplot <- heatmap_selected_proteins(proteome_data = db_execution_phos_protn$normalized_data,
                                                  list_protein = unlist(list_proteins))
        db_execution_phos_protn$protein_heatmap = prot_boxplot$plot
        
        tagList(
          tags$h3("Heatmap selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      } else{
        db_execution_phos_protn$protein_heatmap = NULL
      }
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
      if(input$peptide_diff_barplot_phos_protn){
        ploft_diff_number_pep <- generate_differential_barplots(db_execution_phos_protn$differential_results,
                                                                data_type="peptide")
        db_execution_phos_protn$peptide_differential_barplot = ploft_diff_number_pep$plot
        tagList(
          tags$h3("N° differential phospho-peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos_protn$peptide_differential_barplot = NULL
      }
    })
    
    output$render_peptide_vulcano_phos_protn <- renderUI({
      if(input$peptide_vulcano_phos_protn){
        generate_volcano_plots_peptide <- list()
        for(comp in names(db_execution_phos_protn$formule_contrast)){
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
        tabs_pep_vulcano <- list()
        for (i in seq_along(generate_volcano_plots_peptide)) {
          plot_id <- paste0(names(generate_volcano_plots_peptide)[i], "_pep_phos_protn")
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
          tags$h3("Vulcano Plot differential phospho-peptides"),
          do.call(tabsetPanel, c(list(id = "dynamic_tabs_vulcano_peptide_phos_protn"), tabs_pep_vulcano))
        )
      } else{
        db_execution_phos_protn$peptide_vulcano = NULL
      }
    })
    
    output$render_mds_peptide_diff_phos_protn <- renderUI({
      if(input$mds_diff_peptide_phos_protn){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution_phos_protn$differential_results,
                                                                proteome_data = db_execution_phos_protn$normalized_data,
                                                                type = "peptide")
        db_execution_phos_protn$peptide_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential phospho-peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos_protn$peptide_differential_MDS = NULL
      }
    })
    
    output$render_pca_peptide_diff_phos_protn <- renderUI({
      if(input$pca_diff_peptide_phos_protn){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution_phos_protn$differential_results,
                                                                proteome_data = db_execution_phos_protn$normalized_data,
                                                                type = "peptide")
        db_execution_phos_protn$peptide_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential phospho-peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      } else{
        db_execution_phos_protn$peptide_differential_PCA = NULL
      }
    })
    
  })
  
  ## PhosProTN_with_prot: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn_phos_protn, {
    output$render_enrichement_analysis_phos_protn <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution_phos_protn$enrichmnent_results <- perform_enrichment_analysis(differential_results = db_execution_phos_protn$differential_results,
                                                                        enrichR_custom_DB = T,
                                                                        enrich_filter_DBs=input$DB_enrichment_phos_protn,    
                                                                        overlap_size_enrich_thr=as.double(input$FC_thr_phos_protn),
                                                                        pval_fdr_enrich = input$pval_fdr_phos_protn,
                                                                        pval_enrich_thr=as.double(input$pval_thr_phos_protn),
                                                                        dirOutput=db_execution_phos_protn$dirOutput, 
                                                                        with_background = input$enrich_with_background_phos_protn)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich_phos_protn), regex=","))
        plots_down <- enrichment_figure(enr_df = db_execution_phos_protn$enrichmnent_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution_phos_protn$enrichmnent_results$anno_class), input$DB_enrichment_phos_protn))})
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
          withProgress(message = "Prepraring files to download, please wait!", {
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
            } 
            
            if(input$abundance_plot_phos_protn & !is.null(db_execution_phos_protn$generate_abundance)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/missing_available_abundance_proteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_abundance$proteome_plot, 
                     create.dir = T, width = 7, height = 5)
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/missing_available_abundance_phosphoproteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_abundance$phospho_plot, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.25)
            
            if(input$peptide_distribution_phos_protn & !is.null(db_execution_phos_protn$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_per_protein_proteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_peptide_distribution$proteome_plot, 
                     create.dir = T, width = 7, height = 5)
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_per_protein_phosphoproteomics.pdf"), 
                     plot = db_execution_phos_protn$generate_peptide_distribution$phospho_plot, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.30)
            
            if(input$protein_violin_phos_protn & !is.null(db_execution_phos_protn$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution_phos_protn$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.35)
            
            if(input$peptide_violin_phos_protn & !is.null(db_execution_phos_protn$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution_phos_protn$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.40)
            
            if(input$mds_protein_phos_protn & !is.null(db_execution_phos_protn$protein_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution_phos_protn$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.43)
            
            if(input$mds_peptide_phos_protn & !is.null(db_execution_phos_protn$peptide_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution_phos_protn$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.45)
            
            if(input$pca_protein_phos_protn & !is.null(db_execution_phos_protn$protein_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution_phos_protn$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.47)
            
            if(input$pca_peptide_phos_protn & !is.null(db_execution_phos_protn$peptide_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution_phos_protn$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein_phos_protn & !is.null(db_execution_phos_protn$protein_boxplot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution_phos_protn$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein_phos_protn & !is.null(db_execution_phos_protn$protein_heatmap)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution_phos_protn$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.55)
            
            if(!is.null(db_execution_phos_protn$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution_phos_protn$protein_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.58)
            
            if(!is.null(db_execution_phos_protn$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution_phos_protn$peptide_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.60)
            
            if(!is.null(db_execution_phos_protn$protein_vulcano)){
              dir.create(file.path(paste0(db_execution_phos_protn$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos_protn$protein_vulcano)){
                plotly::save_image(db_execution_phos_protn$protein_vulcano[[comp]], 
                                   file = paste0(db_execution_phos_protn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos_protn$protein_vulcano[[comp]], 
                                        file = paste0(db_execution_phos_protn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_phos_protn$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution_phos_protn$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_phos_protn$peptide_vulcano)){
                plotly::save_image(db_execution_phos_protn$peptide_vulcano[[comp]], 
                                   file = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_phos_protn$peptide_vulcano[[comp]], 
                                        file = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.68)
            
            if(!is.null(db_execution_phos_protn$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution_phos_protn$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.69)
            
            if(!is.null(db_execution_phos_protn$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution_phos_protn$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.70)
            
            if(!is.null(db_execution_phos_protn$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution_phos_protn$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.72)
            
            if(!is.null(db_execution_phos_protn$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution_phos_protn$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution_phos_protn$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.75)
            
            if(length(db_execution_phos_protn$enrichmnent_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich_phos_protn), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution_phos_protn$enrichmnent_results,
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
            
            # Save RData db_execution_phos_protn
            db_results_PhosProTN_with_proteome = reactiveValuesToList(db_execution_phos_protn)
            db_results_PhosProTN_with_proteome <- db_results_PhosProTN_with_proteome[!(unlist(lapply(db_results_PhosProTN_with_proteome, is.null)))]
            save(db_results_PhosProTN_with_proteome, file = paste0(db_results_PhosProTN_with_proteome$dirOutput,"db_results_PhosProTN_with_proteome.RData"))
            
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
    } else if(input$sw_analyzer_interactn == "MQ"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome_interactn", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome_interactn", "Select the EVIDENCE file of the PROTEOMICS..."),
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
        numericInput("min_peptide_protein_interactn", "Minimum peptide per protein", value = 1, min = 1)
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
    }
  })
  
  # InteracTN: Execution pipeline ----
  observeEvent(input$report_proteome_interactn, {
    
    output$protn_results_ui_interactn <- renderUI({
      isolate({
        tryCatch(
          {
            withProgress(message = "Rendering, please wait!", {
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
              file_prot_proteome = if(software=="PD"){input$prot_file_proteome_interactn$name}else{NA}
              file_pep_proteome = input$pep_file_proteome_interactn$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_protn"), showWarnings = FALSE)
              dir_input <- paste(dirOutput_Server, "input_protn", sep = "")
              file.copy(from = input$input_file_proteome_interactn$datapath, to = paste0(dir_input,'/ANNOTATION_',file_input_proteome)) 
              if(software=="PD"){file.copy(from = input$prot_file_proteome_interactn$datapath, to =paste0(dir_input,'/PROT_',file_prot_proteome))} 
              file.copy(from = input$pep_file_proteome_interactn$datapath, to = paste0(dir_input,'/PEP_',file_pep_proteome)) 
              
              # If advance filter
              if(input$advance_filter_interactn){
                NA_allow_condition <- input$NA_allow_condition_interactn
                min_peptide_protein <- input$min_peptide_protein_interactn
              } else{
                NA_allow_condition <- 0
                min_peptide_protein <- 1
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
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col, 
                                                                  filt_absent_value = NA_allow_condition, 
                                                                  min_peptide_protein = min_peptide_protein)
                  } else if(software == "MQ"){
                    db_execution_interactn$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
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
              db_execution_interactn$imputed_data <- impute_intensity(proteome_data = db_execution_interactn$proteome_data)
              db_execution_interactn$normalized_data <- normalization_ProTN(proteome_data = db_execution_interactn$imputed_data)
              if(batch_corr){
                message("Executing batch correction...")
                db_execution_interactn$normalized_data <- batch_correction(proteome_data = db_execution_interactn$normalized_data, 
                                                                 batch_col = str_to_lower(batch_correction_col))
              }
              
              output$c_anno_interactn <- DT::renderDT(db_execution_interactn$proteome_data$c_anno)
              tagList(
                fluidRow(
                  downloadButton("download_proteome_interactn", "Download results (ZIP file)", width = "240px")
                ),
                # html(html = paste0("<p>",msg_read_function,"</p><br>"), id = "messagge_read"),
                # shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE),
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
      if(input$abundance_plot_interactn){
        generate_abundance <- generate_abundance_plot(proteome_data = db_execution_interactn$proteome_data)
        db_execution_interactn$generate_abundance = generate_abundance$plot
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          renderPlot(generate_abundance$plot)
        )
      } else{
        db_execution_interactn$generate_abundance = NULL
      }
    })
    
    output$render_peptide_distribution_interactn <- renderUI({ 
      if(input$peptide_distribution_interactn){
        generate_peptide_distribution <- generate_peptide_distribution_plot(proteome_data = db_execution_interactn$proteome_data)
        db_execution_interactn$generate_peptide_distribution = generate_peptide_distribution$plot
        tagList(
          tags$h3("N° peptides per proteins"),
          renderPlot(generate_peptide_distribution$plot)
        )
      } else{
        db_execution_interactn$generate_peptide_distribution = NULL
      }
    })
    
    output$render_protein_violin_interactn <- renderUI({ 
      if(input$protein_violin_interactn){
        generate_protein_violin <- plot_abundance_distribution(proteome_data = db_execution_interactn$normalized_data,
                                                               type = "protein")
        db_execution_interactn$protein_abundance_distribution = generate_protein_violin$plot
        tagList(
          tags$h3("Distribution protein abundance"),
          renderPlot(generate_protein_violin$plot)
        )
      } else{
        db_execution_interactn$protein_abundance_distribution = NULL
      }
    })
    
    output$render_peptide_violin_interactn <- renderUI({ 
      if(input$peptide_violin_interactn){
        generate_peptide_violin <- plot_abundance_distribution(proteome_data = db_execution_interactn$normalized_data,
                                                               type = "peptide")
        db_execution_interactn$peptide_abundance_distirbution = generate_peptide_violin$plot
        tagList(
          tags$h3("Distribution peptide abundace"),
          renderPlot(generate_peptide_violin$plot)
        )
      } else{
        db_execution_interactn$peptide_abundance_distirbution = NULL
      }
    })
    
    output$render_mds_protein_interactn <- renderUI({ 
      if(input$mds_protein_interactn){
        res_plot <- mds_plot(proteome_data = db_execution_interactn$normalized_data,
                             type = "protein")
        db_execution_interactn$protein_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_interactn$protein_MDS = NULL
      }
    })
    
    output$render_mds_peptide_interactn <- renderUI({ 
      if(input$mds_peptide_interactn){
        res_plot <- mds_plot(proteome_data = db_execution_interactn$normalized_data,
                             type = "peptide")
        db_execution_interactn$peptide_MDS = res_plot$plot
        tagList(
          tags$h3("MDS based on peptides"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_interactn$peptide_MDS = NULL
      }
    })
    
    output$render_pca_protein_interactn <- renderUI({ 
      if(input$pca_protein_interactn){
        res_plot <- pca_plot(proteome_data = db_execution_interactn$normalized_data,
                             type = "protein")
        db_execution_interactn$protein_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on proteins"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_interactn$protein_PCA = NULL
      }
    })
    
    output$render_pca_peptide_interactn <- renderUI({ 
      if(input$pca_peptide_interactn){
        res_plot <- pca_plot(proteome_data = db_execution_interactn$normalized_data,
                             type = "peptide")
        db_execution_interactn$peptide_PCA = res_plot$plot
        tagList(
          tags$h3("PCA based on peptides"),
          renderPlot(res_plot$plot)
        )
      } else{
        db_execution_interactn$peptide_PCA = NULL
      }
    })
    
    output$render_protein_boxplot_interactn <- renderUI({
      if(input$boxplot_protein_interactn){
        req(input$list_proteins_interactn)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_interactn), regex=",")
        prot_boxplot <- plot_selected_proteins(proteome_data = db_execution_interactn$normalized_data,
                                               list_protein = unlist(list_proteins))
        db_execution_interactn$protein_boxplot = prot_boxplot$plot
        
        tagList(
          tags$h3("Boxplot selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      }else{
        db_execution_interactn$protein_boxplot = NULL
      }
    })
    
    output$render_protein_heatmap_interactn <- renderUI({
      if(input$heatmap_protein_interactn){
        req(input$list_proteins_interactn)
        list_proteins <- stri_split(stri_replace_all(regex = " ",replacement = "",str = input$list_proteins_interactn), regex=",")
        prot_boxplot <- heatmap_selected_proteins(proteome_data = db_execution_interactn$normalized_data,
                                                  list_protein = unlist(list_proteins))
        db_execution_interactn$protein_heatmap = prot_boxplot$plot
        
        tagList(
          tags$h3("Heatmap selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      }else{
        db_execution_interactn$protein_heatmap = NULL
      }
    })
  })
  
  ## InteracTN: differential analysis ----
  observeEvent(input$execute_differential_analysis_btn_interactn, {
    output$render_differential_analysis_interactn <- renderUI({
      isolate({
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
      if(input$protein_diff_barplot_interactn){
        ploft_diff_number <- generate_differential_barplots(db_execution_interactn$differential_results,
                                                            data_type="protein")
        db_execution_interactn$protein_differential_barplot = ploft_diff_number$plot
        tagList(
          tags$h3("N° differential proteins"),
          renderPlot(ploft_diff_number$plot)
        )
      }else{
        db_execution_interactn$protein_differential_barplot = NULL
      }
    })
    
    output$render_peptide_diff_barplot_interactn <- renderUI({
      if(input$peptide_diff_barplot_interactn){
        ploft_diff_number_pep <- generate_differential_barplots(db_execution_interactn$differential_results,
                                                                data_type="peptide")
        db_execution_interactn$peptide_differential_barplot = ploft_diff_number_pep$plot
        tagList(
          tags$h3("N° differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution_interactn$peptide_differential_barplot = NULL
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
      if(input$mds_diff_protein_interactn){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                                proteome_data = db_execution_interactn$normalized_data,
                                                                type = "protein")
        db_execution_interactn$protein_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution_interactn$protein_differential_MDS = NULL
      }
    })
    
    output$render_mds_peptide_diff_interactn <- renderUI({
      if(input$mds_diff_peptide_interactn){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                                proteome_data = db_execution_interactn$normalized_data,
                                                                type = "peptide")
        db_execution_interactn$peptide_differential_MDS = ploft_diff_number_pep$plot
        tagList(
          tags$h3("MDS based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution_interactn$peptide_differential_MDS = NULL
      }
    })
    
    output$render_pca_protein_diff_interactn <- renderUI({
      if(input$pca_diff_protein_interactn){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                                proteome_data = db_execution_interactn$normalized_data,
                                                                type = "protein")
        db_execution_interactn$protein_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution_interactn$protein_differential_PCA = NULL
      }
    })
    
    output$render_pca_peptide_diff_interactn <- renderUI({
      if(input$pca_diff_peptide_interactn){
        ploft_diff_number_pep <- pca_differential_analysis_plot(differential_analysis = db_execution_interactn$differential_results,
                                                                proteome_data = db_execution_interactn$normalized_data,
                                                                type = "peptide")
        db_execution_interactn$peptide_differential_PCA = ploft_diff_number_pep$plot
        tagList(
          tags$h3("PCA based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }else{
        db_execution_interactn$peptide_differential_PCA = NULL
      }
    })
    
  })
  
  ## InteracTN: enrichment analysis ----
  observeEvent(input$execute_enrichment_analysis_btn_interactn, {
    output$render_enrichement_analysis_interactn <- renderUI({
      isolate({
        # TODO: gallery of plots
        db_execution_interactn$enrichmnent_results <- perform_enrichment_analysis(differential_results = db_execution_interactn$differential_results,
                                                                        enrichR_custom_DB = T,
                                                                        enrich_filter_DBs=input$DB_enrichment_interactn,    
                                                                        overlap_size_enrich_thr=as.double(input$FC_thr_interactn),
                                                                        pval_fdr_enrich = input$pval_fdr_interactn,
                                                                        pval_enrich_thr=as.double(input$pval_thr_interactn),
                                                                        dirOutput=db_execution_interactn$dirOutput, 
                                                                        with_background = input$enrich_with_background_interactn)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",str = input$terms_enrich_interactn), regex=","))
        plots_down <- enrichment_figure(enr_df = db_execution_interactn$enrichmnent_results,
                                        category = c("down","up"), 
                                        enrich_filter_term = terms_enrich,
                                        save=F)
        
        #LOAD category EnrichR
        dbs_default <- read_tsv("data/dbs_enrichR.txt", col_names = FALSE) %>% as.data.frame()
        dbs_category <- dbs_default %>% split(f = as.factor(.$X2))
        category_db <- lapply(dbs_category, function(x){filter(x, x[,1] %in% intersect(unique(db_execution_interactn$enrichmnent_results$anno_class), input$DB_enrichment_interactn))})
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
          withProgress(message = "Prepraring files to download, please wait!", {
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
            } 
            setProgress(value = 0.25)
            
            if(input$peptide_distribution_interactn & !is.null(db_execution_interactn$generate_peptide_distribution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_per_protein.pdf"), 
                     plot = db_execution_interactn$generate_peptide_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.30)
            
            if(input$protein_violin_interactn & !is.null(db_execution_interactn$protein_abundance_distribution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_abundance_distribution.pdf"), 
                     plot = db_execution_interactn$protein_abundance_distribution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.35)
            
            if(input$peptide_violin_interactn & !is.null(db_execution_interactn$peptide_abundance_distirbution)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_abundance_distribution.pdf"), 
                     plot = db_execution_interactn$peptide_abundance_distirbution, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.40)
            
            if(input$mds_protein_interactn & !is.null(db_execution_interactn$protein_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_MDS.pdf"), 
                     plot = db_execution_interactn$protein_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.43)
            
            if(input$mds_peptide_interactn & !is.null(db_execution_interactn$peptide_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_MDS.pdf"), 
                     plot = db_execution_interactn$peptide_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.45)
            
            if(input$pca_protein_interactn & !is.null(db_execution_interactn$protein_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_PCA.pdf"), 
                     plot = db_execution_interactn$protein_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.47)
            
            if(input$pca_peptide_interactn & !is.null(db_execution_interactn$peptide_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_PCA.pdf"), 
                     plot = db_execution_interactn$peptide_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.50)
            
            # TODO: adapt based on number of protein
            if(input$boxplot_protein_interactn & !is.null(db_execution_interactn$protein_boxplot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_boxplot.pdf"), 
                     plot = db_execution_interactn$protein_boxplot, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.52)
            
            # TODO: adapt based on number of protein
            if(input$heatmap_protein_interactn & !is.null(db_execution_interactn$protein_heatmap)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_heatmap.pdf"), 
                     plot = db_execution_interactn$protein_heatmap, 
                     create.dir = T, width = 8, height = 7)
            } 
            setProgress(value = 0.55)
            
            if(!is.null(db_execution_interactn$protein_differential_barplot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_differential_barplot.pdf"), 
                     plot = db_execution_interactn$protein_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.58)
            
            if(!is.null(db_execution_interactn$peptide_differential_barplot)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_differential_barplot.pdf"), 
                     plot = db_execution_interactn$peptide_differential_barplot, 
                     create.dir = T, width = 8, height = 4)
            } 
            setProgress(value = 0.60)
            
            if(!is.null(db_execution_interactn$protein_vulcano)){
              dir.create(file.path(paste0(db_execution_interactn$dirOutput,"pics/"), "protein_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_interactn$protein_vulcano)){
                plotly::save_image(db_execution_interactn$protein_vulcano[[comp]], 
                                   file = paste0(db_execution_interactn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_interactn$protein_vulcano[[comp]], 
                                        file = paste0(db_execution_interactn$dirOutput,"pics/protein_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.64)
            
            if(!is.null(db_execution_interactn$peptide_vulcano)){
              dir.create(file.path(paste0(db_execution_interactn$dirOutput,"pics/"), "peptide_vulcano"), showWarnings = FALSE)
              for(comp in names(db_execution_interactn$peptide_vulcano)){
                plotly::save_image(db_execution_interactn$peptide_vulcano[[comp]], 
                                   file = paste0(db_execution_interactn$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.png"))
                htmlwidgets::saveWidget(db_execution_interactn$peptide_vulcano[[comp]], 
                                        file = paste0(db_execution_interactn$dirOutput,"pics/peptide_vulcano/",comp,"_protein_vulcano.html"))
              }
            } 
            setProgress(value = 0.68)
            
            if(!is.null(db_execution_interactn$protein_differential_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_differential_MDS.pdf"), 
                     plot = db_execution_interactn$protein_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.69)
            
            if(!is.null(db_execution_interactn$peptide_differential_MDS)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_differential_MDS.pdf"), 
                     plot = db_execution_interactn$peptide_differential_MDS, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.70)
            
            if(!is.null(db_execution_interactn$protein_differential_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/protein_differential_PCA.pdf"), 
                     plot = db_execution_interactn$protein_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.72)
            
            if(!is.null(db_execution_interactn$peptide_differential_PCA)){
              ggsave(filename = paste0(db_execution_interactn$dirOutput,"pics/peptide_differential_PCA.pdf"), 
                     plot = db_execution_interactn$peptide_differential_PCA, 
                     create.dir = T, width = 7, height = 5)
            } 
            setProgress(value = 0.75)
            
            if(length(db_execution_interactn$enrichmnent_results)>0){
              terms_enrich <- unlist(stri_split(stri_replace_all(regex = "\"|;|.",replacement = "",
                                                                 str = input$terms_enrich_interactn), regex=","))
              plots_down <- enrichment_figure(enr_df = db_execution_interactn$enrichmnent_results,
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
            
            # Save RData db_execution_interactn
            db_results_interacTN = reactiveValuesToList(db_execution_interactn)
            db_results_interacTN <- db_results_interacTN[!(unlist(lapply(db_results_interacTN, is.null)))]
            save(db_results_interacTN, file = paste0(db_results_interacTN$dirOutput,"db_results_InteracTN.RData"))
            
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
  
  ##############################################################################
  # ----
  # -- DELETE TEMP FILES WHEN SESSION ENDS -- #
  # session$onSessionEnded(function() {
  #   if (dir.exists(tempdir())){unlink(list.files(tempdir(), full.names = T), recursive = T)}
  # })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 8100))
