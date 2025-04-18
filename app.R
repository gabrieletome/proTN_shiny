# ProTN: an integrative pipeline for complete analysis of proteomics           # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# Issue at: https://github.com/TebaldiLab/ProTN/issues                         #
# PI: Dr. Toma Tebaldi, PhD                                                    #
list.of.packages <- c("shiny","tidyverse","markdown","knitr","shinydashboard",
                      "shinydashboardPlus","shinymaterial","shinyjs","magrittr",
                      "dplyr","stringr","shinyBS","DT","bslib","readr",
                      "plotly","rhandsontable")
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
library(shinyjs)
library(magrittr)
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
        menuItem("Contacts", tabName = "contacts", icon = icon("comment", lib="font-awesome"))
      )
    ),
    #Main body web page
    body=dashboardBody(
      #Load associated file, CSS, JS, logo
      useShinyjs(),
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
                  radioButtons("sw_analyzer", "Software Analyzer", 
                               choiceNames = c("ProteomeDiscoverer", "MaxQuant", "TMT_PD"),
                               choiceValues = c("PD","MQ","TMT_PD"), inline = TRUE),
                ),
                uiOutput("input_proteome"),
                checkboxInput("batch_correction", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui"),
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
                textOutput("messagge_read"),
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
                                 generate_abundance = NULL,
                                 generate_peptide_distribution = NULL,
                                 protein_abundance_distribution = NULL,
                                 peptide_abundance_distirbution = NULL,
                                 protein_MDS = NULL,
                                 peptide_MDS = NULL,
                                 protein_PCA = NULL,
                                 peptide_PCA = NULL,
                                 protein_boxplot = NULL,
                                 protein_heatmap = NULL,
                                 protein_differential_barplot = NULL,
                                 peptide_differential_barplot = NULL,
                                 protein_vulcano = NULL,
                                 peptide_vulcano = NULL,
                                 protein_differential_MDS = NULL,
                                 peptide_differential_MDS = NULL,
                                 protein_differential_PCA = NULL,
                                 peptide_differential_PCA = NULL)
  
  # Optional visibility based on the selection ----
  
  ## PROTN: Visibility of the proteomics files for ProTN ----
  output$input_proteome <- renderUI({
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
    } else if(input$sw_analyzer == "MQ"){
      tagList(
        fluidRow(
          fileInput("input_file_proteome", "Select the SAMPLE_ANNOTATION file of the PROTEOMICS..."),
        ),
        fluidRow(
          fileInput("pep_file_proteome", "Select the EVIDENCE file of the PROTEOMICS..."),
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
        actionButton("execute_differential_analysis_btn", "Run!"),
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
              file_prot_proteome = if(software=="PD"){input$prot_file_proteome$name}else{NA}
              file_pep_proteome = input$pep_file_proteome$name
              
              # Move data in correct folder
              dir.create(file.path(dirOutput_Server, "input_protn"), showWarnings = FALSE)
              dir_input <- paste(dirOutput_Server, "input_protn", sep = "")
              file.copy(from = input$input_file_proteome$datapath, to = paste0(dir_input,'/ANNOTATION_',file_input_proteome)) 
              if(software=="PD"){file.copy(from = input$prot_file_proteome$datapath, to =paste0(dir_input,'/PROT_',file_prot_proteome))} 
              file.copy(from = input$pep_file_proteome$datapath, to = paste0(dir_input,'/PEP_',file_pep_proteome)) 
              
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
              msg_read_function <- c()
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
                                                                  batch_col = batch_correction_col)
                  } else if(software == "MQ"){
                    db_execution$proteome_data <- read_proteomics(software = "MQ",
                                                                  folder = dir_input,
                                                                  peptide_filename = "PEP_",
                                                                  annotation_filename = "ANNOTATION_", 
                                                                  batch_corr_exe = batch_corr, 
                                                                  batch_col = batch_correction_col)
                  }
                },
                message = function(m) {
                  msg_read_function <- c(msg_read_function, m$message)
                  shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE)
                  progress=progress+0.05
                  setProgress(value = progress)
                }
              )
              
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
          plot_id <- names(generate_volcano_plots_protein)[i]
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
          plot_id <- names(generate_volcano_plots_peptide)[i]
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
          # Create an output slot for each plot
          local({
            my_i <- i
            my_plot_id <- plot_id
            output[[my_plot_id]] <- renderPlot({
              plots_down[[names(plots_down)[my_i]]]
            }, height = max(min(20,uniqueN(db_execution$enrichmnent_results[anno_class %in% category_db[[my_plot_id]][,1],
                                                                            "anno_name"])*0.3),3)*85)
          })
          
          tabs[[i]] <- tabPanel(
            title = paste(names(plots_down)[i]),
            plotOutput(plot_id, height = max(min(20,uniqueN(db_execution$enrichmnent_results[anno_class %in% category_db[[names(plots_down)[i]]][,1],
                                                                                             "anno_name"])*0.3),3)*85)
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
    filename = "results_case_study.zip",
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
              readMQ_files = if (input$sw_analyzer == "MQ") {TRUE} else {FALSE},
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
  
  # ----
  # -- DELETE TEMP FILES WHEN SESSION ENDS -- #
  # session$onSessionEnded(function() {
  #   if (dir.exists(tempdir())){unlink(list.files(tempdir(), full.names = T), recursive = T)}
  # })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 8100))
