# ProTN: an integrative pipeline for complete analysis of proteomics           # 
# data from mass spectrometry                                                  #
# Laboratory of RNA and Disease Data Science, University of Trento             #
# Developer: Gabriele Tomè                                                     #
# Issue at: https://github.com/TebaldiLab/ProTN/issues                         #
# PI: Dr. Toma Tebaldi, PhD                                                    #
list.of.packages <- c("shiny","tidyverse","markdown","knitr","shinydashboard",
                      "shinydashboardPlus","shinymaterial","shinyjs","magrittr",
                      "dplyr","stringr","shinyBS","proTN","DT")
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

tmpdir<-stri_replace_all_regex(tempdir(), pattern = "/Rtmp\\w+", replacement = "")
if (!dir.exists(file.path(tmpdir, "ProTN_shiny"))){
  dir.create(file.path(tmpdir, "ProTN_shiny"), showWarnings = FALSE)
}
Sys.setenv(TMPDIR=file.path(tmpdir, "ProTN_shiny"))
unlink(tempdir(), recursive = T)
tempdir(check = T)
#Javascript code to disable click when is running the app
jsCode <- "shinyjs.pageDisable = function(params){
              $('body').css('pointer-events', params);
            };"

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
      tags$link(rel="stylesheet", href="https://fonts.googleapis.com/css?family=El+Messiri"),
      includeCSS("www/css/custom_theme.css"),
      includeCSS("www/css/materialize.css"),
      includeScript("www/js/materialize.js"),
      extendShinyjs(text = jsCode, functions = c("pageDisable")),
      
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
                width = 4,
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
                tags$h3("Select what execute:"),
                checkboxInput("abundance_plot", "% missing values", TRUE),
                checkboxInput("peptide_distribution", "N° peptides per protein", TRUE),
                checkboxInput("protein_violin", "Distribution abundance proteins", FALSE),
                checkboxInput("peptide_violin", "Distribution abundance peptides", FALSE),
                checkboxInput("batch_correction", "Batch Correction", FALSE),
                uiOutput("batch_correction_ui"),
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
                width = 8,
                fluidRow(
                  actionButton("report_proteome", "Generate report"),
                  actionButton("case_study_proteome", "Case Study Example")
                ),
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
                uiOutput("render_enrichement_analysis")
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
                                 differential_results = list(),
                                 enrichmnent_results = list())
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
        textAreaInput("formule_contrast", "Write in each line a different comparison", rows = 4),
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
        #TODO: show parameter
        tags$h3("STRINGdb network:"),
        checkboxInput("stringdb_analysis", "Execute STRINGdb", FALSE)
      )
    } 
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
        actionButton("execute_enrichment_analysis_btn", "Run!")
      )
    } 
  })
  
  # PROTN: Execution pipeline ----
  observeEvent(input$report_proteome, {
    
    output$protn_results_ui <- renderUI({
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
            
            message(software)
            withCallingHandlers(
              {
                shinyjs::html("text", "")
                if(software == "PD"){
                  db_execution$proteome_data <- read_proteomics(software = "PD",
                                                                folder = dir_input,
                                                                peptide_filename = "PEP_",
                                                                annotation_filename = "ANNOTATION_",
                                                                proteinGroup_filename = "PROT_")
                } else if(software == "MQ"){
                  db_execution$proteome_data <- read_proteomics(software = "MQ",
                                                                folder = dir_input,
                                                                peptide_filename = "PEP_",
                                                                annotation_filename = "ANNOTATION_")
                }
              },
              message = function(m) {
                shinyjs::html(id = "messagge_read", html = paste0("<p>",m$message,"</p>"), add = TRUE)
              }
            )
            
            db_execution$data_loaded <- TRUE
            db_execution$imputed_data <- impute_intensity(proteome_data = db_execution$proteome_data)
            db_execution$normalized_data <- normalization_ProTN(proteome_data = db_execution$imputed_data)
            
            
            output$c_anno <- DT::renderDT(db_execution$proteome_data$c_anno)
            DT::DTOutput("c_anno")
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
    
    output$render_abundance_plot <- renderUI({ 
      if(input$abundance_plot){
        generate_abundance <- generate_abundance_plot(proteome_data = db_execution$proteome_data)
        tagList(
          tags$h3("Percentage missing values respect detected abundance"),
          renderPlot(generate_abundance$plot)
        )
      } 
    })
    
    output$render_peptide_distribution <- renderUI({ 
      if(input$peptide_distribution){
        generate_peptide_distribution <- generate_peptide_distribution_plot(proteome_data = db_execution$proteome_data)
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
        tagList(
          tags$h3("PCA based on peptides"),
          renderPlot(res_plot$plot)
        )
      } 
    })
    
    output$render_protein_boxplot <- renderUI({
      if(input$boxplot_protein){
        req(input$list_proteins)
        list_proteins <- stri_split(stri_replace_all(regex = " |\"|;|.",replacement = "",str = input$list_proteins), regex=",")
        prot_boxplot <- plot_selected_proteins(proteome_data = db_execution$normalized_data,
                                               list_protein = unlist(list_proteins))
        
        tagList(
          tags$h3("Boxplot selected proteins"),
          renderPlot(prot_boxplot$plot)
        )
      }
    })
    
    output$render_protein_heatmap <- renderUI({
      if(input$heatmap_protein){
        req(input$list_proteins)
        list_proteins <- stri_split(stri_replace_all(regex = " |\"|;|.",replacement = "",str = input$list_proteins), regex=",")
        prot_boxplot <- heatmap_selected_proteins(proteome_data = db_execution$normalized_data,
                                                  list_protein = unlist(list_proteins))
        
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
      # req(input$formule_contrast)
      isolate({
        formule_diff <- stri_split(input$formule_contrast, regex = "\n")[[1]]
        
        db_execution$formule_contrast <- as.list(formule_diff)
        names(db_execution$formule_contrast) <- stri_replace_all(formule_diff, replacement = "_VS_", regex = "-")
        # db_execution$formule_contrast <- list("A2016vsWT"="A2016T-wt")
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
        tagList(
          tags$h3("N° differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }
    })
    
    output$render_protein_vulcano <- renderUI({
      if(input$protein_vulcano){
        generate_volcano_plots_protein<-generate_volcano_plots(db_execution$differential_results,
                                                               data_type="protein",
                                                               comparison=names(db_execution$formule_contrast)[[1]],
                                                               fc_thr=as.double(input$FC_thr),
                                                               pval_fdr = input$pval_fdr,
                                                               pval_thr=as.double(input$pval_thr))
        tagList(
          tags$h3("Vulcano Plot differential proteins"),
          renderPlotly(generate_volcano_plots_protein[[names(db_execution$formule_contrast)[[1]]]])
        )
      }
    })
    
    output$render_peptide_vulcano <- renderUI({
      if(input$peptide_vulcano){
        generate_volcano_plots_peptide<-generate_volcano_plots(db_execution$differential_results,
                                                               data_type="peptide",
                                                               comparison=names(db_execution$formule_contrast)[[1]],
                                                               fc_thr=as.double(input$FC_thr),
                                                               pval_fdr = input$pval_fdr,
                                                               pval_thr=as.double(input$pval_thr))
        tagList(
          tags$h3("Vulcano Plot differential peptides"),
          renderPlotly(generate_volcano_plots_peptide[[names(db_execution$formule_contrast)[[1]]]])
        )
      }
    })
    
    output$render_mds_protein_diff <- renderUI({
      if(input$mds_diff_protein){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "protein")
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
        tagList(
          tags$h3("MDS based on differential peptides"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }
    })
    
    output$render_pca_protein_diff <- renderUI({
      if(input$pca_diff_protein){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "protein")
        tagList(
          tags$h3("PCA based on differential proteins"),
          renderPlot(ploft_diff_number_pep$plot)
        )
      }
    })
    
    output$render_pca_peptide_diff <- renderUI({
      if(input$pca_diff_peptide){
        ploft_diff_number_pep <- mds_differential_analysis_plot(differential_analysis = db_execution$differential_results,
                                                                proteome_data = db_execution$normalized_data,
                                                                type = "peptide")
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
                                                          dirOutput=db_execution$dirOutput)
        
        terms_enrich <- unlist(stri_split(stri_replace_all(regex = " |\"|;|.",replacement = "",str = input$terms_enrich), regex=","))
        plots_down <- enrichment_figure(enr_df = db_execution$enrichmnent_results,
                                        category = c("down","up"),
                                        save=F)
        navset_card_underline(
          nav_panel("1", renderPlot(plots_down[[1]])),
          nav_panel("2", renderPlot(plots_down[[2]]))
        )
        
      })
    })
  })
  # ----
  # -- DELETE TEMP FILES WHEN SESSION ENDS -- #
  # session$onSessionEnded(function() {
  #   if (dir.exists(tempdir())){unlink(list.files(tempdir(), full.names = T), recursive = T)}
  # })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 8100))
