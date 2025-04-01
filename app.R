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
                  uiOutput("help1")
                ),
                fluidRow(
                  textAreaInput("description_exp", "Brief description", rows = 4),
                  uiOutput("help2")
                ),
                fluidRow(
                  radioButtons("sw_analyzer", "Software Analyzer", 
                               choiceNames = c("ProteomeDiscoverer", "MaxQuant", "TMT_PD"),
                               choiceValues = c("PD","MQ","TMT_PD"), inline = TRUE),
                  uiOutput("help3")
                ),
                uiOutput("input_proteome"),
                tags$h3("Select what execute:"),
                checkboxInput("abundance_plot", "% missing values", TRUE),
                checkboxInput("peptide_distribution", "N° peptides per protein", TRUE),
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
                checkboxInput("differential_analysis", "Execute differential analysis", FALSE), 
                #TODO: all parameter and figure differential analysis
                tags$h3("Enrichment Analysis:"),
                checkboxInput("enrichment_analysis", "Execute enrichment analysis", FALSE), 
                #TODO: show parameter
                tags$h3("STRINGdb network:"),
                checkboxInput("stringdb_analysis", "Execute STRINGdb", FALSE), 
              ),
              column(
                width = 8,
                fluidRow(
                  actionButton("report_proteome", "Generate report"),
                  actionButton("case_study_proteome", "Case Study Example")
                ),
                textOutput("messagge_read"),
                uiOutput("protn_results_ui"),
                uiOutput("render_abundance_plot"),
                uiOutput("render_peptide_distribution"),
                uiOutput("render_mds_protein"),
                uiOutput("render_mds_peptide"),
                uiOutput("render_pca_protein"),
                uiOutput("render_pca_peptide")
                
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
                                 proteome_data = list(),
                                 imputed_data = list(),
                                 normalized_data = list())
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
                shinyjs::html(id = "messagge_read", html = paste0(m$message,"\n"), add = TRUE)
              }
            )
            
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
  })
  
  
  # ----
  # -- DELETE TEMP FILES WHEN SESSION ENDS -- #
  # session$onSessionEnded(function() {
  #   if (dir.exists(tempdir())){unlink(list.files(tempdir(), full.names = T), recursive = T)}
  # })
}

# Run the application
shinyApp(ui = ui, server = server, options = list(port = 8100))
