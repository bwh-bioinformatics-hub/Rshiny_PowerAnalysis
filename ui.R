library(shiny)
library(shinyjs)

fixedPage(
    
    useShinyjs(),
    
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "jumbotron-narrow.css"),
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
        tags$link(rel = "shortcut icon", href = "http://www.humanbraincode.org/img/ico/favicon.png"),
        tags$script("src"="func.js")
    ),
    
    tags$div(class="header",
             tags$div(class="row",
                      tags$div(class="col-sm-2",
                               # tags$a(href="http://www.humanbraincode.org", 
                                      tags$img(src = "powerEQTL.logo.png", width = "180"), style="text-align: center;"),
                      tags$div(class="col-sm-10",
                               tags$p(class="title", style="line-height: 80%;", 
                                      HTML("Power<font color='red'>eQTL</font>")),
                               tags$p(class="title", style="font-size:18px;color:gray;", "An R package and R shiny application for calculating 
                                      sample size and power of bulk tissue and single-cell eQTL analysis"))
             )
    ),

    tabsetPanel(
        tabPanel("Power curves for tissue eQTL", 
                 sidebarLayout(
                     sidebarPanel(
                         
                         textInput(inputId = "myntot", label = "Sample sizes (separated by ',')", 
                                   value = "50, 100, 150, 200, 250, 300"),
                         numericInput(inputId = "sigma1", 
                                      label <- htmltools::doRenderTags(HTML(paste0("Standard deviation of gene expression (σ",tags$sub("y"),")"))),
                                      value = 0.29, min = 0.001, max = 1, step = 0.001),
                         numericInput(inputId = "delta1", label = " ", 
                                      value = 0.29*1.4, min = 0.001, max = 1, step = 0.001),
                         # actionButton("btn", "On hover", icon = icon("info")),
                         numericInput(inputId = "delta2", label = "Mean difference of gene expression (δ2)", 
                                      value = 0.29*1.4, min = 0.001, max = 1, step = 0.001),
                         numericInput(inputId = "alpha1", label = "Family-wise type I error rate (FWER)", 
                                      value = 0.05, min = 0.001, max = 1, step = 0.001),
                         numericInput(inputId = "nTest1", label = "Total number of tests (nTests)", 
                                      value = 10e5, min = 1, step = 1000),
                         
                         radioButtons(inputId = "radio", label = "Model",
                                      choices = c("One-way unbalanced ANOVA", "Simple linear regression")),
                         
                         actionButton("btn_t", "More Options"),
                         downloadButton('export2', "Download Report"),
                         tags$br(),
                         tags$br(),
                         hidden(
                             textInput(inputId = "MAF_t", "MAF range ('min, max')",
                                       value = "0.002, 0.5"),
                             selectInput(inputId = "pos_t", "Legend Position", choices = c("bottomright", "topleft"))
                         ),
                         
                         
                     ),
                     mainPanel(
                         
                         plotOutput(outputId = "tissue", hover = "t_plot_hover", width = "750px", height = "450px"),
                         
                         tags$br(),
                         htmlOutput("description2"),
                         
                         
                         
                     ))),
        tabPanel("Power curves for single-cell eQTL", 
                 sidebarLayout(
                     sidebarPanel(
                         
                         textInput(inputId = "myntot1", label = "Sample sizes (separated by ',')", 
                                   value = "50, 100, 150, 200, 250, 300"),
                         numericInput(inputId = "m", label = "Number of cells from each subject (m)", 
                                      value = 640, min = 0, step = 100),
                         numericInput(inputId = "delta", label = "Slope under alternative hypothesis (β1)", 
                                      value = 0.29*1.5, min = 0.001, max = 1, step = 0.001),
                         numericInput(inputId = "sigma", 
                                      label <- htmltools::doRenderTags(HTML(paste0("Standard deviation of gene expression (σ",tags$sub("y"),")"))),
                                      value = 0.29, min = 0.001, max = 1, step = 0.001),
                         numericInput(inputId = "rho", label = "Intra-class correlation (range [0, 1]; ρ)", 
                                      value = 0.5, min = 0, max = 1, step = 0.001),
                         numericInput(inputId = "alpha", label = "Family-wise type I error rate (FWER)", 
                                      value = 0.05, min = 0.001, max = 1, step = 0.001),
                         numericInput(inputId = "nTest", label = "Total number of tests (nTests)", 
                                      value = 10e5, min = 1, step = 1000),
                         
                         actionButton("btn_sc", "More Options"),
                         downloadButton('export', "Download Report"),
                         tags$br(),
                         tags$br(),
                         hidden(
                             textInput(inputId = "MAF_sc", "MAF range ('min, max')",
                                       value = "0.002, 0.5"),
                             selectInput(inputId = "pos", "Legend Position", choices = c("bottomright", "topleft"))
                         ),
                         
                     ),
                     mainPanel(
                         
                         plotOutput(outputId = "cell", hover = "sc_plot_hover", width = "750px", height = "450px"),
                         
                         htmlOutput("description"),
                         tags$br(),
                         
                         
                     )
                 )
        ),
        
        
        tabPanel("Calculator for tissue eQTL",
                 sidebarLayout(
                     sidebarPanel(
                         sliderInput(inputId = "power_test", label = "Desired power level",
                                     value = 0.8, min = 0.001, max = 1),
                         sliderInput(inputId = "n_test", label = "Number of subjects needed",
                                     value = 1114, min = 10, max = 10000),
                         sliderInput(inputId = "maf_test", label = "Minor allele frequencies (between 0 and 0.5)",
                                     value = 0.02, min = 0.001, max = 0.5),
                         
                         hidden(
                             sliderInput(inputId = "slope_test", label = "Slope of the regression line (β1)",
                                         value = 0.2, min = 0.01, max = 1)
                         ),
                         
                         radioButtons(inputId = "radio_test", label = "Model",
                                      choices = c("One-way unbalanced ANOVA", "Simple linear regression")),
                         
                         actionButton("btn_test", "More Options"),
                         tags$br(),
                         tags$br(),
                         hidden(
                             sliderInput(inputId = "delta1_test", label = "Mean difference of gene expression (δ1)",
                                         value = 0.2, min = 0.01, max = 10),
                             sliderInput(inputId = "delta2_test", label = "Mean difference of gene expression (δ2)",
                                         value = 0.2, min = 0.01, max = 10),
                             sliderInput(inputId = "sigma_test", 
                                         label <- htmltools::doRenderTags(HTML(paste0("Standard deviation of gene expression (σ",tags$sub("y"),")"))),
                                         value = 0.2, min = 0.01, max = 10),
                             sliderInput(inputId = "FWER_test", label = "Family-wise type I error rate (FWER)",
                                         value = 0.05,min = 0.01, max = 1),
                             sliderInput(inputId = "nTest_test", label = "Total number of tests (nTests)",
                                         value = 10e5, min = 0, max = 10e7)
                             
                         ),
                         
                     ),
                     mainPanel(
                         htmlOutput("title"),
                         tableOutput("approx"),
                         # htmlOutput("title1"),
                         # tableOutput("approx1"),
                         htmlOutput("Explanation1")
                     )
                 )
        ),
        tabPanel("Calculator for single-cell eQTL",
                 sidebarLayout(
                     sidebarPanel(
                         
                         sliderInput(inputId = "power_est", label = "Desired power level",
                                     value = 0.8, min = 0.001, max = 1),
                         sliderInput(inputId = "n_est", label = "Number of subjects needed",
                                     value = 507, min = 10, max = 1000),
                         sliderInput(inputId = "maf_est", label = "Minor allele frequencies (between 0 and 0.5)",
                                     value = 0.02, min = 0.001, max = 0.5),
                         sliderInput(inputId = "slope_est", label = "Slope under alternative hypothesis (β1)",
                                     value = 0.2, min = 0.01, max = 1),
                         
                         actionButton("btn_scest", "More Options"),
                         tags$br(),
                         tags$br(),
                         hidden(
                             sliderInput(inputId = "sigma_est",
                                         label <- htmltools::doRenderTags(HTML(paste0("Standard deviation of gene expression (σ",tags$sub("y"),")"))),
                                         value = 0.2, min = 0.01, max = 10),
                             sliderInput(inputId = "FWER_est", label = "Family-wise type I error rate (FWER)",
                                         value = 0.05,min = 0.01, max = 1),
                             sliderInput(inputId = "nTest_est", label = "Total number of tests (nTests)",
                                         value = 10e5, min = 0, max = 10e7),
                             sliderInput(inputId = "rho_est", label = "Intra-class correlation (ρ; range [0, 1])",
                                         value = 0.5, min = 0, max = 1),
                             sliderInput(inputId = "m_est", label = "Number of cells from each subject (m)",
                                         value = 640, min = 500, max = 100000)
                         ),
                         
                     ),
                     mainPanel(
                         htmlOutput("title1"),
                         tableOutput("approx1"),
                         htmlOutput("Explanation")
                     )
                 )
        ),
        tabPanel("About", value = "about",
                 htmlOutput("about")
        )
    ),
    
    
    
    # footer
    tags$div(class="footer",
             tags$div(class="row",
                      tags$div(class="col-sm-10",
                               tags$p(HTML(paste0("&copy; ", tags$a(href="https://bioinformatics.bwh.harvard.edu/", 
                                                                           "Genomics and Bioinformatics Hub "), 
                                                  format(Sys.Date(), "%Y"), "")))),
                      tags$div(class="col-sm-2",
                               tags$a(href="http://www.brighamandwomens.org", 
                                      tags$img(src = "http://www.humanbraincode.org/img/bwh.png", width="50")),
                               tags$a(href="https://hms.harvard.edu", 
                                      tags$img(src = "http://www.humanbraincode.org/img/hms.png", width="45")),
                               tags$a(href="http://www.massgeneral.org", 
                                      tags$img(src = "http://www.humanbraincode.org/img/mgh.png", width="48"))
                      )
             )
    )
)