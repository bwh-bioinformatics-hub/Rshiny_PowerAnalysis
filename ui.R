library(shiny)
library(shinyjs)

fluidPage(
    
    useShinyjs(),
    
    titlePanel("Sample size and power calculator for bulk tissue and single-cell eQTL analysis", ),
    tabsetPanel(type = "tabs",
                tabPanel("Power for single-cell eQTL", 
                         sidebarLayout(
                             sidebarPanel(

                                 textInput(inputId = "myntot1", label = "Sample sizes (separated by ',')", 
                                           value = "50, 100, 150, 200, 250, 300"),
                                 numericInput(inputId = "m", label = "Number of cells from each subject (m)", 
                                              value = 640, min = 0),
                                 numericInput(inputId = "delta", label = "Slope under alternative hypothesis (β)", 
                                              value = 0.29*1.5, min = 0.001, max = 1),
                                 numericInput(inputId = "sigma", 
                                              label = "Standard deviation of gene expression levels in one group of subjects (σ)", 
                                              value = 0.29, min = 0.001, max = 1),
                                 numericInput(inputId = "rho", label = "Intra-class correlation (range [0, 1]; ρ)", 
                                              value = 0.5, min = 0, max = 1),
                                 numericInput(inputId = "alpha", label = "Family-wise type I error rate (FWER)", 
                                              value = 0.05, min = 0.001, max = 1),
                                 numericInput(inputId = "nTest", label = "Total number of tests (nTests)", 
                                              value = 10e5, min = 1),
                                 
                                 actionButton("btn_sc", "More Options"),
                                 tags$br(),
                                 tags$br(),
                                 hidden(
                                     textInput(inputId = "MAF_sc", "MAF range ('min, max')",
                                               value = "0.002, 0.5"),
                                     selectInput(inputId = "pos", "Legend Position", choices = c("bottomright", "topleft"))
                                 )
                             ),
                             mainPanel(
                                 
                                 plotOutput(outputId = "cell", hover = "sc_plot_hover", width = "800px", height = "500px"),

                                 downloadButton('export', "Download Report"),
                                 tags$br(),
                                 img(src='eQTL_demo_SLR.png', align = "center", width = "550px", height = "170px"),
                                 htmlOutput("description"),
                                 
                             )
                         )
                ),
                tabPanel("Power for tissue eQTL", 
                         sidebarLayout(
                             sidebarPanel(

                                 textInput(inputId = "myntot", label = "Sample sizes (separated by ',')", 
                                           value = "50, 100, 150, 200, 250, 300"),
                                 numericInput(inputId = "sigma1", 
                                              label = "Standard deviation of gene expression levels in one group of subjects (σ)", 
                                              value = 0.29, min = 0.001, max = 1),
                                 numericInput(inputId = "delta1", label = " ", 
                                              value = 0.29*1.4, min = 0.001, max = 1),
                                 numericInput(inputId = "delta2", label = "Mean difference of gene expression levels between genotype classes (δ2)", 
                                              value = 0.29*1.4, min = 0.001, max = 1),
                                 numericInput(inputId = "alpha1", label = "Family-wise type I error rate (FWER)", 
                                              value = 0.05, min = 0.001, max = 1),
                                 numericInput(inputId = "nTest1", label = "Total number of tests (nTests)", 
                                              value = 10e5, min = 1),
                                 
                                 radioButtons(inputId = "radio", label = "Model",
                                              choices = c("One-way unbalanced ANOVA", "Simple linear regression")),
                                 
                                 actionButton("btn_t", "More Options"),
                                 hidden(
                                     textInput(inputId = "MAF_t", "MAF range ('min, max')",
                                               value = "0.002, 0.5"),
                                     selectInput(inputId = "pos_t", "Legend Position", choices = c("bottomright", "topleft"))
                                 )

                             ),
                             mainPanel(
                                 
                                 plotOutput(outputId = "tissue", hover = "t_plot_hover", width = "80%", height = "500px"),
                                 downloadButton('export2', "Download Report"),
                                 tags$br(),
                                 uiOutput("demo"),
                                 # img(src='eQTL_demo_ANOVA.png', align = "middle", width = 500, height = 170),
                                 htmlOutput("description2"),
                                 
                                 
                                 
                             ))),
                tabPanel("Sample size estimation",
                         sidebarLayout(
                             sidebarPanel(
                                 
                                 sliderInput(inputId = "power_est", label = "Desired power level",
                                             value = 0.8, min = 0.001, max = 1),
                                 # sliderInput(inputId = "n_est", label = "Number_of_subjects_needed",
                                 #             value = 100, min = 10, max = 1e30),
                                 sliderInput(inputId = "maf_est", label = "Minor allele frequencies (between 0 and 0.5)",
                                              value = 0.02, min = 0.001, max = 0.5),
                                 
                                 sliderInput(inputId = "slope_est", label = "Mean difference of gene expression levels between genotype classes (δ)", 
                                              value = 0.2, min = 0.01, max = 10),
                                 sliderInput(inputId = "slope1_est", label = "Mean difference of gene expression levels between genotype classes (δ2; for ANOVA only)", 
                                             value = 0.2, min = 0.01, max = 10),
                                 sliderInput(inputId = "sigma_est", label = "Standard deviation of gene expression levels in one group of subjects (σ)", 
                                             value = 0.2, min = 0.01, max = 10),
                                 sliderInput(inputId = "FWER_est", label = "Family-wise type I error rate (FWER)", 
                                              value = 0.05,min = 0.01, max = 1),
                                 sliderInput(inputId = "nTest_est", label = "Total number of tests (nTests)", 
                                              value = 10e5, min = 0, max = 10e7),
                                 sliderInput(inputId = "rho_est", label = "Intra-class correlation (ρ; range [0, 1]; for single-cell only)", 
                                              value = 0.5, min = 0, max = 1),
                                 sliderInput(inputId = "m_est", label = "Number of cells from each subject (m; for single-cell only)", 
                                              value = 640, min = 500, max = 100000),

                                ),
                             mainPanel(
                                 htmlOutput("title"),
                                 tableOutput("approx"),
                                 htmlOutput("title1"),
                                 tableOutput("approx1"),
                                 htmlOutput("Explanation")
                                 )
                             )
                         ),
                tabPanel("About",
                         htmlOutput("authors")
                         )
    )
)

